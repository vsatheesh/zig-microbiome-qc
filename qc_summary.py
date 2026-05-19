#!/usr/bin/env python3
"""
qc_summary.py — Python reimplementation of qc_summary.zig

Parses MultiQC and DADA2 output TSVs for microbiome 16S amplicon runs,
joins them by sample ID, classifies samples by sequencing depth and process
metrics, flags quality issues, and emits a summary report.

Usage:
    qc_summary.py --multiqc-general-stats <file> --multiqc-fastqc <file> \
                  [--dada2-stats <file>] [--out-prefix <prefix>]
"""

import sys
import re
from pathlib import Path
from typing import Optional

import click
import pandas as pd

# Spec requires exit 1 for unknown/bad arguments; click defaults to 2.
click.exceptions.UsageError.exit_code = 1

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

PASS_STRONG_DEFAULT = 20_000
PASS_ACCEPTABLE_DEFAULT = 10_000
LOW_DEPTH_DEFAULT = 5_000
LOW_RETENTION_DEFAULT = 0.10
POOR_MERGING_DEFAULT = 0.60
HIGH_CHIMERA_DEFAULT = 0.30
VERY_LOW_FILTERED_DEFAULT = 0.20
HIGH_DUPLICATION_DEFAULT = 50.0
GC_OUTLIER_SD_DEFAULT = 2.0

# ---------------------------------------------------------------------------
# Sample ID normalization
# ---------------------------------------------------------------------------

_FILE_SUFFIXES = [
    "_fastqc.zip", "_fastqc", ".fastq.gz", ".fq.gz", ".fastq", ".fq", ".zip",
]

# Ordered longest → shortest so longest wins
_READ_SUFFIXES = [
    "_R1_001", "_R2_001", ".R1_001", ".R2_001",
    "_R1", "_R2", ".R1", ".R2",
    "_1", "_2",
]

_LANE_RE = re.compile(r"_L\d{3}$")


def normalize_sample_id(name: str) -> str:
    """Strip file-type and read-direction suffixes to get a canonical sample ID.

    Order matches the Zig implementation: file suffixes → read suffixes → lane suffix.
    """
    name = name.strip()
    changed = True
    while changed:
        changed = False
        for sfx in _FILE_SUFFIXES:
            if name.endswith(sfx):
                name = name[: -len(sfx)]
                changed = True
                break
    for sfx in _READ_SUFFIXES:
        if name.endswith(sfx):
            name = name[: -len(sfx)]
            break
    name = _LANE_RE.sub("", name)
    return name


# ---------------------------------------------------------------------------
# FastQC status helpers
# ---------------------------------------------------------------------------

_STATUS_MAP = {"pass": 0, "warn": 1, "fail": 2}
_STATUS_INT_TO_STR = {0: "pass", 1: "warn", 2: "fail"}

_FASTQC_STATUS_COLS = [
    "basic_statistics",
    "per_base_sequence_quality",
    "per_tile_sequence_quality",
    "per_sequence_quality_scores",
    "per_base_sequence_content",
    "per_sequence_gc_content",
    "per_base_n_content",
    "sequence_length_distribution",
    "sequence_duplication_levels",
    "overrepresented_sequences",
    "adapter_content",
]


def _parse_status(val) -> Optional[int]:
    if pd.isna(val):
        return None
    s = str(val).strip().lower()
    return _STATUS_MAP.get(s)


def _status_str(val) -> str:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return "N/A"
    try:
        return _STATUS_INT_TO_STR.get(int(val), "N/A")
    except (TypeError, ValueError):
        return "N/A"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_required_columns(df: pd.DataFrame, path: str, required: list[str]) -> None:
    """Print ALL missing columns, then exit 1 if any are missing."""
    missing = [c for c in required if c not in df.columns]
    if missing:
        for col in missing:
            print(f"ERROR: {path}: missing required column '{col}'", file=sys.stderr)
        raise SystemExit(1)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_general_stats(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    validate_required_columns(df, path, ["Sample"])
    return df


def parse_fastqc_stats(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    validate_required_columns(df, path, ["Sample"])
    return df


def parse_dada2_stats(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)

    # Detect and drop QIIME2 metadata directive rows like "#q2:types"
    first_col = df.columns[0]
    df = df[~df[first_col].str.startswith("#", na=False)].reset_index(drop=True)

    id_col = None
    for candidate in ("sample-id", "sample_id"):
        if candidate in df.columns:
            id_col = candidate
            break
    if id_col is None:
        print(
            f"ERROR: {path}: missing required DADA2 sample id column "
            f"(expected 'sample-id' or 'sample_id')",
            file=sys.stderr,
        )
        raise SystemExit(1)

    required_count_cols = ["input", "filtered", "denoised", "merged", "non-chimeric"]
    missing = [c for c in required_count_cols if c not in df.columns]
    if missing:
        for col in missing:
            print(f"ERROR: {path}: missing required column '{col}'", file=sys.stderr)
        raise SystemExit(1)

    if id_col != "sample-id":
        df = df.rename(columns={id_col: "sample-id"})

    return df


# ---------------------------------------------------------------------------
# MultiQC aggregation (R1 + R2 → one biological sample)
# ---------------------------------------------------------------------------

_GS_FLOAT_COLS = [
    "percent_duplicates",
    "percent_gc",
    "avg_sequence_length",
    "median_sequence_length",
    "percent_fails",
    "total_sequences",
]


def aggregate_multiqc(gs_df: pd.DataFrame, fq_df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise sample IDs, merge FastQC statuses into general-stats rows,
    then aggregate R1+R2 rows into one biological sample per base ID.

    Aggregation rules:
    - total_sequences: sum
    - other numeric metrics: read-count weighted average (arithmetic mean fallback)
    - FastQC statuses: worst observed (fail > warn > pass > missing)
    """
    gs = gs_df.copy()
    gs["base_id"] = gs["Sample"].apply(normalize_sample_id)

    for col in _GS_FLOAT_COLS:
        if col in gs.columns:
            gs[col] = pd.to_numeric(gs[col], errors="coerce")

    # Parse FastQC status columns
    fq = fq_df.copy()
    for col in _FASTQC_STATUS_COLS:
        if col in fq.columns:
            fq[col] = fq[col].apply(_parse_status)

    fq_status_cols = ["Sample"] + [c for c in _FASTQC_STATUS_COLS if c in fq.columns]
    gs = gs.merge(fq[fq_status_cols], on="Sample", how="left")

    rows = []
    for base_id, grp in gs.groupby("base_id", sort=False):
        row: dict = {"sample_id": base_id, "file_count": len(grp)}

        total_seq = grp["total_sequences"].sum() if "total_sequences" in grp.columns else None
        row["total_sequences"] = total_seq if pd.notna(total_seq) else None

        weights = grp["total_sequences"].fillna(0) if "total_sequences" in grp.columns else pd.Series([0.0] * len(grp))

        for col in [c for c in _GS_FLOAT_COLS if c != "total_sequences"]:
            if col not in grp.columns:
                row[col] = None
                continue
            vals = grp[col]
            valid = vals.notna()
            if not valid.any():
                row[col] = None
            else:
                w = weights[valid]
                v = vals[valid]
                row[col] = (v * w).sum() / w.sum() if w.sum() > 0 else v.mean()

        for col in _FASTQC_STATUS_COLS:
            if col not in grp.columns:
                row[col] = None
                continue
            valid_vals = grp[col].dropna()
            row[col] = int(valid_vals.max()) if len(valid_vals) > 0 else None

        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# DADA2 joining and coverage reporting
# ---------------------------------------------------------------------------


def join_dada2(
    samples_df: pd.DataFrame, dada2_df: pd.DataFrame
) -> tuple[pd.DataFrame, dict]:
    """
    Normalise DADA2 sample IDs, join to aggregated MultiQC samples.
    Returns (joined_df, coverage_dict).
    """
    d2 = dada2_df.copy()
    d2["base_id"] = d2["sample-id"].apply(normalize_sample_id)

    count_cols = ["input", "filtered", "denoised", "merged", "non-chimeric"]
    for col in count_cols:
        if col in d2.columns:
            d2[col] = pd.to_numeric(d2[col], errors="coerce")

    # Warn on duplicate DADA2 base IDs, keep first
    dupes_mask = d2.duplicated("base_id", keep=False)
    if dupes_mask.any():
        n_dup = d2.loc[dupes_mask, "sample-id"].nunique()
        print(f"  WARNING: {n_dup} duplicate DADA2 IDs skipped", file=sys.stderr)
    d2 = d2.drop_duplicates("base_id", keep="first")

    joined = samples_df.merge(
        d2[["base_id"] + count_cols],
        left_on="sample_id",
        right_on="base_id",
        how="left",
    ).drop(columns=["base_id"], errors="ignore")

    multiqc_ids = set(samples_df["sample_id"])
    dada2_ids = set(d2["base_id"])

    return joined, {
        "multiqc_only": sorted(multiqc_ids - dada2_ids),
        "dada2_only": sorted(dada2_ids - multiqc_ids),
    }


# ---------------------------------------------------------------------------
# Monotonic validation
# ---------------------------------------------------------------------------


def validate_dada2_monotonic(df: pd.DataFrame) -> None:
    """Emit a per-sample warning if DADA2 counts violate input >= filtered >= ... >= non-chimeric."""
    stages = ["input", "filtered", "denoised", "merged", "non-chimeric"]
    for _, row in df.iterrows():
        vals = {c: row.get(c) for c in stages}
        if any(pd.isna(v) for v in vals.values() if v is not None):
            continue
        if any(vals[c] is None for c in stages):
            continue
        violated = any(vals[a] < vals[b] for a, b in zip(stages, stages[1:]))
        if violated:
            print(
                f"  WARNING: non-monotonic DADA2 counts for '{row['sample_id']}': "
                f"input={int(vals['input'])} filtered={int(vals['filtered'])} "
                f"denoised={int(vals['denoised'])} merged={int(vals['merged'])} "
                f"non_chimeric={int(vals['non-chimeric'])}",
                file=sys.stderr,
            )


# ---------------------------------------------------------------------------
# Derived metrics
# ---------------------------------------------------------------------------


def compute_derived_metrics(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    def _safe_div(a_col: str, b_col: str) -> pd.Series:
        a = pd.to_numeric(df.get(a_col), errors="coerce")
        b = pd.to_numeric(df.get(b_col), errors="coerce")
        return a.where(b.notna() & (b != 0)) / b.replace(0, float("nan"))

    df["retention_rate"] = _safe_div("non-chimeric", "input")
    df["filter_retention"] = _safe_div("filtered", "input")
    df["denoise_retention"] = _safe_div("denoised", "filtered")
    df["merge_efficiency"] = _safe_div("merged", "denoised")

    den = pd.to_numeric(df.get("denoised"), errors="coerce")
    nc = pd.to_numeric(df.get("non-chimeric"), errors="coerce")
    df["chimera_rate"] = (den - nc) / den.replace(0, float("nan"))

    return df


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def classify_samples(
    df: pd.DataFrame,
    pass_strong: int,
    pass_acceptable: int,
    low_depth: int,
) -> pd.DataFrame:
    df = df.copy()

    def _classify(row):
        nc = row.get("non-chimeric")
        if nc is None or (isinstance(nc, float) and pd.isna(nc)):
            return "MULTIQC_ONLY"
        nc = int(nc)
        if nc >= pass_strong:
            return "PASS_STRONG"
        if nc >= pass_acceptable:
            return "PASS_ACCEPTABLE"
        if nc >= low_depth:
            return "LOW_DEPTH"
        return "FAIL_DEPTH"

    df["classification"] = df.apply(_classify, axis=1)
    return df


# ---------------------------------------------------------------------------
# Flag assignment
# ---------------------------------------------------------------------------


def assign_flags(
    df: pd.DataFrame,
    low_retention: float = LOW_RETENTION_DEFAULT,
    poor_merging: float = POOR_MERGING_DEFAULT,
    high_chimera: float = HIGH_CHIMERA_DEFAULT,
    very_low_filtered: float = VERY_LOW_FILTERED_DEFAULT,
    high_duplication: float = HIGH_DUPLICATION_DEFAULT,
    gc_outlier_sd: float = GC_OUTLIER_SD_DEFAULT,
) -> pd.DataFrame:
    df = df.copy()

    def _notna_lt(col: str, thresh: float) -> pd.Series:
        s = pd.to_numeric(df.get(col), errors="coerce")
        return s.notna() & (s < thresh)

    def _notna_gt(col: str, thresh: float) -> pd.Series:
        s = pd.to_numeric(df.get(col), errors="coerce")
        return s.notna() & (s > thresh)

    df["flag_low_retention"] = _notna_lt("retention_rate", low_retention)
    df["flag_poor_merging"] = _notna_lt("merge_efficiency", poor_merging)
    df["flag_high_chimera"] = _notna_gt("chimera_rate", high_chimera)
    df["flag_very_low_filtered"] = _notna_lt("filter_retention", very_low_filtered)

    # LOW_QUALITY: FastQC per_base_sequence_quality or per_sequence_quality_scores == fail (2)
    pbq = pd.to_numeric(df.get("per_base_sequence_quality", pd.Series([None] * len(df))), errors="coerce")
    psq = pd.to_numeric(df.get("per_sequence_quality_scores", pd.Series([None] * len(df))), errors="coerce")
    df["flag_low_quality"] = (pbq == 2) | (psq == 2)

    # HIGH_ADAPTER: FastQC adapter_content == fail (2)
    ac = pd.to_numeric(df.get("adapter_content", pd.Series([None] * len(df))), errors="coerce")
    df["flag_high_adapter"] = ac == 2

    # HIGH_DUPLICATION: percent_duplicates > threshold
    df["flag_high_duplication"] = _notna_gt("percent_duplicates", high_duplication)

    # GC_OUTLIER: outside mean ± gc_outlier_sd * stdev
    gc = pd.to_numeric(df.get("percent_gc", pd.Series([None] * len(df))), errors="coerce")
    gc_valid = gc.dropna()
    if len(gc_valid) >= 2:
        gc_mean = gc_valid.mean()
        gc_std = gc_valid.std(ddof=0)
        df["flag_gc_outlier"] = gc.notna() & ((gc - gc_mean).abs() > gc_outlier_sd * gc_std)
    else:
        df["flag_gc_outlier"] = False

    return df


# ---------------------------------------------------------------------------
# Recommendations
# ---------------------------------------------------------------------------


def assign_recommendations(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    def _recommend(row):
        if row["classification"] == "FAIL_DEPTH" or row.get("flag_low_retention", False):
            return "exclude"
        if row.get("flag_low_quality", False) or row.get("flag_high_adapter", False):
            return "trim/reprocess"
        if (
            row["classification"] == "LOW_DEPTH"
            or row.get("flag_poor_merging", False)
            or row.get("flag_high_chimera", False)
        ):
            return "review manually"
        return "continue"

    df["recommendation"] = df.apply(_recommend, axis=1)
    return df


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------


def _pct(count: int, total: int) -> float:
    return 0.0 if total == 0 else 100.0 * count / total


def _fmt_f(val, decimals: int = 2) -> str:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return "N/A"
    return f"{float(val):.{decimals}f}"


def _fmt_i(val) -> str:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return "N/A"
    return str(int(val))


def _counts(df: pd.DataFrame) -> dict:
    cls = df["classification"].value_counts().to_dict()
    rec = df["recommendation"].value_counts().to_dict() if "recommendation" in df.columns else {}
    return {
        "n_strong": cls.get("PASS_STRONG", 0),
        "n_acceptable": cls.get("PASS_ACCEPTABLE", 0),
        "n_low": cls.get("LOW_DEPTH", 0),
        "n_fail": cls.get("FAIL_DEPTH", 0),
        "n_mqc_only": cls.get("MULTIQC_ONLY", 0),
        "n_low_ret": int(df.get("flag_low_retention", pd.Series([False] * len(df))).sum()),
        "n_poor_merg": int(df.get("flag_poor_merging", pd.Series([False] * len(df))).sum()),
        "n_high_chim": int(df.get("flag_high_chimera", pd.Series([False] * len(df))).sum()),
        "n_vlf": int(df.get("flag_very_low_filtered", pd.Series([False] * len(df))).sum()),
        "n_low_qual": int(df.get("flag_low_quality", pd.Series([False] * len(df))).sum()),
        "n_high_adapt": int(df.get("flag_high_adapter", pd.Series([False] * len(df))).sum()),
        "n_high_dup": int(df.get("flag_high_duplication", pd.Series([False] * len(df))).sum()),
        "n_gc_out": int(df.get("flag_gc_outlier", pd.Series([False] * len(df))).sum()),
        "n_cont": rec.get("continue", 0),
        "n_trim": rec.get("trim/reprocess", 0),
        "n_rev": rec.get("review manually", 0),
        "n_excl": rec.get("exclude", 0),
    }


def build_stdout_report(
    df: pd.DataFrame,
    has_dada2: bool,
    coverage: dict,
    cfg: dict,
) -> str:
    lines = []
    n = len(df)
    n_files = int(df["file_count"].sum())
    c = _counts(df)

    lines.append("\n==========================================================")
    lines.append("         MICROBIOME QC SUMMARY REPORT")
    lines.append("==========================================================\n")
    lines.append(f"Total biological samples: {n}")
    lines.append(f"MultiQC file rows aggregated: {n_files}")

    if not has_dada2:
        lines.append("\nWARNING: MultiQC-only mode. DADA2 stats not provided.")
        lines.append("  Microbiome usability assessment incomplete.")

    if has_dada2:
        n_matched = n - c["n_mqc_only"]
        n_mqc_only = len(coverage["multiqc_only"])
        n_dada2_only = len(coverage["dada2_only"])
        lines.append("\n-- DADA2 SAMPLE COVERAGE --")
        lines.append(f"  Matched (both MultiQC + DADA2): {n_matched}")
        if n_mqc_only > 0:
            lines.append(f"  WARNING: MultiQC only (no DADA2 match): {n_mqc_only}")
        else:
            lines.append("  MultiQC only (no DADA2 match): 0")
        if n_dada2_only > 0:
            lines.append(f"  WARNING: DADA2 only (no MultiQC match): {n_dada2_only}")
        else:
            lines.append("  DADA2 only (no MultiQC match): 0")

    lines.append("\n-- DEPTH CLASSIFICATION --")
    if has_dada2:
        lines.append(f"  PASS_STRONG     (>={cfg['pass_strong']}): {c['n_strong']}  ({_pct(c['n_strong'], n):.1f}%)")
        lines.append(f"  PASS_ACCEPTABLE (>={cfg['pass_acceptable']}): {c['n_acceptable']}  ({_pct(c['n_acceptable'], n):.1f}%)")
        lines.append(f"  LOW_DEPTH       (>={cfg['low_depth']}): {c['n_low']}  ({_pct(c['n_low'], n):.1f}%)")
        lines.append(f"  FAIL_DEPTH      (<{cfg['low_depth']}):  {c['n_fail']}  ({_pct(c['n_fail'], n):.1f}%)")
        lines.append(f"  MULTIQC_ONLY    (no DADA2): {c['n_mqc_only']}  ({_pct(c['n_mqc_only'], n):.1f}%)")
    else:
        lines.append(f"  MULTIQC_ONLY: {n} biological samples")

    lines.append("\n-- PROCESS FLAGS (DADA2) --")
    if has_dada2:
        lines.append(f"  LOW_RETENTION     (<{cfg['low_retention']*100:.0f}%): {c['n_low_ret']}")
        lines.append(f"  POOR_MERGING      (<{cfg['poor_merging']*100:.0f}%): {c['n_poor_merg']}")
        lines.append(f"  HIGH_CHIMERA      (>{cfg['high_chimera']*100:.0f}%): {c['n_high_chim']}")
        lines.append(f"  VERY_LOW_FILTERED (<{cfg['very_low_filtered']*100:.0f}%): {c['n_vlf']}")
    else:
        lines.append("  (not available without DADA2 stats)")

    lines.append("\n-- SEQUENCING FLAGS (MultiQC) --")
    lines.append(f"  LOW_QUALITY:       {c['n_low_qual']}")
    lines.append(f"  HIGH_ADAPTER:      {c['n_high_adapt']}")
    lines.append(f"  HIGH_DUPLICATION:  {c['n_high_dup']}")
    lines.append(f"  GC_OUTLIER (warn): {c['n_gc_out']}")

    lines.append("\n-- RECOMMENDATIONS --")
    lines.append(f"  continue:        {c['n_cont']}  ({_pct(c['n_cont'], n):.1f}%)")
    lines.append(f"  trim/reprocess:  {c['n_trim']}  ({_pct(c['n_trim'], n):.1f}%)")
    lines.append(f"  review manually: {c['n_rev']}  ({_pct(c['n_rev'], n):.1f}%)")
    lines.append(f"  exclude:         {c['n_excl']}  ({_pct(c['n_excl'], n):.1f}%)")

    mqc_only_list = coverage.get("multiqc_only", [])
    if has_dada2 and mqc_only_list:
        lines.append(f"\n-- MultiQC SAMPLES WITHOUT DADA2 ({len(mqc_only_list)}) --")
        for s in mqc_only_list[:10]:
            lines.append(f"  {s}")
        if len(mqc_only_list) > 10:
            lines.append(f"  ... and {len(mqc_only_list) - 10} more")

    lines.append("")
    return "\n".join(lines)


def build_markdown_report(
    df: pd.DataFrame,
    has_dada2: bool,
    coverage: dict,
    cfg: dict,
    gs_path: str,
    fq_path: str,
    dada2_path: Optional[str],
) -> str:
    lines = []
    n = len(df)
    n_files = int(df["file_count"].sum())
    c = _counts(df)

    lines.append("# Microbiome QC Summary Report\n")
    lines.append("**Inputs:**  ")
    lines.append(f"- General stats: `{gs_path}`  ")
    lines.append(f"- FastQC stats:  `{fq_path}`  ")
    if has_dada2:
        lines.append(f"- DADA2 stats:   `{dada2_path}`  \n")
    else:
        lines.append("- DADA2 stats:   _not provided_ (**MultiQC-only mode**)  \n")
        lines.append("> **MultiQC-only mode**: Provide `--dada2-stats` after DADA2 for complete evaluation.\n")

    lines.append(f"## Dataset Overview\n\n**Total biological samples:** {n}  ")
    lines.append(f"**MultiQC file rows aggregated:** {n_files}\n")

    if has_dada2:
        n_matched = n - c["n_mqc_only"]
        n_mqc_only = len(coverage["multiqc_only"])
        n_dada2_only = len(coverage["dada2_only"])
        lines.append("### DADA2 Sample Coverage\n")
        lines.append("| Set | Count |")
        lines.append("|---|---:|")
        lines.append(f"| Matched (MultiQC + DADA2) | {n_matched} |")
        mqc_warn = "**WARNING** " if n_mqc_only > 0 else ""
        lines.append(f"| {mqc_warn}MultiQC only (no DADA2 match) | {n_mqc_only} |")
        d2_warn = "**WARNING** " if n_dada2_only > 0 else ""
        lines.append(f"| {d2_warn}DADA2 only (no MultiQC match) | {n_dada2_only} |\n")

        lines.append("### Depth Classification (non-chimeric reads)\n")
        lines.append("| Classification | Threshold | Count | Pct |")
        lines.append("|---|---|---:|---:|")
        lines.append(f"| PASS_STRONG | >={cfg['pass_strong']} | {c['n_strong']} | {_pct(c['n_strong'], n):.1f}% |")
        lines.append(f"| PASS_ACCEPTABLE | >={cfg['pass_acceptable']} | {c['n_acceptable']} | {_pct(c['n_acceptable'], n):.1f}% |")
        lines.append(f"| LOW_DEPTH | >={cfg['low_depth']} | {c['n_low']} | {_pct(c['n_low'], n):.1f}% |")
        lines.append(f"| FAIL_DEPTH | <{cfg['low_depth']} | {c['n_fail']} | {_pct(c['n_fail'], n):.1f}% |")
        lines.append(f"| MULTIQC_ONLY | no DADA2 match | {c['n_mqc_only']} | {_pct(c['n_mqc_only'], n):.1f}% |\n")

        lines.append("## Process Flags (DADA2)\n")
        lines.append("| Flag | Threshold | N Flagged |")
        lines.append("|---|---|---:|")
        lines.append(f"| LOW_RETENTION | <{cfg['low_retention']*100:.0f}% | {c['n_low_ret']} |")
        lines.append(f"| POOR_MERGING | <{cfg['poor_merging']*100:.0f}% | {c['n_poor_merg']} |")
        lines.append(f"| HIGH_CHIMERA | >{cfg['high_chimera']*100:.0f}% | {c['n_high_chim']} |")
        lines.append(f"| VERY_LOW_FILTERED | <{cfg['very_low_filtered']*100:.0f}% | {c['n_vlf']} |\n")
    else:
        lines.append("_Depth classification unavailable (MultiQC-only mode)._\n")

    lines.append("## Sequencing Flags (MultiQC)\n")
    lines.append("| Flag | N Flagged |")
    lines.append("|---|---:|")
    lines.append(f"| LOW_QUALITY | {c['n_low_qual']} |")
    lines.append(f"| HIGH_ADAPTER | {c['n_high_adapt']} |")
    lines.append(f"| HIGH_DUPLICATION (>{cfg['high_duplication']:.0f}%) | {c['n_high_dup']} |")
    lines.append(f"| GC_OUTLIER (>{cfg['gc_outlier_sd']:.1f}SD, warn) | {c['n_gc_out']} |\n")

    if has_dada2:
        lines.append("## Low/Failing Depth Samples\n")
        low_fail = df[df["classification"].isin(["FAIL_DEPTH", "LOW_DEPTH"])]
        if len(low_fail) > 0:
            lines.append("| Sample | Non-Chimeric | Class | Recommendation |")
            lines.append("|---|---:|---|---|")
            for _, row in low_fail.iterrows():
                lines.append(
                    f"| {row['sample_id']} | {_fmt_i(row.get('non-chimeric'))} "
                    f"| {row['classification']} | {row['recommendation']} |"
                )
        else:
            lines.append("_None below depth thresholds._")
        lines.append("")

        for flag_col, flag_label in [
            ("flag_low_retention", "LOW_RETENTION"),
            ("flag_poor_merging", "POOR_MERGING"),
            ("flag_high_chimera", "HIGH_CHIMERA"),
            ("flag_very_low_filtered", "VERY_LOW_FILTERED"),
        ]:
            lines.append(f"## Process Flag: {flag_label}\n")
            flagged = df[df.get(flag_col, pd.Series([False] * len(df))).astype(bool)]
            if len(flagged) > 0:
                lines.append("| Sample | NonChimeric | Retention | MergeEff | ChimeraRate | FilterRet |")
                lines.append("|---|---:|---:|---:|---:|---:|")
                for _, row in flagged.iterrows():
                    lines.append(
                        f"| {row['sample_id']} | {_fmt_i(row.get('non-chimeric'))} "
                        f"| {_fmt_f(row.get('retention_rate'), 4)} "
                        f"| {_fmt_f(row.get('merge_efficiency'), 4)} "
                        f"| {_fmt_f(row.get('chimera_rate'), 4)} "
                        f"| {_fmt_f(row.get('filter_retention'), 4)} |"
                    )
            else:
                lines.append("_None._")
            lines.append("")

    lines.append("## Sequencing-Flagged Samples\n")
    for flag_col, flag_label in [
        ("flag_low_quality", "LOW_QUALITY"),
        ("flag_high_adapter", "HIGH_ADAPTER"),
        ("flag_high_duplication", "HIGH_DUPLICATION"),
        ("flag_gc_outlier", "GC_OUTLIER"),
    ]:
        lines.append(f"### {flag_label}\n")
        flagged = df[df.get(flag_col, pd.Series([False] * len(df))).astype(bool)]
        if len(flagged) > 0:
            lines.append("| Sample | %GC | %Dup | TotalSeqs | Adapter | BaseQual | SeqQual |")
            lines.append("|---|---:|---:|---:|---|---|---|")
            for _, row in flagged.iterrows():
                lines.append(
                    f"| {row['sample_id']} "
                    f"| {_fmt_f(row.get('percent_gc'), 1)} "
                    f"| {_fmt_f(row.get('percent_duplicates'), 1)} "
                    f"| {_fmt_f(row.get('total_sequences'), 0)} "
                    f"| {_status_str(row.get('adapter_content'))} "
                    f"| {_status_str(row.get('per_base_sequence_quality'))} "
                    f"| {_status_str(row.get('per_sequence_quality_scores'))} |"
                )
        else:
            lines.append("_None._")
        lines.append("")

    # Top/bottom 10
    n_top = min(10, n)
    if has_dada2:
        sorted_df = df.copy()
        sorted_df["_nc_sort"] = pd.to_numeric(sorted_df.get("non-chimeric"), errors="coerce").fillna(0)
        top = sorted_df.sort_values("_nc_sort", ascending=False).head(n_top)
        lines.append(f"## Top {n_top} Samples (Non-Chimeric Reads)\n")
        lines.append("| Rank | Sample | NonChimeric | Retention | Class |")
        lines.append("|---:|---|---:|---:|---|")
        for rank, (_, row) in enumerate(top.iterrows(), 1):
            lines.append(
                f"| {rank} | {row['sample_id']} | {_fmt_i(row.get('non-chimeric'))} "
                f"| {_fmt_f(row.get('retention_rate'), 4)} | {row['classification']} |"
            )
        lines.append("")

        bottom = sorted_df.sort_values("_nc_sort", ascending=True).head(n_top)
        lines.append(f"## Bottom {n_top} Samples (Non-Chimeric Reads)\n")
        lines.append("| Rank | Sample | NonChimeric | Retention | Class |")
        lines.append("|---:|---|---:|---:|---|")
        for rank, (_, row) in enumerate(bottom.iterrows(), 1):
            lines.append(
                f"| {rank} | {row['sample_id']} | {_fmt_i(row.get('non-chimeric'))} "
                f"| {_fmt_f(row.get('retention_rate'), 4)} | {row['classification']} |"
            )
        lines.append("")
    else:
        sorted_df = df.copy()
        sorted_df["_ts_sort"] = pd.to_numeric(sorted_df.get("total_sequences"), errors="coerce").fillna(0)
        top = sorted_df.sort_values("_ts_sort", ascending=False).head(n_top)
        lines.append(f"## Top {n_top} by Total Sequences\n")
        lines.append("| Rank | Sample | TotalSeqs | %GC | %Dup |")
        lines.append("|---:|---|---:|---:|---:|")
        for rank, (_, row) in enumerate(top.iterrows(), 1):
            lines.append(
                f"| {rank} | {row['sample_id']} | {_fmt_f(row.get('total_sequences'), 0)} "
                f"| {_fmt_f(row.get('percent_gc'), 1)} | {_fmt_f(row.get('percent_duplicates'), 1)} |"
            )
        lines.append("")

    lines.append("## Recommended Actions\n")
    lines.append("| Action | Count | Pct |")
    lines.append("|---|---:|---:|")
    lines.append(f"| continue | {c['n_cont']} | {_pct(c['n_cont'], n):.1f}% |")
    lines.append(f"| trim/reprocess | {c['n_trim']} | {_pct(c['n_trim'], n):.1f}% |")
    lines.append(f"| review manually | {c['n_rev']} | {_pct(c['n_rev'], n):.1f}% |")
    lines.append(f"| exclude | {c['n_excl']} | {_pct(c['n_excl'], n):.1f}% |\n")

    lines.append("### Overall Recommendation\n")
    if c["n_excl"] > n // 4:
        lines.append("**ACTION REQUIRED**: >25% samples flagged for exclusion. Review and consider reprocessing.\n")
    elif c["n_trim"] > 0:
        lines.append("**TRIM/REPROCESS**: Sequencing quality issues found. Consider quality trimming.\n")
    elif (c["n_low"] + c["n_fail"] > 0) and has_dada2:
        lines.append("**REVIEW**: Some samples have low depth. Review for exclusion or inclusion with caveats.\n")
    else:
        lines.append("**PROCEED**: Dataset suitable for downstream analysis. Apply standard rarefaction/normalization.\n")

    mqc_only_list = coverage.get("multiqc_only", [])
    if has_dada2 and mqc_only_list:
        lines.append(f"## MultiQC Samples Without DADA2 ({len(mqc_only_list)})\n")
        for s in mqc_only_list:
            lines.append(f"- `{s}`")
        lines.append("")

    dada2_only_list = coverage.get("dada2_only", [])
    if dada2_only_list:
        lines.append(f"## DADA2 Samples Not Matched in MultiQC ({len(dada2_only_list)})\n")
        for s in dada2_only_list:
            lines.append(f"- `{s}`")
        lines.append("")

    lines.append("## Thresholds Used\n")
    lines.append("| Parameter | Value |")
    lines.append("|---|---|")
    lines.append(f"| pass_strong (>=) | {cfg['pass_strong']} |")
    lines.append(f"| pass_acceptable (>=) | {cfg['pass_acceptable']} |")
    lines.append(f"| low_depth (>=) | {cfg['low_depth']} |")
    lines.append(f"| low_retention (<) | {cfg['low_retention']:.2f} |")
    lines.append(f"| poor_merging (<) | {cfg['poor_merging']:.2f} |")
    lines.append(f"| high_chimera (>) | {cfg['high_chimera']:.2f} |")
    lines.append(f"| very_low_filtered (<) | {cfg['very_low_filtered']:.2f} |")
    lines.append(f"| high_duplication (>) | {cfg['high_duplication']:.1f}% |")
    lines.append(f"| gc_outlier_sd | {cfg['gc_outlier_sd']:.1f} |")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# TSV output writers
# ---------------------------------------------------------------------------


def write_samples_tsv(df: pd.DataFrame, has_dada2: bool, prefix: str) -> None:
    cols = ["sample_id", "file_count", "total_sequences", "percent_gc",
            "percent_duplicates", "avg_sequence_length", "median_sequence_length",
            "percent_fails"]
    cols += [c for c in _FASTQC_STATUS_COLS if c in df.columns]
    if has_dada2:
        cols += ["input", "filtered", "denoised", "merged", "non-chimeric",
                 "retention_rate", "filter_retention", "denoise_retention",
                 "merge_efficiency", "chimera_rate"]
    cols += ["classification"]
    if has_dada2:
        cols += ["flag_low_retention", "flag_poor_merging",
                 "flag_high_chimera", "flag_very_low_filtered"]
    cols += ["flag_low_quality", "flag_high_adapter",
             "flag_high_duplication", "flag_gc_outlier", "recommendation"]

    out = df[[c for c in cols if c in df.columns]].copy()

    # Convert FastQC status ints back to pass/warn/fail strings
    for col in _FASTQC_STATUS_COLS:
        if col in out.columns:
            out[col] = out[col].apply(_status_str)

    out.to_csv(f"{prefix}.samples.tsv", sep="\t", index=False)


def write_flags_tsv(df: pd.DataFrame, has_dada2: bool, prefix: str) -> None:
    rows = []

    if has_dada2:
        process_flags = [
            ("flag_low_retention", "LOW_RETENTION", "retention_rate"),
            ("flag_poor_merging", "POOR_MERGING", "merge_efficiency"),
            ("flag_high_chimera", "HIGH_CHIMERA", "chimera_rate"),
            ("flag_very_low_filtered", "VERY_LOW_FILTERED", "filter_retention"),
        ]
    else:
        process_flags = []

    seq_flags = [
        ("flag_low_quality", "LOW_QUALITY", None),
        ("flag_high_adapter", "HIGH_ADAPTER", None),
        ("flag_high_duplication", "HIGH_DUPLICATION", "percent_duplicates"),
        ("flag_gc_outlier", "GC_OUTLIER", "percent_gc"),
    ]

    for _, row in df.iterrows():
        for flag_col, flag_name, val_col in process_flags:
            if row.get(flag_col, False):
                val = _fmt_f(row.get(val_col), 6) if val_col else ""
                rows.append({"sample_id": row["sample_id"], "flag_type": "process",
                             "flag_name": flag_name, "value": val})
        for flag_col, flag_name, val_col in seq_flags:
            if row.get(flag_col, False):
                val = _fmt_f(row.get(val_col), 2) if val_col else ""
                rows.append({"sample_id": row["sample_id"], "flag_type": "sequencing",
                             "flag_name": flag_name, "value": val})

    flags_df = pd.DataFrame(rows, columns=["sample_id", "flag_type", "flag_name", "value"])
    flags_df.to_csv(f"{prefix}.flags.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--multiqc-general-stats", required=True, metavar="<file>",
              help="multiqc_general_stats.txt (tab-separated)")
@click.option("--multiqc-fastqc", required=True, metavar="<file>",
              help="multiqc_fastqc.txt (tab-separated)")
@click.option("--dada2-stats", default=None, metavar="<file>",
              help="DADA2 denoising stats TSV (optional)")
@click.option("--out-prefix", default="microbiome_qc_summary", show_default=True,
              metavar="<prefix>", help="Output file prefix")
@click.option("--pass-strong", default=PASS_STRONG_DEFAULT, type=int, show_default=True)
@click.option("--pass-acceptable", default=PASS_ACCEPTABLE_DEFAULT, type=int, show_default=True)
@click.option("--low-depth", default=LOW_DEPTH_DEFAULT, type=int, show_default=True)
@click.option("--low-retention", default=LOW_RETENTION_DEFAULT, type=float, show_default=True)
@click.option("--poor-merging", default=POOR_MERGING_DEFAULT, type=float, show_default=True)
@click.option("--high-chimera", default=HIGH_CHIMERA_DEFAULT, type=float, show_default=True)
@click.option("--very-low-filtered", default=VERY_LOW_FILTERED_DEFAULT, type=float, show_default=True)
@click.option("--high-duplication", default=HIGH_DUPLICATION_DEFAULT, type=float, show_default=True)
@click.option("--gc-outlier-sd", default=GC_OUTLIER_SD_DEFAULT, type=float, show_default=True)
def main(
    multiqc_general_stats,
    multiqc_fastqc,
    dada2_stats,
    out_prefix,
    pass_strong,
    pass_acceptable,
    low_depth,
    low_retention,
    poor_merging,
    high_chimera,
    very_low_filtered,
    high_duplication,
    gc_outlier_sd,
):
    """Microbiome 16S QC summary tool."""
    cfg = {
        "pass_strong": pass_strong,
        "pass_acceptable": pass_acceptable,
        "low_depth": low_depth,
        "low_retention": low_retention,
        "poor_merging": poor_merging,
        "high_chimera": high_chimera,
        "very_low_filtered": very_low_filtered,
        "high_duplication": high_duplication,
        "gc_outlier_sd": gc_outlier_sd,
    }

    # Parse inputs
    print(f"Parsing: {multiqc_general_stats}", file=sys.stderr)
    gs_df = parse_general_stats(multiqc_general_stats)
    print(f"  -> {len(gs_df)} MultiQC general stats rows loaded", file=sys.stderr)

    print(f"Parsing: {multiqc_fastqc}", file=sys.stderr)
    fq_df = parse_fastqc_stats(multiqc_fastqc)

    # Aggregate R1+R2 rows per biological sample
    samples_df = aggregate_multiqc(gs_df, fq_df)
    print(f"  -> {len(samples_df)} biological samples after aggregation", file=sys.stderr)

    has_dada2 = dada2_stats is not None
    coverage: dict = {"multiqc_only": [], "dada2_only": []}

    if has_dada2:
        print(f"Parsing: {dada2_stats}", file=sys.stderr)
        d2_df = parse_dada2_stats(dada2_stats)
        samples_df, coverage = join_dada2(samples_df, d2_df)

        n_mqc_only = len(coverage["multiqc_only"])
        n_d2_only = len(coverage["dada2_only"])
        if n_mqc_only > 0:
            print(f"  WARNING: {n_mqc_only} MultiQC biological samples had no DADA2 match", file=sys.stderr)
        if n_d2_only > 0:
            print(f"  WARNING: {n_d2_only} DADA2 IDs not matched in MultiQC", file=sys.stderr)

        samples_df = compute_derived_metrics(samples_df)
        validate_dada2_monotonic(samples_df)
    else:
        print("No DADA2 stats -- MultiQC-only mode.", file=sys.stderr)

    samples_df = classify_samples(samples_df, pass_strong, pass_acceptable, low_depth)
    samples_df = assign_flags(
        samples_df,
        low_retention=low_retention,
        poor_merging=poor_merging,
        high_chimera=high_chimera,
        very_low_filtered=very_low_filtered,
        high_duplication=high_duplication,
        gc_outlier_sd=gc_outlier_sd,
    )
    samples_df = assign_recommendations(samples_df)

    # Write outputs
    print(f"Writing outputs: {out_prefix}.*", file=sys.stderr)

    stdout_text = build_stdout_report(samples_df, has_dada2, coverage, cfg)
    print(stdout_text)

    md_text = build_markdown_report(
        samples_df, has_dada2, coverage, cfg,
        multiqc_general_stats, multiqc_fastqc, dada2_stats,
    )
    Path(f"{out_prefix}.report.md").write_text(md_text)

    write_samples_tsv(samples_df, has_dada2, out_prefix)
    write_flags_tsv(samples_df, has_dada2, out_prefix)

    print(
        f"Done:\n  {out_prefix}.report.md\n  {out_prefix}.samples.tsv\n  {out_prefix}.flags.tsv",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
