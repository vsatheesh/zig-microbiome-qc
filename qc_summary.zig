//! microbiome_qc_summary - Decision-oriented QC report for 16S microbiome projects.
//! Zig 0.17.0-dev compatible (std.Io API, unmanaged ArrayList).
//!
//! Usage:
//!   qc_summary --multiqc-general-stats <path> --multiqc-fastqc <path> [options]

const std = @import("std");
const Io = std.Io;
const Allocator = std.mem.Allocator;

// ArrayList is now unmanaged in 0.17.0-dev: use .empty, deinit(alloc), append(alloc, ...)
fn AL(comptime T: type) type {
    return std.ArrayList(T);
}

// ============================================================
// Configuration
// ============================================================

const Config = struct {
    multiqc_general_stats: []const u8 = "",
    multiqc_fastqc: []const u8 = "",
    dada2_stats: []const u8 = "",
    out_prefix: []const u8 = "microbiome_qc_summary",
    pass_strong: u64 = 20_000,
    pass_acceptable: u64 = 10_000,
    low_depth: u64 = 5_000,
    low_retention: f64 = 0.10,
    poor_merging: f64 = 0.60,
    high_chimera: f64 = 0.30,
    very_low_filtered: f64 = 0.20,
    high_dup_pct: f64 = 90.0,
    gc_outlier_sd: f64 = 3.0,
};

// ============================================================
// Data structures
// ============================================================

const Classification = enum {
    pass_strong,
    pass_acceptable,
    low_depth,
    fail_depth,
    multiqc_only,

    fn label(self: Classification) []const u8 {
        return switch (self) {
            .pass_strong => "PASS_STRONG",
            .pass_acceptable => "PASS_ACCEPTABLE",
            .low_depth => "LOW_DEPTH",
            .fail_depth => "FAIL_DEPTH",
            .multiqc_only => "MULTIQC_ONLY",
        };
    }
};

const Recommendation = enum {
    continue_analysis,
    trim_reprocess,
    review_manually,
    exclude,

    fn label(self: Recommendation) []const u8 {
        return switch (self) {
            .continue_analysis => "continue",
            .trim_reprocess => "trim/reprocess",
            .review_manually => "review manually",
            .exclude => "exclude",
        };
    }
};

const NA_F64: f64 = -1.0;
const NA_U64: u64 = std.math.maxInt(u64);

const Sample = struct {
    id: []const u8,
    base_id: []const u8,
    is_r1: bool,
    // MultiQC general stats
    pct_duplicates: f64 = NA_F64,
    pct_gc: f64 = NA_F64,
    avg_seq_len: f64 = NA_F64,
    median_seq_len: f64 = NA_F64,
    pct_fails: f64 = NA_F64,
    total_sequences: f64 = NA_F64,
    // FastQC module statuses: 0=pass 1=warn 2=fail 255=missing
    basic_statistics: u8 = 255,
    per_base_seq_quality: u8 = 255,
    per_tile_seq_quality: u8 = 255,
    per_seq_quality_scores: u8 = 255,
    per_base_seq_content: u8 = 255,
    per_seq_gc_content: u8 = 255,
    per_base_n_content: u8 = 255,
    seq_len_distribution: u8 = 255,
    seq_duplication_levels: u8 = 255,
    overrepresented_seqs: u8 = 255,
    adapter_content: u8 = 255,
    // DADA2 stats (optional)
    dada2_input: u64 = NA_U64,
    dada2_filtered: u64 = NA_U64,
    dada2_denoised: u64 = NA_U64,
    dada2_merged: u64 = NA_U64,
    dada2_non_chimeric: u64 = NA_U64,
    // Derived
    retention_rate: f64 = NA_F64,
    filter_retention: f64 = NA_F64,
    denoise_retention: f64 = NA_F64,
    merge_efficiency: f64 = NA_F64,
    chimera_rate: f64 = NA_F64,
    // Output
    classification: Classification = .multiqc_only,
    flag_low_retention: bool = false,
    flag_poor_merging: bool = false,
    flag_high_chimera: bool = false,
    flag_very_low_filtered: bool = false,
    flag_low_quality: bool = false,
    flag_high_adapter: bool = false,
    flag_high_duplication: bool = false,
    flag_gc_outlier: bool = false,
    recommendation: Recommendation = .review_manually,
};

// ============================================================
// Parsing helpers
// ============================================================

fn parseF64(s: []const u8) f64 {
    const t = std.mem.trim(u8, s, " \t\r\n");
    if (t.len == 0 or std.mem.eql(u8, t, "N/A") or std.mem.eql(u8, t, "NA")) return NA_F64;
    return std.fmt.parseFloat(f64, t) catch NA_F64;
}

fn parseU64fromField(s: []const u8) u64 {
    const t = std.mem.trim(u8, s, " \t\r\n");
    if (t.len == 0 or std.mem.eql(u8, t, "N/A") or std.mem.eql(u8, t, "NA")) return NA_U64;
    if (std.fmt.parseFloat(f64, t)) |f| {
        if (f >= 0.0) return @as(u64, @intFromFloat(f));
        return NA_U64;
    } else |_| {}
    return std.fmt.parseInt(u64, t, 10) catch NA_U64;
}

fn parseStatus(s: []const u8) u8 {
    const t = std.mem.trim(u8, s, " \t\r\n");
    if (std.mem.eql(u8, t, "pass")) return 0;
    if (std.mem.eql(u8, t, "warn")) return 1;
    if (std.mem.eql(u8, t, "fail")) return 2;
    return 255;
}

fn stripReadSuffix(name: []const u8) struct { base: []const u8, is_r1: bool } {
    if (std.mem.endsWith(u8, name, "_R1")) return .{ .base = name[0 .. name.len - 3], .is_r1 = true };
    if (std.mem.endsWith(u8, name, "_R2")) return .{ .base = name[0 .. name.len - 3], .is_r1 = false };
    if (std.mem.endsWith(u8, name, "_1")) return .{ .base = name[0 .. name.len - 2], .is_r1 = true };
    if (std.mem.endsWith(u8, name, "_2")) return .{ .base = name[0 .. name.len - 2], .is_r1 = false };
    return .{ .base = name, .is_r1 = true };
}

fn divSafe(a: f64, b: f64) f64 {
    if (b == 0.0 or b == NA_F64 or a == NA_F64) return NA_F64;
    return a / b;
}

// ============================================================
// TSV parsing
// ============================================================

const Row = AL([]u8);

fn parseTsv(
    io: Io,
    alloc: Allocator,
    path: []const u8,
    headers_out: *AL([]u8),
    rows_out: *AL(Row),
) !void {
    const file = Io.Dir.openFileAbsolute(io, path, .{}) catch |err| {
        std.debug.print("ERROR: Cannot open: {s} ({any})\n", .{ path, err });
        return err;
    };
    defer file.close(io);

    var read_buf: [131072]u8 = undefined;
    var rdr = file.reader(io, &read_buf);
    var first_line = true;

    while (true) {
        const maybe_line = rdr.interface.takeDelimiter('\n') catch |err| switch (err) {
            error.ReadFailed => return error.ReadFailed,
            error.StreamTooLong => continue,
        };
        const raw = maybe_line orelse break;
        const line = std.mem.trim(u8, raw, "\r\n \t");
        if (line.len == 0) continue;

        var fields: Row = .empty;
        var it = std.mem.splitScalar(u8, line, '\t');
        while (it.next()) |field| {
            try fields.append(alloc, try alloc.dupe(u8, field));
        }

        if (first_line) {
            first_line = false;
            for (fields.items) |f| try headers_out.append(alloc, f);
            fields.deinit(alloc);
        } else {
            try rows_out.append(alloc, fields);
        }
    }
}

fn buildColIndex(alloc: Allocator, headers: []const []u8) !std.StringHashMap(usize) {
    var map = std.StringHashMap(usize).init(alloc);
    for (headers, 0..) |h, i| {
        const key = std.mem.trim(u8, h, " \t\r\n");
        try map.put(key, i);
    }
    return map;
}

fn getField(row: []const []u8, col_map: *const std.StringHashMap(usize), name: []const u8) []const u8 {
    if (col_map.get(name)) |idx| {
        if (idx < row.len) return std.mem.trim(u8, row[idx], " \t\r\n");
    }
    return "";
}

// ============================================================
// Main parsing functions
// ============================================================

fn parseGeneralStats(
    io: Io,
    alloc: Allocator,
    path: []const u8,
    sample_map: *std.StringHashMap(usize),
    samples: *AL(Sample),
) !void {
    var headers: AL([]u8) = .empty;
    defer {
        for (headers.items) |h| alloc.free(h);
        headers.deinit(alloc);
    }
    var rows: AL(Row) = .empty;
    defer {
        for (rows.items) |*r| {
            for (r.items) |f| alloc.free(f);
            r.deinit(alloc);
        }
        rows.deinit(alloc);
    }
    try parseTsv(io, alloc, path, &headers, &rows);
    var col = try buildColIndex(alloc, headers.items);
    defer col.deinit();

    for (rows.items) |row| {
        if (row.items.len == 0) continue;
        const raw_id = getField(row.items, &col, "Sample");
        if (raw_id.len == 0) continue;
        const si = stripReadSuffix(raw_id);
        var s = Sample{
            .id = try alloc.dupe(u8, raw_id),
            .base_id = try alloc.dupe(u8, si.base),
            .is_r1 = si.is_r1,
        };
        s.pct_duplicates = parseF64(getField(row.items, &col, "percent_duplicates"));
        s.pct_gc = parseF64(getField(row.items, &col, "percent_gc"));
        s.avg_seq_len = parseF64(getField(row.items, &col, "avg_sequence_length"));
        s.median_seq_len = parseF64(getField(row.items, &col, "median_sequence_length"));
        s.pct_fails = parseF64(getField(row.items, &col, "percent_fails"));
        s.total_sequences = parseF64(getField(row.items, &col, "total_sequences"));
        const idx = samples.items.len;
        try samples.append(alloc, s);
        try sample_map.put(samples.items[idx].id, idx);
    }
}

fn parseFastqcStats(
    io: Io,
    alloc: Allocator,
    path: []const u8,
    sample_map: *std.StringHashMap(usize),
    samples: *AL(Sample),
) !void {
    var headers: AL([]u8) = .empty;
    defer {
        for (headers.items) |h| alloc.free(h);
        headers.deinit(alloc);
    }
    var rows: AL(Row) = .empty;
    defer {
        for (rows.items) |*r| {
            for (r.items) |f| alloc.free(f);
            r.deinit(alloc);
        }
        rows.deinit(alloc);
    }
    try parseTsv(io, alloc, path, &headers, &rows);
    var col = try buildColIndex(alloc, headers.items);
    defer col.deinit();

    for (rows.items) |row| {
        if (row.items.len == 0) continue;
        const raw_id = getField(row.items, &col, "Sample");
        if (raw_id.len == 0) continue;
        const idx = sample_map.get(raw_id) orelse continue;
        var s = &samples.items[idx];
        s.basic_statistics = parseStatus(getField(row.items, &col, "basic_statistics"));
        s.per_base_seq_quality = parseStatus(getField(row.items, &col, "per_base_sequence_quality"));
        s.per_tile_seq_quality = parseStatus(getField(row.items, &col, "per_tile_sequence_quality"));
        s.per_seq_quality_scores = parseStatus(getField(row.items, &col, "per_sequence_quality_scores"));
        s.per_base_seq_content = parseStatus(getField(row.items, &col, "per_base_sequence_content"));
        s.per_seq_gc_content = parseStatus(getField(row.items, &col, "per_sequence_gc_content"));
        s.per_base_n_content = parseStatus(getField(row.items, &col, "per_base_n_content"));
        s.seq_len_distribution = parseStatus(getField(row.items, &col, "sequence_length_distribution"));
        s.seq_duplication_levels = parseStatus(getField(row.items, &col, "sequence_duplication_levels"));
        s.overrepresented_seqs = parseStatus(getField(row.items, &col, "overrepresented_sequences"));
        s.adapter_content = parseStatus(getField(row.items, &col, "adapter_content"));
    }
}

fn parseDada2Stats(
    io: Io,
    alloc: Allocator,
    path: []const u8,
    sample_map: *std.StringHashMap(usize),
    samples: *AL(Sample),
    unmatched: *AL([]u8),
) !void {
    var headers: AL([]u8) = .empty;
    defer {
        for (headers.items) |h| alloc.free(h);
        headers.deinit(alloc);
    }
    var rows: AL(Row) = .empty;
    defer {
        for (rows.items) |*r| {
            for (r.items) |f| alloc.free(f);
            r.deinit(alloc);
        }
        rows.deinit(alloc);
    }
    try parseTsv(io, alloc, path, &headers, &rows);
    var col = try buildColIndex(alloc, headers.items);
    defer col.deinit();

    const id_col: []const u8 = if (col.contains("sample-id")) "sample-id" else "sample_id";

    for (rows.items) |row| {
        if (row.items.len == 0) continue;
        const raw_id = getField(row.items, &col, id_col);
        if (raw_id.len == 0) continue;

        var matched = false;
        if (sample_map.get(raw_id)) |idx| {
            applyDada2(row.items, &col, &samples.items[idx]);
            matched = true;
        } else {
            for (samples.items) |*s| {
                if (std.mem.eql(u8, s.base_id, raw_id)) {
                    applyDada2(row.items, &col, s);
                    matched = true;
                }
            }
        }
        if (!matched) try unmatched.append(alloc, try alloc.dupe(u8, raw_id));
    }
}

fn applyDada2(row: []const []u8, col: *const std.StringHashMap(usize), s: *Sample) void {
    s.dada2_input = parseU64fromField(getField(row, col, "input"));
    s.dada2_filtered = parseU64fromField(getField(row, col, "filtered"));
    s.dada2_denoised = parseU64fromField(getField(row, col, "denoised"));
    s.dada2_merged = parseU64fromField(getField(row, col, "merged"));
    s.dada2_non_chimeric = parseU64fromField(getField(row, col, "non-chimeric"));
}

// ============================================================
// Metrics, classification, flags
// ============================================================

fn computeDerivedMetrics(samples: []Sample) void {
    for (samples) |*s| {
        if (s.dada2_input == NA_U64) continue;
        const inp: f64 = @floatFromInt(s.dada2_input);
        const filt: f64 = if (s.dada2_filtered == NA_U64) NA_F64 else @as(f64, @floatFromInt(s.dada2_filtered));
        const den: f64 = if (s.dada2_denoised == NA_U64) NA_F64 else @as(f64, @floatFromInt(s.dada2_denoised));
        const merg: f64 = if (s.dada2_merged == NA_U64) NA_F64 else @as(f64, @floatFromInt(s.dada2_merged));
        const nc: f64 = if (s.dada2_non_chimeric == NA_U64) NA_F64 else @as(f64, @floatFromInt(s.dada2_non_chimeric));
        s.retention_rate = divSafe(nc, inp);
        s.filter_retention = divSafe(filt, inp);
        s.denoise_retention = divSafe(den, filt);
        s.merge_efficiency = divSafe(merg, den);
        if (merg != NA_F64 and merg > 0.0 and nc != NA_F64) {
            s.chimera_rate = 1.0 - (nc / merg);
        }
    }
}

fn classifySamples(samples: []Sample, cfg: Config) void {
    for (samples) |*s| {
        if (s.dada2_non_chimeric == NA_U64) {
            s.classification = .multiqc_only;
        } else {
            const nc = s.dada2_non_chimeric;
            s.classification = if (nc >= cfg.pass_strong) .pass_strong else if (nc >= cfg.pass_acceptable) .pass_acceptable else if (nc >= cfg.low_depth) .low_depth else .fail_depth;
        }
    }
}

fn computeGcStats(samples: []const Sample) struct { mean: f64, sd: f64 } {
    var sum: f64 = 0.0;
    var count: f64 = 0.0;
    for (samples) |s| {
        if (s.pct_gc != NA_F64) {
            sum += s.pct_gc;
            count += 1.0;
        }
    }
    if (count == 0.0) return .{ .mean = NA_F64, .sd = NA_F64 };
    const mean = sum / count;
    var sq_sum: f64 = 0.0;
    for (samples) |s| {
        if (s.pct_gc != NA_F64) {
            const d = s.pct_gc - mean;
            sq_sum += d * d;
        }
    }
    return .{ .mean = mean, .sd = @sqrt(sq_sum / count) };
}

fn assignFlags(samples: []Sample, cfg: Config) void {
    const gc = computeGcStats(samples);
    const gc_thresh: f64 = if (gc.sd != NA_F64) cfg.gc_outlier_sd * gc.sd else NA_F64;
    for (samples) |*s| {
        if (s.retention_rate != NA_F64) s.flag_low_retention = s.retention_rate < cfg.low_retention;
        if (s.merge_efficiency != NA_F64) s.flag_poor_merging = s.merge_efficiency < cfg.poor_merging;
        if (s.chimera_rate != NA_F64) s.flag_high_chimera = s.chimera_rate > cfg.high_chimera;
        if (s.filter_retention != NA_F64) s.flag_very_low_filtered = s.filter_retention < cfg.very_low_filtered;
        s.flag_low_quality = (s.per_base_seq_quality == 2 or s.per_seq_quality_scores == 2);
        s.flag_high_adapter = (s.adapter_content == 2);
        s.flag_high_duplication = (s.pct_duplicates != NA_F64 and s.pct_duplicates > cfg.high_dup_pct);
        if (gc_thresh != NA_F64 and s.pct_gc != NA_F64 and gc.mean != NA_F64) {
            s.flag_gc_outlier = @abs(s.pct_gc - gc.mean) > gc_thresh;
        }
    }
}

fn assignRecommendations(samples: []Sample) void {
    for (samples) |*s| {
        const severe = (s.classification == .fail_depth) or
            s.flag_very_low_filtered or
            (s.flag_high_chimera and s.flag_low_retention);
        const seq_issue = s.flag_high_adapter or s.flag_low_quality;
        const moderate = (s.classification == .low_depth) or s.flag_poor_merging;
        if (severe) {
            s.recommendation = .exclude;
        } else if (seq_issue) {
            s.recommendation = .trim_reprocess;
        } else if (moderate) {
            s.recommendation = .review_manually;
        } else {
            s.recommendation = .continue_analysis;
        }
    }
}

// ============================================================
// Sorting
// ============================================================

fn cmpNcDesc(_: void, a: Sample, b: Sample) bool {
    const av = if (a.dada2_non_chimeric == NA_U64) 0 else a.dada2_non_chimeric;
    const bv = if (b.dada2_non_chimeric == NA_U64) 0 else b.dada2_non_chimeric;
    return av > bv;
}
fn cmpNcAsc(_: void, a: Sample, b: Sample) bool {
    const av = if (a.dada2_non_chimeric == NA_U64) 0 else a.dada2_non_chimeric;
    const bv = if (b.dada2_non_chimeric == NA_U64) 0 else b.dada2_non_chimeric;
    return av < bv;
}
fn cmpTotalDesc(_: void, a: Sample, b: Sample) bool {
    const av: f64 = if (a.total_sequences == NA_F64) 0.0 else a.total_sequences;
    const bv: f64 = if (b.total_sequences == NA_F64) 0.0 else b.total_sequences;
    return av > bv;
}

// ============================================================
// Formatting helpers
// ============================================================

fn f2s(buf: []u8, val: f64, comptime dec: u8) []const u8 {
    if (val == NA_F64) return "N/A";
    const s = switch (dec) {
        0 => std.fmt.bufPrint(buf, "{d:.0}", .{val}),
        1 => std.fmt.bufPrint(buf, "{d:.1}", .{val}),
        2 => std.fmt.bufPrint(buf, "{d:.2}", .{val}),
        3 => std.fmt.bufPrint(buf, "{d:.3}", .{val}),
        4 => std.fmt.bufPrint(buf, "{d:.4}", .{val}),
        6 => std.fmt.bufPrint(buf, "{d:.6}", .{val}),
        else => std.fmt.bufPrint(buf, "{d:.2}", .{val}),
    };
    return s catch "err";
}

fn u2s(buf: []u8, val: u64) []const u8 {
    if (val == NA_U64) return "N/A";
    return std.fmt.bufPrint(buf, "{d}", .{val}) catch "err";
}

fn statusStr(v: u8) []const u8 {
    return switch (v) {
        0 => "pass",
        1 => "warn",
        2 => "fail",
        else => "N/A",
    };
}

fn boolStr(v: bool) []const u8 {
    return if (v) "true" else "false";
}

fn pct(count: usize, total: usize) f64 {
    if (total == 0) return 0.0;
    return @as(f64, @floatFromInt(count)) / @as(f64, @floatFromInt(total)) * 100.0;
}

// ============================================================
// Count summary
// ============================================================

const Counts = struct {
    n_strong: usize = 0,
    n_acceptable: usize = 0,
    n_low: usize = 0,
    n_fail: usize = 0,
    n_mqc_only: usize = 0,
    n_low_ret: usize = 0,
    n_poor_merg: usize = 0,
    n_high_chim: usize = 0,
    n_vlf: usize = 0,
    n_low_qual: usize = 0,
    n_high_adapt: usize = 0,
    n_high_dup: usize = 0,
    n_gc_out: usize = 0,
    n_cont: usize = 0,
    n_trim: usize = 0,
    n_rev: usize = 0,
    n_excl: usize = 0,
};

fn countSamples(samples: []const Sample) Counts {
    var c = Counts{};
    for (samples) |s| {
        switch (s.classification) {
            .pass_strong => c.n_strong += 1,
            .pass_acceptable => c.n_acceptable += 1,
            .low_depth => c.n_low += 1,
            .fail_depth => c.n_fail += 1,
            .multiqc_only => c.n_mqc_only += 1,
        }
        if (s.flag_low_retention) c.n_low_ret += 1;
        if (s.flag_poor_merging) c.n_poor_merg += 1;
        if (s.flag_high_chimera) c.n_high_chim += 1;
        if (s.flag_very_low_filtered) c.n_vlf += 1;
        if (s.flag_low_quality) c.n_low_qual += 1;
        if (s.flag_high_adapter) c.n_high_adapt += 1;
        if (s.flag_high_duplication) c.n_high_dup += 1;
        if (s.flag_gc_outlier) c.n_gc_out += 1;
        switch (s.recommendation) {
            .continue_analysis => c.n_cont += 1,
            .trim_reprocess => c.n_trim += 1,
            .review_manually => c.n_rev += 1,
            .exclude => c.n_excl += 1,
        }
    }
    return c;
}

// ============================================================
// Writer helper: stderr progress
// ============================================================

fn logInfo(io: Io, comptime fmt: []const u8, args: anytype) void {
    var buf: [512]u8 = undefined;
    var fw = Io.File.stderr().writer(io, &buf);
    fw.interface.print(fmt, args) catch {};
    fw.interface.flush() catch {};
}

// ============================================================
// Output: stdout summary
// ============================================================

fn writeStdout(io: Io, samples: []const Sample, cfg: Config, has_dada2: bool, unmatched: []const []u8) !void {
    var buf: [8192]u8 = undefined;
    var fw = Io.File.stdout().writer(io, &buf);
    const w = &fw.interface;

    const n = samples.len;
    const c = countSamples(samples);

    try w.print("\n==========================================================\n", .{});
    try w.print("         MICROBIOME QC SUMMARY REPORT\n", .{});
    try w.print("==========================================================\n\n", .{});
    try w.print("Total samples (files): {d}\n", .{n});

    if (!has_dada2) {
        try w.print("\nWARNING: MultiQC-only mode. DADA2 stats not provided.\n", .{});
        try w.print("  Microbiome usability assessment incomplete.\n", .{});
    }

    try w.print("\n-- DEPTH CLASSIFICATION --\n", .{});
    if (has_dada2) {
        try w.print("  PASS_STRONG     (>={d}): {d}  ({d:.1}%)\n", .{ cfg.pass_strong, c.n_strong, pct(c.n_strong, n) });
        try w.print("  PASS_ACCEPTABLE (>={d}): {d}  ({d:.1}%)\n", .{ cfg.pass_acceptable, c.n_acceptable, pct(c.n_acceptable, n) });
        try w.print("  LOW_DEPTH       (>={d}): {d}  ({d:.1}%)\n", .{ cfg.low_depth, c.n_low, pct(c.n_low, n) });
        try w.print("  FAIL_DEPTH      (<{d}):  {d}  ({d:.1}%)\n", .{ cfg.low_depth, c.n_fail, pct(c.n_fail, n) });
    } else {
        try w.print("  MULTIQC_ONLY: {d} files\n", .{c.n_mqc_only});
    }

    try w.print("\n-- PROCESS FLAGS (DADA2) --\n", .{});
    if (has_dada2) {
        try w.print("  LOW_RETENTION     (<{d:.0}%): {d}\n", .{ cfg.low_retention * 100.0, c.n_low_ret });
        try w.print("  POOR_MERGING      (<{d:.0}%): {d}\n", .{ cfg.poor_merging * 100.0, c.n_poor_merg });
        try w.print("  HIGH_CHIMERA      (>{d:.0}%): {d}\n", .{ cfg.high_chimera * 100.0, c.n_high_chim });
        try w.print("  VERY_LOW_FILTERED (<{d:.0}%): {d}\n", .{ cfg.very_low_filtered * 100.0, c.n_vlf });
    } else {
        try w.print("  (not available without DADA2 stats)\n", .{});
    }

    try w.print("\n-- SEQUENCING FLAGS (MultiQC) --\n", .{});
    try w.print("  LOW_QUALITY:       {d}\n", .{c.n_low_qual});
    try w.print("  HIGH_ADAPTER:      {d}\n", .{c.n_high_adapt});
    try w.print("  HIGH_DUPLICATION:  {d}\n", .{c.n_high_dup});
    try w.print("  GC_OUTLIER (warn): {d}\n", .{c.n_gc_out});

    try w.print("\n-- RECOMMENDATIONS --\n", .{});
    try w.print("  continue:        {d}  ({d:.1}%)\n", .{ c.n_cont, pct(c.n_cont, n) });
    try w.print("  trim/reprocess:  {d}  ({d:.1}%)\n", .{ c.n_trim, pct(c.n_trim, n) });
    try w.print("  review manually: {d}  ({d:.1}%)\n", .{ c.n_rev, pct(c.n_rev, n) });
    try w.print("  exclude:         {d}  ({d:.1}%)\n", .{ c.n_excl, pct(c.n_excl, n) });

    if (unmatched.len > 0) {
        try w.print("\n-- DADA2 IDs UNMATCHED IN MultiQC ({d}) --\n", .{unmatched.len});
        for (unmatched[0..@min(10, unmatched.len)]) |u| try w.print("  {s}\n", .{u});
        if (unmatched.len > 10) try w.print("  ... and {d} more\n", .{unmatched.len - 10});
    }
    try w.print("\n", .{});
    try w.flush();
}

// ============================================================
// Output: Markdown report
// ============================================================

fn writeMarkdown(
    io: Io,
    alloc: Allocator,
    samples: []Sample,
    cfg: Config,
    has_dada2: bool,
    unmatched: []const []u8,
    gs_path: []const u8,
    fq_path: []const u8,
    d2_path: []const u8,
    prefix: []const u8,
) !void {
    const fname = try std.fmt.allocPrint(alloc, "{s}.report.md", .{prefix});
    defer alloc.free(fname);
    const file = try Io.Dir.createFileAbsolute(io, fname, .{});
    defer file.close(io);
    var file_buf: [65536]u8 = undefined;
    var fw = file.writer(io, &file_buf);
    const w = &fw.interface;

    const n = samples.len;
    const c = countSamples(samples);

    try w.print("# Microbiome QC Summary Report\n\n", .{});
    try w.print("**Inputs:**  \n", .{});
    try w.print("- General stats: `{s}`  \n", .{gs_path});
    try w.print("- FastQC stats:  `{s}`  \n", .{fq_path});
    if (has_dada2) {
        try w.print("- DADA2 stats:   `{s}`  \n\n", .{d2_path});
    } else {
        try w.print("- DADA2 stats:   _not provided_ (**MultiQC-only mode**)  \n\n", .{});
        try w.print("> **MultiQC-only mode**: Provide `--dada2-stats` after DADA2 for complete evaluation.\n\n", .{});
    }

    try w.print("## Dataset Overview\n\n**Total samples (files):** {d}\n\n", .{n});

    if (has_dada2) {
        try w.print("### Depth Classification (non-chimeric reads)\n\n", .{});
        try w.print("| Classification | Threshold | Count | Pct |\n|---|---|---:|---:|\n", .{});
        try w.print("| PASS_STRONG | >={d} | {d} | {d:.1}% |\n", .{ cfg.pass_strong, c.n_strong, pct(c.n_strong, n) });
        try w.print("| PASS_ACCEPTABLE | >={d} | {d} | {d:.1}% |\n", .{ cfg.pass_acceptable, c.n_acceptable, pct(c.n_acceptable, n) });
        try w.print("| LOW_DEPTH | >={d} | {d} | {d:.1}% |\n", .{ cfg.low_depth, c.n_low, pct(c.n_low, n) });
        try w.print("| FAIL_DEPTH | <{d} | {d} | {d:.1}% |\n\n", .{ cfg.low_depth, c.n_fail, pct(c.n_fail, n) });

        try w.print("## Process Flags (DADA2)\n\n", .{});
        try w.print("| Flag | Threshold | N Flagged |\n|---|---|---:|\n", .{});
        try w.print("| LOW_RETENTION | <{d:.0}% | {d} |\n", .{ cfg.low_retention * 100.0, c.n_low_ret });
        try w.print("| POOR_MERGING | <{d:.0}% | {d} |\n", .{ cfg.poor_merging * 100.0, c.n_poor_merg });
        try w.print("| HIGH_CHIMERA | >{d:.0}% | {d} |\n", .{ cfg.high_chimera * 100.0, c.n_high_chim });
        try w.print("| VERY_LOW_FILTERED | <{d:.0}% | {d} |\n\n", .{ cfg.very_low_filtered * 100.0, c.n_vlf });
    } else {
        try w.print("_Depth classification unavailable (MultiQC-only mode)._\n\n", .{});
    }

    try w.print("## Sequencing Flags (MultiQC)\n\n", .{});
    try w.print("| Flag | N Flagged |\n|---|---:|\n", .{});
    try w.print("| LOW_QUALITY | {d} |\n", .{c.n_low_qual});
    try w.print("| HIGH_ADAPTER | {d} |\n", .{c.n_high_adapt});
    try w.print("| HIGH_DUPLICATION (>{d:.0}%) | {d} |\n", .{ cfg.high_dup_pct, c.n_high_dup });
    try w.print("| GC_OUTLIER (>{d:.1}SD, warn) | {d} |\n\n", .{ cfg.gc_outlier_sd, c.n_gc_out });

    if (has_dada2) {
        try w.print("## Low/Failing Depth Samples\n\n", .{});
        var found_any = false;
        for (samples) |s| {
            if (s.classification == .fail_depth or s.classification == .low_depth) {
                if (!found_any) {
                    try w.print("| Sample | Non-Chimeric | Class | Recommendation |\n|---|---:|---|---|\n", .{});
                    found_any = true;
                }
                var b: [32]u8 = undefined;
                try w.print("| {s} | {s} | {s} | {s} |\n", .{ s.id, u2s(&b, s.dada2_non_chimeric), s.classification.label(), s.recommendation.label() });
            }
        }
        if (!found_any) try w.print("_None below depth thresholds._\n", .{});
        try w.print("\n", .{});

        // Process flags detail
        const pflag_info = [_]struct { name: []const u8, field: []const u8 }{
            .{ .name = "LOW_RETENTION", .field = "low_retention" },
            .{ .name = "POOR_MERGING", .field = "poor_merging" },
            .{ .name = "HIGH_CHIMERA", .field = "high_chimera" },
            .{ .name = "VERY_LOW_FILTERED", .field = "very_low_filtered" },
        };
        _ = pflag_info;

        inline for ([_][]const u8{ "LOW_RETENTION", "POOR_MERGING", "HIGH_CHIMERA", "VERY_LOW_FILTERED" }) |fname_flag| {
            try w.print("## Process Flag: {s}\n\n", .{fname_flag});
            var found = false;
            for (samples) |s| {
                const flagged = comptime if (std.mem.eql(u8, fname_flag, "LOW_RETENTION")) true else false;
                _ = flagged;
                const is_flagged = blk: {
                    if (comptime std.mem.eql(u8, fname_flag, "LOW_RETENTION")) break :blk s.flag_low_retention;
                    if (comptime std.mem.eql(u8, fname_flag, "POOR_MERGING")) break :blk s.flag_poor_merging;
                    if (comptime std.mem.eql(u8, fname_flag, "HIGH_CHIMERA")) break :blk s.flag_high_chimera;
                    if (comptime std.mem.eql(u8, fname_flag, "VERY_LOW_FILTERED")) break :blk s.flag_very_low_filtered;
                    break :blk false;
                };
                if (is_flagged) {
                    if (!found) {
                        try w.print("| Sample | NonChimeric | Retention | MergeEff | ChimeraRate | FilterRet |\n|---|---:|---:|---:|---:|---:|\n", .{});
                        found = true;
                    }
                    var b1: [32]u8 = undefined;
                    var b2: [32]u8 = undefined;
                    var b3: [32]u8 = undefined;
                    var b4: [32]u8 = undefined;
                    var b5: [32]u8 = undefined;
                    try w.print("| {s} | {s} | {s} | {s} | {s} | {s} |\n", .{
                        s.id,                          u2s(&b1, s.dada2_non_chimeric),
                        f2s(&b2, s.retention_rate, 4), f2s(&b3, s.merge_efficiency, 4),
                        f2s(&b4, s.chimera_rate, 4),   f2s(&b5, s.filter_retention, 4),
                    });
                }
            }
            if (!found) try w.print("_None._\n", .{});
            try w.print("\n", .{});
        }
    }

    try w.print("## Sequencing-Flagged Samples\n\n", .{});
    inline for ([_][]const u8{ "LOW_QUALITY", "HIGH_ADAPTER", "HIGH_DUPLICATION", "GC_OUTLIER" }) |sflag| {
        try w.print("### {s}\n\n", .{sflag});
        var found = false;
        for (samples) |s| {
            const is_flagged = blk: {
                if (comptime std.mem.eql(u8, sflag, "LOW_QUALITY")) break :blk s.flag_low_quality;
                if (comptime std.mem.eql(u8, sflag, "HIGH_ADAPTER")) break :blk s.flag_high_adapter;
                if (comptime std.mem.eql(u8, sflag, "HIGH_DUPLICATION")) break :blk s.flag_high_duplication;
                if (comptime std.mem.eql(u8, sflag, "GC_OUTLIER")) break :blk s.flag_gc_outlier;
                break :blk false;
            };
            if (is_flagged) {
                if (!found) {
                    try w.print("| Sample | %GC | %Dup | TotalSeqs | Adapter | BaseQual | SeqQual |\n|---|---:|---:|---:|---|---|---|\n", .{});
                    found = true;
                }
                var b1: [32]u8 = undefined;
                var b2: [32]u8 = undefined;
                var b3: [32]u8 = undefined;
                try w.print("| {s} | {s} | {s} | {s} | {s} | {s} | {s} |\n", .{
                    s.id,                                f2s(&b1, s.pct_gc, 1),        f2s(&b2, s.pct_duplicates, 1),
                    f2s(&b3, s.total_sequences, 0),      statusStr(s.adapter_content), statusStr(s.per_base_seq_quality),
                    statusStr(s.per_seq_quality_scores),
                });
            }
        }
        if (!found) try w.print("_None._\n", .{});
        try w.print("\n", .{});
    }

    // Top/bottom tables
    const n_top = @min(10, n);
    var sorted_copy = try alloc.dupe(Sample, samples);
    defer alloc.free(sorted_copy);

    if (has_dada2) {
        std.mem.sort(Sample, sorted_copy, {}, cmpNcDesc);
        try w.print("## Top {d} Samples (Non-Chimeric Reads)\n\n| Rank | Sample | NonChimeric | Retention | Class |\n|---:|---|---:|---:|---|\n", .{n_top});
        for (sorted_copy[0..n_top], 1..) |s, rank| {
            var b1: [32]u8 = undefined;
            var b2: [32]u8 = undefined;
            try w.print("| {d} | {s} | {s} | {s} | {s} |\n", .{ rank, s.id, u2s(&b1, s.dada2_non_chimeric), f2s(&b2, s.retention_rate, 4), s.classification.label() });
        }
        try w.print("\n", .{});
        std.mem.sort(Sample, sorted_copy, {}, cmpNcAsc);
        try w.print("## Bottom {d} Samples (Non-Chimeric Reads)\n\n| Rank | Sample | NonChimeric | Retention | Class |\n|---:|---|---:|---:|---|\n", .{n_top});
        for (sorted_copy[0..n_top], 1..) |s, rank| {
            var b1: [32]u8 = undefined;
            var b2: [32]u8 = undefined;
            try w.print("| {d} | {s} | {s} | {s} | {s} |\n", .{ rank, s.id, u2s(&b1, s.dada2_non_chimeric), f2s(&b2, s.retention_rate, 4), s.classification.label() });
        }
        try w.print("\n", .{});
    } else {
        std.mem.sort(Sample, sorted_copy, {}, cmpTotalDesc);
        try w.print("## Top {d} by Total Sequences\n\n| Rank | Sample | TotalSeqs | %GC | %Dup |\n|---:|---|---:|---:|---:|\n", .{n_top});
        for (sorted_copy[0..n_top], 1..) |s, rank| {
            var b1: [32]u8 = undefined;
            var b2: [32]u8 = undefined;
            var b3: [32]u8 = undefined;
            try w.print("| {d} | {s} | {s} | {s} | {s} |\n", .{ rank, s.id, f2s(&b1, s.total_sequences, 0), f2s(&b2, s.pct_gc, 1), f2s(&b3, s.pct_duplicates, 1) });
        }
        try w.print("\n", .{});
    }

    try w.print("## Recommended Actions\n\n| Action | Count | Pct |\n|---|---:|---:|\n", .{});
    try w.print("| continue | {d} | {d:.1}% |\n", .{ c.n_cont, pct(c.n_cont, n) });
    try w.print("| trim/reprocess | {d} | {d:.1}% |\n", .{ c.n_trim, pct(c.n_trim, n) });
    try w.print("| review manually | {d} | {d:.1}% |\n", .{ c.n_rev, pct(c.n_rev, n) });
    try w.print("| exclude | {d} | {d:.1}% |\n\n", .{ c.n_excl, pct(c.n_excl, n) });

    try w.print("### Overall Recommendation\n\n", .{});
    if (c.n_excl > n / 4) {
        try w.print("**ACTION REQUIRED**: >25% samples flagged for exclusion. Review and consider reprocessing.\n\n", .{});
    } else if (c.n_trim > 0) {
        try w.print("**TRIM/REPROCESS**: Sequencing quality issues found. Consider quality trimming.\n\n", .{});
    } else if ((c.n_low + c.n_fail > 0) and has_dada2) {
        try w.print("**REVIEW**: Some samples have low depth. Review for exclusion or inclusion with caveats.\n\n", .{});
    } else {
        try w.print("**PROCEED**: Dataset suitable for downstream analysis. Apply standard rarefaction/normalization.\n\n", .{});
    }

    if (unmatched.len > 0) {
        try w.print("## DADA2 Samples Not Matched in MultiQC ({d})\n\n", .{unmatched.len});
        for (unmatched) |u| try w.print("- `{s}`\n", .{u});
        try w.print("\n", .{});
    }

    try w.print("## Thresholds Used\n\n| Parameter | Value |\n|---|---|\n", .{});
    try w.print("| pass_strong (>=) | {d} |\n", .{cfg.pass_strong});
    try w.print("| pass_acceptable (>=) | {d} |\n", .{cfg.pass_acceptable});
    try w.print("| low_depth (>=) | {d} |\n", .{cfg.low_depth});
    try w.print("| low_retention (<) | {d:.2} |\n", .{cfg.low_retention});
    try w.print("| poor_merging (<) | {d:.2} |\n", .{cfg.poor_merging});
    try w.print("| high_chimera (>) | {d:.2} |\n", .{cfg.high_chimera});
    try w.print("| very_low_filtered (<) | {d:.2} |\n", .{cfg.very_low_filtered});
    try w.print("| high_dup_pct (>) | {d:.1}% |\n", .{cfg.high_dup_pct});
    try w.print("| gc_outlier_sd | {d:.1} |\n", .{cfg.gc_outlier_sd});

    try w.flush();
}

// ============================================================
// Output: samples TSV
// ============================================================

fn writeSamplesTsv(io: Io, alloc: Allocator, samples: []const Sample, has_dada2: bool, prefix: []const u8) !void {
    const fname = try std.fmt.allocPrint(alloc, "{s}.samples.tsv", .{prefix});
    defer alloc.free(fname);
    const file = try Io.Dir.createFileAbsolute(io, fname, .{});
    defer file.close(io);
    var file_buf: [65536]u8 = undefined;
    var fw = file.writer(io, &file_buf);
    const w = &fw.interface;

    // Header
    try w.print("sample_id\tbase_id\tis_r1\tpct_duplicates\tpct_gc\tavg_seq_len\tmedian_seq_len\tpct_fails\ttotal_sequences\t", .{});
    try w.print("basic_statistics\tper_base_seq_quality\tper_tile_seq_quality\tper_seq_quality_scores\t", .{});
    try w.print("per_base_seq_content\tper_seq_gc_content\tper_base_n_content\tseq_len_distribution\t", .{});
    try w.print("seq_duplication_levels\toverrepresented_seqs\tadapter_content\t", .{});
    if (has_dada2) {
        try w.print("dada2_input\tdada2_filtered\tdada2_denoised\tdada2_merged\tdada2_non_chimeric\t", .{});
        try w.print("retention_rate\tfilter_retention\tdenoise_retention\tmerge_efficiency\tchimera_rate\t", .{});
    }
    try w.print("classification\t", .{});
    if (has_dada2) {
        try w.print("flag_low_retention\tflag_poor_merging\tflag_high_chimera\tflag_very_low_filtered\t", .{});
    }
    try w.print("flag_low_quality\tflag_high_adapter\tflag_high_duplication\tflag_gc_outlier\trecommendation\n", .{});

    for (samples) |s| {
        var b1: [32]u8 = undefined;
        var b2: [32]u8 = undefined;
        var b3: [32]u8 = undefined;
        var b4: [32]u8 = undefined;
        var b5: [32]u8 = undefined;
        var b6: [32]u8 = undefined;
        try w.print("{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t", .{
            s.id,                          s.base_id,                if (s.is_r1) "R1" else "R2",
            f2s(&b1, s.pct_duplicates, 4), f2s(&b2, s.pct_gc, 1),    f2s(&b3, s.avg_seq_len, 1),
            f2s(&b4, s.median_seq_len, 1), f2s(&b5, s.pct_fails, 4), f2s(&b6, s.total_sequences, 0),
        });
        try w.print("{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t", .{
            statusStr(s.basic_statistics),       statusStr(s.per_base_seq_quality),
            statusStr(s.per_tile_seq_quality),   statusStr(s.per_seq_quality_scores),
            statusStr(s.per_base_seq_content),   statusStr(s.per_seq_gc_content),
            statusStr(s.per_base_n_content),     statusStr(s.seq_len_distribution),
            statusStr(s.seq_duplication_levels), statusStr(s.overrepresented_seqs),
            statusStr(s.adapter_content),
        });
        if (has_dada2) {
            var c1: [32]u8 = undefined;
            var c2: [32]u8 = undefined;
            var c3: [32]u8 = undefined;
            var c4: [32]u8 = undefined;
            var c5: [32]u8 = undefined;
            var d1: [32]u8 = undefined;
            var d2: [32]u8 = undefined;
            var d3: [32]u8 = undefined;
            var d4: [32]u8 = undefined;
            var d5: [32]u8 = undefined;
            try w.print("{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t{s}\t", .{
                u2s(&c1, s.dada2_input),         u2s(&c2, s.dada2_filtered),
                u2s(&c3, s.dada2_denoised),      u2s(&c4, s.dada2_merged),
                u2s(&c5, s.dada2_non_chimeric),  f2s(&d1, s.retention_rate, 6),
                f2s(&d2, s.filter_retention, 6), f2s(&d3, s.denoise_retention, 6),
                f2s(&d4, s.merge_efficiency, 6), f2s(&d5, s.chimera_rate, 6),
            });
        }
        try w.print("{s}\t", .{s.classification.label()});
        if (has_dada2) {
            try w.print("{s}\t{s}\t{s}\t{s}\t", .{
                boolStr(s.flag_low_retention), boolStr(s.flag_poor_merging),
                boolStr(s.flag_high_chimera),  boolStr(s.flag_very_low_filtered),
            });
        }
        try w.print("{s}\t{s}\t{s}\t{s}\t{s}\n", .{
            boolStr(s.flag_low_quality),      boolStr(s.flag_high_adapter),
            boolStr(s.flag_high_duplication), boolStr(s.flag_gc_outlier),
            s.recommendation.label(),
        });
    }
    try w.flush();
}

// ============================================================
// Output: flags TSV (melted)
// ============================================================

fn writeFlagsTsv(io: Io, alloc: Allocator, samples: []const Sample, has_dada2: bool, prefix: []const u8) !void {
    const fname = try std.fmt.allocPrint(alloc, "{s}.flags.tsv", .{prefix});
    defer alloc.free(fname);
    const file = try Io.Dir.createFileAbsolute(io, fname, .{});
    defer file.close(io);
    var file_buf: [32768]u8 = undefined;
    var fw = file.writer(io, &file_buf);
    const w = &fw.interface;

    try w.print("sample_id\tflag_type\tflag_name\tvalue\n", .{});
    for (samples) |s| {
        if (has_dada2) {
            if (s.flag_low_retention) {
                var b: [32]u8 = undefined;
                try w.print("{s}\tprocess\tLOW_RETENTION\t{s}\n", .{ s.id, f2s(&b, s.retention_rate, 6) });
            }
            if (s.flag_poor_merging) {
                var b: [32]u8 = undefined;
                try w.print("{s}\tprocess\tPOOR_MERGING\t{s}\n", .{ s.id, f2s(&b, s.merge_efficiency, 6) });
            }
            if (s.flag_high_chimera) {
                var b: [32]u8 = undefined;
                try w.print("{s}\tprocess\tHIGH_CHIMERA\t{s}\n", .{ s.id, f2s(&b, s.chimera_rate, 6) });
            }
            if (s.flag_very_low_filtered) {
                var b: [32]u8 = undefined;
                try w.print("{s}\tprocess\tVERY_LOW_FILTERED\t{s}\n", .{ s.id, f2s(&b, s.filter_retention, 6) });
            }
        }
        if (s.flag_low_quality) try w.print("{s}\tsequencing\tLOW_QUALITY\tbase={s},seq={s}\n", .{ s.id, statusStr(s.per_base_seq_quality), statusStr(s.per_seq_quality_scores) });
        if (s.flag_high_adapter) try w.print("{s}\tsequencing\tHIGH_ADAPTER\t{s}\n", .{ s.id, statusStr(s.adapter_content) });
        if (s.flag_high_duplication) {
            var b: [32]u8 = undefined;
            try w.print("{s}\tsequencing\tHIGH_DUPLICATION\t{s}\n", .{ s.id, f2s(&b, s.pct_duplicates, 2) });
        }
        if (s.flag_gc_outlier) {
            var b: [32]u8 = undefined;
            try w.print("{s}\tsequencing\tGC_OUTLIER\t{s}\n", .{ s.id, f2s(&b, s.pct_gc, 1) });
        }
    }
    try w.flush();
}

// ============================================================
// CLI
// ============================================================

fn printHelp(io: Io) void {
    var buf: [4096]u8 = undefined;
    var fw = Io.File.stderr().writer(io, &buf);
    const w = &fw.interface;
    w.print(
        \\Usage: qc_summary [options]
        \\
        \\Required:
        \\  --multiqc-general-stats <path>   multiqc_general_stats.txt
        \\  --multiqc-fastqc <path>          multiqc_fastqc.txt
        \\
        \\Optional:
        \\  --dada2-stats <path>             DADA2 denoising stats TSV
        \\                                   (columns: sample-id, input, filtered,
        \\                                    denoised, merged, non-chimeric)
        \\  --out-prefix <prefix>            Output prefix (default: microbiome_qc_summary)
        \\
        \\Thresholds:
        \\  --pass-strong <N>      non-chimeric >= N => PASS_STRONG    (default: 20000)
        \\  --pass-acceptable <N>  non-chimeric >= N => PASS_ACCEPTABLE (default: 10000)
        \\  --low-depth <N>        non-chimeric >= N => LOW_DEPTH       (default: 5000)
        \\  --low-retention <f>    flag if retention_rate < f          (default: 0.10)
        \\  --poor-merging <f>     flag if merge_efficiency < f        (default: 0.60)
        \\  --high-chimera <f>     flag if chimera_rate > f            (default: 0.30)
        \\  --very-low-filtered <f> flag if filter_retention < f       (default: 0.20)
        \\  --high-dup-pct <f>     flag if pct_duplicates > f          (default: 90.0)
        \\  --gc-outlier-sd <f>    GC outlier SD threshold             (default: 3.0)
        \\  --help                 Show this help
        \\
    , .{}) catch {};
    w.flush() catch {};
}

fn parseArgs(io: Io, alloc: Allocator, args: []const [:0]const u8) !Config {
    var cfg = Config{};
    var i: usize = 1;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelp(io);
            std.process.exit(0);
        }
        const flags = [_][]const u8{
            "--multiqc-general-stats", "--multiqc-fastqc",  "--dada2-stats",       "--out-prefix",
            "--pass-strong",           "--pass-acceptable", "--low-depth",         "--low-retention",
            "--poor-merging",          "--high-chimera",    "--very-low-filtered", "--high-dup-pct",
            "--gc-outlier-sd",
        };
        var matched = false;
        for (flags) |flag| {
            if (std.mem.eql(u8, arg, flag)) {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("ERROR: {s} requires a value\n", .{flag});
                    std.process.exit(1);
                }
                const val = args[i];
                if (std.mem.eql(u8, flag, "--multiqc-general-stats")) {
                    cfg.multiqc_general_stats = try alloc.dupe(u8, val);
                } else if (std.mem.eql(u8, flag, "--multiqc-fastqc")) {
                    cfg.multiqc_fastqc = try alloc.dupe(u8, val);
                } else if (std.mem.eql(u8, flag, "--dada2-stats")) {
                    cfg.dada2_stats = try alloc.dupe(u8, val);
                } else if (std.mem.eql(u8, flag, "--out-prefix")) {
                    cfg.out_prefix = try alloc.dupe(u8, val);
                } else if (std.mem.eql(u8, flag, "--pass-strong")) {
                    cfg.pass_strong = std.fmt.parseInt(u64, val, 10) catch {
                        std.debug.print("ERROR: bad --pass-strong\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--pass-acceptable")) {
                    cfg.pass_acceptable = std.fmt.parseInt(u64, val, 10) catch {
                        std.debug.print("ERROR: bad --pass-acceptable\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--low-depth")) {
                    cfg.low_depth = std.fmt.parseInt(u64, val, 10) catch {
                        std.debug.print("ERROR: bad --low-depth\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--low-retention")) {
                    cfg.low_retention = std.fmt.parseFloat(f64, val) catch {
                        std.debug.print("ERROR: bad --low-retention\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--poor-merging")) {
                    cfg.poor_merging = std.fmt.parseFloat(f64, val) catch {
                        std.debug.print("ERROR: bad --poor-merging\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--high-chimera")) {
                    cfg.high_chimera = std.fmt.parseFloat(f64, val) catch {
                        std.debug.print("ERROR: bad --high-chimera\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--very-low-filtered")) {
                    cfg.very_low_filtered = std.fmt.parseFloat(f64, val) catch {
                        std.debug.print("ERROR: bad --very-low-filtered\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--high-dup-pct")) {
                    cfg.high_dup_pct = std.fmt.parseFloat(f64, val) catch {
                        std.debug.print("ERROR: bad --high-dup-pct\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--gc-outlier-sd")) {
                    cfg.gc_outlier_sd = std.fmt.parseFloat(f64, val) catch {
                        std.debug.print("ERROR: bad --gc-outlier-sd\n", .{});
                        std.process.exit(1);
                    };
                }
                matched = true;
                break;
            }
        }
        if (!matched) std.debug.print("WARNING: unknown argument ignored: {s}\n", .{arg});
    }

    if (cfg.multiqc_general_stats.len == 0) {
        std.debug.print("ERROR: --multiqc-general-stats is required\n", .{});
        printHelp(io);
        std.process.exit(1);
    }
    if (cfg.multiqc_fastqc.len == 0) {
        std.debug.print("ERROR: --multiqc-fastqc is required\n", .{});
        printHelp(io);
        std.process.exit(1);
    }
    return cfg;
}

fn resolvePath(io: Io, alloc: Allocator, path: []const u8) ![]u8 {
    if (std.fs.path.isAbsolute(path)) return alloc.dupe(u8, path);
    var cwd_buf: [std.fs.max_path_bytes]u8 = undefined;
    const len = std.process.currentPath(io, &cwd_buf) catch return alloc.dupe(u8, path);
    const cwd = cwd_buf[0..len];
    return std.fs.path.join(alloc, &.{ cwd, path });
}

// ============================================================
// Entry point
// ============================================================

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const alloc = init.gpa;
    const arena = init.arena.allocator();

    const args = try init.minimal.args.toSlice(arena);
    const cfg = try parseArgs(io, alloc, args);

    const gs_path = try resolvePath(io, alloc, cfg.multiqc_general_stats);
    defer alloc.free(gs_path);
    const fq_path = try resolvePath(io, alloc, cfg.multiqc_fastqc);
    defer alloc.free(fq_path);
    const d2_path = if (cfg.dada2_stats.len > 0) try resolvePath(io, alloc, cfg.dada2_stats) else try alloc.dupe(u8, "");
    defer alloc.free(d2_path);
    const prefix_path = try resolvePath(io, alloc, cfg.out_prefix);
    defer alloc.free(prefix_path);

    const has_dada2 = d2_path.len > 0;

    var samples: AL(Sample) = .empty;
    defer {
        for (samples.items) |s| {
            alloc.free(s.id);
            alloc.free(s.base_id);
        }
        samples.deinit(alloc);
    }
    var sample_map = std.StringHashMap(usize).init(alloc);
    defer sample_map.deinit();

    logInfo(io, "Parsing: {s}\n", .{gs_path});
    try parseGeneralStats(io, alloc, gs_path, &sample_map, &samples);
    logInfo(io, "  -> {d} samples loaded\n", .{samples.items.len});

    logInfo(io, "Parsing: {s}\n", .{fq_path});
    try parseFastqcStats(io, alloc, fq_path, &sample_map, &samples);

    var unmatched: AL([]u8) = .empty;
    defer {
        for (unmatched.items) |u| alloc.free(u);
        unmatched.deinit(alloc);
    }

    if (has_dada2) {
        logInfo(io, "Parsing: {s}\n", .{d2_path});
        try parseDada2Stats(io, alloc, d2_path, &sample_map, &samples, &unmatched);
        if (unmatched.items.len > 0) {
            logInfo(io, "  WARNING: {d} DADA2 IDs not matched in MultiQC\n", .{unmatched.items.len});
        }
    } else {
        logInfo(io, "No DADA2 stats -- MultiQC-only mode.\n", .{});
    }

    computeDerivedMetrics(samples.items);
    classifySamples(samples.items, cfg);
    assignFlags(samples.items, cfg);
    assignRecommendations(samples.items);

    logInfo(io, "Writing outputs: {s}.*\n", .{prefix_path});
    try writeStdout(io, samples.items, cfg, has_dada2, unmatched.items);
    try writeMarkdown(io, alloc, samples.items, cfg, has_dada2, unmatched.items, gs_path, fq_path, d2_path, prefix_path);
    try writeSamplesTsv(io, alloc, samples.items, has_dada2, prefix_path);
    try writeFlagsTsv(io, alloc, samples.items, has_dada2, prefix_path);

    logInfo(io, "Done:\n  {s}.report.md\n  {s}.samples.tsv\n  {s}.flags.tsv\n", .{ prefix_path, prefix_path, prefix_path });
}
