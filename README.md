# 16S MultiQC + DADA2 QC Summary (Zig)

CLI tool written in Zig to integrate MultiQC and DADA2 outputs for microbiome QC reporting.

## Goals
- Correct sample-level aggregation (not file-level)
- Safe handling of missing data
- Biologically meaningful QC metrics

## Current Status
- Single-file prototype: qc_summary.zig
- Under active refactoring

## Usage (dev)
zig run qc_summary.zig
