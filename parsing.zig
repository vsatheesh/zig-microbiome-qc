//! TSV/CSV parsing, sample-ID normalization, aggregation, and CLI argument parsing.

const std = @import("std");
const Io = std.Io;
const Allocator = std.mem.Allocator;
const types = @import("types.zig");
const val = @import("validation.zig");

pub const AL = types.AL;
pub const Row = types.Row;
const NA_F64 = types.NA_F64;
const NA_U64 = types.NA_U64;
const Config = types.Config;
const ReadSide = types.ReadSide;
const NormalizedSampleId = types.NormalizedSampleId;
const MultiqcFileRow = types.MultiqcFileRow;
const SampleRecord = types.SampleRecord;

// ============================================================
// String / ID helpers
// ============================================================

fn stripSuffix(name: []const u8, suffix: []const u8) ?[]const u8 {
    if (!std.mem.endsWith(u8, name, suffix)) return null;
    return name[0 .. name.len - suffix.len];
}

fn stripKnownFileSuffixes(name: []const u8) []const u8 {
    var base = std.mem.trim(u8, name, " \t\r\n");
    while (true) {
        if (stripSuffix(base, "_fastqc.zip")) |v| { base = v; continue; }
        if (stripSuffix(base, "_fastqc"))     |v| { base = v; continue; }
        if (stripSuffix(base, ".fastq.gz"))   |v| { base = v; continue; }
        if (stripSuffix(base, ".fq.gz"))      |v| { base = v; continue; }
        if (stripSuffix(base, ".fastq"))      |v| { base = v; continue; }
        if (stripSuffix(base, ".fq"))         |v| { base = v; continue; }
        if (stripSuffix(base, ".zip"))        |v| { base = v; continue; }
        break;
    }
    return base;
}

fn isDigit(c: u8) bool {
    return c >= '0' and c <= '9';
}

fn stripLaneSuffix(name: []const u8) []const u8 {
    if (name.len < 5) return name;
    const lane = name[name.len - 5 ..];
    if (lane[0] == '_' and lane[1] == 'L' and isDigit(lane[2]) and isDigit(lane[3]) and isDigit(lane[4])) {
        return name[0 .. name.len - 5];
    }
    return name;
}

pub fn normalizeSampleId(name: []const u8) NormalizedSampleId {
    var base = stripKnownFileSuffixes(name);
    var side: ReadSide = .unknown;
    const patterns = [_]struct { suffix: []const u8, read_side: ReadSide }{
        .{ .suffix = "_R1_001", .read_side = .r1 },
        .{ .suffix = "_R2_001", .read_side = .r2 },
        .{ .suffix = ".R1_001", .read_side = .r1 },
        .{ .suffix = ".R2_001", .read_side = .r2 },
        .{ .suffix = "_R1",     .read_side = .r1 },
        .{ .suffix = "_R2",     .read_side = .r2 },
        .{ .suffix = ".R1",     .read_side = .r1 },
        .{ .suffix = ".R2",     .read_side = .r2 },
        .{ .suffix = "_1",      .read_side = .r1 },
        .{ .suffix = "_2",      .read_side = .r2 },
    };
    for (patterns) |p| {
        if (stripSuffix(base, p.suffix)) |v| {
            base = v;
            side = p.read_side;
            break;
        }
    }
    base = stripLaneSuffix(base);
    return .{ .base = base, .read_side = side };
}

pub fn divSafe(a: f64, b: f64) f64 {
    if (b == 0.0 or b == NA_F64 or a == NA_F64) return NA_F64;
    return a / b;
}

// Strip trailing -R<digits> from manifest sample-ids (e.g. "Sample_1-R2" -> "Sample_1").
pub fn stripManifestReplicateSuffix(id: []const u8) []const u8 {
    var i: usize = id.len;
    while (i > 0 and isDigit(id[i - 1])) i -= 1;
    if (i > 0 and id[i - 1] == 'R') i -= 1;
    if (i > 0 and id[i - 1] == '-') return id[0 .. i - 1];
    return id;
}

// ============================================================
// TSV / CSV low-level parsing
// ============================================================

pub fn parseTsv(
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
        const line = std.mem.trimEnd(u8, raw, "\r\n");
        if (std.mem.trim(u8, line, " \t").len == 0) continue;

        var fields: Row = .empty;
        // Split on tab only; fields are stored verbatim (no trimming) to
        // preserve exact TSV column boundaries.
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

pub fn buildColIndex(alloc: Allocator, headers: []const []u8) !std.StringHashMap(usize) {
    var map = std.StringHashMap(usize).init(alloc);
    for (headers, 0..) |h, i| {
        const key = std.mem.trim(u8, h, " \t\r\n");
        try map.put(key, i);
    }
    return map;
}

pub fn getField(row: []const []u8, col_map: *const std.StringHashMap(usize), name: []const u8) []const u8 {
    if (col_map.get(name)) |idx| {
        if (idx < row.len) return std.mem.trim(u8, row[idx], " \t\r\n");
    }
    return "";
}

pub fn mergeStatus(current: u8, next: u8) u8 {
    if (next == 255) return current;
    if (current == 255) return next;
    return if (next > current) next else current;
}

pub fn findFileRowIndex(
    file_map: *const std.StringHashMap(usize),
    file_rows: []const MultiqcFileRow,
    raw_id: []const u8,
) ?usize {
    if (file_map.get(raw_id)) |idx| return idx;

    const norm = normalizeSampleId(raw_id);
    var found: ?usize = null;
    for (file_rows, 0..) |r, i| {
        if (!std.mem.eql(u8, r.base_id, norm.base)) continue;
        if (norm.read_side != .unknown and r.read_side != norm.read_side) continue;
        if (found != null) return null;
        found = i;
    }
    return found;
}

// ============================================================
// Main parsing functions
// ============================================================

pub fn parseGeneralStats(
    io: Io,
    alloc: Allocator,
    path: []const u8,
    file_map: *std.StringHashMap(usize),
    file_rows: *AL(MultiqcFileRow),
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
    try val.validateRequiredColumns(&col, path, &[_][]const u8{"Sample"});

    for (rows.items, 0..) |row, row_i| {
        if (row.items.len == 0) continue;
        const raw_id = getField(row.items, &col, "Sample");
        if (raw_id.len == 0) {
            std.debug.print("ERROR: {s}: missing Sample value on data row {d}\n", .{ path, row_i + 2 });
            return error.MissingRequiredValue;
        }
        const norm = normalizeSampleId(raw_id);
        var r = MultiqcFileRow{
            .id = try alloc.dupe(u8, raw_id),
            .base_id = try alloc.dupe(u8, norm.base),
            .read_side = norm.read_side,
        };
        r.pct_duplicates = val.parseF64(getField(row.items, &col, "percent_duplicates"));
        r.pct_gc = val.parseF64Validated(getField(row.items, &col, "percent_gc"), raw_id, "percent_gc");
        r.avg_seq_len = val.parseF64(getField(row.items, &col, "avg_sequence_length"));
        r.median_seq_len = val.parseF64(getField(row.items, &col, "median_sequence_length"));
        r.pct_fails = val.parseF64(getField(row.items, &col, "percent_fails"));
        r.total_sequences = val.parseF64Validated(getField(row.items, &col, "total_sequences"), raw_id, "total_sequences");
        const idx = file_rows.items.len;
        try file_rows.append(alloc, r);
        try file_map.put(file_rows.items[idx].id, idx);
    }
}

pub fn parseFastqcStats(
    io: Io,
    alloc: Allocator,
    path: []const u8,
    file_map: *std.StringHashMap(usize),
    file_rows: *AL(MultiqcFileRow),
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
    try val.validateRequiredColumns(&col, path, &[_][]const u8{"Sample"});

    for (rows.items, 0..) |row, row_i| {
        if (row.items.len == 0) continue;
        const raw_id = getField(row.items, &col, "Sample");
        if (raw_id.len == 0) {
            std.debug.print("ERROR: {s}: missing Sample value on data row {d}\n", .{ path, row_i + 2 });
            return error.MissingRequiredValue;
        }
        const idx = findFileRowIndex(file_map, file_rows.items, raw_id) orelse continue;
        var r = &file_rows.items[idx];
        r.basic_statistics       = mergeStatus(r.basic_statistics,       val.parseStatus(getField(row.items, &col, "basic_statistics")));
        r.per_base_seq_quality   = mergeStatus(r.per_base_seq_quality,   val.parseStatus(getField(row.items, &col, "per_base_sequence_quality")));
        r.per_tile_seq_quality   = mergeStatus(r.per_tile_seq_quality,   val.parseStatus(getField(row.items, &col, "per_tile_sequence_quality")));
        r.per_seq_quality_scores = mergeStatus(r.per_seq_quality_scores, val.parseStatus(getField(row.items, &col, "per_sequence_quality_scores")));
        r.per_base_seq_content   = mergeStatus(r.per_base_seq_content,   val.parseStatus(getField(row.items, &col, "per_base_sequence_content")));
        r.per_seq_gc_content     = mergeStatus(r.per_seq_gc_content,     val.parseStatus(getField(row.items, &col, "per_sequence_gc_content")));
        r.per_base_n_content     = mergeStatus(r.per_base_n_content,     val.parseStatus(getField(row.items, &col, "per_base_n_content")));
        r.seq_len_distribution   = mergeStatus(r.seq_len_distribution,   val.parseStatus(getField(row.items, &col, "sequence_length_distribution")));
        r.seq_duplication_levels = mergeStatus(r.seq_duplication_levels, val.parseStatus(getField(row.items, &col, "sequence_duplication_levels")));
        r.overrepresented_seqs   = mergeStatus(r.overrepresented_seqs,   val.parseStatus(getField(row.items, &col, "overrepresented_sequences")));
        r.adapter_content        = mergeStatus(r.adapter_content,        val.parseStatus(getField(row.items, &col, "adapter_content")));
    }
}

// Parse a QIIME2-style manifest CSV (sample-id,absolute-filepath,direction).
// Builds a map from file-stem -> DADA2 sample ID (manifest sample-id with -R<N> stripped).
pub fn parseManifestCsv(
    io: Io,
    alloc: Allocator,
    path: []const u8,
    out: *std.StringHashMap([]u8),
) !void {
    const file = Io.Dir.openFileAbsolute(io, path, .{}) catch |err| {
        std.debug.print("ERROR: Cannot open manifest: {s} ({any})\n", .{ path, err });
        return err;
    };
    defer file.close(io);

    var read_buf: [131072]u8 = undefined;
    var rdr = file.reader(io, &read_buf);
    var first_line = true;
    var sample_id_col: usize = 0;
    var filepath_col: usize = 1;

    while (true) {
        const maybe_line = rdr.interface.takeDelimiter('\n') catch |err| switch (err) {
            error.ReadFailed => return error.ReadFailed,
            error.StreamTooLong => continue,
        };
        const raw = maybe_line orelse break;
        const line = std.mem.trimEnd(u8, raw, "\r\n");
        if (std.mem.trim(u8, line, " \t").len == 0) continue;

        var fields_buf: [8][]const u8 = undefined;
        var n_fields: usize = 0;
        var it = std.mem.splitScalar(u8, line, ',');
        while (it.next()) |f| {
            if (n_fields < fields_buf.len) {
                fields_buf[n_fields] = f;
                n_fields += 1;
            }
        }
        if (n_fields == 0) continue;

        if (first_line) {
            first_line = false;
            for (fields_buf[0..n_fields], 0..) |h, i| {
                const hh = std.mem.trim(u8, h, " \t\r\n");
                if (std.mem.eql(u8, hh, "sample-id") or std.mem.eql(u8, hh, "sample_id")) sample_id_col = i;
                if (std.mem.eql(u8, hh, "absolute-filepath") or std.mem.eql(u8, hh, "filepath") or std.mem.eql(u8, hh, "file-path")) filepath_col = i;
            }
            continue;
        }

        if (n_fields <= sample_id_col or n_fields <= filepath_col) continue;

        const manifest_sample_id = std.mem.trim(u8, fields_buf[sample_id_col], " \t\r\n");
        const filepath = std.mem.trim(u8, fields_buf[filepath_col], " \t\r\n");
        if (manifest_sample_id.len == 0 or filepath.len == 0) continue;

        const basename = std.fs.path.basename(filepath);
        const file_stem = stripKnownFileSuffixes(basename);
        const dada2_id = stripManifestReplicateSuffix(manifest_sample_id);

        if (out.contains(file_stem)) continue;
        const key = try alloc.dupe(u8, file_stem);
        const val_str = try alloc.dupe(u8, dada2_id);
        try out.put(key, val_str);
    }
}

// Replace each file row's base_id with the DADA2 sample ID from the manifest map.
pub fn remapBaseIdsWithManifest(
    alloc: Allocator,
    file_rows: []MultiqcFileRow,
    manifest_map: *const std.StringHashMap([]u8),
) !void {
    for (file_rows) |*row| {
        if (manifest_map.get(row.id)) |mapped_id| {
            alloc.free(row.base_id);
            row.base_id = try alloc.dupe(u8, mapped_id);
        }
    }
}

pub fn applyDada2(row: []const []u8, col: *const std.StringHashMap(usize), s: *SampleRecord, sample_id: []const u8) void {
    s.has_dada2 = true;
    s.dada2_input        = val.parseU64Validated(getField(row, col, "input"),        sample_id, "input");
    s.dada2_filtered     = val.parseU64Validated(getField(row, col, "filtered"),     sample_id, "filtered");
    s.dada2_denoised     = val.parseU64Validated(getField(row, col, "denoised"),     sample_id, "denoised");
    s.dada2_merged       = val.parseU64Validated(getField(row, col, "merged"),       sample_id, "merged");
    s.dada2_non_chimeric = val.parseU64Validated(getField(row, col, "non-chimeric"), sample_id, "non-chimeric");
}

pub fn parseDada2Stats(
    io: Io,
    alloc: Allocator,
    path: []const u8,
    sample_map: *const std.StringHashMap(usize),
    samples: []SampleRecord,
    unmatched_dada2: *AL([]u8),
    duplicate_dada2: *AL([]u8),
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

    const id_col = try val.requireDada2IdColumn(&col, path);
    const required_counts = [_][]const u8{ "input", "filtered", "denoised", "merged", "non-chimeric" };
    try val.validateRequiredColumns(&col, path, &required_counts);

    for (rows.items, 0..) |row, row_i| {
        if (row.items.len == 0) continue;
        const raw_id = getField(row.items, &col, id_col);
        if (raw_id.len == 0) {
            std.debug.print("ERROR: {s}: missing DADA2 sample id on data row {d}\n", .{ path, row_i + 2 });
            return error.MissingRequiredValue;
        }
        // Skip QIIME2 metadata directive rows (e.g. "#q2:types").
        if (std.mem.startsWith(u8, raw_id, "#")) continue;

        // Try exact match first (needed when manifest has already remapped base IDs to
        // DADA2-style IDs like "Sample_1"), then fall back to normalized match so that
        // stripping R1/R2 suffixes still works for the no-manifest code path.
        const idx_opt = sample_map.get(raw_id) orelse blk: {
            const norm = normalizeSampleId(raw_id);
            break :blk sample_map.get(norm.base);
        };

        if (idx_opt) |idx| {
            const s = &samples[idx];
            if (s.has_dada2) {
                try duplicate_dada2.append(alloc, try alloc.dupe(u8, raw_id));
                continue;
            }
            applyDada2(row.items, &col, s, raw_id);
        } else {
            try unmatched_dada2.append(alloc, try alloc.dupe(u8, raw_id));
        }
    }
}

fn applyMultiqcRowToSample(row: MultiqcFileRow, s: *SampleRecord) void {
    s.file_count += 1;
    switch (row.read_side) {
        .r1 => s.seen_r1 = true,
        .r2 => s.seen_r2 = true,
        .unknown => s.seen_unknown_read = true,
    }

    if (row.total_sequences != NA_F64) {
        if (s.total_sequences == NA_F64) s.total_sequences = 0.0;
        s.total_sequences += row.total_sequences;
    }
    const weight = row.total_sequences;
    s.pct_duplicates_agg.add(row.pct_duplicates, weight);
    s.pct_gc_agg.add(row.pct_gc, weight);
    s.avg_seq_len_agg.add(row.avg_seq_len, weight);
    s.median_seq_len_agg.add(row.median_seq_len, weight);
    s.pct_fails_agg.add(row.pct_fails, weight);

    s.basic_statistics       = mergeStatus(s.basic_statistics,       row.basic_statistics);
    s.per_base_seq_quality   = mergeStatus(s.per_base_seq_quality,   row.per_base_seq_quality);
    s.per_tile_seq_quality   = mergeStatus(s.per_tile_seq_quality,   row.per_tile_seq_quality);
    s.per_seq_quality_scores = mergeStatus(s.per_seq_quality_scores, row.per_seq_quality_scores);
    s.per_base_seq_content   = mergeStatus(s.per_base_seq_content,   row.per_base_seq_content);
    s.per_seq_gc_content     = mergeStatus(s.per_seq_gc_content,     row.per_seq_gc_content);
    s.per_base_n_content     = mergeStatus(s.per_base_n_content,     row.per_base_n_content);
    s.seq_len_distribution   = mergeStatus(s.seq_len_distribution,   row.seq_len_distribution);
    s.seq_duplication_levels = mergeStatus(s.seq_duplication_levels, row.seq_duplication_levels);
    s.overrepresented_seqs   = mergeStatus(s.overrepresented_seqs,   row.overrepresented_seqs);
    s.adapter_content        = mergeStatus(s.adapter_content,        row.adapter_content);
}

fn finalizeSampleAggregation(s: *SampleRecord) void {
    s.pct_duplicates = s.pct_duplicates_agg.value();
    s.pct_gc         = s.pct_gc_agg.value();
    s.avg_seq_len    = s.avg_seq_len_agg.value();
    s.median_seq_len = s.median_seq_len_agg.value();
    s.pct_fails      = s.pct_fails_agg.value();
}

fn validateSampleRecordUniqueness(samples: []const SampleRecord, sample_map: *const std.StringHashMap(usize)) void {
    for (samples, 0..) |s, i| {
        std.debug.assert(sample_map.get(s.id).? == i);
    }
}

pub fn aggregateMultiqcRows(
    alloc: Allocator,
    file_rows: []const MultiqcFileRow,
    samples: *AL(SampleRecord),
    sample_map: *std.StringHashMap(usize),
) !void {
    // Aggregation rules:
    // - total_sequences is a count and is summed across all rows for the biological sample.
    // - percent/average fields use read-count weighted averages when total_sequences is available.
    //   If no usable read count exists for a metric, fall back to an arithmetic mean of available values.
    // - FastQC statuses use the worst observed status across rows: fail > warn > pass > missing.
    for (file_rows) |row| {
        var idx: usize = undefined;
        if (sample_map.get(row.base_id)) |existing| {
            idx = existing;
        } else {
            idx = samples.items.len;
            try samples.append(alloc, SampleRecord{ .id = try alloc.dupe(u8, row.base_id) });
            try sample_map.put(samples.items[idx].id, idx);
        }
        applyMultiqcRowToSample(row, &samples.items[idx]);
    }
    for (samples.items) |*s| finalizeSampleAggregation(s);
    validateSampleRecordUniqueness(samples.items, sample_map);
}

pub fn collectMultiqcWithoutDada2(alloc: Allocator, samples: []const SampleRecord, out: *AL([]u8)) !void {
    for (samples) |s| {
        if (!s.has_dada2) try out.append(alloc, try alloc.dupe(u8, s.id));
    }
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
        \\  --manifest <path>                QIIME2 manifest CSV for filename->sample-id mapping
        \\                                   (columns: sample-id, absolute-filepath, direction)
        \\                                   Maps MultiQC filenames to DADA2 sample IDs by
        \\                                   stripping the -R<N> replicate suffix from manifest
        \\                                   sample-ids (e.g. Sample_1-R1 -> Sample_1)
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

pub fn parseArgs(io: Io, args: []const [:0]const u8) !Config {
    var cfg = Config{};
    var i: usize = 1;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelp(io);
            std.process.exit(0);
        }
        const flags = [_][]const u8{
            "--multiqc-general-stats", "--multiqc-fastqc",  "--dada2-stats",       "--manifest",
            "--out-prefix",            "--pass-strong",     "--pass-acceptable",   "--low-depth",
            "--low-retention",         "--poor-merging",    "--high-chimera",      "--very-low-filtered",
            "--high-dup-pct",          "--gc-outlier-sd",
        };
        var matched = false;
        for (flags) |flag| {
            if (std.mem.eql(u8, arg, flag)) {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("ERROR: {s} requires a value\n", .{flag});
                    std.process.exit(1);
                }
                const v = args[i];
                if (std.mem.eql(u8, flag, "--multiqc-general-stats")) {
                    cfg.multiqc_general_stats = v;
                } else if (std.mem.eql(u8, flag, "--multiqc-fastqc")) {
                    cfg.multiqc_fastqc = v;
                } else if (std.mem.eql(u8, flag, "--dada2-stats")) {
                    cfg.dada2_stats = v;
                } else if (std.mem.eql(u8, flag, "--manifest")) {
                    cfg.manifest = v;
                } else if (std.mem.eql(u8, flag, "--out-prefix")) {
                    cfg.out_prefix = v;
                } else if (std.mem.eql(u8, flag, "--pass-strong")) {
                    cfg.pass_strong = std.fmt.parseInt(u64, v, 10) catch {
                        std.debug.print("ERROR: bad --pass-strong\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--pass-acceptable")) {
                    cfg.pass_acceptable = std.fmt.parseInt(u64, v, 10) catch {
                        std.debug.print("ERROR: bad --pass-acceptable\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--low-depth")) {
                    cfg.low_depth = std.fmt.parseInt(u64, v, 10) catch {
                        std.debug.print("ERROR: bad --low-depth\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--low-retention")) {
                    cfg.low_retention = std.fmt.parseFloat(f64, v) catch {
                        std.debug.print("ERROR: bad --low-retention\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--poor-merging")) {
                    cfg.poor_merging = std.fmt.parseFloat(f64, v) catch {
                        std.debug.print("ERROR: bad --poor-merging\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--high-chimera")) {
                    cfg.high_chimera = std.fmt.parseFloat(f64, v) catch {
                        std.debug.print("ERROR: bad --high-chimera\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--very-low-filtered")) {
                    cfg.very_low_filtered = std.fmt.parseFloat(f64, v) catch {
                        std.debug.print("ERROR: bad --very-low-filtered\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--high-dup-pct")) {
                    cfg.high_dup_pct = std.fmt.parseFloat(f64, v) catch {
                        std.debug.print("ERROR: bad --high-dup-pct\n", .{});
                        std.process.exit(1);
                    };
                } else if (std.mem.eql(u8, flag, "--gc-outlier-sd")) {
                    cfg.gc_outlier_sd = std.fmt.parseFloat(f64, v) catch {
                        std.debug.print("ERROR: bad --gc-outlier-sd\n", .{});
                        std.process.exit(1);
                    };
                }
                matched = true;
                break;
            }
        }
        if (!matched) {
            std.debug.print("ERROR: unknown argument: {s}\n", .{arg});
            std.process.exit(1);
        }
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

// ============================================================
// Tests
// ============================================================

test "normalize common MultiQC and FASTQ sample IDs" {
    const r1 = normalizeSampleId("sample_R1_001.fastq.gz");
    try std.testing.expectEqualStrings("sample", r1.base);
    try std.testing.expectEqual(ReadSide.r1, r1.read_side);

    const r2 = normalizeSampleId("sample_L001_R2_001_fastqc.zip");
    try std.testing.expectEqualStrings("sample", r2.base);
    try std.testing.expectEqual(ReadSide.r2, r2.read_side);

    const underscore_pair = normalizeSampleId("sample_1.fq.gz");
    try std.testing.expectEqualStrings("sample", underscore_pair.base);
    try std.testing.expectEqual(ReadSide.r1, underscore_pair.read_side);

    const single = normalizeSampleId("sample.fastq.gz");
    try std.testing.expectEqualStrings("sample", single.base);
    try std.testing.expectEqual(ReadSide.unknown, single.read_side);
}

test "R1 and R2 normalize to the same base sample ID" {
    const norm_r1 = normalizeSampleId("sampleA_R1.fastq.gz");
    const norm_r2 = normalizeSampleId("sampleA_R2.fastq.gz");
    try std.testing.expectEqualStrings(norm_r1.base, norm_r2.base);
    try std.testing.expectEqual(ReadSide.r1, norm_r1.read_side);
    try std.testing.expectEqual(ReadSide.r2, norm_r2.read_side);
}

test "stripManifestReplicateSuffix" {
    try std.testing.expectEqualStrings("Sample_1",  stripManifestReplicateSuffix("Sample_1-R1"));
    try std.testing.expectEqualStrings("Sample_10", stripManifestReplicateSuffix("Sample_10-R2"));
    try std.testing.expectEqualStrings("Sample_1",  stripManifestReplicateSuffix("Sample_1"));
}

test "aggregate paired MultiQC rows into one biological sample" {
    const alloc = std.testing.allocator;
    var rows: AL(MultiqcFileRow) = .empty;
    defer rows.deinit(alloc);

    try rows.append(alloc, MultiqcFileRow{
        .id = "sample_R1",
        .base_id = "sample",
        .read_side = .r1,
        .pct_duplicates = 10.0,
        .pct_gc = 50.0,
        .total_sequences = 100.0,
        .per_base_seq_quality = 0,
    });
    try rows.append(alloc, MultiqcFileRow{
        .id = "sample_R2",
        .base_id = "sample",
        .read_side = .r2,
        .pct_duplicates = 30.0,
        .pct_gc = 60.0,
        .total_sequences = 300.0,
        .per_base_seq_quality = 2,
    });

    var samples: AL(SampleRecord) = .empty;
    defer {
        for (samples.items) |s| alloc.free(s.id);
        samples.deinit(alloc);
    }
    var sample_map = std.StringHashMap(usize).init(alloc);
    defer sample_map.deinit();

    try aggregateMultiqcRows(alloc, rows.items, &samples, &sample_map);
    try std.testing.expectEqual(@as(usize, 1), samples.items.len);

    const s = samples.items[0];
    try std.testing.expectEqualStrings("sample", s.id);
    try std.testing.expectEqual(@as(usize, 2), s.file_count);
    try std.testing.expect(s.seen_r1);
    try std.testing.expect(s.seen_r2);
    try std.testing.expectApproxEqAbs(@as(f64, 400.0), s.total_sequences, 0.0001);
    try std.testing.expectApproxEqAbs(@as(f64, 25.0), s.pct_duplicates, 0.0001);
    try std.testing.expectApproxEqAbs(@as(f64, 57.5), s.pct_gc, 0.0001);
    try std.testing.expectEqual(@as(u8, 2), s.per_base_seq_quality);
}

test "_R1 and _R2 raw IDs produce one merged biological sample" {
    const alloc = std.testing.allocator;

    const raw_ids = [_]struct { id: []const u8, side: ReadSide }{
        .{ .id = "sampleX_R1.fastq.gz", .side = .r1 },
        .{ .id = "sampleX_R2.fastq.gz", .side = .r2 },
    };

    var rows: AL(MultiqcFileRow) = .empty;
    defer rows.deinit(alloc);
    for (raw_ids) |entry| {
        const norm = normalizeSampleId(entry.id);
        try rows.append(alloc, MultiqcFileRow{
            .id = entry.id,
            .base_id = norm.base,
            .read_side = norm.read_side,
            .total_sequences = 1000.0,
        });
    }

    var samples: AL(SampleRecord) = .empty;
    defer {
        for (samples.items) |s| alloc.free(s.id);
        samples.deinit(alloc);
    }
    var sample_map = std.StringHashMap(usize).init(alloc);
    defer sample_map.deinit();

    try aggregateMultiqcRows(alloc, rows.items, &samples, &sample_map);
    try std.testing.expectEqual(@as(usize, 1), samples.items.len);
    try std.testing.expectEqualStrings("sampleX", samples.items[0].id);
    try std.testing.expect(samples.items[0].seen_r1);
    try std.testing.expect(samples.items[0].seen_r2);
    try std.testing.expectApproxEqAbs(@as(f64, 2000.0), samples.items[0].total_sequences, 0.001);
}

test "applyDada2 warns on malformed numeric field and continues processing" {
    const alloc = std.testing.allocator;
    var col = std.StringHashMap(usize).init(alloc);
    defer col.deinit();
    try col.put("input", 0);
    try col.put("filtered", 1);
    try col.put("denoised", 2);
    try col.put("merged", 3);
    try col.put("non-chimeric", 4);

    // Build a mutable row ([]u8 elements) so applyDada2's type matches.
    const strs = [_][]const u8{ "not_a_number", "100", "90", "80", "70" };
    var row: [5][]u8 = undefined;
    for (strs, 0..) |s, i| row[i] = try alloc.dupe(u8, s);
    defer for (row) |f| alloc.free(f);

    var s = SampleRecord{ .id = "sample_bad" };

    // Emits a warning to stderr for "input"; all other fields parse correctly.
    applyDada2(row[0..], &col, &s, "sample_bad");

    try std.testing.expect(s.has_dada2);
    try std.testing.expectEqual(NA_U64, s.dada2_input);
    try std.testing.expectEqual(@as(u64, 100), s.dada2_filtered);
    try std.testing.expectEqual(@as(u64, 90), s.dada2_denoised);
    try std.testing.expectEqual(@as(u64, 80), s.dada2_merged);
    try std.testing.expectEqual(@as(u64, 70), s.dada2_non_chimeric);
}
