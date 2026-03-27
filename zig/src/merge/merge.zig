// Merge engine for bcftools-zig.
//
// Ported from vcfmerge.c. Implements:
//   - Multi-file reading: open N VCF files, iterate in coordinate order
//   - Position matching: find records at the same CHROM:POS across files
//   - Allele merging: create unified REF/ALT allele list
//   - Genotype remapping: remap GT allele indices to unified allele list
//   - Simple INFO merging: take first non-missing value
//
// Performance notes (audit items):
//   #9:  No quadratic copy_string_field; allele building uses direct slices
//   #16-19: Arena allocator per-position avoids per-record heap allocs

const std = @import("std");
const VcfReader = @import("../vcf/reader.zig").VcfReader;
const VcfRecord = @import("../vcf/record.zig").VcfRecord;

/// Merge mode controlling which variant types are collapsed together.
pub const MergeMode = enum {
    none,
    snps,
    indels,
    both,
    all,
};

/// Per-file state: the current record (if any) and reader reference.
const FileState = struct {
    reader: VcfReader,
    rec: VcfRecord,
    has_record: bool,
    sample_offset: u32, // index of this file's first sample in merged output
    n_samples: u32, // number of samples in this file
};

/// Allele mapping entry: maps a file-local allele index to the merged index.
const AlleleMap = struct {
    merged_index: u32,
};

/// Core merge context. Holds N readers and produces merged VCF output.
///
/// Usage:
///   var ctx = try MergeContext.init(allocator, file_paths, .{});
///   defer ctx.deinit();
///   while (try ctx.next()) |merged_line| { ... }
pub const MergeContext = struct {
    allocator: std.mem.Allocator,

    // Per-file state
    files: []FileState,
    n_files: u32,

    // Merged sample list
    all_samples: std.ArrayListUnmanaged([]const u8),
    total_samples: u32,

    // Merged header lines (union of all input headers)
    merged_header_lines: std.ArrayListUnmanaged([]const u8),

    // Options
    merge_mode: MergeMode,
    force_samples: bool,

    // Per-position arena (audit #16-19): reset after each merged record
    arena: std.heap.ArenaAllocator,

    // Reusable scratch buffers
    // Merged alleles for current position
    merged_alleles: std.ArrayListUnmanaged([]const u8),
    // Which files have records at current position
    active_files: std.ArrayListUnmanaged(u32),
    // Per-file allele maps: file_allele_maps[file_idx] maps local allele -> merged allele
    file_allele_maps: std.ArrayListUnmanaged([]AlleleMap),

    // Output buffer
    output_buf: std.ArrayListUnmanaged(u8),

    // Statistics
    n_records_merged: u64,

    pub const Options = struct {
        merge_mode: MergeMode = .both,
        force_samples: bool = false,
    };

    pub fn init(allocator: std.mem.Allocator, file_paths: []const []const u8, opts: Options) !MergeContext {
        if (file_paths.len < 2) return error.TooFewFiles;
        if (file_paths.len > 1024) return error.TooManyFiles;

        const n: u32 = @intCast(file_paths.len);
        const files = try allocator.alloc(FileState, n);
        errdefer allocator.free(files);

        var total_samples: u32 = 0;
        var all_samples: std.ArrayListUnmanaged([]const u8) = .empty;
        errdefer all_samples.deinit(allocator);
        var merged_header_lines: std.ArrayListUnmanaged([]const u8) = .empty;
        errdefer merged_header_lines.deinit(allocator);

        // Track seen header lines to avoid duplicates in merged header
        var seen_headers = std.StringHashMap(void).init(allocator);
        defer seen_headers.deinit();

        // Track sample names for duplicate detection
        var sample_set = std.StringHashMap(u32).init(allocator);
        defer sample_set.deinit();

        var init_count: u32 = 0;
        errdefer {
            for (files[0..init_count]) |*f| {
                f.rec.deinit();
                f.reader.deinit();
            }
        }

        for (file_paths, 0..) |path, idx| {
            var reader = try VcfReader.open(allocator, path);
            errdefer reader.deinit();

            var rec = VcfRecord.init(allocator);
            errdefer rec.deinit();

            // Check for duplicate sample names
            const ns: u32 = @intCast(reader.nSamples());
            for (reader.sample_names.items) |sname| {
                if (sample_set.contains(sname)) {
                    if (!opts.force_samples) return error.DuplicateSampleName;
                }
                try sample_set.put(sname, @intCast(all_samples.items.len));
                try all_samples.append(allocator, sname);
            }

            // Merge header lines (union)
            for (reader.header_lines.items) |hline| {
                if (!seen_headers.contains(hline)) {
                    try seen_headers.put(hline, {});
                    try merged_header_lines.append(allocator, hline);
                }
            }

            // Read first record
            const has = reader.next(&rec) catch false;

            files[idx] = .{
                .reader = reader,
                .rec = rec,
                .has_record = has,
                .sample_offset = total_samples,
                .n_samples = ns,
            };
            init_count += 1;
            total_samples += ns;
        }

        return MergeContext{
            .allocator = allocator,
            .files = files,
            .n_files = n,
            .all_samples = all_samples,
            .total_samples = total_samples,
            .merged_header_lines = merged_header_lines,
            .merge_mode = opts.merge_mode,
            .force_samples = opts.force_samples,
            .arena = std.heap.ArenaAllocator.init(allocator),
            .merged_alleles = .empty,
            .active_files = .empty,
            .file_allele_maps = .empty,
            .output_buf = .empty,
            .n_records_merged = 0,
        };
    }

    pub fn deinit(self: *MergeContext) void {
        for (self.files) |*f| {
            f.rec.deinit();
            f.reader.deinit();
        }
        self.allocator.free(self.files);
        self.all_samples.deinit(self.allocator);
        self.merged_header_lines.deinit(self.allocator);
        self.arena.deinit();
        self.merged_alleles.deinit(self.allocator);
        self.active_files.deinit(self.allocator);
        self.file_allele_maps.deinit(self.allocator);
        self.output_buf.deinit(self.allocator);
        self.* = undefined;
    }

    /// Write the merged VCF header to the given writer.
    pub fn writeHeader(self: *MergeContext, writer: anytype) !void {
        // Write all meta-information lines (##...) except #CHROM
        for (self.merged_header_lines.items) |hline| {
            if (hline.len >= 2 and hline[0] == '#' and hline[1] == '#') {
                try writer.writeAll(hline);
                try writer.writeAll("\n");
            }
        }

        // Write #CHROM line with merged samples
        try writer.writeAll("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if (self.total_samples > 0) {
            try writer.writeAll("\tFORMAT");
            for (self.all_samples.items) |sname| {
                try writer.writeAll("\t");
                try writer.writeAll(sname);
            }
        }
        try writer.writeAll("\n");
    }

    /// Get the next merged record as a VCF text line.
    /// Returns null at EOF (all files exhausted).
    pub fn next(self: *MergeContext) !?[]const u8 {
        // Reset per-position arena (audit #16-19)
        _ = self.arena.reset(.retain_capacity);

        // 1. Find minimum position across all files
        var min_chrom: ?[]const u8 = null;
        var min_rid: i32 = std.math.maxInt(i32);
        var min_pos: u32 = std.math.maxInt(u32);

        for (self.files) |*f| {
            if (!f.has_record) continue;
            const cmp = chromPosCompare(f.rec.rid, f.rec.pos, min_rid, min_pos);
            if (min_chrom == null or cmp == .lt) {
                min_chrom = f.rec.chrom;
                min_rid = f.rec.rid;
                min_pos = f.rec.pos;
            }
        }

        if (min_chrom == null) return null; // EOF

        // 2. Collect all files with records at this position
        self.active_files.clearRetainingCapacity();
        for (self.files, 0..) |*f, idx| {
            if (!f.has_record) continue;
            if (f.rec.rid == min_rid and f.rec.pos == min_pos) {
                try self.active_files.append(self.allocator, @intCast(idx));
            }
        }

        // 3. Merge alleles across active files
        try self.mergeAlleles();

        // 4. Build the output VCF line
        try self.buildOutputLine(min_chrom.?, min_pos);

        // 5. Advance active files to their next record
        for (self.active_files.items) |file_idx| {
            var f = &self.files[file_idx];
            f.has_record = f.reader.next(&f.rec) catch false;
        }

        self.n_records_merged += 1;
        return self.output_buf.items;
    }

    /// Merge alleles from all active files into a unified REF + ALT list.
    /// Also builds per-file allele mappings (local -> merged index).
    fn mergeAlleles(self: *MergeContext) !void {
        const arena_alloc = self.arena.allocator();
        self.merged_alleles.clearRetainingCapacity();
        self.file_allele_maps.clearRetainingCapacity();

        // Ensure we have space for allele maps for all files
        try self.file_allele_maps.resize(self.allocator, self.n_files);
        for (self.file_allele_maps.items) |*m| {
            m.* = &.{};
        }

        // Use a hash map to deduplicate alleles (direct slice comparison, audit #9)
        var allele_index_map = std.StringHashMap(u32).init(arena_alloc);

        for (self.active_files.items) |file_idx| {
            const f = &self.files[file_idx];
            const n_alleles = f.rec.nAllele();

            // Allocate allele map for this file using arena (audit #16-19)
            const amap = try arena_alloc.alloc(AlleleMap, n_alleles);

            for (0..n_alleles) |ai| {
                const allele_str = f.rec.allele(ai);

                if (ai == 0 and self.merged_alleles.items.len == 0) {
                    // First REF allele establishes the reference
                    try self.merged_alleles.append(self.allocator, allele_str);
                    try allele_index_map.put(allele_str, 0);
                    amap[0] = .{ .merged_index = 0 };
                    continue;
                }

                if (ai == 0) {
                    // REF from another file: should match (we map it to index 0)
                    amap[0] = .{ .merged_index = 0 };
                    continue;
                }

                // ALT allele: check if already in merged set
                if (allele_index_map.get(allele_str)) |existing_idx| {
                    amap[ai] = .{ .merged_index = existing_idx };
                } else {
                    const new_idx: u32 = @intCast(self.merged_alleles.items.len);
                    try self.merged_alleles.append(self.allocator, allele_str);
                    try allele_index_map.put(allele_str, new_idx);
                    amap[ai] = .{ .merged_index = new_idx };
                }
            }

            self.file_allele_maps.items[file_idx] = amap;
        }
    }

    /// Build the merged VCF output line for the current position.
    fn buildOutputLine(self: *MergeContext, chrom: []const u8, pos: u32) !void {
        self.output_buf.clearRetainingCapacity();

        // CHROM
        try self.output_buf.appendSlice(self.allocator, chrom);
        try self.output_buf.append(self.allocator, '\t');

        // POS (1-based)
        var pos_buf: [16]u8 = undefined;
        const pos_str = std.fmt.bufPrint(&pos_buf, "{d}", .{pos + 1}) catch unreachable;
        try self.output_buf.appendSlice(self.allocator, pos_str);
        try self.output_buf.append(self.allocator, '\t');

        // ID: take first non-missing ID from active files
        const id = self.pickFirstNonMissing(.id);
        try self.output_buf.appendSlice(self.allocator, id);
        try self.output_buf.append(self.allocator, '\t');

        // REF
        if (self.merged_alleles.items.len > 0) {
            try self.output_buf.appendSlice(self.allocator, self.merged_alleles.items[0]);
        } else {
            try self.output_buf.append(self.allocator, '.');
        }
        try self.output_buf.append(self.allocator, '\t');

        // ALT
        if (self.merged_alleles.items.len > 1) {
            for (self.merged_alleles.items[1..], 0..) |alt, i| {
                if (i > 0) try self.output_buf.append(self.allocator, ',');
                try self.output_buf.appendSlice(self.allocator, alt);
            }
        } else {
            try self.output_buf.append(self.allocator, '.');
        }
        try self.output_buf.append(self.allocator, '\t');

        // QUAL: take best (lowest non-missing) from active files
        const qual = self.pickBestQual();
        try self.output_buf.appendSlice(self.allocator, qual);
        try self.output_buf.append(self.allocator, '\t');

        // FILTER: take PASS if any file has PASS, else first non-missing
        const filter_str = self.pickFilter();
        try self.output_buf.appendSlice(self.allocator, filter_str);
        try self.output_buf.append(self.allocator, '\t');

        // INFO: take first non-missing from active files (simple merge)
        const info = self.pickFirstNonMissing(.info);
        try self.output_buf.appendSlice(self.allocator, info);

        // FORMAT + sample columns (only if there are samples)
        if (self.total_samples > 0) {
            try self.output_buf.append(self.allocator, '\t');
            try self.output_buf.appendSlice(self.allocator, "GT");

            // Build GT for each sample
            try self.buildMergedGenotypes();
        }
    }

    /// Field selector for pickFirstNonMissing.
    const PickField = enum { id, info };

    /// Pick the first non-missing value for a field from active files.
    fn pickFirstNonMissing(self: *MergeContext, field: PickField) []const u8 {
        for (self.active_files.items) |file_idx| {
            const f = &self.files[file_idx];
            const val = switch (field) {
                .id => f.rec.id,
                .info => self.getInfoField(f),
            };
            if (val.len > 0 and !std.mem.eql(u8, val, ".")) {
                return val;
            }
        }
        return ".";
    }

    /// Extract the INFO field (column 7) from a file's current record storage.
    fn getInfoField(self: *MergeContext, f: *const FileState) []const u8 {
        _ = self;
        const storage = f.rec._storage orelse return ".";
        // Find the 7th tab-separated field (0-indexed column 7)
        var col: u32 = 0;
        var start: usize = 0;
        for (storage, 0..) |c, i| {
            if (c == '\t') {
                if (col == 7) return storage[start..i];
                col += 1;
                start = i + 1;
            }
        }
        if (col == 7) return storage[start..];
        return ".";
    }

    /// Pick the best QUAL value (lowest non-missing).
    fn pickBestQual(self: *MergeContext) []const u8 {
        for (self.active_files.items) |file_idx| {
            const f = &self.files[file_idx];
            if (f.rec.qual) |_| {
                // Return the raw text from storage
                const storage = f.rec._storage orelse continue;
                var col: u32 = 0;
                var start: usize = 0;
                for (storage, 0..) |c, i| {
                    if (c == '\t') {
                        if (col == 5) return storage[start..i];
                        col += 1;
                        start = i + 1;
                    }
                }
                if (col == 5) return storage[start..];
            }
        }
        return ".";
    }

    /// Pick filter string: PASS if any active file has PASS, else first non-missing.
    fn pickFilter(self: *MergeContext) []const u8 {
        var first_non_missing: ?[]const u8 = null;
        for (self.active_files.items) |file_idx| {
            const f = &self.files[file_idx];
            if (f.rec.filter.len > 0 and !std.mem.eql(u8, f.rec.filter, ".")) {
                if (std.mem.eql(u8, f.rec.filter, "PASS")) return "PASS";
                if (first_non_missing == null) first_non_missing = f.rec.filter;
            }
        }
        return first_non_missing orelse ".";
    }

    /// Build merged genotype columns for all samples.
    /// For samples in active files: parse and remap GT.
    /// For samples in inactive files: output ./.
    fn buildMergedGenotypes(self: *MergeContext) !void {
        // For each sample (in order), determine which file owns it and
        // whether that file is active at this position.
        var sample_idx: u32 = 0;
        for (self.files, 0..) |*f, file_idx| {
            const is_active = self.isFileActive(@intCast(file_idx));

            for (0..f.n_samples) |local_smpl| {
                try self.output_buf.append(self.allocator, '\t');

                if (!is_active) {
                    // File has no record at this position -> missing GT
                    try self.output_buf.appendSlice(self.allocator, "./.");
                } else {
                    // Parse GT from this file's record and remap alleles
                    try self.remapAndWriteGt(
                        f,
                        @intCast(file_idx),
                        @intCast(local_smpl),
                    );
                }

                sample_idx += 1;
            }
        }
    }

    /// Check if a file index is in the active set.
    fn isFileActive(self: *MergeContext, file_idx: u32) bool {
        for (self.active_files.items) |aidx| {
            if (aidx == file_idx) return true;
        }
        return false;
    }

    /// Parse the GT field from a file's record for a given local sample,
    /// remap allele indices using the allele map, and write to output_buf.
    fn remapAndWriteGt(
        self: *MergeContext,
        f: *const FileState,
        file_idx: u32,
        local_sample: u32,
    ) !void {
        const storage = f.rec._storage orelse {
            try self.output_buf.appendSlice(self.allocator, "./.");
            return;
        };

        // Find the sample column (9 + local_sample)
        const target_col: u32 = 9 + local_sample;
        const gt_text = getColumn(storage, target_col) orelse {
            try self.output_buf.appendSlice(self.allocator, "./.");
            return;
        };

        // The GT field is the first colon-separated subfield in the sample column
        const gt_end = std.mem.indexOfScalar(u8, gt_text, ':') orelse gt_text.len;
        const gt = gt_text[0..gt_end];

        if (gt.len == 0) {
            try self.output_buf.appendSlice(self.allocator, "./.");
            return;
        }

        // Parse and remap the GT: alleles separated by / or |
        const amap = self.file_allele_maps.items[file_idx];
        var i: usize = 0;
        var first = true;
        while (i < gt.len) {
            // Find separator
            var end = i;
            while (end < gt.len and gt[end] != '/' and gt[end] != '|') : (end += 1) {}

            if (!first) {
                // Write the separator that was before this allele
                if (i > 0) {
                    try self.output_buf.append(self.allocator, gt[i - 1]);
                }
            }
            first = false;

            const allele_text = gt[i..end];

            if (std.mem.eql(u8, allele_text, ".")) {
                try self.output_buf.append(self.allocator, '.');
            } else {
                const local_idx = std.fmt.parseInt(u32, allele_text, 10) catch {
                    try self.output_buf.append(self.allocator, '.');
                    i = if (end < gt.len) end + 1 else gt.len;
                    continue;
                };

                if (local_idx < amap.len) {
                    const merged_idx = amap[local_idx].merged_index;
                    var idx_buf: [16]u8 = undefined;
                    const idx_str = std.fmt.bufPrint(&idx_buf, "{d}", .{merged_idx}) catch unreachable;
                    try self.output_buf.appendSlice(self.allocator, idx_str);
                } else {
                    try self.output_buf.append(self.allocator, '.');
                }
            }

            i = if (end < gt.len) end + 1 else gt.len;
        }

        // If we wrote nothing (empty GT), write ./.
        if (first) {
            try self.output_buf.appendSlice(self.allocator, "./.");
        }
    }

    /// Compare two (rid, pos) pairs for coordinate ordering.
    fn chromPosCompare(rid_a: i32, pos_a: u32, rid_b: i32, pos_b: u32) std.math.Order {
        if (rid_a < rid_b) return .lt;
        if (rid_a > rid_b) return .gt;
        if (pos_a < pos_b) return .lt;
        if (pos_a > pos_b) return .gt;
        return .eq;
    }
};

/// Get a tab-separated column from a line by 0-based column index.
fn getColumn(line: []const u8, col: u32) ?[]const u8 {
    var current_col: u32 = 0;
    var start: usize = 0;
    for (line, 0..) |c, i| {
        if (c == '\t') {
            if (current_col == col) return line[start..i];
            current_col += 1;
            start = i + 1;
        }
    }
    if (current_col == col) return line[start..];
    return null;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Helper to create a temp VCF file for testing.
fn writeTempVcf(tmp_dir: std.testing.TmpDir, name: []const u8, content: []const u8) ![]const u8 {
    const f = try tmp_dir.dir.createFile(name, .{});
    try f.writeAll(content);
    f.close();
    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try tmp_dir.dir.realpath(name, &path_buf);
    // Caller needs a stable copy since path_buf is on stack
    return path;
}

test "merge 2 files with non-overlapping samples at same position" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const vcf1 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSampleA\n" ++
        "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n";

    const vcf2 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSampleB\n" ++
        "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\n";

    const f1 = try tmp_dir.dir.createFile("f1.vcf", .{});
    try f1.writeAll(vcf1);
    f1.close();
    const f2 = try tmp_dir.dir.createFile("f2.vcf", .{});
    try f2.writeAll(vcf2);
    f2.close();

    var path_buf1: [std.fs.max_path_bytes]u8 = undefined;
    var path_buf2: [std.fs.max_path_bytes]u8 = undefined;
    const p1 = try tmp_dir.dir.realpath("f1.vcf", &path_buf1);
    const p2 = try tmp_dir.dir.realpath("f2.vcf", &path_buf2);

    // We need stable paths since the slices point into stack buffers
    const path1 = try allocator.dupe(u8, p1);
    defer allocator.free(path1);
    const path2 = try allocator.dupe(u8, p2);
    defer allocator.free(path2);

    const paths: []const []const u8 = &.{ path1, path2 };

    var ctx = try MergeContext.init(allocator, paths, .{});
    defer ctx.deinit();

    // Should have 2 samples total
    try std.testing.expectEqual(@as(u32, 2), ctx.total_samples);

    const line = try ctx.next();
    try std.testing.expect(line != null);

    // Check that merged line contains both samples
    const l = line.?;
    // Should contain "chr1\t100\t"
    try std.testing.expect(std.mem.indexOf(u8, l, "chr1\t100\t") != null);
    // Should contain both GT values (0/1 for SampleA, 0/0 for SampleB)
    // The line format: chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0

    // Verify no more records
    const line2 = try ctx.next();
    try std.testing.expect(line2 == null);
}

test "allele merging: A>T + A>G produces A>T,G" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const vcf1 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" ++
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT\t0/1\n";

    const vcf2 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n" ++
        "chr1\t200\t.\tA\tG\t.\t.\t.\tGT\t0/1\n";

    const f1 = try tmp_dir.dir.createFile("a1.vcf", .{});
    try f1.writeAll(vcf1);
    f1.close();
    const f2 = try tmp_dir.dir.createFile("a2.vcf", .{});
    try f2.writeAll(vcf2);
    f2.close();

    var pb1: [std.fs.max_path_bytes]u8 = undefined;
    var pb2: [std.fs.max_path_bytes]u8 = undefined;
    const rp1 = try tmp_dir.dir.realpath("a1.vcf", &pb1);
    const rp2 = try tmp_dir.dir.realpath("a2.vcf", &pb2);

    const p1 = try allocator.dupe(u8, rp1);
    defer allocator.free(p1);
    const p2 = try allocator.dupe(u8, rp2);
    defer allocator.free(p2);

    const paths: []const []const u8 = &.{ p1, p2 };

    var ctx = try MergeContext.init(allocator, paths, .{});
    defer ctx.deinit();

    const line = try ctx.next();
    try std.testing.expect(line != null);
    const l = line.?;

    // REF=A, ALT=T,G
    try std.testing.expect(std.mem.indexOf(u8, l, "\tA\tT,G\t") != null);
}

test "GT remapping: 0/1 stays 0/1 for first file, becomes 0/2 for second" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const vcf1 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" ++
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT\t0/1\n";

    const vcf2 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n" ++
        "chr1\t200\t.\tA\tG\t.\t.\t.\tGT\t0/1\n";

    const f1 = try tmp_dir.dir.createFile("g1.vcf", .{});
    try f1.writeAll(vcf1);
    f1.close();
    const f2 = try tmp_dir.dir.createFile("g2.vcf", .{});
    try f2.writeAll(vcf2);
    f2.close();

    var pb1: [std.fs.max_path_bytes]u8 = undefined;
    var pb2: [std.fs.max_path_bytes]u8 = undefined;
    const rp1 = try tmp_dir.dir.realpath("g1.vcf", &pb1);
    const rp2 = try tmp_dir.dir.realpath("g2.vcf", &pb2);

    const path1 = try allocator.dupe(u8, rp1);
    defer allocator.free(path1);
    const path2 = try allocator.dupe(u8, rp2);
    defer allocator.free(path2);

    const paths: []const []const u8 = &.{ path1, path2 };

    var ctx = try MergeContext.init(allocator, paths, .{});
    defer ctx.deinit();

    const line = try ctx.next();
    try std.testing.expect(line != null);
    const l = line.?;

    // S1 GT should be 0/1, S2 GT should be 0/2
    // Format: ...GT\t0/1\t0/2
    try std.testing.expect(std.mem.indexOf(u8, l, "GT\t0/1\t0/2") != null);
}

test "non-overlapping positions: records from each file pass through" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const vcf1 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" ++
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\n";

    const vcf2 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n" ++
        "chr1\t200\t.\tG\tC\t.\t.\t.\tGT\t1/1\n";

    const f1 = try tmp_dir.dir.createFile("n1.vcf", .{});
    try f1.writeAll(vcf1);
    f1.close();
    const f2 = try tmp_dir.dir.createFile("n2.vcf", .{});
    try f2.writeAll(vcf2);
    f2.close();

    var pb1: [std.fs.max_path_bytes]u8 = undefined;
    var pb2: [std.fs.max_path_bytes]u8 = undefined;
    const rp1 = try tmp_dir.dir.realpath("n1.vcf", &pb1);
    const rp2 = try tmp_dir.dir.realpath("n2.vcf", &pb2);

    const path1 = try allocator.dupe(u8, rp1);
    defer allocator.free(path1);
    const path2 = try allocator.dupe(u8, rp2);
    defer allocator.free(path2);

    const paths: []const []const u8 = &.{ path1, path2 };

    var ctx = try MergeContext.init(allocator, paths, .{});
    defer ctx.deinit();

    // First record: chr1:100 from file1 only
    const line1 = try ctx.next();
    try std.testing.expect(line1 != null);
    const l1 = line1.?;
    try std.testing.expect(std.mem.indexOf(u8, l1, "chr1\t100\t") != null);
    // S1 has 0/1, S2 has ./. (missing)
    try std.testing.expect(std.mem.indexOf(u8, l1, "0/1\t./.") != null);

    // Second record: chr1:200 from file2 only
    const line2 = try ctx.next();
    try std.testing.expect(line2 != null);
    const l2 = line2.?;
    try std.testing.expect(std.mem.indexOf(u8, l2, "chr1\t200\t") != null);
    // S1 has ./. (missing), S2 has 1/1
    try std.testing.expect(std.mem.indexOf(u8, l2, "./.\t1/1") != null);

    // No more records
    const line3 = try ctx.next();
    try std.testing.expect(line3 == null);
}

test "missing genotypes: samples from file without record get ./." {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // File 1 has 2 samples, file 2 has 1 sample
    // Only file 1 has a record at this position
    const vcf1 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n" ++
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/1\n";

    const vcf2 =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=1000>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS3\n" ++
        "chr1\t999\t.\tG\tC\t.\t.\t.\tGT\t0/1\n";

    const f1 = try tmp_dir.dir.createFile("m1.vcf", .{});
    try f1.writeAll(vcf1);
    f1.close();
    const f2 = try tmp_dir.dir.createFile("m2.vcf", .{});
    try f2.writeAll(vcf2);
    f2.close();

    var pb1: [std.fs.max_path_bytes]u8 = undefined;
    var pb2: [std.fs.max_path_bytes]u8 = undefined;
    const rp1 = try tmp_dir.dir.realpath("m1.vcf", &pb1);
    const rp2 = try tmp_dir.dir.realpath("m2.vcf", &pb2);

    const path1 = try allocator.dupe(u8, rp1);
    defer allocator.free(path1);
    const path2 = try allocator.dupe(u8, rp2);
    defer allocator.free(path2);

    const paths: []const []const u8 = &.{ path1, path2 };

    var ctx = try MergeContext.init(allocator, paths, .{});
    defer ctx.deinit();

    try std.testing.expectEqual(@as(u32, 3), ctx.total_samples);

    // First record at chr1:100 — only file 1 active
    const line1 = try ctx.next();
    try std.testing.expect(line1 != null);
    const l1 = line1.?;
    // S1=0/1, S2=1/1, S3=./.
    try std.testing.expect(std.mem.indexOf(u8, l1, "0/1\t1/1\t./.") != null);

    // Second record at chr1:999 — only file 2 active
    const line2 = try ctx.next();
    try std.testing.expect(line2 != null);
    const l2 = line2.?;
    // S1=./., S2=./., S3=0/1
    try std.testing.expect(std.mem.indexOf(u8, l2, "./.\t./.\t0/1") != null);
}

test "getColumn helper" {
    const line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1";
    try std.testing.expectEqualStrings("chr1", getColumn(line, 0).?);
    try std.testing.expectEqualStrings("100", getColumn(line, 1).?);
    try std.testing.expectEqualStrings("A", getColumn(line, 3).?);
    try std.testing.expectEqualStrings("0/1", getColumn(line, 9).?);
    try std.testing.expect(getColumn(line, 10) == null);
}
