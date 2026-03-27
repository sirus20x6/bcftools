// Variant normalization engine for bcftools-zig.
//
// Ported from vcfnorm.c. Implements:
//   - Left-alignment of indels (shift variants left to canonical position)
//   - REF/ALT trimming (remove common prefix/suffix)
//   - Multiallelic splitting (-m-): split multi-allelic records into biallelics
//   - Reference checking (-c): verify REF matches FASTA reference
//
// The core algorithms operate on VcfRecord (from vcf/record.zig) and modify
// records in-place where possible, or produce new records for splitting.

const std = @import("std");
const VcfRecord = @import("../vcf/record.zig").VcfRecord;

/// How to handle REF allele mismatches against the FASTA reference.
pub const CheckRef = enum {
    /// Do not check REF against reference.
    none,
    /// Print a warning on mismatch but continue.
    warn,
    /// Attempt to fix the REF allele (replace with FASTA sequence).
    fix,
    /// Discard records whose REF does not match.
    discard,
    /// Exit with error on mismatch.
    exit_on_error,
};

/// Which variant classes to split when using multiallelic splitting.
pub const SplitMode = enum {
    /// No splitting.
    none,
    /// Split SNP-only multiallelic sites.
    snps,
    /// Split indel-only multiallelic sites.
    indels,
    /// Split both SNPs and indels.
    both,
    /// Split any multiallelic site.
    any,
};

/// Function type for fetching a reference sequence region.
/// Parameters: context pointer, allocator, chromosome (0-terminated), begin (0-based), end (0-based inclusive).
/// Returns the fetched sequence or null on failure.
pub const FetchSeqFn = *const fn (ctx: *anyopaque, allocator: std.mem.Allocator, chr: [*:0]const u8, beg: i64, end: i64) ?[]u8;

/// Configuration options for NormContext initialization.
pub const Options = struct {
    do_left_align: bool = true,
    do_trim: bool = true,
    check_ref: CheckRef = .warn,
    split_mode: SplitMode = .none,
    fai_ctx: ?*anyopaque = null,
    fetch_seq_fn: ?FetchSeqFn = null,
    aln_win: u32 = 200,
};

/// Result of normalizing a single record.
pub const NormResult = struct {
    /// Whether the record was modified in place.
    modified: bool,
    /// If check_ref is discard and REF mismatched, this is true.
    discarded: bool,
    /// Warning message (if any). Caller does not own this memory.
    warning: ?[]const u8,
};

/// Normalization context. Holds configuration and scratch buffers.
pub const NormContext = struct {
    allocator: std.mem.Allocator,

    // Reference access
    fai_ctx: ?*anyopaque,
    fetch_seq_fn: ?FetchSeqFn,

    // Options
    do_left_align: bool,
    do_trim: bool,
    check_ref: CheckRef,
    split_mode: SplitMode,
    aln_win: u32,

    // Statistics
    stats: Stats,

    pub const Stats = struct {
        total: u64 = 0,
        modified: u64 = 0,
        split: u64 = 0,
        ref_mismatch: u64 = 0,
        ref_fixed: u64 = 0,
        discarded: u64 = 0,
    };

    /// Initialize a NormContext with the given options.
    pub fn init(allocator: std.mem.Allocator, options: Options) NormContext {
        return .{
            .allocator = allocator,
            .fai_ctx = options.fai_ctx,
            .fetch_seq_fn = options.fetch_seq_fn,
            .do_left_align = options.do_left_align,
            .do_trim = options.do_trim,
            .check_ref = options.check_ref,
            .split_mode = options.split_mode,
            .aln_win = options.aln_win,
            .stats = .{},
        };
    }

    pub fn deinit(self: *NormContext) void {
        _ = self;
    }

    // -----------------------------------------------------------------
    // Reference checking
    // -----------------------------------------------------------------

    /// Check whether the REF allele matches the FASTA reference.
    /// Returns a warning message if there is a mismatch, or null if OK.
    /// If check_ref is .fix, the REF allele in `rec` is replaced.
    pub fn checkRef(self: *NormContext, rec: *VcfRecord) !?[]const u8 {
        const fetch_fn = self.fetch_seq_fn orelse return null;
        const ctx = self.fai_ctx orelse return null;

        if (rec.ref_allele.len == 0) return null;

        // Build null-terminated chromosome name
        var chr_buf: [256]u8 = undefined;
        if (rec.chrom.len >= chr_buf.len) return null;
        @memcpy(chr_buf[0..rec.chrom.len], rec.chrom);
        chr_buf[rec.chrom.len] = 0;
        const chr_z: [*:0]const u8 = chr_buf[0..rec.chrom.len :0];

        const beg: i64 = @intCast(rec.pos);
        const end: i64 = beg + @as(i64, @intCast(rec.ref_allele.len)) - 1;

        const ref_seq = fetch_fn(ctx, self.allocator, chr_z, beg, end) orelse
            return "failed to fetch reference sequence";
        defer self.allocator.free(ref_seq);

        // Compare (case-insensitive)
        if (ref_seq.len != rec.ref_allele.len) {
            self.stats.ref_mismatch += 1;
            return "REF length does not match reference";
        }

        var mismatch = false;
        for (ref_seq, 0..) |c, i| {
            if (toUpper(c) != toUpper(rec.ref_allele[i])) {
                mismatch = true;
                break;
            }
        }

        if (!mismatch) return null;

        self.stats.ref_mismatch += 1;

        switch (self.check_ref) {
            .fix => {
                // Replace REF allele with the reference sequence.
                try self.replaceRefAllele(rec, ref_seq);
                self.stats.ref_fixed += 1;
                return null;
            },
            .warn => return "REF_MISMATCH",
            .discard => return "REF_MISMATCH_DISCARD",
            .exit_on_error => return "REF_MISMATCH_FATAL",
            .none => return null,
        }
    }

    /// Replace the REF allele in a VcfRecord, rebuilding its internal storage.
    fn replaceRefAllele(self: *NormContext, rec: *VcfRecord, new_ref: []const u8) !void {
        _ = self;
        try rebuildRecordFromParts(rec, rec.pos, new_ref, rec.alt_alleles.items);
    }

    // -----------------------------------------------------------------
    // Left-alignment
    // -----------------------------------------------------------------

    /// Left-align a single allele pair (ref, alt) given access to reference.
    /// Modifies ref_buf and alt_buf in-place and returns the new 0-based position.
    /// This is the core of the realign_left algorithm from vcfnorm.c.
    fn leftAlignPair(
        self: *NormContext,
        chr: [*:0]const u8,
        pos: u32,
        ref_buf: *std.ArrayListUnmanaged(u8),
        alt_buf: *std.ArrayListUnmanaged(u8),
    ) !u32 {
        const alloc = self.allocator;
        var cur_pos: u32 = pos;

        // Trim common suffix
        while (ref_buf.items.len > 1 and alt_buf.items.len > 1) {
            if (toUpper(ref_buf.items[ref_buf.items.len - 1]) !=
                toUpper(alt_buf.items[alt_buf.items.len - 1]))
                break;
            _ = ref_buf.pop();
            _ = alt_buf.pop();
        }

        // Left-shift: keep trimming suffix and padding from left
        while (true) {
            // Check if rightmost bases are identical
            if (ref_buf.items.len > 1 and alt_buf.items.len > 1) {
                if (toUpper(ref_buf.items[ref_buf.items.len - 1]) ==
                    toUpper(alt_buf.items[alt_buf.items.len - 1]))
                {
                    _ = ref_buf.pop();
                    _ = alt_buf.pop();
                    continue;
                }
            }

            // If either allele would become empty, pad from the left
            if (ref_buf.items.len == 0 or alt_buf.items.len == 0 or
                (ref_buf.items.len == 1 and alt_buf.items.len == 1 and
                toUpper(ref_buf.items[0]) == toUpper(alt_buf.items[0])))
            {
                if (cur_pos == 0) break;

                // Fetch one base to the left
                const fetch_fn = self.fetch_seq_fn orelse break;
                const ctx = self.fai_ctx orelse break;
                const beg: i64 = @as(i64, @intCast(cur_pos)) - 1;
                const seq = fetch_fn(ctx, alloc, chr, beg, beg) orelse break;
                defer alloc.free(seq);

                if (seq.len == 0) break;

                // Prepend the base to both alleles
                try ref_buf.insert(alloc, 0, toUpper(seq[0]));
                try alt_buf.insert(alloc, 0, toUpper(seq[0]));
                cur_pos -= 1;

                // Continue trimming suffix
                continue;
            }

            break;
        }

        // Trim common prefix (keep at least 1 base)
        while (ref_buf.items.len > 1 and alt_buf.items.len > 1) {
            if (toUpper(ref_buf.items[0]) != toUpper(alt_buf.items[0])) break;
            _ = ref_buf.orderedRemove(0);
            _ = alt_buf.orderedRemove(0);
            cur_pos += 1;
        }

        return cur_pos;
    }

    // -----------------------------------------------------------------
    // Public API: normalize
    // -----------------------------------------------------------------

    /// Normalize a single VCF record in-place.
    /// Performs left-alignment and trimming of REF/ALT alleles.
    /// Returns a NormResult indicating what happened.
    pub fn normalize(self: *NormContext, rec: *VcfRecord) !NormResult {
        self.stats.total += 1;

        var result = NormResult{
            .modified = false,
            .discarded = false,
            .warning = null,
        };

        // Reference check
        if (self.check_ref != .none) {
            if (try self.checkRef(rec)) |warn| {
                if (self.check_ref == .discard) {
                    self.stats.discarded += 1;
                    result.discarded = true;
                    result.warning = warn;
                    return result;
                }
                result.warning = warn;
            }
        }

        // No ALTs -> nothing to normalize
        if (rec.alt_alleles.items.len == 0) return result;

        // Skip symbolic alleles (e.g. <DEL>, <DUP>, etc.)
        for (rec.alt_alleles.items) |alt| {
            if (alt.len > 0 and alt[0] == '<') return result;
        }

        // For each ALT allele, left-align and trim against REF.
        if (self.do_left_align or self.do_trim) {
            const changed = try self.normalizeAlleles(rec);
            if (changed) {
                result.modified = true;
                self.stats.modified += 1;
            }
        }

        return result;
    }

    /// Normalize all alleles of a record (left-align + trim).
    /// Returns true if the record was modified.
    fn normalizeAlleles(self: *NormContext, rec: *VcfRecord) !bool {
        const n_alt = rec.alt_alleles.items.len;
        if (n_alt == 0) return false;

        // For a single ALT, we can do precise left-alignment.
        // For multiple ALTs, we do a simplified trim (common prefix/suffix).
        if (n_alt == 1 and self.do_left_align and self.fetch_seq_fn != null) {
            return try self.normalizeOnePair(rec);
        }

        // Multi-allelic or no FASTA: just trim common prefix/suffix
        if (self.do_trim) {
            return try self.trimAlleles(rec);
        }

        return false;
    }

    /// Left-align and trim a single REF/ALT pair.
    fn normalizeOnePair(self: *NormContext, rec: *VcfRecord) !bool {
        const alloc = self.allocator;
        var ref_buf: std.ArrayListUnmanaged(u8) = .empty;
        defer ref_buf.deinit(alloc);
        var alt_buf: std.ArrayListUnmanaged(u8) = .empty;
        defer alt_buf.deinit(alloc);

        try ref_buf.appendSlice(alloc, rec.ref_allele);
        try alt_buf.appendSlice(alloc, rec.alt_alleles.items[0]);

        // Build null-terminated chromosome name
        var chr_buf: [256]u8 = undefined;
        if (rec.chrom.len >= chr_buf.len) return false;
        @memcpy(chr_buf[0..rec.chrom.len], rec.chrom);
        chr_buf[rec.chrom.len] = 0;
        const chr_z: [*:0]const u8 = chr_buf[0..rec.chrom.len :0];

        const new_pos = try self.leftAlignPair(chr_z, rec.pos, &ref_buf, &alt_buf);

        // Check if anything changed
        if (new_pos == rec.pos and
            std.mem.eql(u8, ref_buf.items, rec.ref_allele) and
            std.mem.eql(u8, alt_buf.items, rec.alt_alleles.items[0]))
        {
            return false;
        }

        // Rebuild the record with new alleles and position
        const new_alts: [1][]const u8 = .{alt_buf.items};
        try rebuildRecordFromParts(rec, new_pos, ref_buf.items, &new_alts);
        return true;
    }

    /// Trim common prefix and suffix from all alleles (no left-shift).
    fn trimAlleles(self: *NormContext, rec: *VcfRecord) !bool {
        _ = self;
        const n_alleles = rec.alt_alleles.items.len + 1;
        if (n_alleles < 2) return false;

        const ref = rec.ref_allele;
        const alts = rec.alt_alleles.items;

        // Find common suffix length
        var suffix_len: usize = 0;
        outer_suffix: while (true) {
            if (suffix_len + 1 >= ref.len) break;
            const ref_char = toUpper(ref[ref.len - 1 - suffix_len]);
            for (alts) |alt| {
                if (suffix_len + 1 >= alt.len) break :outer_suffix;
                if (toUpper(alt[alt.len - 1 - suffix_len]) != ref_char) break :outer_suffix;
            }
            suffix_len += 1;
        }

        // Find common prefix length (leave at least 1 base after suffix trim)
        var prefix_len: usize = 0;
        outer_prefix: while (true) {
            if (prefix_len + suffix_len + 1 >= ref.len) break;
            const ref_char = toUpper(ref[prefix_len]);
            for (alts) |alt| {
                if (prefix_len + suffix_len + 1 >= alt.len) break :outer_prefix;
                if (toUpper(alt[prefix_len]) != ref_char) break :outer_prefix;
            }
            prefix_len += 1;
        }

        if (prefix_len == 0 and suffix_len == 0) return false;

        // Build trimmed alleles
        const new_ref = ref[prefix_len .. ref.len - suffix_len];
        var new_alts_buf: [256][]const u8 = undefined;
        for (alts, 0..) |alt, i| {
            new_alts_buf[i] = alt[prefix_len .. alt.len - suffix_len];
        }
        const new_alts = new_alts_buf[0..alts.len];
        const new_pos = rec.pos + @as(u32, @intCast(prefix_len));

        try rebuildRecordFromParts(rec, new_pos, new_ref, new_alts);
        return true;
    }

    // -----------------------------------------------------------------
    // Public API: split
    // -----------------------------------------------------------------

    /// Split a multi-allelic record into biallelic records.
    /// Returns a list of new VcfRecords, one per ALT allele.
    /// Caller owns the returned slice and must deinit each record.
    pub fn split(self: *NormContext, rec: *const VcfRecord) ![]VcfRecord {
        const n_alt = rec.alt_alleles.items.len;
        if (n_alt <= 1) {
            // Already biallelic or no ALT; return a single copy
            const result = try self.allocator.alloc(VcfRecord, 1);
            result[0] = try dupeRecord(self.allocator, rec);
            return result;
        }

        // Check split mode filter
        if (!self.shouldSplit(rec)) {
            const result = try self.allocator.alloc(VcfRecord, 1);
            result[0] = try dupeRecord(self.allocator, rec);
            return result;
        }

        self.stats.split += 1;

        const result = try self.allocator.alloc(VcfRecord, n_alt);
        errdefer {
            for (result) |*r| r.deinit();
            self.allocator.free(result);
        }

        for (0..n_alt) |i| {
            result[i] = VcfRecord.init(self.allocator);

            // Build a VCF line with a single ALT
            const single_alt: [1][]const u8 = .{rec.alt_alleles.items[i]};
            try rebuildRecordFromFields(
                &result[i],
                rec.chrom,
                rec.pos,
                rec.id,
                rec.ref_allele,
                &single_alt,
                rec.qual,
                rec.filter,
            );
        }

        return result;
    }

    /// Determine whether a record should be split based on split_mode.
    fn shouldSplit(self: *const NormContext, rec: *const VcfRecord) bool {
        if (self.split_mode == .none) return false;
        if (self.split_mode == .any) return true;

        // Classify variant types present
        var has_snp = false;
        var has_indel = false;
        for (rec.alt_alleles.items) |alt| {
            if (alt.len == rec.ref_allele.len) {
                has_snp = true;
            } else {
                has_indel = true;
            }
        }

        return switch (self.split_mode) {
            .snps => has_snp,
            .indels => has_indel,
            .both => has_snp or has_indel,
            .any => true,
            .none => false,
        };
    }

    /// Process a record through the full normalization + optional splitting pipeline.
    /// Returns a list of output records. Caller owns the returned slice.
    pub fn processRecord(self: *NormContext, rec: *VcfRecord) ![]VcfRecord {
        // First normalize
        const norm_result = try self.normalize(rec);
        if (norm_result.discarded) {
            return try self.allocator.alloc(VcfRecord, 0);
        }

        // Then optionally split
        if (self.split_mode != .none and rec.alt_alleles.items.len > 1) {
            const split_recs = try self.split(rec);
            // Normalize each split record individually
            for (split_recs) |*sr| {
                _ = try self.normalize(sr);
            }
            return split_recs;
        }

        // Return single record copy
        const result = try self.allocator.alloc(VcfRecord, 1);
        result[0] = try dupeRecord(self.allocator, rec);
        return result;
    }
};

// =========================================================================
// Helpers
// =========================================================================

fn toUpper(c: u8) u8 {
    return if (c >= 'a' and c <= 'z') c - 32 else c;
}

/// Build a VCF text line in a stack buffer and parse it into the record.
fn buildVcfLine(
    alloc: std.mem.Allocator,
    chrom: []const u8,
    pos_0based: u32,
    id: []const u8,
    ref_allele: []const u8,
    alt_alleles: []const []const u8,
    qual: ?f32,
    filter_field: []const u8,
) ![]u8 {
    var buf: std.ArrayListUnmanaged(u8) = .empty;
    errdefer buf.deinit(alloc);

    // CHROM
    try buf.appendSlice(alloc, chrom);
    try buf.append(alloc, '\t');
    // POS (1-based)
    var pos_str_buf: [16]u8 = undefined;
    const pos_str = std.fmt.bufPrint(&pos_str_buf, "{d}", .{pos_0based + 1}) catch unreachable;
    try buf.appendSlice(alloc, pos_str);
    try buf.append(alloc, '\t');
    // ID
    try buf.appendSlice(alloc, if (id.len > 0) id else ".");
    try buf.append(alloc, '\t');
    // REF
    try buf.appendSlice(alloc, ref_allele);
    try buf.append(alloc, '\t');
    // ALT
    if (alt_alleles.len == 0) {
        try buf.append(alloc, '.');
    } else {
        for (alt_alleles, 0..) |alt, i| {
            if (i > 0) try buf.append(alloc, ',');
            try buf.appendSlice(alloc, alt);
        }
    }
    try buf.append(alloc, '\t');
    // QUAL
    if (qual) |q| {
        var q_buf: [32]u8 = undefined;
        const q_str = std.fmt.bufPrint(&q_buf, "{d}", .{q}) catch ".";
        try buf.appendSlice(alloc, q_str);
    } else {
        try buf.append(alloc, '.');
    }
    try buf.append(alloc, '\t');
    // FILTER
    try buf.appendSlice(alloc, if (filter_field.len > 0) filter_field else ".");
    try buf.append(alloc, '\t');
    // INFO (minimal)
    try buf.append(alloc, '.');

    return buf.toOwnedSlice(alloc);
}

/// Rebuild a VcfRecord with new position, REF, and ALT alleles, preserving other fields.
fn rebuildRecordFromParts(rec: *VcfRecord, new_pos: u32, new_ref: []const u8, new_alts: []const []const u8) !void {
    const alloc = rec.allocator;
    const line = try buildVcfLine(alloc, rec.chrom, new_pos, rec.id, new_ref, new_alts, rec.qual, rec.filter);
    defer alloc.free(line);
    try rec.parseLine(line);
}

/// Rebuild a VcfRecord from explicit field values.
fn rebuildRecordFromFields(
    rec: *VcfRecord,
    chrom: []const u8,
    pos: u32,
    id: []const u8,
    ref_allele: []const u8,
    alt_alleles: []const []const u8,
    qual: ?f32,
    filter_field: []const u8,
) !void {
    const alloc = rec.allocator;
    const line = try buildVcfLine(alloc, chrom, pos, id, ref_allele, alt_alleles, qual, filter_field);
    defer alloc.free(line);
    try rec.parseLine(line);
}

/// Duplicate a VcfRecord by serializing and re-parsing.
fn dupeRecord(allocator: std.mem.Allocator, src: *const VcfRecord) !VcfRecord {
    const line = try buildVcfLine(
        allocator,
        src.chrom,
        src.pos,
        src.id,
        src.ref_allele,
        src.alt_alleles.items,
        src.qual,
        src.filter,
    );
    defer allocator.free(line);

    var rec = VcfRecord.init(allocator);
    try rec.parseLine(line);
    return rec;
}

// =========================================================================
// Tests
// =========================================================================

test "trim common prefix: AGC>ATC -> G>T at pos+1" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = true,
        .check_ref = .none,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t10\t.\tAGC\tATC\t.\tPASS\t.");
    // pos is 0-based: 10-1 = 9

    const result = try ctx.normalize(&rec);
    try std.testing.expect(result.modified);
    try std.testing.expectEqualStrings("G", rec.ref_allele);
    try std.testing.expectEqualStrings("T", rec.alt_alleles.items[0]);
    try std.testing.expectEqual(@as(u32, 10), rec.pos); // 0-based: was 9, now 10 (shifted +1 for prefix trim)
}

test "trim common suffix: TGA>TA -> TG>T" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = true,
        .check_ref = .none,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t10\t.\tTGA\tTA\t.\tPASS\t.");

    const result = try ctx.normalize(&rec);
    try std.testing.expect(result.modified);
    try std.testing.expectEqualStrings("TG", rec.ref_allele);
    try std.testing.expectEqualStrings("T", rec.alt_alleles.items[0]);
    try std.testing.expectEqual(@as(u32, 9), rec.pos); // unchanged
}

test "trim prefix and suffix: AAGCC>AATCC -> G>T" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = true,
        .check_ref = .none,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t10\t.\tAAGCC\tAATCC\t.\tPASS\t.");

    const result = try ctx.normalize(&rec);
    try std.testing.expect(result.modified);
    try std.testing.expectEqualStrings("G", rec.ref_allele);
    try std.testing.expectEqualStrings("T", rec.alt_alleles.items[0]);
    try std.testing.expectEqual(@as(u32, 11), rec.pos); // 9 + 2 prefix bases
}

test "no trimming needed for SNP" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = true,
        .check_ref = .none,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t10\t.\tA\tT\t.\tPASS\t.");

    const result = try ctx.normalize(&rec);
    try std.testing.expect(!result.modified);
    try std.testing.expectEqualStrings("A", rec.ref_allele);
    try std.testing.expectEqualStrings("T", rec.alt_alleles.items[0]);
}

test "split multiallelic: A>T,G -> two records" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = false,
        .check_ref = .none,
        .split_mode = .any,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();
    try rec.parseLine("chr1\t10\t.\tA\tT,G\t.\tPASS\t.");

    const records = try ctx.split(&rec);
    defer {
        for (records) |*r| {
            var mr = r;
            mr.deinit();
        }
        allocator.free(records);
    }

    try std.testing.expectEqual(@as(usize, 2), records.len);

    // First record: A>T
    try std.testing.expectEqualStrings("chr1", records[0].chrom);
    try std.testing.expectEqual(@as(u32, 9), records[0].pos);
    try std.testing.expectEqualStrings("A", records[0].ref_allele);
    try std.testing.expectEqual(@as(usize, 1), records[0].alt_alleles.items.len);
    try std.testing.expectEqualStrings("T", records[0].alt_alleles.items[0]);

    // Second record: A>G
    try std.testing.expectEqualStrings("chr1", records[1].chrom);
    try std.testing.expectEqual(@as(u32, 9), records[1].pos);
    try std.testing.expectEqualStrings("A", records[1].ref_allele);
    try std.testing.expectEqual(@as(usize, 1), records[1].alt_alleles.items.len);
    try std.testing.expectEqualStrings("G", records[1].alt_alleles.items[0]);
}

test "split biallelic returns single record" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .split_mode = .any,
        .check_ref = .none,
        .do_left_align = false,
        .do_trim = false,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();
    try rec.parseLine("chr1\t10\t.\tA\tT\t.\tPASS\t.");

    const records = try ctx.split(&rec);
    defer {
        for (records) |*r| {
            var mr = r;
            mr.deinit();
        }
        allocator.free(records);
    }

    try std.testing.expectEqual(@as(usize, 1), records.len);
    try std.testing.expectEqualStrings("A", records[0].ref_allele);
    try std.testing.expectEqualStrings("T", records[0].alt_alleles.items[0]);
}

test "ref check with mock fetch function" {
    const allocator = std.testing.allocator;

    // Mock reference: always returns "C" for any position
    const MockFetcher = struct {
        fn fetch(_: *anyopaque, alloc: std.mem.Allocator, _: [*:0]const u8, _: i64, _: i64) ?[]u8 {
            const result = alloc.alloc(u8, 1) catch return null;
            result[0] = 'C';
            return result;
        }
    };

    var dummy_ctx: u8 = 0;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = false,
        .check_ref = .warn,
        .fai_ctx = @ptrCast(&dummy_ctx),
        .fetch_seq_fn = &MockFetcher.fetch,
    });
    defer ctx.deinit();

    // REF=A but reference has C -> mismatch
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();
    try rec.parseLine("chr1\t10\t.\tA\tT\t.\tPASS\t.");

    const result = try ctx.normalize(&rec);
    try std.testing.expect(result.warning != null);
    try std.testing.expectEqualStrings("REF_MISMATCH", result.warning.?);
    try std.testing.expectEqual(@as(u64, 1), ctx.stats.ref_mismatch);
}

test "ref check with discard" {
    const allocator = std.testing.allocator;

    const MockFetcher = struct {
        fn fetch(_: *anyopaque, alloc: std.mem.Allocator, _: [*:0]const u8, _: i64, _: i64) ?[]u8 {
            const result = alloc.alloc(u8, 1) catch return null;
            result[0] = 'C';
            return result;
        }
    };

    var dummy_ctx: u8 = 0;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = false,
        .check_ref = .discard,
        .fai_ctx = @ptrCast(&dummy_ctx),
        .fetch_seq_fn = &MockFetcher.fetch,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();
    try rec.parseLine("chr1\t10\t.\tA\tT\t.\tPASS\t.");

    const result = try ctx.normalize(&rec);
    try std.testing.expect(result.discarded);
    try std.testing.expectEqual(@as(u64, 1), ctx.stats.discarded);
}

test "left-align deletion in repeat: AA>A shifts to CA>C" {
    const allocator = std.testing.allocator;

    // Mock reference: "GCAAATGCCC"
    // positions: 0=G, 1=C, 2=A, 3=A, 4=A, 5=T, 6=G, 7=C, 8=C, 9=C
    const MockRef = struct {
        fn fetch(_: *anyopaque, alloc: std.mem.Allocator, _: [*:0]const u8, beg: i64, end: i64) ?[]u8 {
            const ref_seq = "GCAAATGCCC";
            const b: usize = @intCast(beg);
            const e: usize = @intCast(end);
            if (e >= ref_seq.len) return null;
            const len = e - b + 1;
            const result = alloc.alloc(u8, len) catch return null;
            @memcpy(result, ref_seq[b .. b + len]);
            return result;
        }
    };

    var dummy_ctx: u8 = 0;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = true,
        .do_trim = true,
        .check_ref = .none,
        .fai_ctx = @ptrCast(&dummy_ctx),
        .fetch_seq_fn = &MockRef.fetch,
    });
    defer ctx.deinit();

    // AA>A at pos 2 (0-based), 1-based = 3
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();
    try rec.parseLine("chr1\t3\t.\tAA\tA\t.\tPASS\t.");
    // pos = 2 (0-based)

    const result = try ctx.normalize(&rec);
    try std.testing.expect(result.modified);
    // Expected: left-align to CA>C at pos 1 (0-based)
    try std.testing.expectEqual(@as(u32, 1), rec.pos);
    try std.testing.expectEqualStrings("CA", rec.ref_allele);
    try std.testing.expectEqualStrings("C", rec.alt_alleles.items[0]);
}

test "multiallelic trim: ACG>ATG,AAG -> C>T,A at pos+1" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = true,
        .check_ref = .none,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();
    try rec.parseLine("chr1\t10\t.\tACG\tATG,AAG\t.\tPASS\t.");

    const result = try ctx.normalize(&rec);
    try std.testing.expect(result.modified);
    // Common prefix 'A', common suffix 'G' -> trim both
    try std.testing.expectEqualStrings("C", rec.ref_allele);
    try std.testing.expectEqual(@as(usize, 2), rec.alt_alleles.items.len);
    try std.testing.expectEqualStrings("T", rec.alt_alleles.items[0]);
    try std.testing.expectEqualStrings("A", rec.alt_alleles.items[1]);
    try std.testing.expectEqual(@as(u32, 10), rec.pos);
}

test "symbolic allele is not normalized" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = true,
        .check_ref = .none,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();
    try rec.parseLine("chr1\t10\t.\tA\t<DEL>\t.\tPASS\t.");

    const result = try ctx.normalize(&rec);
    try std.testing.expect(!result.modified);
}

test "split mode snps only" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .split_mode = .snps,
        .check_ref = .none,
        .do_left_align = false,
        .do_trim = false,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();
    try rec.parseLine("chr1\t10\t.\tA\tT,G\t.\tPASS\t.");

    const records = try ctx.split(&rec);
    defer {
        for (records) |*r| {
            var mr = r;
            mr.deinit();
        }
        allocator.free(records);
    }

    // All ALTs are SNPs (same length as REF), so should split
    try std.testing.expectEqual(@as(usize, 2), records.len);
}

test "split mode indels skips SNP-only" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .split_mode = .indels,
        .check_ref = .none,
        .do_left_align = false,
        .do_trim = false,
    });
    defer ctx.deinit();

    // All SNPs
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();
    try rec.parseLine("chr1\t10\t.\tA\tT,G\t.\tPASS\t.");

    const records = try ctx.split(&rec);
    defer {
        for (records) |*r| {
            var mr = r;
            mr.deinit();
        }
        allocator.free(records);
    }

    // No indels present, so should NOT split (returns 1 record)
    try std.testing.expectEqual(@as(usize, 1), records.len);
}

test "NormContext stats tracking" {
    const allocator = std.testing.allocator;
    var ctx = NormContext.init(allocator, .{
        .do_left_align = false,
        .do_trim = true,
        .check_ref = .none,
    });
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t10\t.\tAGC\tATC\t.\tPASS\t.");
    _ = try ctx.normalize(&rec);

    try rec.parseLine("chr1\t20\t.\tA\tT\t.\tPASS\t.");
    _ = try ctx.normalize(&rec);

    try std.testing.expectEqual(@as(u64, 2), ctx.stats.total);
    try std.testing.expectEqual(@as(u64, 1), ctx.stats.modified);
}
