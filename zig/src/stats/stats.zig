const std = @import("std");
const vcf_record = @import("../vcf/record.zig");
const VcfRecord = vcf_record.VcfRecord;

/// Per-sample statistics — stored as array-of-structs for cache locality
/// (addresses performance audit item #13: the C code uses struct-of-arrays
/// layout with separate smpl_hets[], smpl_homRR[], etc. arrays that scatter
/// accesses across memory; this AoS layout keeps all per-sample counters
/// in a single cache line).
pub const SampleStats = struct {
    n_variants: u64 = 0,
    n_het: u64 = 0,
    n_hom_alt: u64 = 0,
    n_hom_ref: u64 = 0,
    n_missing: u64 = 0,
    n_singletons: u64 = 0,
    n_ts: u64 = 0,
    n_tv: u64 = 0,
};

/// Indel length histogram offset: maps indel length [-50..50] to array index [0..100].
const indel_offset: i32 = 50;
const indel_hist_len: usize = 101;

/// Quality histogram: bins 0..255, each bin covers QUAL values [bin, bin+1).
const qual_hist_len: usize = 256;

/// Core statistics context for the `stats` subcommand.
///
/// Accumulates variant-level and per-sample counts from a stream of VCF records.
/// The design mirrors the C `stats_t` struct but uses array-of-structs layout
/// for per-sample data (audit #13).
pub const StatsContext = struct {
    allocator: std.mem.Allocator,

    // ---- Variant-level counters ----
    n_snps: u64 = 0,
    n_indels: u64 = 0,
    n_mnps: u64 = 0,
    n_other: u64 = 0,
    n_records: u64 = 0,

    // ---- Ts/Tv ----
    n_transitions: u64 = 0,
    n_transversions: u64 = 0,

    // ---- Indel length histogram: index = length + indel_offset ----
    indel_lengths: [indel_hist_len]u64 = [_]u64{0} ** indel_hist_len,

    // ---- Per-sample stats (AoS, audit #13) ----
    sample_stats: []SampleStats,
    sample_names: []const []const u8,
    n_samples: u32,

    // ---- Quality histogram ----
    qual_hist: [qual_hist_len]u64 = [_]u64{0} ** qual_hist_len,

    pub fn init(allocator: std.mem.Allocator, n_samples: u32, sample_names: []const []const u8) !StatsContext {
        const stats = try allocator.alloc(SampleStats, n_samples);
        @memset(stats, SampleStats{});

        return StatsContext{
            .allocator = allocator,
            .sample_stats = stats,
            .sample_names = sample_names,
            .n_samples = n_samples,
        };
    }

    pub fn deinit(self: *StatsContext) void {
        self.allocator.free(self.sample_stats);
        self.* = undefined;
    }

    /// Process a single VCF record, updating all counters.
    ///
    /// `raw_line` is the full VCF text line (tab-separated) needed for parsing
    /// FORMAT + sample genotype columns that VcfRecord.parseLine currently skips.
    pub fn addRecord(self: *StatsContext, rec: *const VcfRecord, raw_line: ?[]const u8) void {
        self.n_records += 1;

        // Quality histogram
        if (rec.qual) |q| {
            const bin: usize = @intFromFloat(@min(@max(q, 0.0), @as(f32, @floatFromInt(qual_hist_len - 1))));
            self.qual_hist[bin] += 1;
        }

        // Classify each ALT allele
        const ref = rec.ref_allele;
        for (rec.alt_alleles.items) |alt| {
            if (alt.len == 0 or std.mem.eql(u8, alt, ".") or std.mem.eql(u8, alt, "*")) continue;

            const vtype = classifyVariant(ref, alt);
            switch (vtype) {
                .snp => {
                    self.n_snps += 1;
                    if (isTransition(ref[0], alt[0])) {
                        self.n_transitions += 1;
                    } else {
                        self.n_transversions += 1;
                    }
                },
                .indel => {
                    self.n_indels += 1;
                    const ilen = indelLength(ref, alt);
                    const idx = @as(usize, @intCast(std.math.clamp(ilen + indel_offset, 0, @as(i32, @intCast(indel_hist_len - 1)))));
                    self.indel_lengths[idx] += 1;
                },
                .mnp => {
                    self.n_mnps += 1;
                },
                .other => {
                    self.n_other += 1;
                },
            }
        }

        // Per-sample genotype stats
        if (self.n_samples > 0) {
            if (raw_line) |line| {
                self.processSampleGenotypes(rec, line);
            }
        }
    }

    /// Parse FORMAT + sample columns from the raw VCF line and update per-sample stats.
    fn processSampleGenotypes(self: *StatsContext, rec: *const VcfRecord, raw_line: []const u8) void {
        // Find the FORMAT column (column 8, 0-indexed) and sample columns (9+).
        // We need to skip the first 8 columns (CHROM through INFO).
        var col: u32 = 0;
        var pos: usize = 0;
        var format_start: usize = 0;
        var format_end: usize = 0;
        var samples_start: usize = 0;

        for (raw_line, 0..) |c, i| {
            if (c == '\t') {
                col += 1;
                if (col == 8) {
                    format_start = i + 1;
                } else if (col == 9) {
                    format_end = i;
                    samples_start = i + 1;
                    pos = i + 1;
                    break;
                }
            }
        }
        if (col < 9) return; // no FORMAT or samples columns

        // Find GT position within FORMAT (it should be the first sub-field)
        const format_field = raw_line[format_start..format_end];
        const gt_index = findFormatField(format_field, "GT") orelse return;

        // Determine variant type for first ALT (for per-sample ts/tv)
        const ref = rec.ref_allele;
        const first_alt_type: ?VariantType = if (rec.alt_alleles.items.len > 0)
            classifyVariant(ref, rec.alt_alleles.items[0])
        else
            null;

        const first_alt_is_ts: bool = if (first_alt_type != null and first_alt_type.? == .snp and rec.alt_alleles.items.len > 0)
            isTransition(ref[0], rec.alt_alleles.items[0][0])
        else
            false;

        // Count allele occurrences across all samples (for singleton detection)
        var allele_counts_buf: [256]u32 = [_]u32{0} ** 256;
        var sample_alleles_buf: [256][2]i16 = undefined;
        const max_samples = @min(self.n_samples, 256);

        // First pass: parse genotypes and count alleles
        var sample_idx: u32 = 0;
        var field_pos: usize = samples_start;

        while (sample_idx < max_samples and field_pos < raw_line.len) {
            const sample_end = findNextTab(raw_line, field_pos) orelse raw_line.len;
            const sample_data = raw_line[field_pos..sample_end];

            const gt_str = extractSubField(sample_data, gt_index);
            const parsed = parseGenotype(gt_str);

            sample_alleles_buf[sample_idx] = .{ parsed.allele1, parsed.allele2 };

            if (parsed.allele1 >= 0) {
                const a1: usize = @intCast(parsed.allele1);
                if (a1 < allele_counts_buf.len) allele_counts_buf[a1] += 1;
            }
            if (parsed.allele2 >= 0) {
                const a2: usize = @intCast(parsed.allele2);
                if (a2 < allele_counts_buf.len) allele_counts_buf[a2] += 1;
            }

            sample_idx += 1;
            field_pos = if (sample_end < raw_line.len) sample_end + 1 else raw_line.len;
        }

        // Second pass: classify genotypes and detect singletons
        for (0..sample_idx) |si| {
            const a1 = sample_alleles_buf[si][0];
            const a2 = sample_alleles_buf[si][1];
            var ss = &self.sample_stats[si];

            if (a1 < 0 or a2 < 0) {
                // Missing genotype
                ss.n_missing += 1;
                continue;
            }

            if (a1 == 0 and a2 == 0) {
                // Hom ref
                ss.n_hom_ref += 1;
                continue;
            }

            ss.n_variants += 1;

            if (a1 == a2) {
                // Hom alt
                ss.n_hom_alt += 1;
            } else if (a1 == 0 or a2 == 0) {
                // Het (one ref, one alt)
                ss.n_het += 1;
            } else {
                // Het (two different non-ref alleles)
                ss.n_het += 1;
            }

            // Per-sample ts/tv (based on first alt allele presence)
            if (first_alt_type != null and first_alt_type.? == .snp) {
                const has_alt1 = (a1 == 1 or a2 == 1);
                if (has_alt1) {
                    if (first_alt_is_ts) {
                        ss.n_ts += 1;
                    } else {
                        ss.n_tv += 1;
                    }
                }
            }

            // Singleton detection: allele appears exactly once across all samples
            const alt_allele = if (a1 != 0) a1 else a2;
            if (alt_allele > 0 and alt_allele < @as(i16, @intCast(allele_counts_buf.len))) {
                if (allele_counts_buf[@intCast(alt_allele)] == 1) {
                    ss.n_singletons += 1;
                }
            }
        }
    }

    /// Write the stats report in bcftools stats text format.
    pub fn writeReport(self: *const StatsContext, writer: *std.Io.Writer) !void {
        // Summary numbers
        try writer.writeAll("# SN, Summary numbers.\n");
        try writer.writeAll("# Use 'grep ^SN | cut -f 2-' to extract this part.\n");
        try writer.print("SN\t0\tnumber of records:\t{d}\n", .{self.n_records});
        try writer.print("SN\t0\tnumber of SNPs:\t{d}\n", .{self.n_snps});
        try writer.print("SN\t0\tnumber of indels:\t{d}\n", .{self.n_indels});
        try writer.print("SN\t0\tnumber of MNPs:\t{d}\n", .{self.n_mnps});
        try writer.print("SN\t0\tnumber of others:\t{d}\n", .{self.n_other});

        // Ts/Tv ratio
        const tstv: f64 = if (self.n_transversions > 0)
            @as(f64, @floatFromInt(self.n_transitions)) / @as(f64, @floatFromInt(self.n_transversions))
        else
            0.0;
        try writer.print("SN\t0\tts/tv:\t{d:.2}\n", .{tstv});

        // Indel length distribution
        try writer.writeAll("# IDD, InDel distribution.\n");
        try writer.writeAll("# Use 'grep ^IDD | cut -f 2-' to extract this part.\n");
        for (0..indel_hist_len) |i| {
            const count = self.indel_lengths[i];
            if (count > 0) {
                const length: i32 = @as(i32, @intCast(i)) - indel_offset;
                try writer.print("IDD\t0\t{d}\t{d}\n", .{ length, count });
            }
        }

        // Quality distribution
        try writer.writeAll("# QUAL, Variant quality distribution.\n");
        try writer.writeAll("# Use 'grep ^QUAL | cut -f 2-' to extract this part.\n");
        for (0..qual_hist_len) |i| {
            const count = self.qual_hist[i];
            if (count > 0) {
                try writer.print("QUAL\t0\t{d}\t{d}\n", .{ i, count });
            }
        }

        // Per-sample counts
        if (self.n_samples > 0) {
            try writer.writeAll("# PSC, Per-sample counts.\n");
            try writer.writeAll("# Use 'grep ^PSC | cut -f 2-' to extract this part.\n");
            try writer.writeAll("# PSC\t[2]id\t[3]sample\t[4]nVariants\t[5]nHet\t[6]nHomAlt\t[7]nHomRef\t[8]nMissing\t[9]nSingletons\t[10]ts/tv\n");

            for (0..self.n_samples) |i| {
                const ss = &self.sample_stats[i];
                const sample_name: []const u8 = if (i < self.sample_names.len)
                    self.sample_names[i]
                else
                    "unknown";

                const sample_tstv: f64 = if (ss.n_tv > 0)
                    @as(f64, @floatFromInt(ss.n_ts)) / @as(f64, @floatFromInt(ss.n_tv))
                else
                    0.0;

                try writer.print("PSC\t{d}\t{s}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d:.2}\n", .{
                    i,
                    sample_name,
                    ss.n_variants,
                    ss.n_het,
                    ss.n_hom_alt,
                    ss.n_hom_ref,
                    ss.n_missing,
                    ss.n_singletons,
                    sample_tstv,
                });
            }
        }
    }
};

// =========================================================================
// Variant classification helpers
// =========================================================================

pub const VariantType = enum {
    snp,
    indel,
    mnp,
    other,
};

/// Classify a variant by comparing REF and ALT allele strings.
pub fn classifyVariant(ref: []const u8, alt: []const u8) VariantType {
    if (ref.len == 0 or alt.len == 0) return .other;

    // Symbolic alleles (e.g., <DEL>, <INS>)
    if (alt[0] == '<') return .other;

    // Breakend notation
    if (alt[0] == ']' or alt[0] == '[') return .other;
    if (alt.len > 1 and (alt[alt.len - 1] == ']' or alt[alt.len - 1] == '[')) return .other;

    if (ref.len == 1 and alt.len == 1) {
        return .snp;
    } else if (ref.len == alt.len) {
        return .mnp;
    } else {
        return .indel;
    }
}

/// Compute indel length: positive = insertion, negative = deletion.
/// Uses simple length difference (like the C code for VCF-normalized indels).
pub fn indelLength(ref: []const u8, alt: []const u8) i32 {
    return @as(i32, @intCast(alt.len)) - @as(i32, @intCast(ref.len));
}

/// Test if a single-base substitution is a transition (purine<->purine or
/// pyrimidine<->pyrimidine).
pub fn isTransition(ref_base: u8, alt_base: u8) bool {
    const r = std.ascii.toUpper(ref_base);
    const a = std.ascii.toUpper(alt_base);

    return (r == 'A' and a == 'G') or
        (r == 'G' and a == 'A') or
        (r == 'C' and a == 'T') or
        (r == 'T' and a == 'C');
}

// =========================================================================
// Genotype parsing helpers
// =========================================================================

const ParsedGT = struct {
    allele1: i16, // -1 = missing
    allele2: i16, // -1 = missing
};

/// Parse a diploid GT string like "0/1", "0|1", "1/1", "./.", "0", etc.
fn parseGenotype(gt_str: ?[]const u8) ParsedGT {
    const s = gt_str orelse return .{ .allele1 = -1, .allele2 = -1 };
    if (s.len == 0) return .{ .allele1 = -1, .allele2 = -1 };

    // Find separator (/ or |)
    var sep_pos: ?usize = null;
    for (s, 0..) |c, i| {
        if (c == '/' or c == '|') {
            sep_pos = i;
            break;
        }
    }

    if (sep_pos) |sp| {
        const a1 = parseAlleleIndex(s[0..sp]);
        const a2 = parseAlleleIndex(s[sp + 1 ..]);
        return .{ .allele1 = a1, .allele2 = a2 };
    } else {
        // Haploid
        const a1 = parseAlleleIndex(s);
        return .{ .allele1 = a1, .allele2 = a1 };
    }
}

/// Parse a single allele index from a GT sub-field. Returns -1 for missing (".").
fn parseAlleleIndex(s: []const u8) i16 {
    if (s.len == 0 or s[0] == '.') return -1;
    return std.fmt.parseInt(i16, s, 10) catch -1;
}

/// Find the 0-based index of a sub-field name within a colon-separated FORMAT string.
fn findFormatField(format: []const u8, name: []const u8) ?u32 {
    var idx: u32 = 0;
    var start: usize = 0;
    for (format, 0..) |c, i| {
        if (c == ':') {
            if (std.mem.eql(u8, format[start..i], name)) return idx;
            idx += 1;
            start = i + 1;
        }
    }
    // Check last field
    if (std.mem.eql(u8, format[start..], name)) return idx;
    return null;
}

/// Extract the nth colon-separated sub-field from a sample data string.
fn extractSubField(sample_data: []const u8, field_idx: u32) ?[]const u8 {
    var idx: u32 = 0;
    var start: usize = 0;
    for (sample_data, 0..) |c, i| {
        if (c == ':') {
            if (idx == field_idx) return sample_data[start..i];
            idx += 1;
            start = i + 1;
        }
    }
    if (idx == field_idx) return sample_data[start..];
    return null;
}

/// Find the position of the next tab character starting from `start`.
fn findNextTab(data: []const u8, start: usize) ?usize {
    for (data[start..], start..) |c, i| {
        if (c == '\t') return i;
    }
    return null;
}

// =========================================================================
// Tests
// =========================================================================

test "classifyVariant" {
    const testing = std.testing;
    try testing.expectEqual(VariantType.snp, classifyVariant("A", "G"));
    try testing.expectEqual(VariantType.snp, classifyVariant("C", "T"));
    try testing.expectEqual(VariantType.indel, classifyVariant("A", "AT"));
    try testing.expectEqual(VariantType.indel, classifyVariant("ATG", "A"));
    try testing.expectEqual(VariantType.mnp, classifyVariant("AT", "GC"));
    try testing.expectEqual(VariantType.other, classifyVariant("A", "<DEL>"));
    try testing.expectEqual(VariantType.other, classifyVariant("", "A"));
}

test "isTransition" {
    const testing = std.testing;
    try testing.expect(isTransition('A', 'G'));
    try testing.expect(isTransition('G', 'A'));
    try testing.expect(isTransition('C', 'T'));
    try testing.expect(isTransition('T', 'C'));
    try testing.expect(!isTransition('A', 'T'));
    try testing.expect(!isTransition('A', 'C'));
    try testing.expect(!isTransition('G', 'T'));
    try testing.expect(!isTransition('G', 'C'));
}

test "indelLength" {
    const testing = std.testing;
    try testing.expectEqual(@as(i32, 1), indelLength("A", "AT"));
    try testing.expectEqual(@as(i32, -2), indelLength("ATG", "A"));
    try testing.expectEqual(@as(i32, 3), indelLength("A", "ATCG"));
}

test "parseGenotype" {
    const testing = std.testing;

    const gt1 = parseGenotype("0/1");
    try testing.expectEqual(@as(i16, 0), gt1.allele1);
    try testing.expectEqual(@as(i16, 1), gt1.allele2);

    const gt2 = parseGenotype("1|1");
    try testing.expectEqual(@as(i16, 1), gt2.allele1);
    try testing.expectEqual(@as(i16, 1), gt2.allele2);

    const gt3 = parseGenotype("./.");
    try testing.expectEqual(@as(i16, -1), gt3.allele1);
    try testing.expectEqual(@as(i16, -1), gt3.allele2);

    const gt4 = parseGenotype("0");
    try testing.expectEqual(@as(i16, 0), gt4.allele1);
    try testing.expectEqual(@as(i16, 0), gt4.allele2);
}

test "findFormatField" {
    const testing = std.testing;
    try testing.expectEqual(@as(?u32, 0), findFormatField("GT:DP:GQ", "GT"));
    try testing.expectEqual(@as(?u32, 1), findFormatField("GT:DP:GQ", "DP"));
    try testing.expectEqual(@as(?u32, 2), findFormatField("GT:DP:GQ", "GQ"));
    try testing.expectEqual(@as(?u32, null), findFormatField("GT:DP:GQ", "AD"));
}

test "extractSubField" {
    const testing = std.testing;
    try testing.expectEqualStrings("0/1", extractSubField("0/1:30:99", 0).?);
    try testing.expectEqualStrings("30", extractSubField("0/1:30:99", 1).?);
    try testing.expectEqualStrings("99", extractSubField("0/1:30:99", 2).?);
    try testing.expectEqual(@as(?[]const u8, null), extractSubField("0/1:30:99", 3));
}

test "SNP counting" {
    const allocator = std.testing.allocator;
    const empty_names: []const []const u8 = &.{};
    var ctx = try StatsContext.init(allocator, 0, empty_names);
    defer ctx.deinit();

    // Parse 3 SNPs
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t100\t.\tA\tG\t30\tPASS\t.");
    ctx.addRecord(&rec, null);
    try rec.parseLine("chr1\t200\t.\tC\tT\t40\tPASS\t.");
    ctx.addRecord(&rec, null);
    try rec.parseLine("chr1\t300\t.\tG\tA\t50\tPASS\t.");
    ctx.addRecord(&rec, null);

    try std.testing.expectEqual(@as(u64, 3), ctx.n_snps);
    try std.testing.expectEqual(@as(u64, 3), ctx.n_records);
    try std.testing.expectEqual(@as(u64, 0), ctx.n_indels);
}

test "Ts/Tv ratio" {
    const allocator = std.testing.allocator;
    const empty_names: []const []const u8 = &.{};
    var ctx = try StatsContext.init(allocator, 0, empty_names);
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    // A>G is transition
    try rec.parseLine("chr1\t100\t.\tA\tG\t.\tPASS\t.");
    ctx.addRecord(&rec, null);
    // A>T is transversion
    try rec.parseLine("chr1\t200\t.\tA\tT\t.\tPASS\t.");
    ctx.addRecord(&rec, null);

    try std.testing.expectEqual(@as(u64, 1), ctx.n_transitions);
    try std.testing.expectEqual(@as(u64, 1), ctx.n_transversions);
    // ratio should be 1.0
}

test "Indel length distribution" {
    const allocator = std.testing.allocator;
    const empty_names: []const []const u8 = &.{};
    var ctx = try StatsContext.init(allocator, 0, empty_names);
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    // Insertion of length 1: A -> AT
    try rec.parseLine("chr1\t100\t.\tA\tAT\t.\tPASS\t.");
    ctx.addRecord(&rec, null);

    // Deletion of length 2: ATG -> A
    try rec.parseLine("chr1\t200\t.\tATG\tA\t.\tPASS\t.");
    ctx.addRecord(&rec, null);

    try std.testing.expectEqual(@as(u64, 2), ctx.n_indels);
    // Length +1 -> index 51
    try std.testing.expectEqual(@as(u64, 1), ctx.indel_lengths[51]);
    // Length -2 -> index 48
    try std.testing.expectEqual(@as(u64, 1), ctx.indel_lengths[48]);
}

test "Per-sample het/hom counting" {
    const allocator = std.testing.allocator;
    const names = [_][]const u8{ "SAMPLE1", "SAMPLE2" };
    const names_slice: []const []const u8 = &names;
    var ctx = try StatsContext.init(allocator, 2, names_slice);
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    const line1 = "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\t1/1";
    try rec.parseLine(line1);
    ctx.addRecord(&rec, line1);

    // SAMPLE1: 0/1 = het
    try std.testing.expectEqual(@as(u64, 1), ctx.sample_stats[0].n_het);
    try std.testing.expectEqual(@as(u64, 0), ctx.sample_stats[0].n_hom_alt);
    // SAMPLE2: 1/1 = hom alt
    try std.testing.expectEqual(@as(u64, 0), ctx.sample_stats[1].n_het);
    try std.testing.expectEqual(@as(u64, 1), ctx.sample_stats[1].n_hom_alt);
}

test "Quality histogram" {
    const allocator = std.testing.allocator;
    const empty_names: []const []const u8 = &.{};
    var ctx = try StatsContext.init(allocator, 0, empty_names);
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t100\t.\tA\tG\t30\tPASS\t.");
    ctx.addRecord(&rec, null);
    try rec.parseLine("chr1\t200\t.\tC\tT\t30\tPASS\t.");
    ctx.addRecord(&rec, null);
    try rec.parseLine("chr1\t300\t.\tG\tA\t50\tPASS\t.");
    ctx.addRecord(&rec, null);

    try std.testing.expectEqual(@as(u64, 2), ctx.qual_hist[30]);
    try std.testing.expectEqual(@as(u64, 1), ctx.qual_hist[50]);
}

test "writeReport produces expected output" {
    const allocator = std.testing.allocator;
    const empty_names: []const []const u8 = &.{};
    var ctx = try StatsContext.init(allocator, 0, empty_names);
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t100\t.\tA\tG\t30\tPASS\t.");
    ctx.addRecord(&rec, null);

    var aw: std.Io.Writer.Allocating = .init(allocator);
    defer aw.deinit();
    try ctx.writeReport(&aw.writer);

    const al = aw.toArrayList();
    defer allocator.free(al.allocatedSlice());
    const output = al.items;
    // Should contain SN lines
    try std.testing.expect(std.mem.indexOf(u8, output, "SN\t0\tnumber of SNPs:\t1") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "SN\t0\tnumber of indels:\t0") != null);
}

test "singleton detection" {
    const allocator = std.testing.allocator;
    const names = [_][]const u8{ "S1", "S2", "S3" };
    const names_slice: []const []const u8 = &names;
    var ctx = try StatsContext.init(allocator, 3, names_slice);
    defer ctx.deinit();

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    // S1 is 0/1, S2 is 0/0, S3 is 0/0: alt allele 1 appears once -> singleton for S1
    const line1 = "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\t0/0\t0/0";
    try rec.parseLine(line1);
    ctx.addRecord(&rec, line1);

    try std.testing.expectEqual(@as(u64, 1), ctx.sample_stats[0].n_singletons);
    try std.testing.expectEqual(@as(u64, 0), ctx.sample_stats[1].n_singletons);
}
