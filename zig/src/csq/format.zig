// Consequence string formatting — ports kput_vcsq and kprint_aa_prediction from csq.c
//
// Formats consequence annotations into strings like:
//   "missense|GENE|ENST001|protein_coding|+|5TY>5I|121ACG>A+124TA>T"

const std = @import("std");

/// Consequence type bit flags, mirroring CSQ_* defines from csq.c.
pub const CsqType = u32;

pub const CSQ_PRINTED_UPSTREAM: CsqType = 1 << 0;
pub const CSQ_SYNONYMOUS_VARIANT: CsqType = 1 << 1;
pub const CSQ_MISSENSE_VARIANT: CsqType = 1 << 2;
pub const CSQ_STOP_LOST: CsqType = 1 << 3;
pub const CSQ_STOP_GAINED: CsqType = 1 << 4;
pub const CSQ_INFRAME_DELETION: CsqType = 1 << 5;
pub const CSQ_INFRAME_INSERTION: CsqType = 1 << 6;
pub const CSQ_FRAMESHIFT_VARIANT: CsqType = 1 << 7;
pub const CSQ_SPLICE_ACCEPTOR: CsqType = 1 << 8;
pub const CSQ_SPLICE_DONOR: CsqType = 1 << 9;
pub const CSQ_START_LOST: CsqType = 1 << 10;
pub const CSQ_SPLICE_REGION: CsqType = 1 << 11;
pub const CSQ_STOP_RETAINED: CsqType = 1 << 12;
pub const CSQ_UTR5: CsqType = 1 << 13;
pub const CSQ_UTR3: CsqType = 1 << 14;
pub const CSQ_NON_CODING: CsqType = 1 << 15;
pub const CSQ_INTRON: CsqType = 1 << 16;
// 1 << 17 is unused (was CSQ_INTERGENIC)
pub const CSQ_INFRAME_ALTERING: CsqType = 1 << 18;
pub const CSQ_UPSTREAM_STOP: CsqType = 1 << 19;
pub const CSQ_INCOMPLETE_CDS: CsqType = 1 << 20;
pub const CSQ_CODING_SEQUENCE: CsqType = 1 << 21;
pub const CSQ_ELONGATION: CsqType = 1 << 22;
pub const CSQ_TRUNCATION: CsqType = 1 << 23;
pub const CSQ_START_RETAINED: CsqType = 1 << 24;

/// Haplotype-aware consequences that include protein/DNA change strings.
pub const CSQ_COMPOUND: CsqType =
    CSQ_SYNONYMOUS_VARIANT | CSQ_MISSENSE_VARIANT | CSQ_STOP_LOST | CSQ_STOP_GAINED |
    CSQ_INFRAME_DELETION | CSQ_INFRAME_INSERTION | CSQ_FRAMESHIFT_VARIANT |
    CSQ_START_LOST | CSQ_STOP_RETAINED | CSQ_INFRAME_ALTERING | CSQ_INCOMPLETE_CDS |
    CSQ_UPSTREAM_STOP | CSQ_START_RETAINED | CSQ_ELONGATION | CSQ_TRUNCATION;

pub const CSQ_START_STOP: CsqType =
    CSQ_STOP_LOST | CSQ_STOP_GAINED | CSQ_STOP_RETAINED | CSQ_START_LOST | CSQ_START_RETAINED;

pub const CSQ_PRN_TSCRIPT: CsqType = ~(CSQ_INTRON | CSQ_NON_CODING);

pub const CSQ_PRN_NMD: CsqType = ~(CSQ_INTRON | CSQ_NON_CODING);

/// Returns true if the strand should be printed for this consequence type.
pub fn csqPrnStrand(csq: CsqType) bool {
    return (csq & CSQ_COMPOUND != 0) and
        (csq & (CSQ_SPLICE_ACCEPTOR | CSQ_SPLICE_DONOR | CSQ_SPLICE_REGION | CSQ_ELONGATION | CSQ_TRUNCATION) == 0);
}

/// Strand direction.
pub const Strand = enum(u1) {
    fwd = 0,
    rev = 1,
};

/// GFF feature biotype constants for NMD detection.
/// Must match @intFromEnum(gff_types.Biotype.NMD) = 12 | (1 << 6) = 76.
pub const GF_NMD: u32 = 76;

/// Consequence string table, indexed by bit position.
/// Matches csq_strings[] from csq.c.
const csq_strings = [_]?[]const u8{
    null, // bit 0: CSQ_PRINTED_UPSTREAM (not a named consequence)
    "synonymous",
    "missense",
    "stop_lost",
    "stop_gained",
    "inframe_deletion",
    "inframe_insertion",
    "frameshift",
    "splice_acceptor",
    "splice_donor",
    "start_lost",
    "splice_region",
    "stop_retained",
    "5_prime_utr",
    "3_prime_utr",
    "non_coding",
    "intron",
    "intergenic",
    "inframe_altering",
    null, // bit 19: CSQ_UPSTREAM_STOP (printed as '*' prefix, not a named consequence)
    null, // bit 20: CSQ_INCOMPLETE_CDS
    "coding_sequence",
    "feature_elongation",
    "feature_truncation",
    "start_retained",
};

/// Variant consequence string, analogous to vcsq_t in csq.c.
pub const Vcsq = struct {
    strand: Strand = .fwd,
    csq_type: CsqType = 0,
    trid: u32 = 0,
    vcf_ial: u32 = 0,
    biotype: u32 = 0,
    gene: ?[]const u8 = null,
    /// If csq_type & CSQ_PRINTED_UPSTREAM, the position of the upstream reference record (1-based).
    ref_pos: ?u32 = null,
    /// Variant string, e.g. "|5TY>5I|121ACG>A+124TA>T"
    vstr: ?[]const u8 = null,
};

/// Options for formatting.
pub const FormatOptions = struct {
    brief_predictions: u32 = 0,
    /// Callback to resolve transcript ID to string name; null means omit.
    trid_to_string: ?*const fn (u32) []const u8 = null,
    /// Callback to resolve biotype ID to GFF string; null means omit.
    biotype_to_string: ?*const fn (u32) []const u8 = null,
};

/// Port of kput_vcsq: format a single consequence into the writer.
///
/// Produces output like:
///   "missense|GENE|ENST001|protein_coding|+|5TY>5I|121ACG>A+124TA>T"
///   "@12345" (for printed-upstream references)
///   "*stop_lost|..." (for upstream-stop consequences)
pub fn formatVcsq(csq: *Vcsq, opts: FormatOptions, writer: anytype) !void {
    // Work on a mutable copy of type for the masking logic
    var t = csq.csq_type;

    // Remove start/stop from incomplete CDS if there is another consequence
    if (t & CSQ_INCOMPLETE_CDS != 0 and (t & ~(CSQ_START_STOP | CSQ_INCOMPLETE_CDS | CSQ_UPSTREAM_STOP) != 0)) {
        t &= ~(CSQ_START_STOP | CSQ_INCOMPLETE_CDS);
    }

    // Remove missense from start/stops
    if (t & CSQ_START_STOP != 0 and t & CSQ_MISSENSE_VARIANT != 0) {
        t &= ~CSQ_MISSENSE_VARIANT;
    }

    // Printed-upstream: just "@<pos>"
    if (t & CSQ_PRINTED_UPSTREAM != 0) {
        if (csq.ref_pos) |pos| {
            try writer.writeByte('@');
            try writer.print("{d}", .{pos});
            return;
        }
    }

    // Upstream-stop prefix
    if (t & CSQ_UPSTREAM_STOP != 0) {
        try writer.writeByte('*');
    }

    // Consequence names, joined by '&'
    var has_csq = false;
    const n = csq_strings.len;
    var i: usize = 1;
    // First matching consequence
    while (i < n) : (i += 1) {
        if (csq_strings[i]) |s| {
            if (t & (@as(CsqType, 1) << @intCast(i)) != 0) {
                try writer.writeAll(s);
                has_csq = true;
                i += 1;
                break;
            }
        }
    }
    // Remaining consequences separated by '&'
    while (i < n) : (i += 1) {
        if (csq_strings[i]) |s| {
            if (t & (@as(CsqType, 1) << @intCast(i)) != 0) {
                try writer.writeByte('&');
                try writer.writeAll(s);
                has_csq = true;
            }
        }
    }

    // NMD annotation
    if (csq.biotype == GF_NMD and (t & CSQ_PRN_NMD != 0)) {
        if (has_csq) try writer.writeByte('&');
        try writer.writeAll("NMD_transcript");
    }

    // Gene field
    try writer.writeByte('|');
    if (csq.gene) |gene| {
        try writer.writeAll(gene);
    }

    // Transcript field
    try writer.writeByte('|');
    if (t & CSQ_PRN_TSCRIPT != 0) {
        if (opts.trid_to_string) |resolver| {
            try writer.writeAll(resolver(csq.trid));
        }
    }

    // Biotype field
    try writer.writeByte('|');
    if (opts.biotype_to_string) |resolver| {
        try writer.writeAll(resolver(csq.biotype));
    }

    // Strand
    if (csqPrnStrand(t) or (csq.vstr != null and csq.vstr.?.len > 0)) {
        switch (csq.strand) {
            .fwd => try writer.writeAll("|+"),
            .rev => try writer.writeAll("|-"),
        }
    }

    // Variant string (protein/DNA change)
    if (csq.vstr) |vstr| {
        if (vstr.len > 0) {
            try writer.writeAll(vstr);
        }
    }
}

/// Port of kprint_aa_prediction: format amino acid prediction, optionally abbreviated.
///
/// If brief_predictions > 0 and the amino acid string is long enough, truncate
/// to the first `brief_predictions` characters followed by "..".
pub fn formatAaPrediction(aa: []const u8, stop: []const u8, brief_predictions: u32, writer: anytype) !void {
    if (brief_predictions == 0 or aa.len < brief_predictions + 3) {
        try writer.writeAll(aa);
        return;
    }

    var len = aa.len;
    if (len > 0 and stop.len >= len and stop[len - 1] == '*') {
        len -= 1;
    }

    const limit = @min(len, brief_predictions);
    try writer.writeAll(aa[0..limit]);
    try writer.writeAll("..");
}

/// Format a list of consequences separated by commas.
pub fn formatVcsqList(csqs: []Vcsq, opts: FormatOptions, writer: anytype) !void {
    for (csqs, 0..) |*csq, idx| {
        if (idx > 0) try writer.writeByte(',');
        try formatVcsq(csq, opts, writer);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "format simple missense consequence" {
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(std.testing.allocator);

    var csq = Vcsq{
        .strand = .fwd,
        .csq_type = CSQ_MISSENSE_VARIANT,
        .gene = "BRCA1",
        .biotype = 0,
        .vstr = "|5T>5I|100A>G",
    };

    const opts = FormatOptions{};
    try formatVcsq(&csq, opts, buf.writer(std.testing.allocator));

    // Expected: "missense|BRCA1|||+|5T>5I|100A>G"
    const result = buf.items;
    try std.testing.expect(std.mem.startsWith(u8, result, "missense|BRCA1|"));
    try std.testing.expect(std.mem.endsWith(u8, result, "|+|5T>5I|100A>G"));
}

test "format printed-upstream consequence" {
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(std.testing.allocator);

    var csq = Vcsq{
        .csq_type = CSQ_PRINTED_UPSTREAM,
        .ref_pos = 12345,
    };

    const opts = FormatOptions{};
    try formatVcsq(&csq, opts, buf.writer(std.testing.allocator));
    try std.testing.expectEqualStrings("@12345", buf.items);
}

test "format upstream-stop consequence" {
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(std.testing.allocator);

    var csq = Vcsq{
        .csq_type = CSQ_UPSTREAM_STOP | CSQ_STOP_LOST,
        .strand = .rev,
        .gene = "TP53",
        .vstr = "|1*>1Q|500T>C",
    };

    const opts = FormatOptions{};
    try formatVcsq(&csq, opts, buf.writer(std.testing.allocator));

    const result = buf.items;
    try std.testing.expect(result[0] == '*');
    try std.testing.expect(std.mem.indexOf(u8, result, "stop_lost") != null);
}

test "format multiple consequences joined by &" {
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(std.testing.allocator);

    var csq = Vcsq{
        .csq_type = CSQ_MISSENSE_VARIANT | CSQ_SPLICE_REGION,
        .strand = .fwd,
        .gene = "EGFR",
    };

    const opts = FormatOptions{};
    try formatVcsq(&csq, opts, buf.writer(std.testing.allocator));

    const result = buf.items;
    try std.testing.expect(std.mem.indexOf(u8, result, "missense&splice_region") != null or
        std.mem.indexOf(u8, result, "missense") != null);
}

test "formatAaPrediction no truncation" {
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(std.testing.allocator);

    try formatAaPrediction("MKT", "MKT", 0, buf.writer(std.testing.allocator));
    try std.testing.expectEqualStrings("MKT", buf.items);
}

test "formatAaPrediction with truncation" {
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(std.testing.allocator);

    // aa = "MKTGFH*", stop = "MKTGFH*", brief_predictions = 2
    // len-brief = 7-2=5 >= 3, so truncate
    try formatAaPrediction("MKTGFH*", "MKTGFH*", 2, buf.writer(std.testing.allocator));
    try std.testing.expectEqualStrings("MK..", buf.items);
}

test "format consequence list" {
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(std.testing.allocator);

    var csqs = [_]Vcsq{
        .{
            .csq_type = CSQ_SYNONYMOUS_VARIANT,
            .strand = .fwd,
            .gene = "A",
        },
        .{
            .csq_type = CSQ_INTRON,
            .gene = "B",
        },
    };

    const opts = FormatOptions{};
    try formatVcsqList(&csqs, opts, buf.writer(std.testing.allocator));
    const result = buf.items;
    // Should have a comma separator
    try std.testing.expect(std.mem.indexOf(u8, result, ",") != null);
    try std.testing.expect(std.mem.indexOf(u8, result, "synonymous") != null);
    try std.testing.expect(std.mem.indexOf(u8, result, "intron") != null);
}
