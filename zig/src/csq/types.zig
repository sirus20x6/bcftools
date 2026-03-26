const std = @import("std");
const Allocator = std.mem.Allocator;
const ArrayList = std.ArrayList;

// ---------------------------------------------------------------------------
// Opaque placeholder for htslib's bcf1_t -- will be replaced by a real
// wrapper once we have Zig bindings for htslib.
// ---------------------------------------------------------------------------
pub const VcfRecord = opaque {};

// Opaque placeholder for gf_tscript_t (from the GFF module).
pub const GfTscript = opaque {};

// ---------------------------------------------------------------------------
// CsqType -- packed bitmask of consequence flags (mirrors the C #define CSQ_*)
// ---------------------------------------------------------------------------
pub const CsqType = packed struct(u32) {
    printed_upstream: bool = false, // bit 0
    synonymous_variant: bool = false, // bit 1
    missense_variant: bool = false, // bit 2
    stop_lost: bool = false, // bit 3
    stop_gained: bool = false, // bit 4
    inframe_deletion: bool = false, // bit 5
    inframe_insertion: bool = false, // bit 6
    frameshift_variant: bool = false, // bit 7
    splice_acceptor: bool = false, // bit 8
    splice_donor: bool = false, // bit 9
    start_lost: bool = false, // bit 10
    splice_region: bool = false, // bit 11
    stop_retained: bool = false, // bit 12
    utr5: bool = false, // bit 13
    utr3: bool = false, // bit 14
    non_coding: bool = false, // bit 15
    intron: bool = false, // bit 16
    _reserved17: bool = false, // bit 17 (unused)
    inframe_altering: bool = false, // bit 18
    upstream_stop: bool = false, // bit 19
    incomplete_cds: bool = false, // bit 20
    coding_sequence: bool = false, // bit 21
    elongation: bool = false, // bit 22
    truncation: bool = false, // bit 23
    start_retained: bool = false, // bit 24
    _pad: u7 = 0, // bits 25-31

    // The consequence string table, indexed by bit position.
    // Entry 0 is reserved for printed_upstream which has no display string.
    pub const strings = [_]?[]const u8{
        null, // 0  printed_upstream
        "synonymous", // 1
        "missense", // 2
        "stop_lost", // 3
        "stop_gained", // 4
        "inframe_deletion", // 5
        "inframe_insertion", // 6
        "frameshift", // 7
        "splice_acceptor", // 8
        "splice_donor", // 9
        "start_lost", // 10
        "splice_region", // 11
        "stop_retained", // 12
        "5_prime_utr", // 13
        "3_prime_utr", // 14
        "non_coding", // 15
        "intron", // 16
        "intergenic", // 17
        "inframe_altering", // 18
        null, // 19  upstream_stop (prefix, no own string)
        null, // 20  incomplete_cds
        "coding_sequence", // 21
        "feature_elongation", // 22
        "feature_truncation", // 23
        "start_retained", // 24
    };

    const compound_mask: u32 = (1 << 1) | // synonymous_variant
        (1 << 2) | // missense_variant
        (1 << 3) | // stop_lost
        (1 << 4) | // stop_gained
        (1 << 5) | // inframe_deletion
        (1 << 6) | // inframe_insertion
        (1 << 7) | // frameshift_variant
        (1 << 10) | // start_lost
        (1 << 12) | // stop_retained
        (1 << 18) | // inframe_altering
        (1 << 19) | // upstream_stop
        (1 << 20) | // incomplete_cds
        (1 << 22) | // elongation
        (1 << 23) | // truncation
        (1 << 24); // start_retained

    const start_stop_mask: u32 = (1 << 3) | // stop_lost
        (1 << 4) | // stop_gained
        (1 << 10) | // start_lost
        (1 << 12) | // stop_retained
        (1 << 24); // start_retained

    const splice_mask: u32 = (1 << 8) | // splice_acceptor
        (1 << 9) | // splice_donor
        (1 << 11); // splice_region

    const prn_tscript_mask: u32 = ~(@as(u32, (1 << 16) | (1 << 15))); // ~(INTRON|NON_CODING)
    const prn_biotype_mask: u32 = 1 << 15; // NON_CODING

    /// Returns true when any compound-relevant consequence bit is set.
    pub fn isCompound(self: CsqType) bool {
        return (self.toInt() & compound_mask) != 0;
    }

    /// Returns true when any start/stop consequence bit is set.
    pub fn isStartStop(self: CsqType) bool {
        return (self.toInt() & start_stop_mask) != 0;
    }

    /// Mirrors CSQ_PRN_STRAND: true when the consequence is compound AND
    /// none of splice_acceptor, splice_donor, splice_region, elongation, or
    /// truncation are set.
    pub fn prnStrand(self: CsqType) bool {
        const raw = self.toInt();
        const splice_elong_trunc = splice_mask | (1 << 22) | (1 << 23);
        return (raw & compound_mask) != 0 and (raw & splice_elong_trunc) == 0;
    }

    /// Whether this consequence should include the transcript id.
    pub fn prnTscript(self: CsqType) bool {
        return (self.toInt() & prn_tscript_mask) != 0;
    }

    /// Whether the biotype should be printed.
    pub fn prnBiotype(self: CsqType) bool {
        return (self.toInt() & prn_biotype_mask) != 0;
    }

    /// Reinterpret the packed struct as a plain u32.
    pub fn toInt(self: CsqType) u32 {
        return @bitCast(self);
    }

    /// Construct a CsqType from a raw u32 bitmask.
    pub fn fromInt(raw: u32) CsqType {
        return @bitCast(raw);
    }

    /// Write all set consequence strings (separated by `&`) into `writer`.
    /// This is the Zig equivalent of iterating csq_strings[] for set bits.
    pub fn format(self: CsqType, writer: anytype) !void {
        const raw = self.toInt();
        var first = true;
        for (0..strings.len) |i| {
            if (raw & (@as(u32, 1) << @intCast(i)) != 0) {
                if (strings[i]) |s| {
                    if (!first) try writer.writeByte('&');
                    try writer.writeAll(s);
                    first = false;
                }
            }
        }
    }

    /// Convenience: return the formatted string as an owned slice.
    pub fn formatAlloc(self: CsqType, allocator: Allocator) ![]u8 {
        var buf = ArrayList(u8).empty;
        defer buf.deinit(allocator);
        try self.format(buf.writer(allocator));
        return buf.toOwnedSlice(allocator);
    }
};

// ---------------------------------------------------------------------------
// Phase -- how to handle phasing
// ---------------------------------------------------------------------------
pub const Phase = enum(u3) {
    require = 0,
    merge = 1,
    as_is = 2,
    skip = 3,
    non_ref = 4,
    drop_gt = 5,
};

// ---------------------------------------------------------------------------
// FilterMode
// ---------------------------------------------------------------------------
pub const FilterMode = enum(u2) {
    none = 0,
    include = 1,
    exclude = 2,
};

/// Number of padding bases added to each end of reference sequences.
pub const n_ref_pad: u32 = 10;

// ---------------------------------------------------------------------------
// Vcsq -- consequence annotation for a single allele (was vcsq_t)
// ---------------------------------------------------------------------------
pub const Vcsq = struct {
    strand: bool = false,
    csq_type: CsqType = .{},
    trid: u32 = 0,
    vcf_ial: u32 = 0,
    biotype: u32 = 0,
    gene: ?[]const u8 = null,
    ref: ?*VcfRecord = null,
    /// Variant string, e.g. "5TY>5I|121ACG>A+124TA>T".  Replaces kstring_t.
    vstr: ArrayList(u8) = .empty,

    pub fn deinit(self: *Vcsq, allocator: Allocator) void {
        self.vstr.deinit(allocator);
    }
};

// ---------------------------------------------------------------------------
// Vrec -- a VCF record held in the processing buffer (was vrec_t)
// ---------------------------------------------------------------------------
pub const Vrec = struct {
    line: ?*VcfRecord = null,
    /// Per-sample consequence bitmask (first/second haplotype interleaved).
    fmt_bm: ?[]u32 = null,
    /// Number of u32s per sample in fmt_bm (max 15, stored in 4 bits in C).
    nfmt: u4 = 0,
    /// Consequence annotations attached to this record.
    vcsq: ArrayList(Vcsq) = .empty,

    pub fn deinit(self: *Vrec, allocator: Allocator) void {
        for (self.vcsq.items) |*v| v.deinit(allocator);
        self.vcsq.deinit(allocator);
    }
};

// ---------------------------------------------------------------------------
// Csq -- a top-level consequence tied to a haplotype (was csq_t)
// ---------------------------------------------------------------------------
pub const Csq = struct {
    pos: u32 = 0,
    vrec: ?*Vrec = null,
    idx: i32 = 0,
    type_info: Vcsq = .{},
    /// For CSQ_PRINTED_UPSTREAM: the 1-based position of the upstream reference record.
    ref_pos: ?u32 = null,

    pub fn deinit(self: *Csq, allocator: Allocator) void {
        self.type_info.deinit(allocator);
    }
};

// ---------------------------------------------------------------------------
// Vbuf -- buffer of VCF records that share the same genomic position
// ---------------------------------------------------------------------------
pub const Vbuf = struct {
    vrec: ArrayList(*Vrec) = .empty,
    keep_until: u32 = 0,

    pub fn deinit(self: *Vbuf, allocator: Allocator) void {
        self.vrec.deinit(allocator);
    }
};

// ---------------------------------------------------------------------------
// HapNode -- node in the per-transcript haplotype tree (was hap_node_t)
//
// The C version uses a `type` field (HAP_ROOT / HAP_CDS / HAP_SSS) plus a
// `seq` pointer that is only valid for HAP_CDS.  We model this as a tagged
// union inside the struct so the invariant is enforced at the type level.
// ---------------------------------------------------------------------------
pub const HapNodeType = enum(u2) {
    root = 1,
    cds = 0,
    sss = 2,
};

pub const HapNode = struct {
    /// The payload differs by node type.  Only `cds` carries a `seq` slice.
    payload: union(HapNodeType) {
        root: void,
        cds: struct {
            /// CDS segment from the parent node to this node.
            seq: ?[]u8 = null,
        },
        sss: void,
    },

    /// Variant description "ref>alt".
    var_str: ?[]const u8 = null,
    /// This node's consequence flags.
    csq: CsqType = .{},
    /// Alt length minus ref length: <0 del, >0 ins, 0 substitution.
    dlen: i32 = 0,
    /// Variant's VCF position (0-based, inclusive).
    rbeg: u32 = 0,
    /// Variant's reference length; alt length = rlen + dlen.
    rlen: i32 = 0,
    /// Position on the spliced reference transcript (0-based, exclusive of N_REF_PAD).
    sbeg: u32 = 0,
    /// Index of the exon this variant overlaps.
    icds: u32 = 0,

    /// Children in the haplotype tree.
    children: ArrayList(*HapNode) = .empty,
    /// Previous coding node.
    prev: ?*HapNode = null,

    /// The VCF record that created this node (type-erased pointer for identity comparison).
    rec: ?*const anyopaque = null,
    /// Current VCF record during traversal (type-erased pointer for identity comparison).
    cur_rec: ?*const anyopaque = null,
    /// Which VCF ALT allele generated this node.
    vcf_ial: i32 = 0,
    /// Number of haplotypes ending at this node.
    nend: u32 = 0,

    /// Mapping from allele index to the currently active child.
    cur_child: ArrayList(i32) = .empty,
    /// List of haplotype consequences, broken by position.
    csq_list: ArrayList(Csq) = .empty,

    pub fn init(node_type: HapNodeType) HapNode {
        return .{
            .payload = switch (node_type) {
                .root => .{ .root = {} },
                .cds => .{ .cds = .{} },
                .sss => .{ .sss = {} },
            },
        };
    }

    pub fn deinit(self: *HapNode, allocator: Allocator) void {
        for (self.csq_list.items) |*c| c.deinit(allocator);
        self.csq_list.deinit(allocator);
        self.cur_child.deinit(allocator);
        self.children.deinit(allocator);
    }

    /// Convenience: get the node type from the tagged union.
    pub fn nodeType(self: HapNode) HapNodeType {
        return self.payload;
    }
};

// ---------------------------------------------------------------------------
// Tscript -- per-transcript haplotype context (was tscript_t)
// ---------------------------------------------------------------------------
pub const Tscript = struct {
    /// Reference sequence, padded with n_ref_pad bases on both ends.
    ref_seq: ?[]u8 = null,
    /// Spliced reference sequence, padded with n_ref_pad bases on both ends.
    sref: ?[]u8 = null,
    /// Root of the haplotype tree.
    root: ?*HapNode = null,
    /// Pointers to haplotype leaf nodes (two per sample).
    hap: ArrayList(*HapNode) = .empty,
    /// Length of sref including 2*n_ref_pad.
    nsref: i32 = 0,

    pub fn deinit(self: *Tscript, allocator: Allocator) void {
        self.hap.deinit(allocator);
    }
};

// ---------------------------------------------------------------------------
// Hstack -- single frame on the haplotype-tree traversal stack
// ---------------------------------------------------------------------------
pub const Hstack = struct {
    node: ?*HapNode = null,
    ichild: i32 = 0,
    dlen: i32 = 0,
    slen: usize = 0,
};

// ---------------------------------------------------------------------------
// Hap -- full traversal context for a transcript's haplotype tree
// ---------------------------------------------------------------------------
pub const Hap = struct {
    stack: ArrayList(Hstack) = .empty,
    tr: ?*GfTscript = null,
    /// Spliced haplotype sequence (ref strand).
    sseq: ArrayList(u8) = .empty,
    /// Variable part of translated haplotype transcript (coding strand).
    tseq: ArrayList(u8) = .empty,
    /// Variable part of translated reference transcript (coding strand).
    tref: ArrayList(u8) = .empty,
    /// Stop/start codons found in tseq.
    tseq_stop: ArrayList(u8) = .empty,
    /// Stop/start codons found in tref.
    tref_stop: ArrayList(u8) = .empty,
    /// Stack's sbeg, for when the first node's type is sss.
    sbeg: u32 = 0,
    upstream_stop: i32 = 0,

    pub fn deinit(self: *Hap, allocator: Allocator) void {
        self.stack.deinit(allocator);
        self.sseq.deinit(allocator);
        self.tseq.deinit(allocator);
        self.tref.deinit(allocator);
        self.tseq_stop.deinit(allocator);
        self.tref_stop.deinit(allocator);
    }
};

// ===========================================================================
// Tests
// ===========================================================================
test "CsqType bit layout matches C defines" {
    // Each named flag must map to the expected bit position.
    const t = comptime CsqType{ .printed_upstream = true };
    try std.testing.expectEqual(@as(u32, 1 << 0), t.toInt());

    try std.testing.expectEqual(@as(u32, 1 << 1), (CsqType{ .synonymous_variant = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 2), (CsqType{ .missense_variant = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 3), (CsqType{ .stop_lost = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 4), (CsqType{ .stop_gained = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 5), (CsqType{ .inframe_deletion = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 6), (CsqType{ .inframe_insertion = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 7), (CsqType{ .frameshift_variant = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 8), (CsqType{ .splice_acceptor = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 9), (CsqType{ .splice_donor = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 10), (CsqType{ .start_lost = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 11), (CsqType{ .splice_region = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 12), (CsqType{ .stop_retained = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 13), (CsqType{ .utr5 = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 14), (CsqType{ .utr3 = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 15), (CsqType{ .non_coding = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 16), (CsqType{ .intron = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 18), (CsqType{ .inframe_altering = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 19), (CsqType{ .upstream_stop = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 20), (CsqType{ .incomplete_cds = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 21), (CsqType{ .coding_sequence = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 22), (CsqType{ .elongation = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 23), (CsqType{ .truncation = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 24), (CsqType{ .start_retained = true }).toInt());
}

test "CsqType roundtrip fromInt/toInt" {
    const raw: u32 = (1 << 2) | (1 << 7) | (1 << 18);
    const csq = CsqType.fromInt(raw);
    try std.testing.expect(csq.missense_variant);
    try std.testing.expect(csq.frameshift_variant);
    try std.testing.expect(csq.inframe_altering);
    try std.testing.expect(!csq.stop_gained);
    try std.testing.expectEqual(raw, csq.toInt());
}

test "CsqType.isCompound" {
    const syn = CsqType{ .synonymous_variant = true };
    try std.testing.expect(syn.isCompound());

    const intron = CsqType{ .intron = true };
    try std.testing.expect(!intron.isCompound());

    const splice = CsqType{ .splice_acceptor = true };
    try std.testing.expect(!splice.isCompound());
}

test "CsqType.isStartStop" {
    const sl = CsqType{ .stop_lost = true };
    try std.testing.expect(sl.isStartStop());

    const sr = CsqType{ .start_retained = true };
    try std.testing.expect(sr.isStartStop());

    const miss = CsqType{ .missense_variant = true };
    try std.testing.expect(!miss.isStartStop());
}

test "CsqType.prnStrand" {
    // A compound consequence without splice/elongation/truncation => true.
    const miss = CsqType{ .missense_variant = true };
    try std.testing.expect(miss.prnStrand());

    // Splice donor alone is not compound => false.
    const sd = CsqType{ .splice_donor = true };
    try std.testing.expect(!sd.prnStrand());

    // Compound + splice_acceptor => false (splice mask blocks it).
    const both = CsqType{ .missense_variant = true, .splice_acceptor = true };
    try std.testing.expect(!both.prnStrand());

    // Elongation alone is compound but blocked by its own mask => false.
    const elong = CsqType{ .elongation = true };
    try std.testing.expect(!elong.prnStrand());
}

test "CsqType.format produces expected strings" {
    const allocator = std.testing.allocator;

    const csq = CsqType{ .missense_variant = true, .splice_region = true };
    const s = try csq.formatAlloc(allocator);
    defer allocator.free(s);
    try std.testing.expectEqualStrings("missense&splice_region", s);
}

test "CsqType default is zero" {
    const csq = CsqType{};
    try std.testing.expectEqual(@as(u32, 0), csq.toInt());
}

test "HapNodeType values match C defines" {
    try std.testing.expectEqual(@as(u2, 0), @intFromEnum(HapNodeType.cds));
    try std.testing.expectEqual(@as(u2, 1), @intFromEnum(HapNodeType.root));
    try std.testing.expectEqual(@as(u2, 2), @intFromEnum(HapNodeType.sss));
}

test "HapNode tagged union discriminates correctly" {
    const allocator = std.testing.allocator;

    var root = HapNode.init(.root);
    defer root.deinit(allocator);
    try std.testing.expect(root.payload == .root);

    var cds_node = HapNode.init(.cds);
    defer cds_node.deinit(allocator);
    try std.testing.expect(cds_node.payload == .cds);

    var sss = HapNode.init(.sss);
    defer sss.deinit(allocator);
    try std.testing.expect(sss.payload == .sss);
}

test "Phase enum values match C defines" {
    try std.testing.expectEqual(@as(u3, 0), @intFromEnum(Phase.require));
    try std.testing.expectEqual(@as(u3, 1), @intFromEnum(Phase.merge));
    try std.testing.expectEqual(@as(u3, 2), @intFromEnum(Phase.as_is));
    try std.testing.expectEqual(@as(u3, 3), @intFromEnum(Phase.skip));
    try std.testing.expectEqual(@as(u3, 4), @intFromEnum(Phase.non_ref));
    try std.testing.expectEqual(@as(u3, 5), @intFromEnum(Phase.drop_gt));
}

test "FilterMode enum values match C defines" {
    try std.testing.expectEqual(@as(u2, 1), @intFromEnum(FilterMode.include));
    try std.testing.expectEqual(@as(u2, 2), @intFromEnum(FilterMode.exclude));
}
