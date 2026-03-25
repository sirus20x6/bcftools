const std = @import("std");

// ---------------------------------------------------------------------------
// Splice constants
// ---------------------------------------------------------------------------

pub const n_splice_donor: u32 = 2;
pub const n_splice_region_exon: u32 = 3;
pub const n_splice_region_intron: u32 = 8;

// ---------------------------------------------------------------------------
// Small enums
// ---------------------------------------------------------------------------

pub const Strand = enum(u2) {
    reverse = 0,
    forward = 1,
    unknown = 2,
};

pub const Trim = enum(u2) {
    none = 0,
    prime5 = 1,
    prime3 = 2,
};

pub const CdsPhase = enum(u2) {
    phase0 = 0,
    phase1 = 1,
    phase2 = 2,
    unknown = 3,
};

pub const UtrType = enum {
    prime3,
    prime5,
};

// ---------------------------------------------------------------------------
// Biotype
// ---------------------------------------------------------------------------

const coding_bit: u5 = 6;

pub const Biotype = enum(u32) {
    // Non-coding biotypes (1..48)
    mt_rRNA = 1,
    mt_tRNA = 2,
    lincRNA = 3,
    miRNA = 4,
    misc_RNA = 5,
    rRNA = 6,
    snRNA = 7,
    snoRNA = 8,
    processed_transcript = 9,
    antisense = 10,
    macro_lncRNA = 11,
    ribozyme = 12,
    sRNA = 13,
    scRNA = 14,
    scaRNA = 15,
    sense_intronic = 16,
    sense_overlapping = 17,
    pseudogene = 18,
    processed_pseudogene = 19,
    artifact = 20,
    IG_pseudogene = 21,
    IG_C_pseudogene = 22,
    IG_J_pseudogene = 23,
    IG_V_pseudogene = 24,
    TR_V_pseudogene = 25,
    TR_J_pseudogene = 26,
    mt_tRNA_pseudogene = 27,
    misc_RNA_pseudogene = 28,
    miRNA_pseudogene = 29,
    ribozyme_pseudogene = 30,
    retained_intron = 31,
    retrotransposed = 32,
    tRNA_pseudogene = 33,
    transcribed_processed_pseudogene = 34,
    transcribed_unprocessed_pseudogene = 35,
    transcribed_unitary_pseudogene = 36,
    translated_unprocessed_pseudogene = 37,
    translated_processed_pseudogene = 38,
    known_ncrna = 39,
    unitary_pseudogene = 40,
    unprocessed_pseudogene = 41,
    lrg_gene = 42,
    three_prime_overlapping_ncRNA = 43,
    disrupted_domain = 44,
    vaultRNA = 45,
    bidirectional_promoter_lncRNA = 46,
    ambiguous_orf = 47,
    lncRNA = 48,

    // Coding biotypes (n | (1 << 6))
    protein_coding = 1 | (1 << coding_bit), // 65
    polymorphic_pseudogene = 2 | (1 << coding_bit),
    IG_C = 3 | (1 << coding_bit),
    IG_D = 4 | (1 << coding_bit),
    IG_J = 5 | (1 << coding_bit),
    IG_LV = 6 | (1 << coding_bit),
    IG_V = 7 | (1 << coding_bit),
    TR_C = 8 | (1 << coding_bit),
    TR_D = 9 | (1 << coding_bit),
    TR_J = 10 | (1 << coding_bit),
    TR_V = 11 | (1 << coding_bit),
    NMD = 12 | (1 << coding_bit),
    non_stop_decay = 13 | (1 << coding_bit),

    // Special types ((1 << 7) + n)
    CDS = (1 << (coding_bit + 1)) + 1, // 129
    exon = (1 << (coding_bit + 1)) + 2,
    UTR3 = (1 << (coding_bit + 1)) + 3,
    UTR5 = (1 << (coding_bit + 1)) + 4,

    /// Returns true if this biotype represents a coding gene/transcript.
    pub fn isCoding(self: Biotype) bool {
        return (@intFromEnum(self) & (1 << coding_bit)) != 0;
    }

    /// Returns true for special structural types (CDS, exon, UTR3, UTR5).
    pub fn isSpecial(self: Biotype) bool {
        return @intFromEnum(self) >= (1 << (coding_bit + 1));
    }

    /// Returns the canonical GFF biotype string for this value.
    pub fn toGffString(self: Biotype) []const u8 {
        return switch (self) {
            .mt_rRNA => "Mt_rRNA",
            .mt_tRNA => "Mt_tRNA",
            .lincRNA => "lincRNA",
            .miRNA => "miRNA",
            .misc_RNA => "misc_RNA",
            .rRNA => "rRNA",
            .snRNA => "snRNA",
            .snoRNA => "snoRNA",
            .processed_transcript => "processed_transcript",
            .antisense => "antisense",
            .macro_lncRNA => "macro_lncRNA",
            .ribozyme => "ribozyme",
            .sRNA => "sRNA",
            .scRNA => "scRNA",
            .scaRNA => "scaRNA",
            .sense_intronic => "sense_intronic",
            .sense_overlapping => "sense_overlapping",
            .pseudogene => "pseudogene",
            .processed_pseudogene => "processed_pseudogene",
            .artifact => "artifact",
            .IG_pseudogene => "IG_pseudogene",
            .IG_C_pseudogene => "IG_C_pseudogene",
            .IG_J_pseudogene => "IG_J_pseudogene",
            .IG_V_pseudogene => "IG_V_pseudogene",
            .TR_V_pseudogene => "TR_V_pseudogene",
            .TR_J_pseudogene => "TR_J_pseudogene",
            .mt_tRNA_pseudogene => "Mt_tRNA_pseudogene",
            .misc_RNA_pseudogene => "misc_RNA_pseudogene",
            .miRNA_pseudogene => "miRNA_pseudogene",
            .ribozyme_pseudogene => "ribozyme",
            .retained_intron => "retained_intron",
            .retrotransposed => "retrotransposed",
            .tRNA_pseudogene => "tRNA_pseudogene",
            .transcribed_processed_pseudogene => "transcribed_processed_pseudogene",
            .transcribed_unprocessed_pseudogene => "transcribed_unprocessed_pseudogene",
            .transcribed_unitary_pseudogene => "transcribed_unitary_pseudogene",
            .translated_unprocessed_pseudogene => "translated_unprocessed_pseudogene",
            .translated_processed_pseudogene => "translated_processed_pseudogene",
            .known_ncrna => "known_ncrna",
            .unitary_pseudogene => "unitary_pseudogene",
            .unprocessed_pseudogene => "unprocessed_pseudogene",
            .lrg_gene => "LRG_gene",
            .three_prime_overlapping_ncRNA => "3prime_overlapping_ncRNA",
            .disrupted_domain => "disrupted_domain",
            .vaultRNA => "vaultRNA",
            .bidirectional_promoter_lncRNA => "bidirectional_promoter_lncRNA",
            .ambiguous_orf => "ambiguous_orf",
            .lncRNA => "lncRNA",
            .protein_coding => "protein_coding",
            .polymorphic_pseudogene => "polymorphic_pseudogene",
            .IG_C => "IG_C_gene",
            .IG_D => "IG_D_gene",
            .IG_J => "IG_J_gene",
            .IG_LV => "IG_LV_gene",
            .IG_V => "IG_V_gene",
            .TR_C => "TR_C_gene",
            .TR_D => "TR_D_gene",
            .TR_J => "TR_J_gene",
            .TR_V => "TR_V_gene",
            .NMD => "nonsense_mediated_decay",
            .non_stop_decay => "non_stop_decay",
            .CDS => "CDS",
            .exon => "exon",
            .UTR3 => "three_prime_UTR",
            .UTR5 => "five_prime_UTR",
        };
    }
};

/// Comptime string map for parsing GFF biotype strings into `Biotype` values.
pub const biotype_map = std.StaticStringMap(Biotype).initComptime(.{
    // Non-coding
    .{ "Mt_rRNA", .mt_rRNA },
    .{ "Mt_tRNA", .mt_tRNA },
    .{ "lincRNA", .lincRNA },
    .{ "miRNA", .miRNA },
    .{ "misc_RNA", .misc_RNA },
    .{ "rRNA", .rRNA },
    .{ "snRNA", .snRNA },
    .{ "snoRNA", .snoRNA },
    .{ "processed_transcript", .processed_transcript },
    .{ "antisense", .antisense },
    .{ "macro_lncRNA", .macro_lncRNA },
    .{ "ribozyme", .ribozyme },
    .{ "sRNA", .sRNA },
    .{ "scRNA", .scRNA },
    .{ "scaRNA", .scaRNA },
    .{ "sense_intronic", .sense_intronic },
    .{ "sense_overlapping", .sense_overlapping },
    .{ "pseudogene", .pseudogene },
    .{ "processed_pseudogene", .processed_pseudogene },
    .{ "artifact", .artifact },
    .{ "IG_pseudogene", .IG_pseudogene },
    .{ "IG_C_pseudogene", .IG_C_pseudogene },
    .{ "IG_J_pseudogene", .IG_J_pseudogene },
    .{ "IG_V_pseudogene", .IG_V_pseudogene },
    .{ "TR_V_pseudogene", .TR_V_pseudogene },
    .{ "TR_J_pseudogene", .TR_J_pseudogene },
    .{ "Mt_tRNA_pseudogene", .mt_tRNA_pseudogene },
    .{ "misc_RNA_pseudogene", .misc_RNA_pseudogene },
    .{ "miRNA_pseudogene", .miRNA_pseudogene },
    .{ "retained_intron", .retained_intron },
    .{ "retrotransposed", .retrotransposed },
    .{ "tRNA_pseudogene", .tRNA_pseudogene },
    .{ "transcribed_processed_pseudogene", .transcribed_processed_pseudogene },
    .{ "transcribed_unprocessed_pseudogene", .transcribed_unprocessed_pseudogene },
    .{ "transcribed_unitary_pseudogene", .transcribed_unitary_pseudogene },
    .{ "translated_unprocessed_pseudogene", .translated_unprocessed_pseudogene },
    .{ "translated_processed_pseudogene", .translated_processed_pseudogene },
    .{ "known_ncrna", .known_ncrna },
    .{ "unitary_pseudogene", .unitary_pseudogene },
    .{ "unprocessed_pseudogene", .unprocessed_pseudogene },
    .{ "LRG_gene", .lrg_gene },
    .{ "3prime_overlapping_ncRNA", .three_prime_overlapping_ncRNA },
    .{ "disrupted_domain", .disrupted_domain },
    .{ "vaultRNA", .vaultRNA },
    .{ "bidirectional_promoter_lncRNA", .bidirectional_promoter_lncRNA },
    .{ "ambiguous_orf", .ambiguous_orf },
    .{ "lncRNA", .lncRNA },
    // Coding
    .{ "protein_coding", .protein_coding },
    .{ "polymorphic_pseudogene", .polymorphic_pseudogene },
    .{ "IG_C_gene", .IG_C },
    .{ "IG_D_gene", .IG_D },
    .{ "IG_J_gene", .IG_J },
    .{ "IG_LV_gene", .IG_LV },
    .{ "IG_V_gene", .IG_V },
    .{ "TR_C_gene", .TR_C },
    .{ "TR_D_gene", .TR_D },
    .{ "TR_J_gene", .TR_J },
    .{ "TR_V_gene", .TR_V },
    .{ "nonsense_mediated_decay", .NMD },
    .{ "NMD", .NMD },
    .{ "non_stop_decay", .non_stop_decay },
    // Special
    .{ "CDS", .CDS },
    .{ "exon", .exon },
    .{ "three_prime_UTR", .UTR3 },
    .{ "five_prime_UTR", .UTR5 },
});

// ---------------------------------------------------------------------------
// Genomic feature structures
// ---------------------------------------------------------------------------

/// A single CDS (coding sequence) segment within a transcript.
pub const CdsEntry = struct {
    tr: *Transcript,
    beg: u32,
    pos: u32,
    len: u32,
    icds: u30,
    phase: CdsPhase,
};

/// A gene record.
pub const Gene = struct {
    name: ?[*:0]u8,
    iseq: u32,
    id: u32,
    beg: u32,
    end: u32,
    strand: Strand,
    used: bool,
};

/// An exon within a transcript.
pub const Exon = struct {
    beg: u32,
    end: u32,
    tr: *Transcript,
};

/// A UTR (untranslated region) within a transcript.
pub const Utr = struct {
    which: UtrType,
    beg: u32,
    end: u32,
    tr: *Transcript,
};

/// A transcript record, replacing the C `gf_tscript_t` struct.
pub const Transcript = struct {
    id: u32,
    beg: u32,
    end: u32,
    strand: Strand,
    used: bool,
    cds: std.ArrayList(*CdsEntry),
    trim: Trim,
    biotype: Biotype,
    gene: ?*Gene,
    aux: ?*anyopaque,
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator) Transcript {
        return .{
            .id = 0,
            .beg = 0,
            .end = 0,
            .strand = .unknown,
            .used = false,
            .cds = .empty,
            .trim = .none,
            .biotype = .processed_transcript,
            .gene = null,
            .aux = null,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Transcript) void {
        self.cds.deinit(self.allocator);
    }
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "biotype_map lookups" {
    const testing = std.testing;

    try testing.expectEqual(Biotype.protein_coding, biotype_map.get("protein_coding").?);
    try testing.expectEqual(Biotype.lncRNA, biotype_map.get("lncRNA").?);
    try testing.expectEqual(Biotype.miRNA, biotype_map.get("miRNA").?);
    try testing.expectEqual(Biotype.pseudogene, biotype_map.get("pseudogene").?);
    try testing.expectEqual(Biotype.NMD, biotype_map.get("nonsense_mediated_decay").?);
    try testing.expectEqual(Biotype.NMD, biotype_map.get("NMD").?);
    try testing.expectEqual(Biotype.CDS, biotype_map.get("CDS").?);
    try testing.expectEqual(Biotype.UTR3, biotype_map.get("three_prime_UTR").?);
    try testing.expectEqual(Biotype.UTR5, biotype_map.get("five_prime_UTR").?);
    try testing.expectEqual(Biotype.IG_C, biotype_map.get("IG_C_gene").?);
    try testing.expectEqual(Biotype.TR_V, biotype_map.get("TR_V_gene").?);
    try testing.expectEqual(Biotype.lrg_gene, biotype_map.get("LRG_gene").?);
    try testing.expectEqual(Biotype.three_prime_overlapping_ncRNA, biotype_map.get("3prime_overlapping_ncRNA").?);

    // Unknown string returns null
    try testing.expect(biotype_map.get("not_a_biotype") == null);
}

test "biotype numeric values match C defines" {
    const testing = std.testing;

    try testing.expectEqual(@as(u32, 1), @intFromEnum(Biotype.mt_rRNA));
    try testing.expectEqual(@as(u32, 48), @intFromEnum(Biotype.lncRNA));
    try testing.expectEqual(@as(u32, 65), @intFromEnum(Biotype.protein_coding));
    try testing.expectEqual(@as(u32, 129), @intFromEnum(Biotype.CDS));
    try testing.expectEqual(@as(u32, 130), @intFromEnum(Biotype.exon));
    try testing.expectEqual(@as(u32, 131), @intFromEnum(Biotype.UTR3));
    try testing.expectEqual(@as(u32, 132), @intFromEnum(Biotype.UTR5));
}

test "isCoding" {
    const testing = std.testing;

    try testing.expect(Biotype.protein_coding.isCoding());
    try testing.expect(Biotype.IG_C.isCoding());
    try testing.expect(Biotype.NMD.isCoding());
    try testing.expect(!Biotype.lncRNA.isCoding());
    try testing.expect(!Biotype.pseudogene.isCoding());
    try testing.expect(!Biotype.miRNA.isCoding());

    // Special types (CDS/exon/UTR) have bit 7 set but NOT bit 6,
    // so they are NOT considered coding by GF_is_coding().
    try testing.expect(!Biotype.CDS.isCoding());
}

test "isSpecial" {
    const testing = std.testing;

    try testing.expect(Biotype.CDS.isSpecial());
    try testing.expect(Biotype.exon.isSpecial());
    try testing.expect(Biotype.UTR3.isSpecial());
    try testing.expect(Biotype.UTR5.isSpecial());
    try testing.expect(!Biotype.protein_coding.isSpecial());
    try testing.expect(!Biotype.lncRNA.isSpecial());
}

test "toGffString round-trip" {
    const testing = std.testing;

    // For every biotype that has a canonical GFF string, the map should
    // resolve that string back to the same biotype.
    const cases = [_]Biotype{
        .protein_coding, .lncRNA,  .miRNA,     .CDS,
        .exon,           .UTR3,    .UTR5,      .NMD,
        .pseudogene,     .rRNA,    .IG_C,      .TR_V,
        .vaultRNA,       .snoRNA,  .lrg_gene,
    };
    for (cases[0..]) |bt| {
        const s = bt.toGffString();
        const resolved = biotype_map.get(s);
        try testing.expect(resolved != null);
        try testing.expectEqual(bt, resolved.?);
    }
}
