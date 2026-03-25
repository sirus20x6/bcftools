// translate.zig — Genetic code tables and nucleotide/codon translation logic.
// Ported from bcftools csq.c (lines 463-589).

const std = @import("std");

/// A genetic code table mapping codons to amino acids and stop/start signals.
pub const GeneticCode = struct {
    id: i32,
    name: []const u8,
    code: [64]u8,
    stop: [64]u8,
};

/// All NCBI genetic code tables.
pub const gencode_tables = [_]GeneticCode{
    .{ .id = 0, .name = "Standard simplified", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF".*, .stop = "--------------M---------------------------------*-*-----*-------".* },
    .{ .id = 1, .name = "Standard", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF".*, .stop = "--------------M---------------M-----------------*-*-----*-----M-".* },
    .{ .id = 2, .name = "Vertebrate Mitochondrial", .code = "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF".*, .stop = "--------*-*-MMMM------------------------------M-*-*-------------".* },
    .{ .id = 3, .name = "Yeast Mitochondrial", .code = "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF".*, .stop = "------------M-M-------------------------------M-*-*-------------".* },
    .{ .id = 4, .name = "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF".*, .stop = "------------MMMM--------------M---------------M-*-*---------M-M-".* },
    .{ .id = 5, .name = "Invertebrate Mitochondrial", .code = "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF".*, .stop = "------------MMMM------------------------------M-*-*-----------M-".* },
    .{ .id = 6, .name = "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF".*, .stop = "--------------M-----------------------------------------*-------".* },
    .{ .id = 9, .name = "Echinoderm Mitochondrial; Flatworm Mitochondrial", .code = "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF".*, .stop = "--------------M-------------------------------M-*-*-------------".* },
    .{ .id = 10, .name = "Euplotid Nuclear", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF".*, .stop = "--------------M---------------------------------*-*-------------".* },
    .{ .id = 11, .name = "Bacterial, Archaeal and Plant Plastid", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF".*, .stop = "------------MMMM--------------M---------------M-*-*-----*-----M-".* },
    .{ .id = 12, .name = "Alternative Yeast Nuclear", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF".*, .stop = "--------------M---------------M-----------------*-*-----*-------".* },
    .{ .id = 13, .name = "Ascidian Mitochondrial", .code = "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF".*, .stop = "------------M-M-------------------------------M-*-*-----------M-".* },
    .{ .id = 14, .name = "Alternative Flatworm Mitochondrial", .code = "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF".*, .stop = "--------------M-----------------------------------*-------------".* },
    .{ .id = 15, .name = "Blepharisma Nuclear", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF".*, .stop = "--------------M---------------------------------*-------*-------".* },
    .{ .id = 16, .name = "Chlorophycean Mitochondrial", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF".*, .stop = "--------------M---------------------------------*-------*-------".* },
    .{ .id = 21, .name = "Trematode Mitochondrial", .code = "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF".*, .stop = "--------------M-------------------------------M-*-*-------------".* },
    .{ .id = 22, .name = "Scenedesmus obliquus Mitochondrial Code", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF".*, .stop = "--------------M---------------------------------*---*---*-------".* },
    .{ .id = 23, .name = "Thraustochytrium mitochondrial code", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF".*, .stop = "--------------MM------------------------------M-*-*-----*---*---".* },
    .{ .id = 24, .name = "Pterobranchia Mitochondrial", .code = "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF".*, .stop = "--------------M---------------M---------------M-*-*-----------M-".* },
    .{ .id = 25, .name = "Candidate Division SR1 and Gracilibacteria", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF".*, .stop = "--------------M-------------------------------M-*-*-----------M-".* },
    .{ .id = 26, .name = "Pachysolen tannophilus Nuclear Code", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLALEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF".*, .stop = "--------------M---------------M-----------------*-*-----*-------".* },
    .{ .id = 27, .name = "Karyorelict Nuclear", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSSWCWCLFLF".*, .stop = "--------------M-----------------------------------------*-------".* },
    .{ .id = 28, .name = "Condylostoma Nuclear", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSSWCWCLFLF".*, .stop = "--------------M---------------------------------*-*-----*-------".* },
    .{ .id = 29, .name = "Mesodinium Nuclear", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYYYSSSS*CWCLFLF".*, .stop = "--------------M-----------------------------------------*-------".* },
    .{ .id = 30, .name = "Peritrich Nuclear", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVEYEYSSSS*CWCLFLF".*, .stop = "--------------M-----------------------------------------*-------".* },
    .{ .id = 31, .name = "Blastocrithidia Nuclear", .code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVEYEYSSSSWCWCLFLF".*, .stop = "--------------M---------------------------------*-*-------------".* },
    .{ .id = 33, .name = "Cephalodiscidae Mitochondrial UAA-Tyr", .code = "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF".*, .stop = "--------------M---------------M---------------M---*-----------M-".* },
};

/// Nucleotide to 2-bit encoding: A=0, C=1, G=2, T=3, other=4.
pub const nt4 = blk: {
    var table: [256]u8 = [_]u8{4} ** 256;
    table['A'] = 0;
    table['a'] = 0;
    table['C'] = 1;
    table['c'] = 1;
    table['G'] = 2;
    table['g'] = 2;
    table['T'] = 3;
    table['t'] = 3;
    break :blk table;
};

/// Complement nucleotide to 2-bit encoding: A->3(T), C->2(G), G->1(C), T->0(A), other=4.
pub const cnt4 = blk: {
    var table: [256]u8 = [_]u8{4} ** 256;
    table['A'] = 3;
    table['a'] = 3;
    table['C'] = 2;
    table['c'] = 2;
    table['G'] = 1;
    table['g'] = 1;
    table['T'] = 0;
    table['t'] = 0;
    break :blk table;
};

/// Compute codon index from three 2-bit encoded nucleotides.
/// Matches the C macro: `(a << 4) | (b << 2) | c`
pub inline fn codonIndex(a: u8, b: u8, c: u8) u8 {
    return (a << 4) | (b << 2) | c;
}

/// Encode a 3-base DNA sequence to a codon index (forward strand).
/// Returns null if any nucleotide is invalid.
pub inline fn dnaIndex(seq: *const [3]u8) ?u8 {
    const a = nt4[seq[0]];
    const b = nt4[seq[1]];
    const c = nt4[seq[2]];
    if (a > 3 or b > 3 or c > 3) return null;
    return codonIndex(a, b, c);
}

/// Encode a 3-base DNA sequence to a codon index (reverse complement).
/// The complement is read in reverse order: complement(seq[2]), complement(seq[1]), complement(seq[0]).
/// Returns null if any nucleotide is invalid.
pub inline fn cdnaIndex(seq: *const [3]u8) ?u8 {
    const a = cnt4[seq[2]];
    const b = cnt4[seq[1]];
    const c = cnt4[seq[0]];
    if (a > 3 or b > 3 or c > 3) return null;
    return codonIndex(a, b, c);
}

/// Translate a forward-strand codon to its amino acid character.
/// Returns null if the codon contains an invalid nucleotide.
pub inline fn dna2aa(code: *const GeneticCode, seq: *const [3]u8) ?u8 {
    const idx = dnaIndex(seq) orelse return null;
    return code.code[idx];
}

/// Translate a reverse-complement codon to its amino acid character.
/// Returns null if the codon contains an invalid nucleotide.
pub inline fn cdna2aa(code: *const GeneticCode, seq: *const [3]u8) ?u8 {
    const idx = cdnaIndex(seq) orelse return null;
    return code.code[idx];
}

/// Look up the stop/start annotation for a forward-strand codon.
/// Returns the character from the stop table ('M' for start, '*' for stop, '-' otherwise).
/// Returns null if the codon contains an invalid nucleotide.
pub inline fn dna2stop(code: *const GeneticCode, seq: *const [3]u8) ?u8 {
    const idx = dnaIndex(seq) orelse return null;
    return code.stop[idx];
}

/// Look up the stop/start annotation for a reverse-complement codon.
/// Returns the character from the stop table ('M' for start, '*' for stop, '-' otherwise).
/// Returns null if the codon contains an invalid nucleotide.
pub inline fn cdna2stop(code: *const GeneticCode, seq: *const [3]u8) ?u8 {
    const idx = cdnaIndex(seq) orelse return null;
    return code.stop[idx];
}

/// Find a genetic code table by its NCBI ID.
/// Returns a pointer to the table entry, or null if the ID is not found.
pub fn findGeneticCode(id: i32) ?*const GeneticCode {
    for (&gencode_tables) |*entry| {
        if (entry.id == id) return entry;
    }
    return null;
}

/// Convert ncsq2 to nfmt: `(ncsq2 >> 1) + 1`.
pub inline fn ncsq2ToNfmt(ncsq2: u32) u32 {
    return (ncsq2 >> 1) + 1;
}

/// Convert icsq2 to bit index: `icsq2 >> 1`.
pub inline fn icsq2ToBit(icsq2: u32) u32 {
    return icsq2 >> 1;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "ATG translates to M (methionine/start) using standard code" {
    const code = findGeneticCode(1) orelse unreachable;
    const seq: [3]u8 = "ATG".*;
    const aa = dna2aa(code, &seq) orelse unreachable;
    try std.testing.expectEqual(@as(u8, 'M'), aa);
}

test "TAA translates to * (stop) using standard code" {
    const code = findGeneticCode(1) orelse unreachable;
    const seq: [3]u8 = "TAA".*;
    const aa = dna2aa(code, &seq) orelse unreachable;
    try std.testing.expectEqual(@as(u8, '*'), aa);
}

test "reverse complement of CAT translates to M (ATG complement)" {
    const code = findGeneticCode(1) orelse unreachable;
    // CAT reverse-complemented: complement(T)=A, complement(A)=T, complement(C)=G -> ATG -> M
    const seq: [3]u8 = "CAT".*;
    const aa = cdna2aa(code, &seq) orelse unreachable;
    try std.testing.expectEqual(@as(u8, 'M'), aa);
}

test "all 28 genetic code table IDs can be found" {
    const expected_ids = [_]i32{
        0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16,
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33,
    };
    // Verify we have exactly 27 tables (matching the C source).
    try std.testing.expectEqual(@as(usize, 27), gencode_tables.len);
    // Verify the expected_ids list matches the table count.
    try std.testing.expectEqual(gencode_tables.len, expected_ids.len);
    for (expected_ids) |id| {
        const found = findGeneticCode(id);
        try std.testing.expect(found != null);
        try std.testing.expectEqual(id, found.?.id);
    }
}

test "invalid nucleotide returns null" {
    const code = findGeneticCode(1) orelse unreachable;
    const seq: [3]u8 = "ANT".*;
    try std.testing.expect(dna2aa(code, &seq) == null);
    try std.testing.expect(cdna2aa(code, &seq) == null);
}

test "nt4 encoding" {
    try std.testing.expectEqual(@as(u8, 0), nt4['A']);
    try std.testing.expectEqual(@as(u8, 0), nt4['a']);
    try std.testing.expectEqual(@as(u8, 1), nt4['C']);
    try std.testing.expectEqual(@as(u8, 2), nt4['G']);
    try std.testing.expectEqual(@as(u8, 3), nt4['T']);
    try std.testing.expectEqual(@as(u8, 4), nt4['N']);
    try std.testing.expectEqual(@as(u8, 4), nt4[0]);
}

test "cnt4 complement encoding" {
    try std.testing.expectEqual(@as(u8, 3), cnt4['A']); // A -> T
    try std.testing.expectEqual(@as(u8, 2), cnt4['C']); // C -> G
    try std.testing.expectEqual(@as(u8, 1), cnt4['G']); // G -> C
    try std.testing.expectEqual(@as(u8, 0), cnt4['T']); // T -> A
}

test "codonIndex matches C macro" {
    // _codon_idx(a,b,c) = (a<<4 | b<<2 | c)
    try std.testing.expectEqual(@as(u8, 0), codonIndex(0, 0, 0)); // AAA
    try std.testing.expectEqual(@as(u8, 63), codonIndex(3, 3, 3)); // TTT
    try std.testing.expectEqual(@as(u8, 14), codonIndex(0, 3, 2)); // ATG = (0<<4)|(3<<2)|2 = 14
}

test "dna2stop identifies start and stop codons" {
    const code = findGeneticCode(1) orelse unreachable;
    const atg: [3]u8 = "ATG".*;
    const taa: [3]u8 = "TAA".*;
    const ggg: [3]u8 = "GGG".*;
    try std.testing.expectEqual(@as(u8, 'M'), dna2stop(code, &atg).?);
    try std.testing.expectEqual(@as(u8, '*'), dna2stop(code, &taa).?);
    try std.testing.expectEqual(@as(u8, '-'), dna2stop(code, &ggg).?);
}

test "ncsq2ToNfmt and icsq2ToBit" {
    try std.testing.expectEqual(@as(u32, 1), ncsq2ToNfmt(0));
    try std.testing.expectEqual(@as(u32, 2), ncsq2ToNfmt(2));
    try std.testing.expectEqual(@as(u32, 3), ncsq2ToNfmt(4));
    try std.testing.expectEqual(@as(u32, 0), icsq2ToBit(0));
    try std.testing.expectEqual(@as(u32, 1), icsq2ToBit(2));
    try std.testing.expectEqual(@as(u32, 2), icsq2ToBit(4));
}

test "lowercase sequences work" {
    const code = findGeneticCode(1) orelse unreachable;
    const seq: [3]u8 = "atg".*;
    try std.testing.expectEqual(@as(u8, 'M'), dna2aa(code, &seq).?);
}
