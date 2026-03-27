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
// SIMD-accelerated nucleotide encoding
// ---------------------------------------------------------------------------

/// Encode 16 nucleotides to nt4 values using SIMD.
/// Returns a vector of u8 values: 0=A, 1=C, 2=G, 3=T, 4=other.
pub fn encodeNt4x16(input: [16]u8) @Vector(16, u8) {
    const vec: @Vector(16, u8) = input;
    // Clear bit 5 (0x20) to handle both upper and lowercase
    const upper = vec & @as(@Vector(16, u8), @splat(@as(u8, 0xDF)));

    const is_a = upper == @as(@Vector(16, u8), @splat(@as(u8, 'A')));
    const is_c = upper == @as(@Vector(16, u8), @splat(@as(u8, 'C')));
    const is_g = upper == @as(@Vector(16, u8), @splat(@as(u8, 'G')));
    const is_t = upper == @as(@Vector(16, u8), @splat(@as(u8, 'T')));

    const zero: @Vector(16, u8) = @splat(0);
    const one: @Vector(16, u8) = @splat(1);
    const two: @Vector(16, u8) = @splat(2);
    const three: @Vector(16, u8) = @splat(3);
    const four: @Vector(16, u8) = @splat(4);

    var result = four; // default: invalid
    result = @select(u8, is_t, three, result);
    result = @select(u8, is_g, two, result);
    result = @select(u8, is_c, one, result);
    result = @select(u8, is_a, zero, result);

    return result;
}

/// Batch translate a DNA sequence to amino acids using SIMD-accelerated encoding.
/// Input: DNA sequence (codons read from seq, length must be a multiple of 3).
/// Output: amino acid characters written to `out`.
/// Returns number of amino acids written.
pub fn batchTranslate(code: *const GeneticCode, seq: []const u8, out: []u8) usize {
    var aa_idx: usize = 0;
    var i: usize = 0;

    // Process 48 nucleotides at a time (16 complete codons) using SIMD encoding
    while (i + 48 <= seq.len and aa_idx + 16 <= out.len) {
        // Encode 48 nucleotides in 3 SIMD passes of 16
        const enc0 = encodeNt4x16(seq[i..][0..16].*);
        const enc1 = encodeNt4x16(seq[i + 16 ..][0..16].*);
        const enc2 = encodeNt4x16(seq[i + 32 ..][0..16].*);

        // Store encoded values for codon assembly
        var encoded: [48]u8 = undefined;
        inline for (0..16) |k| {
            encoded[k] = enc0[k];
        }
        inline for (0..16) |k| {
            encoded[16 + k] = enc1[k];
        }
        inline for (0..16) |k| {
            encoded[32 + k] = enc2[k];
        }

        // Assemble 16 codons and translate
        for (0..16) |j| {
            const a = encoded[j * 3];
            const b = encoded[j * 3 + 1];
            const c = encoded[j * 3 + 2];
            if (a > 3 or b > 3 or c > 3) {
                out[aa_idx] = 'X'; // invalid nucleotide
            } else {
                out[aa_idx] = code.code[codonIndex(a, b, c)];
            }
            aa_idx += 1;
        }
        i += 48;
    }

    // Scalar fallback for remaining codons
    while (i + 3 <= seq.len and aa_idx < out.len) {
        const aa = dna2aa(code, seq[i..][0..3]);
        out[aa_idx] = aa orelse 'X';
        aa_idx += 1;
        i += 3;
    }

    return aa_idx;
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

test "encodeNt4x16 encodes uppercase nucleotides" {
    const input: [16]u8 = "ACGTACGTACGTACGT".*;
    const result = encodeNt4x16(input);
    const expected = [16]u8{ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
    for (0..16) |i| {
        try std.testing.expectEqual(expected[i], result[i]);
    }
}

test "encodeNt4x16 encodes lowercase nucleotides" {
    const input: [16]u8 = "acgtacgtacgtacgt".*;
    const result = encodeNt4x16(input);
    const expected = [16]u8{ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
    for (0..16) |i| {
        try std.testing.expectEqual(expected[i], result[i]);
    }
}

test "encodeNt4x16 marks invalid bases as 4" {
    const input: [16]u8 = "ACNXACGTACGT1234".*;
    const result = encodeNt4x16(input);
    // N=4, X=4, 1=4, 2=4, 3=4, 4=4
    try std.testing.expectEqual(@as(u8, 0), result[0]); // A
    try std.testing.expectEqual(@as(u8, 1), result[1]); // C
    try std.testing.expectEqual(@as(u8, 4), result[2]); // N
    try std.testing.expectEqual(@as(u8, 4), result[3]); // X
    try std.testing.expectEqual(@as(u8, 4), result[12]); // 1
    try std.testing.expectEqual(@as(u8, 4), result[13]); // 2
}

test "encodeNt4x16 matches scalar nt4 for all valid bases" {
    const bases = "AaCcGgTtNn.@ACGT";
    const input: [16]u8 = bases[0..16].*;
    const result = encodeNt4x16(input);
    for (0..16) |i| {
        try std.testing.expectEqual(nt4[bases[i]], result[i]);
    }
}

test "batchTranslate matches scalar dna2aa" {
    const code = findGeneticCode(1) orelse unreachable;
    // 54 bases = 18 codons (16 via SIMD + 2 scalar fallback)
    const seq = "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG";
    var out: [20]u8 = undefined;
    const n = batchTranslate(code, seq, &out);
    try std.testing.expectEqual(@as(usize, 18), n);

    // Verify each codon matches scalar translation
    var i: usize = 0;
    while (i + 3 <= seq.len) : (i += 3) {
        const expected = dna2aa(code, seq[i..][0..3]) orelse 'X';
        try std.testing.expectEqual(expected, out[i / 3]);
    }
}

test "batchTranslate with mixed case and invalid bases" {
    const code = findGeneticCode(1) orelse unreachable;
    // 48 bases = 16 codons (all via SIMD), includes lowercase
    const seq = "atgATGatgATGatgATGatgATGatgATGatgATGatgATGatgATGNNN";
    var out: [20]u8 = undefined;
    const n = batchTranslate(code, seq, &out);
    // 51 bases => 17 codons
    try std.testing.expectEqual(@as(usize, 17), n);
    // First 16 should be M (ATG -> Met)
    for (0..16) |i| {
        try std.testing.expectEqual(@as(u8, 'M'), out[i]);
    }
    // Last one has NNN -> invalid -> 'X'
    try std.testing.expectEqual(@as(u8, 'X'), out[16]);
}

test "batchTranslate with all 64 codons" {
    const code = findGeneticCode(1) orelse unreachable;
    const bases = "ACGT";
    // Generate all 64 codons = 192 bases
    var seq: [192]u8 = undefined;
    var idx: usize = 0;
    for (0..4) |a| {
        for (0..4) |b| {
            for (0..4) |c| {
                seq[idx] = bases[a];
                seq[idx + 1] = bases[b];
                seq[idx + 2] = bases[c];
                idx += 3;
            }
        }
    }

    var out: [64]u8 = undefined;
    const n = batchTranslate(code, &seq, &out);
    try std.testing.expectEqual(@as(usize, 64), n);

    // Verify against scalar
    for (0..64) |i| {
        const expected = dna2aa(code, seq[i * 3 ..][0..3]) orelse 'X';
        try std.testing.expectEqual(expected, out[i]);
    }
}
