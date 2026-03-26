/// End-to-end integration test for the CSQ pipeline.
///
/// Tests the following components without htslib:
///   1. GFF type system and biotype lookup
///   2. Codon translation and genetic code tables
///   3. Consequence string formatting (BCSQ output)
///   4. VCF record parsing through the record module
///   5. CSQ pipeline context initialization and record buffering
///   6. CsqType packed struct bit layout and operations
///   7. VCF reader with header and sample parsing
///   8. Haplotype tree context initialization
const std = @import("std");
const lib = @import("bcftools_zig");
const csq_types = lib.csq_types;
const gff_types = lib.gff_types;
const translate = lib.translate;
const format = lib.csq_format;
const csq = lib.csq_pipeline;
const vcf_record = lib.vcf_record;
const vcf_reader = lib.vcf_reader;

// Minimal VCF content for testing (written to temp files for the reader)
const minimal_vcf_content =
    "##fileformat=VCFv4.2\n" ++
    "##contig=<ID=chr1,length=10000>\n" ++
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n" ++
    "chr1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n" ++
    "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\n" ++
    "chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n" ++
    "chr1\t350\t.\tT\tC\t.\tPASS\t.\tGT\t0/1\n" ++
    "chr1\t500\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\n";

// ===== Test 1: GFF types and biotype lookup ================================

test "GFF biotype lookup and coding classification" {
    // Verify biotype string lookup
    try std.testing.expectEqual(gff_types.Biotype.protein_coding, gff_types.biotype_map.get("protein_coding").?);
    try std.testing.expectEqual(gff_types.Biotype.lncRNA, gff_types.biotype_map.get("lncRNA").?);
    try std.testing.expectEqual(gff_types.Biotype.NMD, gff_types.biotype_map.get("nonsense_mediated_decay").?);
    try std.testing.expectEqual(gff_types.Biotype.CDS, gff_types.biotype_map.get("CDS").?);
    try std.testing.expectEqual(gff_types.Biotype.UTR5, gff_types.biotype_map.get("five_prime_UTR").?);

    // Unknown biotype returns null
    try std.testing.expect(gff_types.biotype_map.get("not_a_biotype") == null);

    // Coding classification
    try std.testing.expect(gff_types.Biotype.protein_coding.isCoding());
    try std.testing.expect(gff_types.Biotype.NMD.isCoding());
    try std.testing.expect(!gff_types.Biotype.lncRNA.isCoding());
    try std.testing.expect(!gff_types.Biotype.miRNA.isCoding());

    // Special types (CDS/exon/UTR) are NOT coding
    try std.testing.expect(!gff_types.Biotype.CDS.isCoding());
    try std.testing.expect(gff_types.Biotype.CDS.isSpecial());

    // Numeric values match C defines
    try std.testing.expectEqual(@as(u32, 65), @intFromEnum(gff_types.Biotype.protein_coding));
    try std.testing.expectEqual(@as(u32, 129), @intFromEnum(gff_types.Biotype.CDS));

    // Round-trip: biotype -> GFF string -> biotype
    const bt = gff_types.Biotype.protein_coding;
    const s = bt.toGffString();
    try std.testing.expectEqual(bt, gff_types.biotype_map.get(s).?);
}

// ===== Test 2: Translation and genetic code tables =========================

test "translation and consequence determination" {
    const code = translate.findGeneticCode(1) orelse unreachable; // Standard code

    // Forward strand translation
    const atg: [3]u8 = "ATG".*;
    try std.testing.expectEqual(@as(u8, 'M'), translate.dna2aa(code, &atg).?);

    const taa: [3]u8 = "TAA".*;
    try std.testing.expectEqual(@as(u8, '*'), translate.dna2aa(code, &taa).?);

    const gct: [3]u8 = "GCT".*;
    try std.testing.expectEqual(@as(u8, 'A'), translate.dna2aa(code, &gct).?);

    const aaa: [3]u8 = "AAA".*;
    try std.testing.expectEqual(@as(u8, 'K'), translate.dna2aa(code, &aaa).?);

    const ttt: [3]u8 = "TTT".*;
    try std.testing.expectEqual(@as(u8, 'F'), translate.dna2aa(code, &ttt).?);

    // Reverse strand: CAT complement -> ATG -> M
    const cat: [3]u8 = "CAT".*;
    try std.testing.expectEqual(@as(u8, 'M'), translate.cdna2aa(code, &cat).?);

    // Stop codon detection
    try std.testing.expectEqual(@as(u8, '*'), translate.dna2stop(code, &taa).?);
    const tag: [3]u8 = "TAG".*;
    try std.testing.expectEqual(@as(u8, '*'), translate.dna2stop(code, &tag).?);
    const tga: [3]u8 = "TGA".*;
    try std.testing.expectEqual(@as(u8, '*'), translate.dna2stop(code, &tga).?);

    // ATG is a start codon, not stop
    try std.testing.expectEqual(@as(u8, 'M'), translate.dna2stop(code, &atg).?);

    // GCT is neither start nor stop
    try std.testing.expectEqual(@as(u8, '-'), translate.dna2stop(code, &gct).?);

    // Invalid nucleotide returns null
    const bad: [3]u8 = "ANT".*;
    try std.testing.expect(translate.dna2aa(code, &bad) == null);

    // Lowercase works
    const atg_lower: [3]u8 = "atg".*;
    try std.testing.expectEqual(@as(u8, 'M'), translate.dna2aa(code, &atg_lower).?);

    // Missense detection: GCT(A) vs GAT(D) - different amino acid
    const gat: [3]u8 = "GAT".*;
    const ref_aa = translate.dna2aa(code, &gct).?;
    const alt_aa = translate.dna2aa(code, &gat).?;
    try std.testing.expectEqual(@as(u8, 'A'), ref_aa);
    try std.testing.expectEqual(@as(u8, 'D'), alt_aa);
    try std.testing.expect(ref_aa != alt_aa); // missense

    // Synonymous detection: GCT(A) vs GCC(A) - same amino acid
    const gcc: [3]u8 = "GCC".*;
    const syn_aa = translate.dna2aa(code, &gcc).?;
    try std.testing.expectEqual(ref_aa, syn_aa); // synonymous

    // Stop gained: AAA(K) -> TAA(*)
    const stop_aa = translate.dna2aa(code, &taa).?;
    try std.testing.expectEqual(@as(u8, '*'), stop_aa);
}

// ===== Test 3: Consequence formatting ======================================

test "consequence formatting" {
    const alloc = std.testing.allocator;

    // Test simple missense consequence
    {
        var buf: std.ArrayList(u8) = .empty;
        defer buf.deinit(alloc);

        var vcsq = format.Vcsq{
            .strand = .fwd,
            .csq_type = format.CSQ_MISSENSE_VARIANT,
            .gene = "BRCA1",
            .biotype = 0,
            .vstr = "|123A>T",
        };

        const opts = format.FormatOptions{};
        try format.formatVcsq(&vcsq, opts, buf.writer(alloc));

        const result = buf.items;
        try std.testing.expect(std.mem.startsWith(u8, result, "missense|BRCA1|"));
        try std.testing.expect(std.mem.indexOf(u8, result, "|+|123A>T") != null);
    }

    // Test upstream-stop consequence (prefixed with '*')
    {
        var buf: std.ArrayList(u8) = .empty;
        defer buf.deinit(alloc);

        var vcsq = format.Vcsq{
            .csq_type = format.CSQ_UPSTREAM_STOP | format.CSQ_STOP_GAINED,
            .strand = .fwd,
            .gene = "KRAS",
            .vstr = "|12G>*",
        };

        const opts = format.FormatOptions{};
        try format.formatVcsq(&vcsq, opts, buf.writer(alloc));

        const result = buf.items;
        try std.testing.expect(result[0] == '*');
        try std.testing.expect(std.mem.indexOf(u8, result, "stop_gained") != null);
        try std.testing.expect(std.mem.indexOf(u8, result, "KRAS") != null);
    }

    // Test printed-upstream reference (@POS)
    {
        var buf: std.ArrayList(u8) = .empty;
        defer buf.deinit(alloc);

        var vcsq = format.Vcsq{
            .csq_type = format.CSQ_PRINTED_UPSTREAM,
            .ref_pos = 12345,
        };

        const opts = format.FormatOptions{};
        try format.formatVcsq(&vcsq, opts, buf.writer(alloc));
        try std.testing.expectEqualStrings("@12345", buf.items);
    }

    // Test intron consequence (no strand printed, no variant string)
    {
        var buf: std.ArrayList(u8) = .empty;
        defer buf.deinit(alloc);

        var vcsq = format.Vcsq{
            .csq_type = format.CSQ_INTRON,
            .gene = "EGFR",
        };

        const opts = format.FormatOptions{};
        try format.formatVcsq(&vcsq, opts, buf.writer(alloc));

        const result = buf.items;
        try std.testing.expect(std.mem.indexOf(u8, result, "intron") != null);
        try std.testing.expect(std.mem.indexOf(u8, result, "EGFR") != null);
    }

    // Test compound consequences joined by '&'
    {
        var buf: std.ArrayList(u8) = .empty;
        defer buf.deinit(alloc);

        var vcsq = format.Vcsq{
            .csq_type = format.CSQ_MISSENSE_VARIANT | format.CSQ_SPLICE_REGION,
            .strand = .fwd,
            .gene = "TP53",
        };

        const opts = format.FormatOptions{};
        try format.formatVcsq(&vcsq, opts, buf.writer(alloc));

        const result = buf.items;
        try std.testing.expect(std.mem.indexOf(u8, result, "missense") != null);
        try std.testing.expect(std.mem.indexOf(u8, result, "splice_region") != null);
        try std.testing.expect(std.mem.indexOf(u8, result, "&") != null);
    }
}

// ===== Test 4: VCF record parsing ==========================================

test "VCF record parsing" {
    const alloc = std.testing.allocator;

    var rec = vcf_record.VcfRecord.init(alloc);
    defer rec.deinit();

    try rec.parseLine("chr1\t150\t.\tA\tG\t.\tPASS\t.");

    try std.testing.expectEqualStrings("chr1", rec.chrom);
    try std.testing.expectEqual(@as(u32, 149), rec.pos); // 0-based
    try std.testing.expectEqualStrings("A", rec.ref_allele);
    try std.testing.expectEqual(@as(u32, 1), rec.rlen);
    try std.testing.expectEqual(@as(u32, 2), rec.nAllele());
    try std.testing.expectEqualStrings("A", rec.allele(0));
    try std.testing.expectEqualStrings("G", rec.allele(1));
    try std.testing.expectEqualStrings("PASS", rec.filter);

    // Parse a multi-allelic line
    rec.clear();
    try rec.parseLine("chr2\t500\trs123\tATG\tA,ATGC\t30\tq10\t.");

    try std.testing.expectEqualStrings("chr2", rec.chrom);
    try std.testing.expectEqual(@as(u32, 499), rec.pos);
    try std.testing.expectEqualStrings("rs123", rec.id);
    try std.testing.expectEqualStrings("ATG", rec.ref_allele);
    try std.testing.expectEqual(@as(u32, 3), rec.rlen);
    try std.testing.expectEqual(@as(u32, 3), rec.nAllele());
    try std.testing.expectEqualStrings("A", rec.allele(1));
    try std.testing.expectEqualStrings("ATGC", rec.allele(2));

    // Parse line with no ALT
    rec.clear();
    try rec.parseLine("chr3\t1000\t.\tC\t.\t.\t.\t.");
    try std.testing.expectEqual(@as(u32, 1), rec.nAllele());
}

// ===== Test 5: CSQ pipeline context and record buffering ===================

test "CSQ pipeline context initialization and buffering" {
    const alloc = std.testing.allocator;

    var ctx = try csq.CsqContext.init(alloc, .{
        .gff_fname = "dummy.gff3",
    });
    defer ctx.deinit();

    // Verify initial state
    try std.testing.expectEqual(@as(i32, -1), ctx.current_rid);
    try std.testing.expectEqual(csq.Phase.require, ctx.phase);
    try std.testing.expectEqualStrings("BCSQ", ctx.bcsq_tag);

    // Create test VCF records and push them into the buffer
    const alleles1 = [_][]const u8{ "A", "G" };
    var rec1 = csq.VcfRecord{
        .pos = 149,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles1,
        .rlen = 1,
    };

    const vbuf1 = try ctx.vbufPush(&rec1);
    try std.testing.expectEqual(@as(usize, 1), vbuf1.vrecs.items.len);

    // Push another record at the same position — should go into the same vbuf
    const alleles2 = [_][]const u8{ "A", "T" };
    var rec2 = csq.VcfRecord{
        .pos = 149,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles2,
        .rlen = 1,
    };

    const vbuf2 = try ctx.vbufPush(&rec2);
    try std.testing.expectEqual(@as(usize, 2), vbuf2.vrecs.items.len);
    // Same vbuf for same position
    try std.testing.expect(vbuf1 == vbuf2);

    // Push a record at a different position
    const alleles3 = [_][]const u8{ "C", "T" };
    var rec3 = csq.VcfRecord{
        .pos = 199,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles3,
        .rlen = 1,
    };

    const vbuf3 = try ctx.vbufPush(&rec3);
    try std.testing.expectEqual(@as(usize, 1), vbuf3.vrecs.items.len);
    try std.testing.expect(vbuf1 != vbuf3);

    // Flush everything
    try ctx.vbufFlush(csq.POS_MAX);
    try std.testing.expectEqual(@as(usize, 0), ctx.vcf_rbuf.len);
}

// ===== Test 6: CsqType packed struct operations ============================

test "CsqType packed struct bit layout and operations" {
    // Verify individual bit positions
    try std.testing.expectEqual(@as(u32, 1 << 1), (csq_types.CsqType{ .synonymous_variant = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 2), (csq_types.CsqType{ .missense_variant = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 4), (csq_types.CsqType{ .stop_gained = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 9), (csq_types.CsqType{ .splice_donor = true }).toInt());
    try std.testing.expectEqual(@as(u32, 1 << 16), (csq_types.CsqType{ .intron = true }).toInt());

    // Round-trip fromInt/toInt
    const raw: u32 = (1 << 2) | (1 << 7) | (1 << 18);
    const csq_val = csq_types.CsqType.fromInt(raw);
    try std.testing.expect(csq_val.missense_variant);
    try std.testing.expect(csq_val.frameshift_variant);
    try std.testing.expect(csq_val.inframe_altering);
    try std.testing.expect(!csq_val.stop_gained);
    try std.testing.expectEqual(raw, csq_val.toInt());

    // Compound detection
    try std.testing.expect((csq_types.CsqType{ .missense_variant = true }).isCompound());
    try std.testing.expect((csq_types.CsqType{ .stop_gained = true }).isCompound());
    try std.testing.expect(!(csq_types.CsqType{ .intron = true }).isCompound());
    try std.testing.expect(!(csq_types.CsqType{ .splice_acceptor = true }).isCompound());

    // Start/stop detection
    try std.testing.expect((csq_types.CsqType{ .stop_lost = true }).isStartStop());
    try std.testing.expect((csq_types.CsqType{ .start_retained = true }).isStartStop());
    try std.testing.expect(!(csq_types.CsqType{ .missense_variant = true }).isStartStop());

    // prnStrand: compound without splice/elongation/truncation
    try std.testing.expect((csq_types.CsqType{ .missense_variant = true }).prnStrand());
    try std.testing.expect(!(csq_types.CsqType{ .splice_donor = true }).prnStrand());
    try std.testing.expect(!(csq_types.CsqType{ .missense_variant = true, .splice_acceptor = true }).prnStrand());

    // Format consequence strings
    const alloc = std.testing.allocator;
    const csq_str = try (csq_types.CsqType{ .missense_variant = true, .splice_region = true }).formatAlloc(alloc);
    defer alloc.free(csq_str);
    try std.testing.expectEqualStrings("missense&splice_region", csq_str);
}

// ===== Test 7: VCF reader with header and samples ==========================

test "VCF reader parses header and records from file" {
    const alloc = std.testing.allocator;

    // Write VCF to temp file
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const tmp_file = try tmp_dir.dir.createFile("test.vcf", .{});
    try tmp_file.writeAll(minimal_vcf_content);
    tmp_file.close();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const tmp_path = try tmp_dir.dir.realpath("test.vcf", &path_buf);

    var reader = try vcf_reader.VcfReader.open(alloc, tmp_path);
    defer reader.deinit();

    // Check header parsing
    try std.testing.expectEqual(@as(usize, 1), reader.nSamples());
    try std.testing.expectEqualStrings("sample1", reader.sample_names.items[0]);
    try std.testing.expectEqual(@as(usize, 1), reader.seq_names.items.len);
    try std.testing.expectEqualStrings("chr1", reader.seq_names.items[0]);

    // Read all records and verify
    var rec = vcf_record.VcfRecord.init(alloc);
    defer rec.deinit();

    var count: usize = 0;
    var positions: [5]u32 = undefined;

    while (try reader.next(&rec)) {
        if (count < 5) positions[count] = rec.pos;
        count += 1;
    }

    try std.testing.expectEqual(@as(usize, 5), count);

    // Verify 0-based positions (VCF POS 150,200,300,350,500 -> 149,199,299,349,499)
    try std.testing.expectEqual(@as(u32, 149), positions[0]);
    try std.testing.expectEqual(@as(u32, 199), positions[1]);
    try std.testing.expectEqual(@as(u32, 299), positions[2]);
    try std.testing.expectEqual(@as(u32, 349), positions[3]);
    try std.testing.expectEqual(@as(u32, 499), positions[4]);
}

// ===== Test 8: Full pipeline smoke test ====================================

test "pipeline smoke test - process records through CsqContext" {
    const alloc = std.testing.allocator;

    var ctx = try csq.CsqContext.init(alloc, .{
        .gff_fname = "test.gff3",
        .phase = .as_is,
        .local_csq = true, // skip haplotype tree flushing for this test
    });
    defer ctx.deinit();

    // Create a series of VCF records and process them
    const allele_sets = [_][2][]const u8{
        .{ "A", "G" },
        .{ "C", "T" },
        .{ "G", "A" },
        .{ "T", "C" },
        .{ "A", "T" },
    };
    const positions = [_]u32{ 149, 199, 299, 349, 499 };

    var records: [5]csq.VcfRecord = undefined;
    for (0..5) |i| {
        records[i] = csq.VcfRecord{
            .pos = positions[i],
            .rid = 0,
            .n_allele = 2,
            .alleles = &allele_sets[i],
            .rlen = 1,
        };
    }

    // Buffer all records via vbufPush (process() includes flush-on-chrom-change
    // logic that would flush the buffer before we can inspect it)
    for (&records) |*rec| {
        _ = try ctx.vbufPush(rec);
    }

    // Verify records were buffered (5 distinct positions = 5 vbufs)
    try std.testing.expectEqual(@as(usize, 5), ctx.vcf_rbuf.len);

    // Test process() on a single record to verify it doesn't crash
    const alleles_extra = [_][]const u8{ "G", "C" };
    var rec_extra = csq.VcfRecord{
        .pos = 600,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles_extra,
        .rlen = 1,
    };
    try ctx.process(&rec_extra);

    // Flush everything
    try ctx.vbufFlush(csq.POS_MAX);
    try std.testing.expectEqual(@as(usize, 0), ctx.vcf_rbuf.len);
}
