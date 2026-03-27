// bench.zig -- Micro-benchmarks for bcftools-zig hot paths.
//
// Measures throughput of the key optimized paths in the Zig port:
//   1. Codon translation (translate.zig)
//   2. CDS translation (haplotype.zig)
//   3. RegionIndex overlap queries (region.zig)
//   4. Consequence formatting (format.zig)
//   5. VCF record buffering (csq.zig)
//   6. GFF parsing (gff.zig)
//
// Build & run:
//   zig build bench

const std = @import("std");
const lib = @import("bcftools_zig");

const translate = lib.translate;
const haplotype = lib.haplotype;
const region = lib.core_region;
const format = lib.csq_format;
const csq = lib.csq_pipeline;
const gff = lib.gff;

// ---------------------------------------------------------------------------
// Output helper — uses std.debug.print (stderr) for reliable unbuffered output
// ---------------------------------------------------------------------------

fn print(comptime fmt: []const u8, args: anytype) void {
    std.debug.print(fmt, args);
}

// ---------------------------------------------------------------------------
// Timing helpers
// ---------------------------------------------------------------------------

fn nanoTimestamp() i128 {
    return std.time.nanoTimestamp();
}

fn elapsed_ms(start: i128, end: i128) f64 {
    const ns: u64 = @intCast(end - start);
    return @as(f64, @floatFromInt(ns)) / 1_000_000.0;
}

fn ops_per_sec(iterations: usize, ms: f64) f64 {
    if (ms <= 0.0) return 0.0;
    return @as(f64, @floatFromInt(iterations)) / (ms / 1000.0);
}

fn printResult(name: []const u8, iterations: usize, ms: f64) void {
    const rate = ops_per_sec(iterations, ms);
    print("  {s}: {d:.1} ms  ({d:.0} ops/sec, n={d})\n", .{ name, ms, rate, iterations });
}

// ---------------------------------------------------------------------------
// PRNG for generating random sequences
// ---------------------------------------------------------------------------

const Rng = std.Random.Xoshiro256;

fn makeRng() Rng {
    return Rng.init(0xDEADBEEF_CAFEBABE);
}

// ---------------------------------------------------------------------------
// 1. Codon Translation
// ---------------------------------------------------------------------------

fn benchCodonTranslation() void {
    print("\n[1] Codon Translation (dna2aa)\n", .{});

    const code = translate.findGeneticCode(1) orelse {
        print("  ERROR: genetic code table 1 not found\n", .{});
        return;
    };

    const bases = "ACGT";
    const n_codons: usize = 1_000_000;

    // Pre-generate random codons
    var rng = makeRng();
    var codons: [1_000_000][3]u8 = undefined;
    for (&codons) |*c| {
        c[0] = bases[rng.random().uintLessThan(usize, 4)];
        c[1] = bases[rng.random().uintLessThan(usize, 4)];
        c[2] = bases[rng.random().uintLessThan(usize, 4)];
    }

    // Benchmark forward translation
    var sink: u32 = 0; // prevent dead-code elimination
    const start_fwd = nanoTimestamp();
    for (&codons) |*c| {
        sink +%= translate.dna2aa(code, c) orelse 0;
    }
    const end_fwd = nanoTimestamp();
    printResult("dna2aa (forward)", n_codons, elapsed_ms(start_fwd, end_fwd));

    // Benchmark reverse-complement translation
    const start_rev = nanoTimestamp();
    for (&codons) |*c| {
        sink +%= translate.cdna2aa(code, c) orelse 0;
    }
    const end_rev = nanoTimestamp();
    printResult("cdna2aa (revcomp)", n_codons, elapsed_ms(start_rev, end_rev));

    // Prevent optimizer from removing sink
    if (sink == 0xFFFFFFFF) print("", .{});
}

// ---------------------------------------------------------------------------
// 2. CDS Translation
// ---------------------------------------------------------------------------

fn benchCdsTranslation(allocator: std.mem.Allocator) !void {
    print("\n[2] CDS Translation (cdsTranslate)\n", .{});

    const code = translate.findGeneticCode(1) orelse {
        print("  ERROR: genetic code table 1 not found\n", .{});
        return;
    };

    // Build a realistic ~1002-base CDS (334 codons) starting with ATG
    const n_ref_pad = 10;
    const cds_len: usize = 1002;
    const sref_len = cds_len + 2 * n_ref_pad;

    var sref_buf: [sref_len]u8 = undefined;
    // Fill padding with reference bases
    for (&sref_buf) |*b| b.* = 'A';

    // Fill CDS region with a realistic coding sequence
    var rng = makeRng();
    const bases = "ACGT";
    // Start codon
    sref_buf[n_ref_pad + 0] = 'A';
    sref_buf[n_ref_pad + 1] = 'T';
    sref_buf[n_ref_pad + 2] = 'G';
    // Random middle
    for (sref_buf[n_ref_pad + 3 .. n_ref_pad + cds_len - 3]) |*b| {
        b.* = bases[rng.random().uintLessThan(usize, 4)];
    }
    // Stop codon
    sref_buf[n_ref_pad + cds_len - 3] = 'T';
    sref_buf[n_ref_pad + cds_len - 2] = 'A';
    sref_buf[n_ref_pad + cds_len - 1] = 'A';

    const seq = sref_buf[n_ref_pad .. n_ref_pad + cds_len];

    var result: std.ArrayList(u8) = .empty;
    defer result.deinit(allocator);
    var result_stop: std.ArrayList(u8) = .empty;
    defer result_stop.deinit(allocator);

    const iterations: usize = 10_000;

    // Benchmark forward strand
    const start_fwd = nanoTimestamp();
    for (0..iterations) |_| {
        try haplotype.cdsTranslate(
            allocator,
            &sref_buf,
            sref_len,
            seq,
            cds_len,
            0, // seq_beg
            0, // ref_beg
            @intCast(cds_len), // ref_end
            .forward,
            &result,
            &result_stop,
            0, // fill
            code,
        );
    }
    const end_fwd = nanoTimestamp();
    printResult("cdsTranslate (fwd, 1002bp)", iterations, elapsed_ms(start_fwd, end_fwd));

    // Benchmark reverse strand
    const start_rev = nanoTimestamp();
    for (0..iterations) |_| {
        try haplotype.cdsTranslate(
            allocator,
            &sref_buf,
            sref_len,
            seq,
            cds_len,
            0,
            0,
            @intCast(cds_len),
            .reverse,
            &result,
            &result_stop,
            0,
            code,
        );
    }
    const end_rev = nanoTimestamp();
    printResult("cdsTranslate (rev, 1002bp)", iterations, elapsed_ms(start_rev, end_rev));
}

// ---------------------------------------------------------------------------
// 3. RegionIndex Overlap Queries
// ---------------------------------------------------------------------------

fn benchRegionOverlap(allocator: std.mem.Allocator) !void {
    print("\n[3] RegionIndex Overlap Queries\n", .{});

    const RIdx = region.RegionIndex(u32);
    var idx = RIdx.init(allocator);
    defer idx.deinit();

    // Build 25 chromosome names
    const n_chroms: usize = 25;
    var chrom_names: [n_chroms][]const u8 = undefined;
    var name_bufs: [n_chroms][8]u8 = undefined;
    for (0..n_chroms) |i| {
        const name = std.fmt.bufPrint(&name_bufs[i], "chr{d}", .{i + 1}) catch unreachable;
        chrom_names[i] = name;
    }

    // Insert 50,000 intervals: 2,000 per chromosome
    const intervals_per_chrom: usize = 2_000;
    var rng = makeRng();
    const genome_size: u32 = 250_000_000; // 250 Mb

    for (0..n_chroms) |ci| {
        for (0..intervals_per_chrom) |j| {
            const beg = rng.random().uintLessThan(u32, genome_size);
            const span = rng.random().uintLessThan(u32, 50_000) + 100; // 100 bp to 50 kb
            const end = beg +| span;
            try idx.insert(chrom_names[ci], beg, end, @intCast(ci * intervals_per_chrom + j));
        }
    }

    // Force sort (first overlap call will sort, but let's be explicit)
    idx.sort();

    // Perform 1 million random overlap queries
    const n_queries: usize = 1_000_000;
    var hits: usize = 0;

    const start = nanoTimestamp();
    for (0..n_queries) |_| {
        const ci = rng.random().uintLessThan(usize, n_chroms);
        const qbeg = rng.random().uintLessThan(u32, genome_size);
        const qend = qbeg +| rng.random().uintLessThan(u32, 1000);
        var it = idx.overlap(chrom_names[ci], qbeg, qend);
        while (it.next()) |_| {
            hits += 1;
        }
    }
    const end = nanoTimestamp();
    const ms = elapsed_ms(start, end);
    printResult("overlap queries", n_queries, ms);
    print("    {d} intervals across {d} chroms, {d} total hits\n", .{
        n_chroms * intervals_per_chrom,
        n_chroms,
        hits,
    });
}

// ---------------------------------------------------------------------------
// 4. Consequence Formatting
// ---------------------------------------------------------------------------

fn benchFormatting(allocator: std.mem.Allocator) !void {
    print("\n[4] Consequence Formatting (formatVcsq)\n", .{});

    // Create diverse Vcsq entries
    const templates = [_]format.Vcsq{
        // Missense with variant string
        .{
            .strand = .fwd,
            .csq_type = format.CSQ_MISSENSE_VARIANT,
            .gene = "BRCA1",
            .vstr = "|5T>5I|100A>G",
        },
        // Synonymous
        .{
            .strand = .fwd,
            .csq_type = format.CSQ_SYNONYMOUS_VARIANT,
            .gene = "TP53",
            .vstr = "|10L|300C>T",
        },
        // Frameshift
        .{
            .strand = .rev,
            .csq_type = format.CSQ_FRAMESHIFT_VARIANT,
            .gene = "EGFR",
            .vstr = "|15fs|450delA",
        },
        // Compound: missense + splice_region
        .{
            .strand = .fwd,
            .csq_type = format.CSQ_MISSENSE_VARIANT | format.CSQ_SPLICE_REGION,
            .gene = "KRAS",
        },
        // Stop gained
        .{
            .strand = .fwd,
            .csq_type = format.CSQ_STOP_GAINED,
            .gene = "APC",
            .vstr = "|100R>*|298C>T",
        },
        // Upstream stop
        .{
            .strand = .rev,
            .csq_type = format.CSQ_UPSTREAM_STOP | format.CSQ_STOP_LOST,
            .gene = "MYC",
            .vstr = "|1*>1Q|500T>C",
        },
        // Printed upstream (reference)
        .{
            .csq_type = format.CSQ_PRINTED_UPSTREAM,
            .ref_pos = 12345,
        },
        // Intron
        .{
            .csq_type = format.CSQ_INTRON,
            .gene = "BRAF",
        },
        // UTR5
        .{
            .csq_type = format.CSQ_UTR5,
            .gene = "PIK3CA",
        },
        // Inframe deletion
        .{
            .strand = .fwd,
            .csq_type = format.CSQ_INFRAME_DELETION,
            .gene = "PTEN",
            .vstr = "|50delRKL|148delAGAAAACTT",
        },
    };

    const opts = format.FormatOptions{};
    const iterations: usize = 100_000;

    // Pre-allocate output buffer
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);
    try buf.ensureTotalCapacity(allocator, 256);

    const start = nanoTimestamp();
    for (0..iterations) |i| {
        buf.clearRetainingCapacity();
        var entry = templates[i % templates.len];
        try format.formatVcsq(&entry, opts, buf.writer(allocator));
    }
    const end = nanoTimestamp();
    printResult("formatVcsq (mixed types)", iterations, elapsed_ms(start, end));
}

// ---------------------------------------------------------------------------
// 5. VCF Record Buffering
// ---------------------------------------------------------------------------

fn benchVbufPush(allocator: std.mem.Allocator) !void {
    print("\n[5] VCF Record Buffering (vbufPush)\n", .{});

    var ctx = try csq.CsqContext.init(allocator, .{
        .gff_fname = "",
    });
    defer ctx.deinit();

    const iterations: usize = 100_000;

    // Create records: mix of same-position and different-position
    var records: [100_000]csq.VcfRecord = undefined;
    var rng = makeRng();

    var pos: u32 = 1000;
    for (0..iterations) |i| {
        // ~30% of records share position with previous
        if (i > 0 and rng.random().uintLessThan(u32, 100) < 30) {
            // same position
        } else {
            pos += rng.random().uintLessThan(u32, 100) + 1;
        }
        records[i] = .{
            .pos = pos,
            .rid = 0,
            .n_allele = 2,
            .alleles = &[_][]const u8{ "A", "T" },
            .rlen = 1,
            .chr = "chr1",
        };
    }

    const start = nanoTimestamp();
    for (0..iterations) |i| {
        _ = try ctx.vbufPush(&records[i]);
    }
    const end = nanoTimestamp();
    printResult("vbufPush (mixed positions)", iterations, elapsed_ms(start, end));
}

// ---------------------------------------------------------------------------
// 6. GFF Parsing
// ---------------------------------------------------------------------------

fn benchGffParsing(allocator: std.mem.Allocator) !void {
    print("\n[6] GFF Parsing (parseContent)\n", .{});

    // Generate synthetic GFF3 content: ~10,000 genes, each with 1 mRNA + 3 CDS + 3 exons
    const n_genes: usize = 10_000;

    var content: std.ArrayList(u8) = .empty;
    defer content.deinit(allocator);

    // GFF3 header
    try content.appendSlice(allocator, "##gff-version 3\n");

    const w = content.writer(allocator);

    var pos: u32 = 1000;
    for (0..n_genes) |gi| {
        const gene_beg = pos;
        const gene_end = pos + 10_000;

        // Gene line
        try w.print("chr1\tensembl\tgene\t{d}\t{d}\t.\t+\t.\tID=gene{d};Name=GENE{d};biotype=protein_coding\n", .{ gene_beg, gene_end, gi, gi });

        // mRNA line
        try w.print("chr1\tensembl\tmRNA\t{d}\t{d}\t.\t+\t.\tID=tr{d};Parent=gene{d};biotype=protein_coding\n", .{ gene_beg, gene_end, gi, gi });

        // 3 CDS entries per gene
        var cds_pos = gene_beg + 100;
        for (0..3) |_| {
            const cds_end = cds_pos + 300;
            try w.print("chr1\tensembl\tCDS\t{d}\t{d}\t.\t+\t0\tParent=tr{d}\n", .{ cds_pos, cds_end, gi });
            cds_pos = cds_end + 500;
        }

        // 3 exon entries matching CDS
        cds_pos = gene_beg + 100;
        for (0..3) |_| {
            const exon_end = cds_pos + 300;
            try w.print("chr1\tensembl\texon\t{d}\t{d}\t.\t+\t.\tParent=tr{d}\n", .{ cds_pos, exon_end, gi });
            cds_pos = exon_end + 500;
        }

        pos = gene_end + 5000;
    }

    print("    Generated {d} bytes of GFF3 ({d} genes)\n", .{ content.items.len, n_genes });

    // Parse the content
    var parser = gff.GffParser.init(allocator);
    defer parser.deinit();

    const start = nanoTimestamp();
    try parser.parseContent(content.items);
    const end = nanoTimestamp();
    const ms = elapsed_ms(start, end);

    print("  parseContent: {d:.1} ms  ({d:.0} genes/sec)\n", .{ ms, ops_per_sec(n_genes, ms) });
    print("    Transcripts registered: {d}\n", .{parser.transcripts.count()});
}

// ---------------------------------------------------------------------------
// 7. SIMD Nucleotide Encoding (encodeNt4x16)
// ---------------------------------------------------------------------------

fn benchSimdEncoding() void {
    print("\n[7] SIMD Nucleotide Encoding (encodeNt4x16 vs scalar nt4)\n", .{});

    // Generate a random DNA sequence
    var rng = makeRng();
    const bases = "ACGTacgt";
    var seq: [1024]u8 = undefined;
    for (&seq) |*b| b.* = bases[rng.random().uintLessThan(usize, 8)];

    const iterations: usize = 1_000_000;

    // Benchmark scalar nt4 lookup
    var sink: u32 = 0;
    const start_scalar = nanoTimestamp();
    for (0..iterations) |_| {
        for (seq) |b| {
            sink +%= translate.nt4[b];
        }
    }
    const end_scalar = nanoTimestamp();
    printResult("scalar nt4 (1024 bases)", iterations, elapsed_ms(start_scalar, end_scalar));

    // Benchmark SIMD encodeNt4x16
    const start_simd = nanoTimestamp();
    for (0..iterations) |_| {
        var i: usize = 0;
        while (i + 16 <= seq.len) : (i += 16) {
            const result = translate.encodeNt4x16(seq[i..][0..16].*);
            sink +%= @reduce(.Add, result);
        }
    }
    const end_simd = nanoTimestamp();
    printResult("SIMD encodeNt4x16 (1024 bases)", iterations, elapsed_ms(start_simd, end_simd));

    if (sink == 0xFFFFFFFF) print("", .{});
}

// ---------------------------------------------------------------------------
// 8. SIMD Batch Translation (batchTranslate vs scalar dna2aa)
// ---------------------------------------------------------------------------

fn benchBatchTranslate() void {
    print("\n[8] Batch Translation (batchTranslate vs scalar dna2aa)\n", .{});

    const code = translate.findGeneticCode(1) orelse {
        print("  ERROR: genetic code table 1 not found\n", .{});
        return;
    };

    var rng = makeRng();
    const bases = "ACGT";
    // 960 bases = 320 codons (divisible by 48 for clean SIMD path)
    var seq: [960]u8 = undefined;
    for (&seq) |*b| b.* = bases[rng.random().uintLessThan(usize, 4)];

    var out: [320]u8 = undefined;
    const iterations: usize = 100_000;

    // Benchmark scalar translation
    var sink: u32 = 0;
    const start_scalar = nanoTimestamp();
    for (0..iterations) |_| {
        var i: usize = 0;
        var j: usize = 0;
        while (i + 3 <= seq.len) : (i += 3) {
            out[j] = translate.dna2aa(code, seq[i..][0..3]) orelse 'X';
            j += 1;
        }
        sink +%= out[0];
    }
    const end_scalar = nanoTimestamp();
    printResult("scalar dna2aa (960 bases)", iterations, elapsed_ms(start_scalar, end_scalar));

    // Benchmark SIMD batch translation
    const start_simd = nanoTimestamp();
    for (0..iterations) |_| {
        _ = translate.batchTranslate(code, &seq, &out);
        sink +%= out[0];
    }
    const end_simd = nanoTimestamp();
    printResult("SIMD batchTranslate (960 bases)", iterations, elapsed_ms(start_simd, end_simd));

    if (sink == 0xFFFFFFFF) print("", .{});
}

// ---------------------------------------------------------------------------
// 9. SIMD Protein Comparison (simdFirstMismatch)
// ---------------------------------------------------------------------------

fn benchProteinComparison() void {
    print("\n[9] Protein Comparison (SIMD simdFirstMismatch vs scalar)\n", .{});

    var rng = makeRng();
    const aas = "ACDEFGHIKLMNPQRSTVWY*";

    // Generate a 500-AA protein pair that differs at position 450
    var prot_a: [500]u8 = undefined;
    var prot_b: [500]u8 = undefined;
    for (0..500) |i| {
        prot_a[i] = aas[rng.random().uintLessThan(usize, aas.len)];
        prot_b[i] = prot_a[i];
    }
    // Introduce a mismatch near the end
    prot_b[450] = if (prot_a[450] == 'X') 'Y' else 'X';

    const iterations: usize = 5_000_000;

    // Benchmark scalar comparison
    var sink: usize = 0;
    const start_scalar = nanoTimestamp();
    for (0..iterations) |iter| {
        // Vary the mismatch position to defeat branch prediction
        prot_b[450] = prot_a[450]; // restore
        const mismatch_pos = 400 + (iter % 100);
        const saved = prot_b[mismatch_pos];
        prot_b[mismatch_pos] = if (prot_a[mismatch_pos] == 'X') 'Y' else 'X';
        for (0..500) |i| {
            if (prot_a[i] != prot_b[i]) {
                sink +%= i;
                break;
            }
        }
        prot_b[mismatch_pos] = saved;
    }
    const end_scalar = nanoTimestamp();
    printResult("scalar first-mismatch (500 AA)", iterations, elapsed_ms(start_scalar, end_scalar));

    // Benchmark SIMD comparison
    const start_simd = nanoTimestamp();
    for (0..iterations) |iter| {
        prot_b[450] = prot_a[450]; // restore
        const mismatch_pos = 400 + (iter % 100);
        const saved = prot_b[mismatch_pos];
        prot_b[mismatch_pos] = if (prot_a[mismatch_pos] == 'X') 'Y' else 'X';
        if (haplotype.simdFirstMismatch(&prot_a, &prot_b)) |pos| {
            sink +%= pos;
        }
        prot_b[mismatch_pos] = saved;
    }
    const end_simd = nanoTimestamp();
    printResult("SIMD simdFirstMismatch (500 AA)", iterations, elapsed_ms(start_simd, end_simd));

    if (sink == 0xFFFFFFFF) print("", .{});
}

// ---------------------------------------------------------------------------
// 10. SIMD Uppercase Conversion
// ---------------------------------------------------------------------------

fn benchUppercase() void {
    print("\n[10] Uppercase Conversion (SIMD vs scalar)\n", .{});

    var rng = makeRng();
    const iterations: usize = 1_000_000;

    // Generate a mixed-case sequence of 1024 bytes
    var seq_template: [1024]u8 = undefined;
    for (&seq_template) |*b| {
        const base: u8 = @intCast('a' + rng.random().uintLessThan(u8, 26));
        // 50% lowercase, 50% uppercase
        b.* = if (rng.random().uintLessThan(u32, 2) == 0) base else base - 32;
    }

    // Benchmark scalar uppercase
    var seq_buf: [1024]u8 = undefined;
    var sink: u32 = 0;
    const start_scalar = nanoTimestamp();
    for (0..iterations) |_| {
        @memcpy(&seq_buf, &seq_template);
        for (&seq_buf) |*c| {
            if (c.* >= 'a' and c.* <= 'z') c.* -= 32;
        }
        sink +%= seq_buf[0];
    }
    const end_scalar = nanoTimestamp();
    printResult("scalar uppercase (1024 bytes)", iterations, elapsed_ms(start_scalar, end_scalar));

    // Benchmark SIMD uppercase (using the CsqContext.uppercaseInPlace is private,
    // so we replicate the SIMD logic here for benchmarking)
    const start_simd = nanoTimestamp();
    for (0..iterations) |_| {
        @memcpy(&seq_buf, &seq_template);
        var i: usize = 0;
        while (i + 16 <= seq_buf.len) {
            var chunk: @Vector(16, u8) = seq_buf[i..][0..16].*;
            const lower_a: @Vector(16, u8) = @splat('a');
            const lower_z: @Vector(16, u8) = @splat('z');
            const ge_a: @Vector(16, bool) = chunk >= lower_a;
            const le_z: @Vector(16, bool) = chunk <= lower_z;
            const is_lower = @select(bool, ge_a, le_z, @as(@Vector(16, bool), @splat(false)));
            const adjustment = @select(u8, is_lower, @as(@Vector(16, u8), @splat(32)), @as(@Vector(16, u8), @splat(0)));
            chunk -= adjustment;
            seq_buf[i..][0..16].* = chunk;
            i += 16;
        }
        while (i < seq_buf.len) {
            if (seq_buf[i] >= 'a' and seq_buf[i] <= 'z') seq_buf[i] -= 32;
            i += 1;
        }
        sink +%= seq_buf[0];
    }
    const end_simd = nanoTimestamp();
    printResult("SIMD uppercase (1024 bytes)", iterations, elapsed_ms(start_simd, end_simd));

    if (sink == 0xFFFFFFFF) print("", .{});
}

// ---------------------------------------------------------------------------
// 11. Parallel hapFinalize Thread Pool
// ---------------------------------------------------------------------------

fn benchParallelHapFinalize(allocator: std.mem.Allocator) !void {
    print("\n[11] Parallel hapFinalize (thread pool dispatch)\n", .{});

    const code = translate.findGeneticCode(1) orelse {
        print("  ERROR: genetic code table 1 not found\n", .{});
        return;
    };

    // Build synthetic transcripts with haplotype trees.
    // Each transcript has a root node with several CDS child nodes ending in
    // leaf nodes (nend > 0). hapFinalize traverses the tree, splices the ref,
    // translates codons, and determines consequences -- this is the expensive
    // per-transcript work we want to parallelize.
    const gff_t = lib.gff_types;
    const csq_t = lib.csq_types;
    const n_transcripts: usize = 200;
    const cds_per_transcript: usize = 5;
    const cds_len: u32 = 300; // 100 codons per CDS
    const n_ref_pad: u32 = @intCast(csq_t.n_ref_pad);

    // Pre-build transcript and tree structures
    var transcripts: [n_transcripts]*gff_t.Transcript = undefined;
    var tscripts: [n_transcripts]*csq_t.Tscript = undefined;
    var genes: [n_transcripts]gff_t.Gene = undefined;
    var cds_entries: [n_transcripts][cds_per_transcript]*gff_t.CdsEntry = undefined;

    // RNG for sequence generation
    var rng = makeRng();
    const bases = "ACGT";

    for (0..n_transcripts) |ti| {
        // Gene
        genes[ti] = .{
            .name = null,
            .iseq = 0,
            .id = @intCast(ti),
            .beg = 0,
            .end = cds_per_transcript * (cds_len + 200) + 2 * n_ref_pad,
            .strand = .forward,
            .used = true,
        };

        // Transcript
        const tr = try allocator.create(gff_t.Transcript);
        tr.* = gff_t.Transcript.init(allocator);
        tr.id = @intCast(ti);
        tr.beg = 0;
        tr.end = genes[ti].end;
        tr.strand = .forward;
        tr.biotype = .protein_coding;
        tr.gene = &genes[ti];
        transcripts[ti] = tr;

        // CDS entries
        var cds_offset: u32 = 0;
        for (0..cds_per_transcript) |ci| {
            const entry = try allocator.create(gff_t.CdsEntry);
            entry.* = .{
                .tr = tr,
                .beg = cds_offset,
                .pos = cds_offset,
                .len = cds_len,
                .icds = @intCast(ci),
                .phase = .phase0,
            };
            try tr.cds.append(allocator, entry);
            cds_entries[ti][ci] = entry;
            cds_offset += cds_len + 200; // gap between CDS
        }

        // Tscript auxiliary data
        const tscript = try allocator.create(csq_t.Tscript);
        tscript.* = .{};
        tscripts[ti] = tscript;
        tr.aux = tscript;

        // Build ref_seq (entire transcript region with padding)
        const ref_len = tr.end + 2 * n_ref_pad;
        const ref_seq = try allocator.alloc(u8, ref_len);
        // Start codon ATG at first CDS
        for (ref_seq) |*b| b.* = bases[rng.random().uintLessThan(usize, 4)];
        // Ensure start codon
        if (n_ref_pad + 0 < ref_seq.len) ref_seq[n_ref_pad] = 'A';
        if (n_ref_pad + 1 < ref_seq.len) ref_seq[n_ref_pad + 1] = 'T';
        if (n_ref_pad + 2 < ref_seq.len) ref_seq[n_ref_pad + 2] = 'G';
        tscript.ref_seq = ref_seq;

        // Build haplotype tree: root -> CDS children -> leaf
        const root = try allocator.create(csq_t.HapNode);
        root.* = csq_t.HapNode.init(.root);
        tscript.root = root;

        // Create a chain: root has one CDS child per exon, last is a leaf
        var parent = root;
        for (0..cds_per_transcript) |ci| {
            const child = try allocator.create(csq_t.HapNode);
            // Build a CDS node with a single substitution
            const seq = try allocator.alloc(u8, cds_len);
            for (seq) |*b| b.* = bases[rng.random().uintLessThan(usize, 4)];
            // ATG at start
            if (ci == 0) {
                seq[0] = 'A';
                seq[1] = 'T';
                seq[2] = 'G';
            }
            child.* = csq_t.HapNode.init(.cds);
            child.payload = .{ .cds = .{ .seq = seq } };
            child.sbeg = @intCast(ci * cds_len);
            child.icds = @intCast(ci);
            child.dlen = 0;
            child.rlen = 1;
            child.rbeg = @intCast(ci * (cds_len + 200) + 50);
            child.rec_pos = child.rbeg;

            // Last child is a leaf
            if (ci == cds_per_transcript - 1) {
                child.nend = 1;
            }

            try parent.children.append(allocator, child);
            parent = child;
        }
    }

    defer {
        for (0..n_transcripts) |ti| {
            // Free tree nodes
            const tscript = tscripts[ti];
            if (tscript.root) |root| {
                var stack: [32]*csq_t.HapNode = undefined;
                var sp: usize = 0;
                stack[sp] = root;
                sp += 1;
                while (sp > 0) {
                    sp -= 1;
                    const node = stack[sp];
                    for (node.children.items) |ch| {
                        if (sp < stack.len) {
                            stack[sp] = ch;
                            sp += 1;
                        }
                    }
                    if (node.payload == .cds) {
                        if (node.payload.cds.seq) |s| allocator.free(s);
                    }
                    node.deinit(allocator);
                    allocator.destroy(node);
                }
            }
            if (tscript.sref) |s| allocator.free(s);
            if (tscript.ref_seq) |r| allocator.free(r);
            tscript.deinit(allocator);
            allocator.destroy(tscript);
            for (cds_entries[ti][0..cds_per_transcript]) |entry| allocator.destroy(entry);
            transcripts[ti].cds.deinit(allocator);
            allocator.destroy(transcripts[ti]);
        }
    }

    // ---- Benchmark sequential hapFinalize ----
    const iterations: usize = 50;

    const start_seq = nanoTimestamp();
    for (0..iterations) |_| {
        // Reset sref so hapFinalize rebuilds it each time
        for (0..n_transcripts) |ti| {
            if (tscripts[ti].sref) |s| allocator.free(s);
            tscripts[ti].sref = null;
            tscripts[ti].nsref = 0;
        }
        for (0..n_transcripts) |ti| {
            var ctx = haplotype.HapContext.init(allocator, code);
            defer ctx.deinit();
            ctx.tr = transcripts[ti];
            haplotype.hapFinalize(&ctx) catch {};
        }
    }
    const end_seq = nanoTimestamp();
    const ms_seq = elapsed_ms(start_seq, end_seq);
    printResult("sequential hapFinalize", n_transcripts * iterations, ms_seq);

    // ---- Benchmark parallel hapFinalize using std.Thread.Pool ----
    const thread_counts = [_]u32{ 2, 4, 8 };
    for (thread_counts) |n_threads| {
        var pool: std.Thread.Pool = undefined;
        pool.init(.{
            .allocator = allocator,
            .n_jobs = n_threads,
        }) catch {
            print("  {d} threads: failed to create pool\n", .{n_threads});
            continue;
        };
        defer pool.deinit();

        const start_par = nanoTimestamp();
        for (0..iterations) |_| {
            for (0..n_transcripts) |ti| {
                if (tscripts[ti].sref) |s| allocator.free(s);
                tscripts[ti].sref = null;
                tscripts[ti].nsref = 0;
            }

            var wg = std.Thread.WaitGroup{};
            for (0..n_transcripts) |ti| {
                pool.spawnWg(&wg, struct {
                    fn work(alloc: std.mem.Allocator, gc: *const translate.GeneticCode, tr: *gff_t.Transcript) void {
                        var ctx = haplotype.HapContext.init(alloc, gc);
                        defer ctx.deinit();
                        ctx.tr = tr;
                        haplotype.hapFinalize(&ctx) catch {};
                    }
                }.work, .{ allocator, code, transcripts[ti] });
            }
            wg.wait();
        }
        const end_par = nanoTimestamp();
        const ms_par = elapsed_ms(start_par, end_par);
        const speedup = ms_seq / ms_par;
        print("  {d} threads: {d:.1} ms  ({d:.2}x speedup)\n", .{ n_threads, ms_par, speedup });
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    print("bcftools-zig micro-benchmarks\n", .{});
    print("=============================\n", .{});

    benchCodonTranslation();
    try benchCdsTranslation(allocator);
    try benchRegionOverlap(allocator);
    try benchFormatting(allocator);
    try benchVbufPush(allocator);
    try benchGffParsing(allocator);
    benchSimdEncoding();
    benchBatchTranslate();
    benchProteinComparison();
    benchUppercase();
    try benchParallelHapFinalize(allocator);

    print("\nDone.\n", .{});
}
