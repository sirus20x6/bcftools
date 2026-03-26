const std = @import("std");
const lib = @import("bcftools_zig");
const csq_mod = lib.csq_pipeline;
const gff_mod = lib.gff;
const VcfReader = lib.vcf_reader.VcfReader;
const VcfRecord = lib.vcf_record.VcfRecord;
const CsqContext = csq_mod.CsqContext;
const GffParser = gff_mod.GffParser;
const Phase = csq_mod.Phase;

const usage_text =
    \\
    \\Usage: bcftools-zig <command> [options]
    \\
    \\Commands:
    \\  csq        Haplotype-aware consequence caller
    \\
    \\Options:
    \\  --help     Show this help message
    \\  --version  Show version information
    \\
;

const csq_usage_text =
    \\
    \\About: Haplotype-aware consequence caller.
    \\Usage: bcftools-zig csq [OPTIONS] in.vcf.gz
    \\
    \\Required options:
    \\  -f, --fasta-ref FILE        Reference file in FASTA format
    \\  -g, --gff-annot FILE        GFF3 annotation file
    \\
    \\CSQ options:
    \\  -B, --trim-protein-seq INT  Brief protein predictions, show max INT aa [0: unlimited]
    \\  -c, --custom-tag STRING     Use this tag instead of BCSQ [BCSQ]
    \\  -l, --local-csq             Localized predictions, do not walk haplotypes
    \\  -n, --ncsq INT              Maximum number of per-haplotype consequences [15]
    \\  -p, --phase a|m|r|R|s       How to handle unphased heterozygous genotypes [r]
    \\      --force                  Run even under sub-optimal conditions
    \\
    \\General options:
    \\  -e, --exclude EXPR          Exclude sites for which the expression is true
    \\  -i, --include EXPR          Include only sites for which the expression is true
    \\  -o, --output FILE           Write output to FILE [standard output]
    \\  -O, --output-type b|u|z|v   b: compressed BCF, u: uncompressed BCF,
    \\                              z: compressed VCF, v: uncompressed VCF [v]
    \\  -r, --regions REGION        Restrict to comma-separated list of regions
    \\  -R, --regions-file FILE     Restrict to regions listed in FILE
    \\  -s, --samples LIST          Samples to include
    \\  -S, --samples-file FILE     Samples to include from FILE
    \\  -t, --targets REGION        Similar to -r but streams rather than index-jumps
    \\  -T, --targets-file FILE     Similar to -R but streams rather than index-jumps
    \\      --threads INT            Use multithreading with INT worker threads [0]
    \\      --write-index            Automatically index the output file
    \\
;

const stdout_file = std.fs.File{ .handle = std.posix.STDOUT_FILENO };
const stderr_file = std.fs.File{ .handle = std.posix.STDERR_FILENO };

// -------------------------------------------------------------------------
// CSQ subcommand options
// -------------------------------------------------------------------------

const CsqOptions = struct {
    fasta_fname: ?[]const u8 = null,
    gff_fname: ?[]const u8 = null,
    input_fname: ?[]const u8 = null,
    output_fname: ?[]const u8 = null,
    output_type: u8 = 'v',
    phase: Phase = .require,
    local_csq: bool = false,
    force: bool = false,
    ncsq: u32 = 15,
    custom_tag: []const u8 = "BCSQ",
    brief_predictions: u32 = 0,
    show_help: bool = false,
};

const CsqArgError = error{
    MissingArgValue,
    InvalidPhase,
    InvalidNcsq,
    InvalidOutputType,
    InvalidBriefPredictions,
    UnknownOption,
};

fn parseCsqArgs(args_iter: *std.process.ArgIterator) CsqArgError!CsqOptions {
    var opts = CsqOptions{};
    while (args_iter.next()) |arg| {
        if (std.mem.eql(u8, arg, "-f") or std.mem.eql(u8, arg, "--fasta-ref")) {
            opts.fasta_fname = args_iter.next() orelse return error.MissingArgValue;
        } else if (std.mem.eql(u8, arg, "-g") or std.mem.eql(u8, arg, "--gff-annot")) {
            opts.gff_fname = args_iter.next() orelse return error.MissingArgValue;
        } else if (std.mem.eql(u8, arg, "-p") or std.mem.eql(u8, arg, "--phase")) {
            const val = args_iter.next() orelse return error.MissingArgValue;
            if (val.len == 0) return error.InvalidPhase;
            opts.phase = switch (val[0]) {
                'a' => .as_is,
                'm' => .merge,
                'r' => .require,
                'R' => .non_ref,
                's' => .skip,
                else => return error.InvalidPhase,
            };
        } else if (std.mem.eql(u8, arg, "-l") or std.mem.eql(u8, arg, "--local-csq")) {
            opts.local_csq = true;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            opts.output_fname = args_iter.next() orelse return error.MissingArgValue;
        } else if (std.mem.eql(u8, arg, "-O") or std.mem.eql(u8, arg, "--output-type")) {
            const val = args_iter.next() orelse return error.MissingArgValue;
            if (val.len == 0) return error.InvalidOutputType;
            switch (val[0]) {
                'v', 'z', 'b', 'u' => opts.output_type = val[0],
                else => return error.InvalidOutputType,
            }
        } else if (std.mem.eql(u8, arg, "-n") or std.mem.eql(u8, arg, "--ncsq")) {
            const val = args_iter.next() orelse return error.MissingArgValue;
            opts.ncsq = std.fmt.parseInt(u32, val, 10) catch return error.InvalidNcsq;
            if (opts.ncsq == 0) return error.InvalidNcsq;
        } else if (std.mem.eql(u8, arg, "-c") or std.mem.eql(u8, arg, "--custom-tag")) {
            opts.custom_tag = args_iter.next() orelse return error.MissingArgValue;
        } else if (std.mem.eql(u8, arg, "-B") or std.mem.eql(u8, arg, "--trim-protein-seq")) {
            const val = args_iter.next() orelse return error.MissingArgValue;
            opts.brief_predictions = std.fmt.parseInt(u32, val, 10) catch return error.InvalidBriefPredictions;
        } else if (std.mem.eql(u8, arg, "--force")) {
            opts.force = true;
        } else if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            opts.show_help = true;
        } else if (arg.len > 0 and arg[0] == '-') {
            // Skip unknown options that take a value (consume next arg to avoid
            // treating it as a positional).  Single-dash flags like -e, -i, -r,
            // -R, -s, -S, -t, -T and long options like --threads, --write-index
            // are not yet implemented but should not cause a hard error.
            _ = args_iter.next();
        } else {
            // Positional argument: input VCF file
            opts.input_fname = arg;
        }
    }
    return opts;
}

// -------------------------------------------------------------------------
// Build a csq_mod.VcfRecord (pipeline type) from a vcf.VcfRecord
// -------------------------------------------------------------------------

/// Scratch space for allele slices when bridging VcfRecord types.
/// This avoids per-record allocation.
var allele_scratch: [256][]const u8 = undefined;

fn buildCsqRecord(rec: *const VcfRecord) csq_mod.VcfRecord {
    const n: u32 = rec.nAllele();
    const count = @min(n, allele_scratch.len);
    if (count > 0) {
        allele_scratch[0] = rec.ref_allele;
    }
    for (0..@min(rec.alt_alleles.items.len, allele_scratch.len - 1)) |i| {
        allele_scratch[i + 1] = rec.alt_alleles.items[i];
    }
    return .{
        .pos = rec.pos,
        .rid = rec.rid,
        .n_allele = count,
        .alleles = allele_scratch[0..count],
        .rlen = rec.rlen,
    };
}

// -------------------------------------------------------------------------
// Write a VCF record line, appending BCSQ to INFO if consequences exist
// -------------------------------------------------------------------------

fn writeVcfRecord(
    out: std.fs.File,
    rec: *const VcfRecord,
    csq_ctx: *CsqContext,
) !void {
    // The record's _storage holds the full original tab-separated line.
    const line = rec._storage orelse return;

    // If the CSQ context has formatted output, we need to inject BCSQ into INFO.
    // For now the pipeline stubs produce no output, so we pass through as-is.
    // When the pipeline is complete, csq_ctx.output will contain the BCSQ string
    // after vbufFlush processes the record.
    _ = csq_ctx;

    try out.writeAll(line);
    try out.writeAll("\n");
}

// -------------------------------------------------------------------------
// CSQ subcommand entry point
// -------------------------------------------------------------------------

fn runCsq(args_iter: *std.process.ArgIterator) !void {
    const opts = parseCsqArgs(args_iter) catch |err| {
        const msg = switch (err) {
            error.MissingArgValue => "Error: missing value for option\n",
            error.InvalidPhase => "Error: invalid phase value, expected one of: a, m, r, R, s\n",
            error.InvalidNcsq => "Error: expected positive integer with --ncsq\n",
            error.InvalidOutputType => "Error: invalid output type, expected one of: v, z, b, u\n",
            error.InvalidBriefPredictions => "Error: expected non-negative integer with --trim-protein-seq\n",
            error.UnknownOption => "Error: unknown option\n",
        };
        stderr_file.writeAll(msg) catch {};
        stderr_file.writeAll(csq_usage_text) catch {};
        std.process.exit(1);
    };

    if (opts.show_help) {
        stdout_file.writeAll(csq_usage_text) catch {};
        return;
    }

    // Validate required arguments
    if (opts.fasta_fname == null) {
        stderr_file.writeAll("Error: missing the --fasta-ref option\n") catch {};
        std.process.exit(1);
    }
    if (opts.gff_fname == null) {
        stderr_file.writeAll("Error: missing the --gff-annot option\n") catch {};
        std.process.exit(1);
    }
    if (opts.input_fname == null) {
        stderr_file.writeAll("Error: no input VCF file specified\n") catch {};
        stderr_file.writeAll(csq_usage_text) catch {};
        std.process.exit(1);
    }

    // Only text VCF output is supported for now (no htslib BCF writing).
    if (opts.output_type != 'v') {
        stderr_file.writeAll("Error: only uncompressed VCF output (-O v) is currently supported\n") catch {};
        std.process.exit(1);
    }

    const allocator = std.heap.page_allocator;

    // ---- Parse GFF annotations ----
    var gff = GffParser.init(allocator);
    defer gff.deinit();
    gff.verbosity = if (opts.force) 0 else 1;
    gff.force = opts.force;

    gff.parse(opts.gff_fname.?) catch |err| {
        std.debug.print("Error: failed to parse GFF file '{s}': {}\n", .{ opts.gff_fname.?, err });
        std.process.exit(1);
    };

    // ---- Open VCF input ----
    var reader = VcfReader.open(allocator, opts.input_fname.?) catch |err| {
        std.debug.print("Error: failed to open VCF file '{s}': {}\n", .{ opts.input_fname.?, err });
        std.process.exit(1);
    };
    defer reader.deinit();

    // ---- Open output ----
    var out_file: std.fs.File = undefined;
    var out_file_needs_close = false;
    if (opts.output_fname) |fname| {
        out_file = std.fs.cwd().createFile(fname, .{}) catch |err| {
            std.debug.print("Error: failed to open output file '{s}': {}\n", .{ fname, err });
            std.process.exit(1);
        };
        out_file_needs_close = true;
    } else {
        out_file = stdout_file;
    }
    defer if (out_file_needs_close) out_file.close();

    // ---- Write VCF header ----
    // First, write all existing header lines except the #CHROM line.
    // Insert the BCSQ INFO definition before the #CHROM line.
    var chrom_line: ?[]const u8 = null;
    for (reader.header_lines.items) |hline| {
        if (hline.len > 0 and hline[0] == '#' and (hline.len < 2 or hline[1] != '#')) {
            // This is the #CHROM line; save it for after we inject BCSQ header.
            chrom_line = hline;
        } else {
            out_file.writeAll(hline) catch |err| {
                std.debug.print("Error: failed to write header: {}\n", .{err});
                std.process.exit(1);
            };
            out_file.writeAll("\n") catch {};
        }
    }

    // Inject BCSQ INFO header line
    out_file.writeAll("##INFO=<ID=") catch {};
    out_file.writeAll(opts.custom_tag) catch {};
    out_file.writeAll(",Number=.,Type=String,Description=\"Haplotype-aware consequence annotation from BCFtools/csq\">\n") catch {};

    // Inject BCSQ FORMAT header line (for per-sample bitmask)
    out_file.writeAll("##FORMAT=<ID=") catch {};
    out_file.writeAll(opts.custom_tag) catch {};
    out_file.writeAll(",Number=.,Type=Integer,Description=\"Bitmask of indexes to consequence types listed in the INFO/") catch {};
    out_file.writeAll(opts.custom_tag) catch {};
    out_file.writeAll(" tag\">\n") catch {};

    // Write the #CHROM line
    if (chrom_line) |cl| {
        out_file.writeAll(cl) catch {};
        out_file.writeAll("\n") catch {};
    }

    // ---- Initialize CSQ context ----
    var csq_ctx = CsqContext.init(allocator, .{
        .gff_fname = opts.gff_fname.?,
        .fasta_fname = opts.fasta_fname.?,
        .phase = opts.phase,
        .local_csq = opts.local_csq,
        .verbosity = 1,
        .force = opts.force,
        .bcsq_tag = opts.custom_tag,
        .ncsq2_max = opts.ncsq * 2,
        .brief_predictions = opts.brief_predictions,
    }) catch |err| {
        std.debug.print("Error: failed to initialize CSQ context: {}\n", .{err});
        std.process.exit(1);
    };
    defer csq_ctx.deinit();

    // ---- Main processing loop ----
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    var n_records: u64 = 0;
    var n_errors: u64 = 0;

    while (true) {
        const has_record = reader.next(&rec) catch |err| {
            n_errors += 1;
            if (n_errors <= 10) {
                std.debug.print("Warning: failed to parse VCF record: {}\n", .{err});
            }
            continue;
        };
        if (!has_record) break;

        n_records += 1;

        // Build a pipeline-compatible record and feed it through CsqContext
        const csq_rec = buildCsqRecord(&rec);
        csq_ctx.process(&csq_rec) catch |err| {
            n_errors += 1;
            if (n_errors <= 10) {
                std.debug.print("Warning: CSQ processing error at {s}:{d}: {}\n", .{ rec.chrom, rec.pos + 1, err });
            }
        };

        // Write the record to output (pass-through for now; BCSQ injection
        // will be wired once vbufFlush emits formatted consequences).
        writeVcfRecord(out_file, &rec, &csq_ctx) catch |err| {
            std.debug.print("Error: failed to write record: {}\n", .{err});
            std.process.exit(1);
        };
    }

    // ---- Flush remaining buffered records ----
    csq_ctx.flush() catch |err| {
        std.debug.print("Warning: error flushing CSQ buffer: {}\n", .{err});
    };

    if (n_errors > 0) {
        std.debug.print("Processed {d} records with {d} warnings/errors\n", .{ n_records, n_errors });
    }
}

// -------------------------------------------------------------------------
// Main entry point
// -------------------------------------------------------------------------

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    var args = try std.process.argsWithAllocator(allocator);
    defer args.deinit();

    // Skip program name
    _ = args.skip();

    const command = args.next() orelse {
        stderr_file.writeAll(usage_text) catch {};
        std.process.exit(1);
    };

    if (std.mem.eql(u8, command, "--help") or std.mem.eql(u8, command, "-h")) {
        try stdout_file.writeAll(usage_text);
        return;
    }

    if (std.mem.eql(u8, command, "--version")) {
        try stdout_file.writeAll("bcftools-zig 0.1.0\n");
        return;
    }

    if (std.mem.eql(u8, command, "csq")) {
        try runCsq(&args);
        return;
    }

    std.debug.print("Unknown command: {s}\n", .{command});
    stderr_file.writeAll(usage_text) catch {};
    std.process.exit(1);
}
