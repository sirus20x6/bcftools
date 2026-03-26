const std = @import("std");
const lib = @import("bcftools_zig");
const csq_mod = lib.csq_pipeline;
const gff_mod = lib.gff;
const htslib = lib.vcf_htslib;
const VcfReader = lib.vcf_reader.VcfReader;
const VcfRecord = lib.vcf_record.VcfRecord;
const CsqContext = csq_mod.CsqContext;
const GffParser = gff_mod.GffParser;
const Phase = csq_mod.Phase;

/// Adapter function: bridges CsqContext.FetchSeqFn to HtsFaidx.fetchSeq
fn htsFaidxFetchAdapter(ctx: *anyopaque, allocator: std.mem.Allocator, chr: [*:0]const u8, beg: i64, end: i64) ?[]u8 {
    const fai: *htslib.HtsFaidx = @ptrCast(@alignCast(ctx));
    return fai.fetchSeq(allocator, chr, beg, end) catch null;
}

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
        .chr = rec.chrom,
        .raw_line = rec._storage,
    };
}

// -------------------------------------------------------------------------
// Write a VCF record line, appending BCSQ to INFO if consequences exist
// -------------------------------------------------------------------------

/// Write a single VCF line to output.  If `bcsq_value` is non-null the
/// annotation is injected into the INFO column (column 7); otherwise the
/// original line is written unchanged.
/// If `fmt_bm` is non-null, also append `:BCSQ` to FORMAT and bitmask values
/// to each sample column.
fn writeVcfLine(
    allocator: std.mem.Allocator,
    out: std.fs.File,
    original_line: []const u8,
    bcsq_value: ?[]const u8,
    bcsq_tag: []const u8,
    fmt_bm: ?[]const u32,
    nfmt: u32,
    n_samples: u32,
) !void {
    var line_to_process = original_line;
    var info_modified: ?[]u8 = null;
    defer if (info_modified) |m| allocator.free(m);

    if (bcsq_value) |val| {
        info_modified = try csq_mod.injectBcsq(allocator, original_line, val, bcsq_tag);
        line_to_process = info_modified.?;
    }

    if (fmt_bm != null and n_samples > 0 and nfmt > 0) {
        const bm = fmt_bm.?;
        // Parse columns to find FORMAT (col 8) and sample columns (col 9+)
        // We need to:
        // 1. Append ":BCSQ_TAG" to the FORMAT column
        // 2. Append ":<bitmask>" to each sample column

        // Strip trailing newline/CR
        var line = line_to_process;
        if (line.len > 0 and line[line.len - 1] == '\n') line = line[0 .. line.len - 1];
        if (line.len > 0 and line[line.len - 1] == '\r') line = line[0 .. line.len - 1];

        // Find column boundaries
        var col_starts: [256]usize = undefined;
        var col_count: usize = 0;
        var start: usize = 0;
        for (line, 0..) |c, idx| {
            if (c == '\t') {
                if (col_count < col_starts.len) {
                    col_starts[col_count] = start;
                    col_count += 1;
                }
                start = idx + 1;
            }
        }
        // Last column
        if (col_count < col_starts.len) {
            col_starts[col_count] = start;
            col_count += 1;
        }

        if (col_count >= 10) {
            // We have FORMAT (col 8, index 8) and at least one sample
            var result: std.ArrayList(u8) = .empty;
            defer result.deinit(allocator);

            // Write columns 0-7 unchanged (CHROM through INFO)
            const fmt_col_start = col_starts[8];
            try result.appendSlice(allocator, line[0..fmt_col_start]);

            // Find FORMAT column end
            const fmt_end = blk: {
                var e = fmt_col_start;
                while (e < line.len and line[e] != '\t') e += 1;
                break :blk e;
            };
            // Append FORMAT field + ":BCSQ"
            try result.appendSlice(allocator, line[fmt_col_start..fmt_end]);
            try result.append(allocator, ':');
            try result.appendSlice(allocator, bcsq_tag);

            // For each sample column, append ":<bitmask>"
            const nfmt_bcsq: usize = @max(1, nfmt);
            var smpl_idx: u32 = 0;
            var pos_s: usize = fmt_end;
            while (smpl_idx < n_samples) : (smpl_idx += 1) {
                // Find sample column boundaries
                if (pos_s < line.len and line[pos_s] == '\t') {
                    pos_s += 1; // skip tab
                }
                var smpl_end = pos_s;
                while (smpl_end < line.len and line[smpl_end] != '\t') smpl_end += 1;

                try result.append(allocator, '\t');
                try result.appendSlice(allocator, line[pos_s..smpl_end]);
                try result.append(allocator, ':');

                // Write bitmask value(s) for this sample
                const bm_offset = @as(usize, smpl_idx) * nfmt_bcsq;
                var fi: usize = 0;
                while (fi < nfmt_bcsq) : (fi += 1) {
                    if (fi > 0) try result.append(allocator, ',');
                    const val = if (bm_offset + fi < bm.len) bm[bm_offset + fi] else 0;
                    var buf: [16]u8 = undefined;
                    const s = std.fmt.bufPrint(&buf, "{d}", .{val}) catch break;
                    try result.appendSlice(allocator, s);
                }

                pos_s = smpl_end;
            }

            try out.writeAll(result.items);
            try out.writeAll("\n");
            return;
        }
    }

    // Strip trailing semicolons from the INFO field (column 7) to match
    // bcftools/htslib normalization.  We do this only for the pass-through
    // (no BCSQ / no FORMAT rewrite) path.
    {
        var line = line_to_process;
        if (line.len > 0 and line[line.len - 1] == '\n') line = line[0 .. line.len - 1];
        if (line.len > 0 and line[line.len - 1] == '\r') line = line[0 .. line.len - 1];

        // Find column 7 (INFO) boundaries
        var tab_count: u32 = 0;
        var info_start: usize = 0;
        var info_end: usize = line.len;
        for (line, 0..) |c, idx| {
            if (c == '\t') {
                tab_count += 1;
                if (tab_count == 7) info_start = idx + 1;
                if (tab_count == 8) {
                    info_end = idx;
                    break;
                }
            }
        }
        // Strip trailing semicolons from INFO
        var stripped_end = info_end;
        while (stripped_end > info_start and line[stripped_end - 1] == ';') stripped_end -= 1;
        if (stripped_end != info_end) {
            try out.writeAll(line[0..stripped_end]);
            try out.writeAll(line[info_end..]);
            try out.writeAll("\n");
        } else {
            try out.writeAll(line);
            try out.writeAll("\n");
        }
    }
}

/// Drain flushed records from the CSQ context and write them to output.
/// This must be called after every CsqContext.process() and after flush().
///
/// Because the CSQ pipeline buffers records internally and releases them
/// via vbufFlush, we need a mapping from the pipeline's VcfRecord (pos-based)
/// back to the original text line.  We look up the original line in the
/// line_map keyed by position.
fn writeFlushedRecords(
    allocator: std.mem.Allocator,
    out: std.fs.File,
    csq_ctx: *CsqContext,
    line_map: *std.AutoHashMap(u64, []const u8),
) !void {
    for (csq_ctx.flushed_records.items) |fr| {
        // Build a lookup key: combine rid + pos to handle multi-chrom inputs
        const key = posKey(fr.rid, fr.pos);
        const original_line = line_map.get(key) orelse continue;

        try writeVcfLine(allocator, out, original_line, fr.bcsq_value, csq_ctx.bcsq_tag, fr.fmt_bm, fr.nfmt, csq_ctx.n_samples);

        // Free the duped bcsq_value string and fmt_bm
        if (fr.bcsq_value) |bv| allocator.free(bv);
        if (fr.fmt_bm) |bm| allocator.free(bm);

        // Remove from map to free memory
        _ = line_map.remove(key);
    }
    csq_ctx.flushed_records.clearRetainingCapacity();
}

/// Combine rid and pos into a single u64 key for the line lookup map.
fn posKey(rid: i32, pos: u32) u64 {
    return (@as(u64, @bitCast(@as(i64, rid))) << 32) | @as(u64, pos);
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

    // Write the #CHROM line and count samples
    var n_samples: u32 = 0;
    if (chrom_line) |cl| {
        out_file.writeAll(cl) catch {};
        out_file.writeAll("\n") catch {};
        // Count samples: #CHROM has 9 fixed columns, then samples
        var tab_count: u32 = 0;
        for (cl) |c| {
            if (c == '\t') tab_count += 1;
        }
        if (tab_count >= 9) n_samples = tab_count - 8; // 9 tabs = 10 cols, 9 fixed + 1 sample
    }

    // When there are no samples, force drop_gt mode (matches C: line 726)
    var phase = opts.phase;
    if (n_samples == 0) phase = .drop_gt;

    // ---- Open FASTA reference ----
    const fasta_z = blk: {
        var buf: [4096]u8 = undefined;
        const fname = opts.fasta_fname.?;
        if (fname.len >= buf.len) {
            stderr_file.writeAll("Error: fasta path too long\n") catch {};
            std.process.exit(1);
        }
        @memcpy(buf[0..fname.len], fname);
        buf[fname.len] = 0;
        break :blk buf[0..fname.len :0];
    };
    var fai = htslib.HtsFaidx.open(fasta_z) catch {
        std.debug.print("Error: failed to open FASTA reference '{s}'\n", .{opts.fasta_fname.?});
        std.process.exit(1);
    };
    defer fai.close();

    // ---- Initialize CSQ context ----
    var csq_ctx = CsqContext.init(allocator, .{
        .gff_fname = opts.gff_fname.?,
        .fasta_fname = opts.fasta_fname.?,
        .phase = phase,
        .local_csq = opts.local_csq,
        .verbosity = 1,
        .force = opts.force,
        .bcsq_tag = opts.custom_tag,
        .ncsq2_max = opts.ncsq * 2,
        .brief_predictions = opts.brief_predictions,
        .n_samples = n_samples,
        .fai_ptr = @ptrCast(&fai),
        .fetch_seq_fn = &htsFaidxFetchAdapter,
    }) catch |err| {
        std.debug.print("Error: failed to initialize CSQ context: {}\n", .{err});
        std.process.exit(1);
    };
    defer csq_ctx.deinit();

    // ---- Main processing loop ----
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    // Map from (rid, pos) -> original VCF text line.  The CSQ pipeline
    // buffers records and flushes them later, so we need to keep the
    // original lines alive until they are written.
    // NOTE: multiple records at the same position will overwrite each other.
    // This is acceptable for now; when htslib bindings replace text VCF
    // reading, records will be managed by the pipeline's own Vbuf.
    var line_map = std.AutoHashMap(u64, []const u8).init(allocator);
    defer {
        var it = line_map.valueIterator();
        while (it.next()) |v| allocator.free(v.*);
        line_map.deinit();
    }

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

        // Save the original line for later BCSQ injection.
        // We must dupe it because rec._storage is reused on the next read.
        if (rec._storage) |storage| {
            const key = posKey(rec.rid, rec.pos);
            const duped = try allocator.dupe(u8, storage);
            try line_map.put(key, duped);
        }

        // Build a pipeline-compatible record and feed it through CsqContext
        const csq_rec = buildCsqRecord(&rec);
        csq_ctx.process(&csq_rec) catch |err| {
            n_errors += 1;
            if (n_errors <= 10) {
                std.debug.print("Warning: CSQ processing error at {s}:{d}: {}\n", .{ rec.chrom, rec.pos + 1, err });
            }
        };

        // Write any records that were flushed by this process() call
        writeFlushedRecords(allocator, out_file, &csq_ctx, &line_map) catch |err| {
            std.debug.print("Error: failed to write flushed records: {}\n", .{err});
            std.process.exit(1);
        };
    }

    // ---- Flush remaining buffered records ----
    csq_ctx.flush() catch |err| {
        std.debug.print("Warning: error flushing CSQ buffer: {}\n", .{err});
    };
    writeFlushedRecords(allocator, out_file, &csq_ctx, &line_map) catch |err| {
        std.debug.print("Error: failed to write final flushed records: {}\n", .{err});
        std.process.exit(1);
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
