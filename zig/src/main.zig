const std = @import("std");

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

const stdout = std.fs.File{ .handle = std.posix.STDOUT_FILENO };
const stderr = std.fs.File{ .handle = std.posix.STDERR_FILENO };

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    var args = try std.process.argsWithAllocator(allocator);
    defer args.deinit();

    // Skip program name
    _ = args.skip();

    const command = args.next() orelse {
        stderr.writeAll(usage_text) catch {};
        std.process.exit(1);
    };

    if (std.mem.eql(u8, command, "--help") or std.mem.eql(u8, command, "-h")) {
        try stdout.writeAll(usage_text);
        return;
    }

    if (std.mem.eql(u8, command, "--version")) {
        try stdout.writeAll("bcftools-zig 0.1.0\n");
        return;
    }

    if (std.mem.eql(u8, command, "csq")) {
        // TODO: Parse csq-specific arguments and run
        try stderr.writeAll(csq_usage_text);
        std.process.exit(1);
    }

    std.debug.print("Unknown command: {s}\n", .{command});
    stderr.writeAll(usage_text) catch {};
    std.process.exit(1);
}
