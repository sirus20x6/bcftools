const std = @import("std");

// Module re-exports
pub const csq_types = @import("csq/types.zig");
pub const gff_types = @import("gff/types.zig");
pub const gff = @import("gff/gff.zig");
pub const translate = @import("csq/translate.zig");
pub const haplotype = @import("csq/haplotype.zig");
pub const core_allocator = @import("core/allocator.zig");
pub const core_buffer = @import("core/buffer.zig");
pub const core_region = @import("core/region.zig");
pub const core_hash = @import("core/hash.zig");
pub const csq_format = @import("csq/format.zig");
pub const csq_pipeline = @import("csq/csq.zig");
pub const vcf_record = @import("vcf/record.zig");
pub const vcf_reader = @import("vcf/reader.zig");
pub const vcf_htslib = @import("vcf/htslib.zig");

const CsqContext = csq_pipeline.CsqContext;
const VcfRecord = csq_pipeline.VcfRecord;
const htslib = vcf_htslib;

// Version
pub const version = "0.1.0";

// -------------------------------------------------------------------------
// C-compatible API
//
// These exports allow C code to create, drive, and destroy a CsqContext
// through an opaque pointer.  bcftools_csq_process bridges htslib's
// bcf1_t to our VcfRecord by unpacking alleles and resolving the
// chromosome name via the bcf_hdr_t stored during init.
// -------------------------------------------------------------------------

var g_allocator = std.heap.c_allocator;

/// Initialize a CSQ context for consequence calling.
/// @param gff_path   Path to GFF3 annotation file (null-terminated).
/// @param fasta_path Path to FASTA reference file (null-terminated).
/// @param hdr_ptr    Pointer to a bcf_hdr_t (used for rid -> chrom name lookup).
///                   May be null if chromosome resolution is not needed.
/// Returns an opaque pointer on success, null on failure.
export fn bcftools_csq_init(
    gff_path: [*:0]const u8,
    fasta_path: [*:0]const u8,
    hdr_ptr: ?*anyopaque,
) ?*anyopaque {
    const ctx = g_allocator.create(CsqContext) catch return null;
    ctx.* = CsqContext.init(g_allocator, .{
        .gff_fname = std.mem.span(gff_path),
        .fasta_fname = std.mem.span(fasta_path),
        .hdr_ptr = hdr_ptr,
    }) catch {
        g_allocator.destroy(ctx);
        return null;
    };
    return @ptrCast(ctx);
}

/// Process a single VCF record through the consequence caller.
/// Returns 0 on success, -1 on error.
///
/// Bridges htslib's bcf1_t to our VcfRecord by wrapping it in an HtsRecord,
/// extracting alleles and chromosome name, and feeding it through the
/// Zig-native CsqContext.process() pipeline.
export fn bcftools_csq_process(
    ctx_ptr: ?*anyopaque,
    rec_ptr: ?*anyopaque,
) c_int {
    const ctx: *CsqContext = @ptrCast(@alignCast(ctx_ptr orelse return -1));

    // Wrap the raw bcf1_t pointer in our HtsRecord/HtsHeader wrappers
    // which handle field access (including C bitfields) correctly.
    const hts_rec = htslib.HtsRecord{ .raw = @ptrCast(@alignCast(rec_ptr orelse return -1)) };

    // Build allele slices from the HtsRecord's unpacked allele array
    var allele_buf: [256][]const u8 = undefined;
    const n_allele = hts_rec.nAllele();
    const count = @min(n_allele, allele_buf.len);
    for (0..count) |i| {
        allele_buf[i] = hts_rec.allele(i);
    }

    // Resolve chromosome name from rid using the stored bcf_hdr_t pointer
    const chr: []const u8 = if (ctx.hdr_ptr) |hdr| blk: {
        const hdr_wrap = htslib.HtsHeader{ .raw = @ptrCast(@alignCast(hdr)) };
        break :blk hdr_wrap.seqName(hts_rec.rid());
    } else "unknown";

    const rec = VcfRecord{
        .pos = hts_rec.pos(),
        .rid = hts_rec.rid(),
        .n_allele = count,
        .alleles = allele_buf[0..count],
        .rlen = hts_rec.rlen(),
        .chr = chr,
    };

    ctx.process(&rec) catch return -1;
    return 0;
}

/// Flush remaining buffered records.
export fn bcftools_csq_flush(ctx_ptr: ?*anyopaque) void {
    const ctx: *CsqContext = @ptrCast(@alignCast(ctx_ptr orelse return));
    ctx.flush() catch {};
}

/// Destroy a CSQ context and free all resources.
export fn bcftools_csq_destroy(ctx_ptr: ?*anyopaque) void {
    const ctx: *CsqContext = @ptrCast(@alignCast(ctx_ptr orelse return));
    ctx.deinit();
    g_allocator.destroy(ctx);
}

// Pull in tests from all modules
comptime {
    _ = csq_types;
    _ = gff_types;
    _ = gff;
    _ = translate;
    _ = haplotype;
    _ = core_buffer;
    _ = core_region;
    _ = csq_format;
    _ = csq_pipeline;
    _ = vcf_record;
    _ = vcf_reader;
}
