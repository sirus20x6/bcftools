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

// Version
pub const version = "0.1.0";

// -------------------------------------------------------------------------
// C-compatible API
//
// These exports allow C code to create, drive, and destroy a CsqContext
// through an opaque pointer.  bcftools_csq_process is a documented stub:
// bridging htslib's bcf1_t to our VcfRecord requires htslib type
// definitions that are not available in this compilation unit.
// -------------------------------------------------------------------------

var g_allocator = std.heap.c_allocator;

/// Initialize a CSQ context for consequence calling.
/// Returns an opaque pointer on success, null on failure.
export fn bcftools_csq_init(
    gff_path: [*:0]const u8,
    fasta_path: [*:0]const u8,
) ?*anyopaque {
    const ctx = g_allocator.create(CsqContext) catch return null;
    ctx.* = CsqContext.init(g_allocator, .{
        .gff_fname = std.mem.span(gff_path),
        .fasta_fname = std.mem.span(fasta_path),
    }) catch {
        g_allocator.destroy(ctx);
        return null;
    };
    return @ptrCast(ctx);
}

/// Process a single VCF record through the consequence caller.
/// Returns 0 on success, -1 on error.
///
/// NOTE: This is currently a stub that always returns 0.  Fully bridging
/// htslib's bcf1_t to our VcfRecord requires htslib type definitions and
/// field accessors that are not yet wired into this compilation unit.
/// Use the Zig-native CsqContext.process() API for full functionality.
export fn bcftools_csq_process(
    ctx_ptr: ?*anyopaque,
    rec: ?*anyopaque,
) c_int {
    const ctx: *CsqContext = @ptrCast(@alignCast(ctx_ptr orelse return -1));
    _ = rec; // TODO: bridge htslib bcf1_t to our VcfRecord
    _ = ctx;
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
