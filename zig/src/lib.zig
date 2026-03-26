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

// Version
pub const version = "0.1.0";

// C-compatible API (Phase 8 - stubs for now)

/// Initialize a CSQ context for consequence calling.
/// Returns null on failure.
export fn bcftools_csq_init(
    _gff_path: [*:0]const u8,
    _fasta_path: [*:0]const u8,
) ?*anyopaque {
    _ = _gff_path;
    _ = _fasta_path;
    // TODO: Implement in Phase 7
    return null;
}

/// Process a single VCF record through the consequence caller.
/// Returns 0 on success, -1 on error.
export fn bcftools_csq_process(
    _ctx: ?*anyopaque,
    _rec: ?*anyopaque,
) c_int {
    _ = _ctx;
    _ = _rec;
    // TODO: Implement in Phase 7
    return -1;
}

/// Flush remaining buffered records.
export fn bcftools_csq_flush(_ctx: ?*anyopaque) void {
    _ = _ctx;
    // TODO: Implement in Phase 7
}

/// Destroy a CSQ context and free all resources.
export fn bcftools_csq_destroy(_ctx: ?*anyopaque) void {
    _ = _ctx;
    // TODO: Implement in Phase 7
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
