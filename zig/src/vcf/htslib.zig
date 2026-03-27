const std = @import("std");

pub const c = @cImport({
    @cInclude("htslib/hts.h");
    @cInclude("htslib/vcf.h");
    @cInclude("htslib/synced_bcf_reader.h");
    @cInclude("htslib/faidx.h");
    @cInclude("htslib/kstring.h");
    @cInclude("vcf/bcf_compat.h");
});

pub const HtsHeader = struct {
    raw: *c.bcf_hdr_t,

    pub fn nSamples(self: *const HtsHeader) u32 {
        return @intCast(c.bcf_hdr_nsamples(self.raw));
    }

    pub fn seqName(self: *const HtsHeader, rid: i32) []const u8 {
        const name = c.bcf_hdr_id2name(self.raw, rid);
        if (name == null) return "";
        return std.mem.span(name);
    }

    pub fn appendVersion(self: *HtsHeader, argc: c_int, argv: [*]const [*:0]const u8, cmd: [*:0]const u8) void {
        _ = c.bcf_hdr_append_version(self.raw, argc, argv, cmd);
    }
};

pub const HtsRecord = struct {
    raw: *c.bcf1_t,

    pub fn pos(self: *const HtsRecord) u32 {
        return @intCast(c.bcf_compat_pos(self.raw));
    }

    pub fn rid(self: *const HtsRecord) i32 {
        return @intCast(c.bcf_compat_rid(self.raw));
    }

    pub fn nAllele(self: *const HtsRecord) u32 {
        return @intCast(c.bcf_compat_n_allele(self.raw));
    }

    pub fn rlen(self: *const HtsRecord) u32 {
        return @intCast(c.bcf_compat_rlen(self.raw));
    }

    pub fn allele(self: *const HtsRecord, idx: usize) []const u8 {
        // Need to unpack first
        _ = c.bcf_unpack(self.raw, c.BCF_UN_STR);
        const alleles = c.bcf_compat_alleles(self.raw);
        if (alleles == null) return "";
        return std.mem.span(alleles[idx]);
    }

    pub fn seqname(self: *const HtsRecord, hdr: *const HtsHeader) []const u8 {
        return hdr.seqName(@intCast(c.bcf_compat_rid(self.raw)));
    }
};

pub const HtsVcfReader = struct {
    sr: *c.bcf_srs_t,
    hdr: HtsHeader,

    pub fn open(fname: [*:0]const u8) !HtsVcfReader {
        const sr = c.bcf_sr_init() orelse return error.HtsOpenFailed;
        if (c.bcf_sr_add_reader(sr, fname) != 1) {
            c.bcf_sr_destroy(sr);
            return error.HtsOpenFailed;
        }
        return .{
            .sr = sr,
            .hdr = .{ .raw = sr.readers[0].header },
        };
    }

    pub fn next(self: *HtsVcfReader) ?HtsRecord {
        if (c.bcf_sr_next_line(self.sr) == 0) return null;
        const rec = c.bcf_sr_get_line(self.sr, 0);
        if (rec == null) return null;
        return .{ .raw = rec.? };
    }

    pub fn close(self: *HtsVcfReader) void {
        c.bcf_sr_destroy(self.sr);
    }

    pub fn setRegions(self: *HtsVcfReader, regions: [*:0]const u8, is_file: bool) !void {
        if (c.bcf_sr_set_regions(self.sr, regions, if (is_file) @as(c_int, 1) else 0) != 0)
            return error.HtsRegionFailed;
    }

    pub fn setTargets(self: *HtsVcfReader, targets: [*:0]const u8, is_file: bool) !void {
        if (c.bcf_sr_set_targets(self.sr, targets, if (is_file) @as(c_int, 1) else 0, 0) != 0)
            return error.HtsTargetFailed;
    }
};

pub const HtsVcfWriter = struct {
    fp: *c.htsFile,
    hdr: *c.bcf_hdr_t,

    pub fn open(fname: [*:0]const u8, mode: [*:0]const u8) !HtsVcfWriter {
        const fp = c.hts_open(fname, mode) orelse return error.HtsOpenFailed;
        return .{ .fp = fp, .hdr = undefined };
    }

    pub fn writeHeader(self: *HtsVcfWriter, hdr: *c.bcf_hdr_t) !void {
        self.hdr = hdr;
        if (c.bcf_hdr_write(self.fp, hdr) != 0) return error.HtsWriteFailed;
    }

    pub fn writeRecord(self: *HtsVcfWriter, rec: *c.bcf1_t) !void {
        if (c.bcf_write(self.fp, self.hdr, rec) != 0) return error.HtsWriteFailed;
    }

    pub fn close(self: *HtsVcfWriter) void {
        _ = c.hts_close(self.fp);
    }
};

pub const HtsFaidx = struct {
    fai: *c.faidx_t,

    pub fn open(fname: [*:0]const u8) !HtsFaidx {
        const fai = c.fai_load(fname) orelse return error.FaidxOpenFailed;
        return .{ .fai = fai };
    }

    pub fn close(self: *HtsFaidx) void {
        c.fai_destroy(self.fai);
    }

    pub fn hasSeq(self: *const HtsFaidx, seq: [*:0]const u8) bool {
        return c.faidx_has_seq(self.fai, seq) != 0;
    }

    /// Fetch a sequence region. Caller owns the returned slice.
    pub fn fetchSeq(self: *const HtsFaidx, alloc: std.mem.Allocator, seq: [*:0]const u8, beg: i64, end: i64) ![]u8 {
        var len: i64 = 0;
        const raw = c.faidx_fetch_seq64(self.fai, seq, beg, end, &len);
        if (raw == null) return error.FaidxFetchFailed;
        defer std.c.free(raw);
        const slice = raw[0..@intCast(len)];
        const result = try alloc.alloc(u8, slice.len);
        @memcpy(result, slice);
        return result;
    }
};

pub const BcfUpdateError = error{HtsUpdateFailed};

pub fn bcfUpdateInfoString(hdr: *c.bcf_hdr_t, rec: *c.bcf1_t, tag: [*:0]const u8, val: [*:0]const u8) !void {
    if (c.bcf_update_info_string(hdr, rec, tag, val) != 0)
        return error.HtsUpdateFailed;
}

pub fn bcfUpdateFormatInt32(hdr: *c.bcf_hdr_t, rec: *c.bcf1_t, tag: [*:0]const u8, vals: []const i32) !void {
    if (c.bcf_update_format_int32(hdr, rec, tag, vals.ptr, @intCast(vals.len)) != 0)
        return error.HtsUpdateFailed;
}

pub fn bcfGetGenotypes(hdr: *c.bcf_hdr_t, rec: *c.bcf1_t, alloc: std.mem.Allocator) ![]i32 {
    var gt: [*c]i32 = null;
    var ngt: c_int = 0;
    const ret = c.bcf_get_genotypes(hdr, rec, @ptrCast(&gt), &ngt);
    if (ret <= 0) return error.HtsUpdateFailed;
    defer std.c.free(gt);
    const slice = gt[0..@intCast(ret)];
    const result = try alloc.alloc(i32, slice.len);
    @memcpy(result, slice);
    return result;
}
