const std = @import("std");

/// A VCF record with parsed fields.
/// Can be populated from parsed VCF lines or (eventually) from htslib via @cImport.
pub const VcfRecord = struct {
    allocator: std.mem.Allocator,

    // Core fields
    chrom: []const u8,
    pos: u32, // 0-based position
    id: []const u8,
    ref_allele: []const u8,
    alt_alleles: std.ArrayListUnmanaged([]const u8),
    qual: ?f32,
    filter: []const u8,

    // Derived
    rid: i32, // chromosome index (-1 = unset)
    rlen: u32, // reference allele length

    // Genotype cache (audit fix #10)
    gt_cache: ?[]i32,
    gt_cache_valid: bool,

    // Internal storage: we dupe strings from the parsed line so the caller
    // can reuse or free the input buffer freely.
    _storage: ?[]u8,

    pub fn init(allocator: std.mem.Allocator) VcfRecord {
        return .{
            .allocator = allocator,
            .chrom = &.{},
            .pos = 0,
            .id = &.{},
            .ref_allele = &.{},
            .alt_alleles = .empty,
            .qual = null,
            .filter = &.{},
            .rid = -1,
            .rlen = 0,
            .gt_cache = null,
            .gt_cache_valid = false,
            ._storage = null,
        };
    }

    pub fn deinit(self: *VcfRecord) void {
        self.alt_alleles.deinit(self.allocator);
        if (self.gt_cache) |cache| {
            self.allocator.free(cache);
        }
        if (self._storage) |s| {
            self.allocator.free(s);
        }
        self.* = undefined;
    }

    /// Total number of alleles (including ref).
    pub fn nAllele(self: *const VcfRecord) u32 {
        // ref + alts
        return @as(u32, @intCast(self.alt_alleles.items.len)) + 1;
    }

    /// Get allele by index (0 = ref).
    pub fn allele(self: *const VcfRecord, idx: usize) []const u8 {
        if (idx == 0) return self.ref_allele;
        return self.alt_alleles.items[idx - 1];
    }

    /// Parse a VCF line into this record (for testing without htslib).
    ///
    /// Fields parsed:
    ///   0: CHROM
    ///   1: POS  (1-based in VCF → stored 0-based)
    ///   2: ID
    ///   3: REF
    ///   4: ALT  (comma-separated)
    ///   5: QUAL
    ///   6: FILTER
    ///   7+: INFO, FORMAT, samples — skipped for now
    pub fn parseLine(self: *VcfRecord, line: []const u8) !void {
        self.clear();

        // Strip trailing newline / carriage-return if present.
        var trimmed = line;
        if (trimmed.len > 0 and trimmed[trimmed.len - 1] == '\n') trimmed = trimmed[0 .. trimmed.len - 1];
        if (trimmed.len > 0 and trimmed[trimmed.len - 1] == '\r') trimmed = trimmed[0 .. trimmed.len - 1];

        // Dupe the whole line so our slices remain valid after the caller
        // reuses or frees the input buffer.
        const storage = try self.allocator.dupe(u8, trimmed);
        self._storage = storage;

        // Split on tabs — we need at least 7 fields (indices 0..6).
        var col: usize = 0;
        var start: usize = 0;
        var fields: [8][]const u8 = undefined;

        for (storage, 0..) |c, i| {
            if (c == '\t') {
                if (col < 8) {
                    fields[col] = storage[start..i];
                }
                col += 1;
                start = i + 1;
                if (col >= 8) break;
            }
        }
        // Last field (or if fewer than 8 tabs).
        if (col < 8) {
            fields[col] = storage[start..];
            col += 1;
        }

        if (col < 7) return error.TooFewFields;

        // 0: CHROM
        self.chrom = fields[0];

        // 1: POS (1-based → 0-based)
        const pos_1based = std.fmt.parseInt(u32, fields[1], 10) catch return error.InvalidPos;
        if (pos_1based == 0) return error.InvalidPos;
        self.pos = pos_1based - 1;

        // 2: ID
        self.id = fields[2];

        // 3: REF
        self.ref_allele = fields[3];
        self.rlen = @intCast(fields[3].len);

        // 4: ALT — split on commas
        const alt_field = fields[4];
        if (alt_field.len > 0 and !std.mem.eql(u8, alt_field, ".")) {
            var alt_start: usize = 0;
            for (alt_field, 0..) |c, i| {
                if (c == ',') {
                    try self.alt_alleles.append(self.allocator, alt_field[alt_start..i]);
                    alt_start = i + 1;
                }
            }
            try self.alt_alleles.append(self.allocator, alt_field[alt_start..]);
        }

        // 5: QUAL
        if (std.mem.eql(u8, fields[5], ".")) {
            self.qual = null;
        } else {
            self.qual = std.fmt.parseFloat(f32, fields[5]) catch null;
        }

        // 6: FILTER
        self.filter = fields[6];
    }

    /// Clear record for reuse.
    pub fn clear(self: *VcfRecord) void {
        self.alt_alleles.clearRetainingCapacity();
        self.chrom = &.{};
        self.pos = 0;
        self.id = &.{};
        self.ref_allele = &.{};
        self.qual = null;
        self.filter = &.{};
        self.rid = -1;
        self.rlen = 0;
        self.invalidateGtCache();
        if (self._storage) |s| {
            self.allocator.free(s);
            self._storage = null;
        }
    }

    /// Invalidate genotype cache (call when record changes).
    pub fn invalidateGtCache(self: *VcfRecord) void {
        self.gt_cache_valid = false;
    }
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "parse VCF line" {
    const allocator = std.testing.allocator;
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t1001\t.\tA\tC,G\t.\tPASS\t.");

    try std.testing.expectEqualStrings("chr1", rec.chrom);
    try std.testing.expectEqual(@as(u32, 1000), rec.pos);
    try std.testing.expectEqualStrings(".", rec.id);
    try std.testing.expectEqualStrings("A", rec.ref_allele);
    try std.testing.expectEqual(@as(u32, 1), rec.rlen);
    try std.testing.expectEqual(@as(?f32, null), rec.qual);
    try std.testing.expectEqualStrings("PASS", rec.filter);
}

test "nAllele and allele accessors" {
    const allocator = std.testing.allocator;
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t1001\t.\tA\tC,G\t.\tPASS\t.");

    try std.testing.expectEqual(@as(u32, 3), rec.nAllele());
    try std.testing.expectEqualStrings("A", rec.allele(0));
    try std.testing.expectEqualStrings("C", rec.allele(1));
    try std.testing.expectEqualStrings("G", rec.allele(2));
}

test "parse line with no ALT" {
    const allocator = std.testing.allocator;
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr2\t500\trs123\tATG\t.\t30\t.\t.");

    try std.testing.expectEqualStrings("chr2", rec.chrom);
    try std.testing.expectEqual(@as(u32, 499), rec.pos);
    try std.testing.expectEqualStrings("rs123", rec.id);
    try std.testing.expectEqualStrings("ATG", rec.ref_allele);
    try std.testing.expectEqual(@as(u32, 3), rec.rlen);
    try std.testing.expectEqual(@as(u32, 1), rec.nAllele());
    try std.testing.expectEqual(@as(?f32, 30.0), rec.qual);
}

test "clear and reuse record" {
    const allocator = std.testing.allocator;
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    try rec.parseLine("chr1\t100\t.\tA\tT\t.\tPASS\t.");
    try std.testing.expectEqual(@as(u32, 2), rec.nAllele());

    try rec.parseLine("chrX\t200\t.\tG\tA,C,T\t50\tq10\t.");
    try std.testing.expectEqual(@as(u32, 4), rec.nAllele());
    try std.testing.expectEqualStrings("chrX", rec.chrom);
    try std.testing.expectEqual(@as(u32, 199), rec.pos);
}
