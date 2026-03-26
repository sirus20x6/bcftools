const std = @import("std");
const VcfRecord = @import("record.zig").VcfRecord;

/// Simple text VCF reader (not bgzf/bcf — those need htslib).
/// Reads a plain-text .vcf file, parsing header metadata and data records.
/// The entire file is read into memory on open.
pub const VcfReader = struct {
    allocator: std.mem.Allocator,
    content: []u8,
    remaining: []const u8,
    header_lines: std.ArrayListUnmanaged([]const u8),
    sample_names: std.ArrayListUnmanaged([]const u8),
    seq_names: std.ArrayListUnmanaged([]const u8),
    seq_map: std.StringHashMapUnmanaged(i32),

    // Heap-allocated storage for duped strings (header lines, sample/seq names).
    _header_storage: std.ArrayListUnmanaged([]u8),

    pub fn open(allocator: std.mem.Allocator, path: []const u8) !VcfReader {
        const file = try std.fs.cwd().openFile(path, .{});
        defer file.close();

        const content = try file.readToEndAlloc(allocator, 256 * 1024 * 1024);
        errdefer allocator.free(content);

        var self = VcfReader{
            .allocator = allocator,
            .content = content,
            .remaining = content,
            .header_lines = .empty,
            .sample_names = .empty,
            .seq_names = .empty,
            .seq_map = .empty,
            ._header_storage = .empty,
        };

        try self.readHeader();
        return self;
    }

    pub fn deinit(self: *VcfReader) void {
        // Free duped header line storage.
        for (self._header_storage.items) |s| {
            self.allocator.free(s);
        }
        self._header_storage.deinit(self.allocator);
        self.header_lines.deinit(self.allocator);

        // Free duped sample names.
        for (self.sample_names.items) |s| {
            self.allocator.free(s);
        }
        self.sample_names.deinit(self.allocator);

        // Free duped seq names.
        for (self.seq_names.items) |s| {
            self.allocator.free(s);
        }
        self.seq_names.deinit(self.allocator);

        self.seq_map.deinit(self.allocator);
        self.allocator.free(self.content);
        self.* = undefined;
    }

    /// Read the next record. Returns false at EOF.
    pub fn next(self: *VcfReader, rec: *VcfRecord) !bool {
        while (true) {
            const line = self.nextLine() orelse return false;

            // Skip empty lines.
            const trimmed = std.mem.trimRight(u8, line, "\r");
            if (trimmed.len == 0) continue;

            // Skip any stray comment lines in the body.
            if (trimmed[0] == '#') continue;

            try rec.parseLine(trimmed);
            self.resolveRid(rec);
            return true;
        }
    }

    /// Get number of samples.
    pub fn nSamples(self: *const VcfReader) usize {
        return self.sample_names.items.len;
    }

    // ------------------------------------------------------------------
    // Internal helpers
    // ------------------------------------------------------------------

    /// Return the next line from the remaining content, advancing the cursor.
    fn nextLine(self: *VcfReader) ?[]const u8 {
        if (self.remaining.len == 0) return null;
        if (std.mem.indexOfScalar(u8, self.remaining, '\n')) |idx| {
            const line = self.remaining[0..idx];
            self.remaining = self.remaining[idx + 1 ..];
            return line;
        }
        // Last line without trailing newline.
        const line = self.remaining;
        self.remaining = self.remaining[self.remaining.len..];
        return line;
    }

    fn readHeader(self: *VcfReader) !void {
        while (true) {
            const raw_line = self.nextLine() orelse return;
            const trimmed = std.mem.trimRight(u8, raw_line, "\r");
            if (trimmed.len == 0) continue;

            if (trimmed.len >= 2 and trimmed[0] == '#' and trimmed[1] == '#') {
                // Meta-information header line.
                const duped = try self.allocator.dupe(u8, trimmed);
                try self._header_storage.append(self.allocator, duped);
                try self.header_lines.append(self.allocator, duped);

                // Extract contig names from ##contig=<ID=...> lines.
                self.parseContigLine(trimmed) catch {};
            } else if (trimmed[0] == '#') {
                // #CHROM header line — extract sample names.
                const duped = try self.allocator.dupe(u8, trimmed);
                try self._header_storage.append(self.allocator, duped);
                try self.header_lines.append(self.allocator, duped);
                try self.parseChromLine(trimmed);
            } else {
                // First data line — rewind so next() can read it.
                // We can simply point remaining back to include this line.
                // Since raw_line is a slice of content, and remaining was advanced
                // past it, we reset remaining to start at raw_line.
                self.remaining = raw_line.ptr[0 .. raw_line.len + 1 + self.remaining.len];
                return;
            }
        }
    }

    fn parseContigLine(self: *VcfReader, line: []const u8) !void {
        // Expected format: ##contig=<ID=chr1,...>
        const prefix = "##contig=<ID=";
        if (!std.mem.startsWith(u8, line, prefix)) return;

        const rest = line[prefix.len..];
        // Find end of ID (comma or '>').
        var end: usize = 0;
        while (end < rest.len and rest[end] != ',' and rest[end] != '>') : (end += 1) {}
        if (end == 0) return;

        const name = rest[0..end];
        const idx: i32 = @intCast(self.seq_names.items.len);

        const duped = try self.allocator.dupe(u8, name);
        try self.seq_names.append(self.allocator, duped);
        try self.seq_map.put(self.allocator, duped, idx);
    }

    fn parseChromLine(self: *VcfReader, line: []const u8) !void {
        // #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO[\tFORMAT\tSample1\t...]
        // Skip the first 9 (or 8) mandatory column names, then collect sample names.
        var col: usize = 0;
        var start: usize = 0;
        for (line, 0..) |c, i| {
            if (c == '\t') {
                col += 1;
                start = i + 1;
                if (col >= 9) break;
            }
        }
        // If fewer than 9 columns, there are no samples.
        if (col < 9) return;

        // Remaining tab-separated tokens are sample names.
        var s = start;
        for (line[start..], start..) |c, i| {
            if (c == '\t') {
                if (i > s) {
                    const duped = try self.allocator.dupe(u8, line[s..i]);
                    try self.sample_names.append(self.allocator, duped);
                }
                s = i + 1;
            }
        }
        // Last sample.
        if (s < line.len) {
            const duped = try self.allocator.dupe(u8, line[s..]);
            try self.sample_names.append(self.allocator, duped);
        }
    }

    fn resolveRid(self: *VcfReader, rec: *VcfRecord) void {
        if (self.seq_map.get(rec.chrom)) |rid| {
            rec.rid = rid;
        } else {
            rec.rid = -1;
        }
    }
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "read VCF from temp file" {
    const allocator = std.testing.allocator;

    const vcf_text =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=249250621>\n" ++
        "##contig=<ID=chr2,length=243199373>\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSampleA\tSampleB\n" ++
        "chr1\t1001\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/1\t1/2\n" ++
        "chr2\t5000\trs456\tTG\tT\t30\t.\t.\tGT\t0/0\t0/1\n";

    // Write to a temp file.
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const tmp_file = try tmp_dir.dir.createFile("test.vcf", .{});
    try tmp_file.writeAll(vcf_text);
    tmp_file.close();

    // Re-open via VcfReader using the real path.
    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const tmp_path = try tmp_dir.dir.realpath("test.vcf", &path_buf);

    var reader = try VcfReader.open(allocator, tmp_path);
    defer reader.deinit();

    // Header checks.
    try std.testing.expectEqual(@as(usize, 2), reader.nSamples());
    try std.testing.expectEqualStrings("SampleA", reader.sample_names.items[0]);
    try std.testing.expectEqualStrings("SampleB", reader.sample_names.items[1]);
    try std.testing.expectEqual(@as(usize, 2), reader.seq_names.items.len);

    // Record 1.
    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    const ok1 = try reader.next(&rec);
    try std.testing.expect(ok1);
    try std.testing.expectEqualStrings("chr1", rec.chrom);
    try std.testing.expectEqual(@as(u32, 1000), rec.pos);
    try std.testing.expectEqual(@as(u32, 3), rec.nAllele());
    try std.testing.expectEqual(@as(i32, 0), rec.rid); // chr1 -> index 0

    // Record 2.
    const ok2 = try reader.next(&rec);
    try std.testing.expect(ok2);
    try std.testing.expectEqualStrings("chr2", rec.chrom);
    try std.testing.expectEqual(@as(u32, 4999), rec.pos);
    try std.testing.expectEqual(@as(u32, 2), rec.nAllele());
    try std.testing.expectEqualStrings("TG", rec.ref_allele);
    try std.testing.expectEqual(@as(i32, 1), rec.rid); // chr2 -> index 1

    // EOF.
    const ok3 = try reader.next(&rec);
    try std.testing.expect(!ok3);
}

test "VCF with no samples" {
    const allocator = std.testing.allocator;

    const vcf_text =
        "##fileformat=VCFv4.2\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" ++
        "chr1\t100\t.\tA\tT\t.\t.\t.\n";

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const tmp_file = try tmp_dir.dir.createFile("nosamp.vcf", .{});
    try tmp_file.writeAll(vcf_text);
    tmp_file.close();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const tmp_path = try tmp_dir.dir.realpath("nosamp.vcf", &path_buf);

    var reader = try VcfReader.open(allocator, tmp_path);
    defer reader.deinit();

    try std.testing.expectEqual(@as(usize, 0), reader.nSamples());

    var rec = VcfRecord.init(allocator);
    defer rec.deinit();

    const ok = try reader.next(&rec);
    try std.testing.expect(ok);
    try std.testing.expectEqual(@as(u32, 99), rec.pos);
    try std.testing.expectEqual(@as(u32, 2), rec.nAllele());
}
