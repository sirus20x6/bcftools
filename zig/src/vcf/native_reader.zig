const std = @import("std");
const simd = @import("../core/simd.zig");

/// A high-performance pure-Zig VCF reader that does not depend on htslib.
///
/// Design principles:
///   - Buffered I/O with 64KB read chunks
///   - Zero-copy parsing: field accessors return slices into the read buffer
///   - Lazy field parsing: tab positions are only computed when a field is accessed
///   - SIMD-accelerated tab finding using @Vector(16, u8)
///
/// The record data returned by field accessors is valid only until the next
/// call to `next()`.
pub const NativeReader = struct {
    allocator: std.mem.Allocator,
    source: Source,
    buffer: []u8,
    buf_len: usize,
    buf_pos: usize,
    eof: bool,

    // Header
    header: Header,

    // Current record (zero-copy slice into buffer)
    current_line: []const u8,

    // Field position cache (lazy-computed tab positions)
    // field_offsets[i] stores the byte offset of the (i+1)-th tab character,
    // so field i spans from (field_offsets[i-1]+1) .. field_offsets[i].
    // field_offsets[0] is the offset of the first tab (end of CHROM).
    field_offsets: [MAX_FIELDS]u32,
    fields_found: u8,
    fields_parsed: bool,

    const BUFFER_SIZE: usize = 64 * 1024;
    const MAX_FIELDS: usize = 10;

    const Source = union(enum) {
        file: std.fs.File,
        memory: struct {
            data: []const u8,
            pos: usize,
        },
    };

    /// Open a VCF file by path.
    pub fn open(allocator: std.mem.Allocator, path: []const u8) !NativeReader {
        const file = try std.fs.cwd().openFile(path, .{});
        errdefer file.close();

        const buffer = try allocator.alloc(u8, BUFFER_SIZE);
        errdefer allocator.free(buffer);

        var self = NativeReader{
            .allocator = allocator,
            .source = .{ .file = file },
            .buffer = buffer,
            .buf_len = 0,
            .buf_pos = 0,
            .eof = false,
            .header = Header.init(allocator),
            .current_line = &.{},
            .field_offsets = undefined,
            .fields_found = 0,
            .fields_parsed = false,
        };

        try self.readHeader();
        return self;
    }

    /// Create a NativeReader from an in-memory buffer (useful for testing).
    pub fn fromMemory(allocator: std.mem.Allocator, data: []const u8) !NativeReader {
        // We still need the read buffer for line assembly
        const buffer = try allocator.alloc(u8, BUFFER_SIZE);
        errdefer allocator.free(buffer);

        var self = NativeReader{
            .allocator = allocator,
            .source = .{ .memory = .{ .data = data, .pos = 0 } },
            .buffer = buffer,
            .buf_len = 0,
            .buf_pos = 0,
            .eof = false,
            .header = Header.init(allocator),
            .current_line = &.{},
            .field_offsets = undefined,
            .fields_found = 0,
            .fields_parsed = false,
        };

        try self.readHeader();
        return self;
    }

    pub fn deinit(self: *NativeReader) void {
        switch (self.source) {
            .file => |f| f.close(),
            .memory => {},
        }
        self.header.deinit();
        self.allocator.free(self.buffer);
        self.* = undefined;
    }

    /// Read the next record. Returns false at EOF.
    /// The record data (accessible via field methods) is valid until the next
    /// call to next().
    pub fn next(self: *NativeReader) !bool {
        self.fields_parsed = false;
        self.fields_found = 0;

        while (true) {
            const line = (try self.readLine()) orelse return false;

            // Skip empty lines
            if (line.len == 0) continue;

            // Skip comment lines (shouldn't appear after header, but be safe)
            if (line[0] == '#') continue;

            self.current_line = line;
            return true;
        }
    }

    // -----------------------------------------------------------------------
    // Zero-copy field accessors
    // -----------------------------------------------------------------------

    /// CHROM field (column 0)
    pub fn chrom(self: *NativeReader) []const u8 {
        return self.getField(0);
    }

    /// POS field (column 1), 0-based. VCF stores 1-based, we subtract 1.
    pub fn pos(self: *NativeReader) u32 {
        const field = self.getField(1);
        const val = std.fmt.parseInt(u32, field, 10) catch return 0;
        if (val == 0) return 0;
        return val - 1;
    }

    /// ID field (column 2)
    pub fn id(self: *NativeReader) []const u8 {
        return self.getField(2);
    }

    /// REF allele (column 3)
    pub fn refAllele(self: *NativeReader) []const u8 {
        return self.getField(3);
    }

    /// Iterator over ALT alleles (column 4), split by comma.
    pub fn altAlleles(self: *NativeReader) AltIterator {
        const field = self.getField(4);
        // "." means no alt alleles
        if (field.len == 1 and field[0] == '.') {
            return .{ .data = &.{}, .pos = 0 };
        }
        return .{ .data = field, .pos = 0 };
    }

    /// QUAL field (column 5). Returns null if ".".
    pub fn qual(self: *NativeReader) ?f32 {
        const field = self.getField(5);
        if (field.len == 1 and field[0] == '.') return null;
        return std.fmt.parseFloat(f32, field) catch null;
    }

    /// FILTER field (column 6)
    pub fn filter(self: *NativeReader) []const u8 {
        return self.getField(6);
    }

    /// INFO field (column 7)
    pub fn info(self: *NativeReader) []const u8 {
        return self.getField(7);
    }

    /// FORMAT field (column 8). May be empty if no samples.
    pub fn format(self: *NativeReader) []const u8 {
        return self.getField(8);
    }

    // -----------------------------------------------------------------------
    // Lazy tab-position parsing with SIMD acceleration
    // -----------------------------------------------------------------------

    fn getField(self: *NativeReader, col: usize) []const u8 {
        if (!self.fields_parsed) self.parseFields();
        if (col >= self.fields_found) return &.{};

        const start: usize = if (col == 0) 0 else @as(usize, self.field_offsets[col - 1]) + 1;
        const end: usize = @as(usize, self.field_offsets[col]);
        return self.current_line[start..end];
    }

    fn parseFields(self: *NativeReader) void {
        if (self.fields_parsed) return;

        var col: u8 = 0;
        var i: usize = 0;
        const line = self.current_line;

        // SIMD: scan 16 bytes at a time looking for tabs
        while (i + 16 <= line.len and col < MAX_FIELDS) {
            const chunk: @Vector(16, u8) = line[i..][0..16].*;
            const tab_splat: @Vector(16, u8) = @splat(@as(u8, '\t'));
            const matches = chunk == tab_splat;
            const mask: u16 = @bitCast(matches);

            if (mask != 0) {
                // Process all tab positions in this chunk
                var remaining_mask = mask;
                while (remaining_mask != 0 and col < MAX_FIELDS) {
                    const bit_pos = @ctz(remaining_mask);
                    self.field_offsets[col] = @intCast(i + bit_pos);
                    col += 1;
                    remaining_mask &= remaining_mask - 1; // clear lowest set bit
                }
                // Advance past the last tab we found in this chunk
                if (col < MAX_FIELDS) {
                    i += 16;
                } else {
                    break;
                }
            } else {
                i += 16;
            }
        }

        // Scalar fallback for remaining bytes
        while (i < line.len and col < MAX_FIELDS) {
            if (line[i] == '\t') {
                self.field_offsets[col] = @intCast(i);
                col += 1;
            }
            i += 1;
        }

        // The last field ends at the end of the line
        if (col < MAX_FIELDS) {
            self.field_offsets[col] = @intCast(line.len);
            col += 1;
        }

        self.fields_found = col;
        self.fields_parsed = true;
    }

    // -----------------------------------------------------------------------
    // Buffered line reading
    // -----------------------------------------------------------------------

    /// Read the next newline-delimited line from the source.
    /// Returns a slice into self.buffer, or null at EOF.
    fn readLine(self: *NativeReader) !?[]const u8 {
        // Strategy: find newline in current buffer. If not found, compact
        // unconsumed data to front of buffer, read more, try again.
        var line_start = self.buf_pos;

        while (true) {
            // Search for newline in buffered data
            const search_region = self.buffer[line_start..self.buf_len];
            if (simd.findByte(search_region, '\n')) |rel_idx| {
                const newline_pos = line_start + rel_idx;
                var end = newline_pos;
                // Strip \r before \n
                if (end > self.buf_pos and self.buffer[end - 1] == '\r') {
                    end -= 1;
                }
                const line = self.buffer[self.buf_pos..end];
                self.buf_pos = newline_pos + 1;
                return line;
            }

            // No newline found. If we've hit EOF, return remaining data.
            if (self.eof) {
                if (self.buf_pos < self.buf_len) {
                    const line = self.buffer[self.buf_pos..self.buf_len];
                    self.buf_pos = self.buf_len;
                    // Strip trailing \r
                    const trimmed = if (line.len > 0 and line[line.len - 1] == '\r')
                        line[0 .. line.len - 1]
                    else
                        line;
                    if (trimmed.len == 0) return null;
                    return trimmed;
                }
                return null;
            }

            // Compact: move unconsumed data to front of buffer
            const unconsumed = self.buf_len - self.buf_pos;
            if (unconsumed > 0 and self.buf_pos > 0) {
                std.mem.copyForwards(u8, self.buffer[0..unconsumed], self.buffer[self.buf_pos..self.buf_len]);
            }
            self.buf_len = unconsumed;
            self.buf_pos = 0;
            line_start = unconsumed;

            // Read more data
            const space = self.buffer[self.buf_len..];
            if (space.len == 0) {
                // Line is longer than buffer -- this shouldn't happen with
                // reasonable VCF files and a 64KB buffer. Return what we have.
                const line = self.buffer[0..self.buf_len];
                self.buf_len = 0;
                self.buf_pos = 0;
                return line;
            }

            const n = try self.readFromSource(space);
            if (n == 0) {
                self.eof = true;
                // Loop back to handle remaining data
                continue;
            }
            self.buf_len += n;
        }
    }

    fn readFromSource(self: *NativeReader, dest: []u8) !usize {
        switch (self.source) {
            .file => |f| {
                return f.read(dest);
            },
            .memory => |*m| {
                const remaining = m.data[m.pos..];
                const n = @min(remaining.len, dest.len);
                @memcpy(dest[0..n], remaining[0..n]);
                m.pos += n;
                return n;
            },
        }
    }

    // -----------------------------------------------------------------------
    // Header parsing
    // -----------------------------------------------------------------------

    fn readHeader(self: *NativeReader) !void {
        while (true) {
            const line = (try self.readLine()) orelse return;

            if (line.len == 0) continue;

            if (line.len >= 2 and line[0] == '#' and line[1] == '#') {
                // Meta-information line
                const duped = try self.allocator.dupe(u8, line);
                try self.header._storage.append(self.allocator, duped);
                try self.header.raw_lines.append(self.allocator, duped);

                // Parse ##contig=<ID=...>
                self.header.parseContigLine(line) catch {};

                // Parse ##INFO=<ID=...>
                self.header.parseFieldLine(line, "##INFO=<", &self.header.info_fields) catch {};

                // Parse ##FORMAT=<ID=...>
                self.header.parseFieldLine(line, "##FORMAT=<", &self.header.format_fields) catch {};
            } else if (line[0] == '#') {
                // #CHROM header line
                const duped = try self.allocator.dupe(u8, line);
                try self.header._storage.append(self.allocator, duped);
                try self.header.raw_lines.append(self.allocator, duped);
                try self.header.parseChromLine(line);
            } else {
                // First data line -- we need to "unread" it.
                // Since readLine advanced buf_pos past the newline,
                // we rewind buf_pos to the start of this line.
                // The line is a slice into self.buffer, so:
                const line_start = @intFromPtr(line.ptr) - @intFromPtr(self.buffer.ptr);
                self.buf_pos = line_start;
                return;
            }
        }
    }
};

/// Parsed VCF header information.
pub const Header = struct {
    allocator: std.mem.Allocator,
    sample_names: std.ArrayListUnmanaged([]const u8),
    contig_map: std.StringHashMapUnmanaged(i32),
    contig_names: std.ArrayListUnmanaged([]const u8),
    info_fields: std.StringHashMapUnmanaged(HeaderField),
    format_fields: std.StringHashMapUnmanaged(HeaderField),
    raw_lines: std.ArrayListUnmanaged([]const u8),

    // All duped strings for cleanup
    _storage: std.ArrayListUnmanaged([]u8),

    pub fn init(allocator: std.mem.Allocator) Header {
        return .{
            .allocator = allocator,
            .sample_names = .empty,
            .contig_map = .empty,
            .contig_names = .empty,
            .info_fields = .empty,
            .format_fields = .empty,
            .raw_lines = .empty,
            ._storage = .empty,
        };
    }

    pub fn deinit(self: *Header) void {
        for (self._storage.items) |s| self.allocator.free(s);
        self._storage.deinit(self.allocator);

        for (self.sample_names.items) |s| self.allocator.free(s);
        self.sample_names.deinit(self.allocator);

        for (self.contig_names.items) |s| self.allocator.free(s);
        self.contig_names.deinit(self.allocator);

        // raw_lines storage was already tracked in _storage
        self.raw_lines.deinit(self.allocator);
        self.contig_map.deinit(self.allocator);
        self.info_fields.deinit(self.allocator);
        self.format_fields.deinit(self.allocator);
    }

    pub fn nSamples(self: *const Header) usize {
        return self.sample_names.items.len;
    }

    fn parseContigLine(self: *Header, line: []const u8) !void {
        const prefix = "##contig=<ID=";
        if (!std.mem.startsWith(u8, line, prefix)) return;
        const rest = line[prefix.len..];
        var end: usize = 0;
        while (end < rest.len and rest[end] != ',' and rest[end] != '>') : (end += 1) {}
        if (end == 0) return;

        const name = rest[0..end];
        const idx: i32 = @intCast(self.contig_names.items.len);
        const duped = try self.allocator.dupe(u8, name);
        try self.contig_names.append(self.allocator, duped);
        try self.contig_map.put(self.allocator, duped, idx);
    }

    fn parseFieldLine(
        self: *Header,
        line: []const u8,
        prefix: []const u8,
        map: *std.StringHashMapUnmanaged(HeaderField),
    ) !void {
        if (!std.mem.startsWith(u8, line, prefix)) return;
        const rest = line[prefix.len..];

        // Extract ID=...
        if (!std.mem.startsWith(u8, rest, "ID=")) return;
        const id_start = 3; // skip "ID="
        var id_end: usize = id_start;
        while (id_end < rest.len and rest[id_end] != ',' and rest[id_end] != '>') : (id_end += 1) {}
        const id_str = rest[id_start..id_end];
        if (id_str.len == 0) return;

        // Extract Number=...
        var number: HeaderNumber = .{ .fixed = 1 };
        if (std.mem.indexOf(u8, rest, "Number=")) |num_pos| {
            const num_start = num_pos + 7;
            var num_end = num_start;
            while (num_end < rest.len and rest[num_end] != ',' and rest[num_end] != '>') : (num_end += 1) {}
            const num_str = rest[num_start..num_end];
            if (std.mem.eql(u8, num_str, ".")) {
                number = .variable;
            } else if (std.mem.eql(u8, num_str, "A")) {
                number = .per_alt;
            } else if (std.mem.eql(u8, num_str, "R")) {
                number = .per_allele;
            } else if (std.mem.eql(u8, num_str, "G")) {
                number = .per_genotype;
            } else {
                number = .{ .fixed = std.fmt.parseInt(u32, num_str, 10) catch 1 };
            }
        }

        // Extract Type=...
        var field_type: HeaderFieldType = .string;
        if (std.mem.indexOf(u8, rest, "Type=")) |type_pos| {
            const type_start = type_pos + 5;
            var type_end = type_start;
            while (type_end < rest.len and rest[type_end] != ',' and rest[type_end] != '>') : (type_end += 1) {}
            const type_str = rest[type_start..type_end];
            if (std.mem.eql(u8, type_str, "Integer")) {
                field_type = .integer;
            } else if (std.mem.eql(u8, type_str, "Float")) {
                field_type = .float;
            } else if (std.mem.eql(u8, type_str, "Flag")) {
                field_type = .flag;
            } else if (std.mem.eql(u8, type_str, "Character")) {
                field_type = .character;
            }
        }

        const duped_id = try self.allocator.dupe(u8, id_str);
        errdefer self.allocator.free(duped_id);
        // Track duped_id for cleanup
        const gop = try map.getOrPut(self.allocator, duped_id);
        if (gop.found_existing) {
            self.allocator.free(duped_id);
        } else {
            try self._storage.append(self.allocator, duped_id);
        }
        gop.value_ptr.* = .{
            .number = number,
            .field_type = field_type,
        };
    }

    fn parseChromLine(self: *Header, line: []const u8) !void {
        // #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO[\tFORMAT\tSample1\t...]
        var col: usize = 0;
        var start: usize = 0;
        for (line, 0..) |c, i| {
            if (c == '\t') {
                col += 1;
                start = i + 1;
                if (col >= 9) break;
            }
        }
        if (col < 9) return; // no samples

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
        if (s < line.len) {
            const duped = try self.allocator.dupe(u8, line[s..]);
            try self.sample_names.append(self.allocator, duped);
        }
    }
};

/// Metadata about a header field (INFO or FORMAT).
pub const HeaderField = struct {
    number: HeaderNumber,
    field_type: HeaderFieldType,
};

pub const HeaderNumber = union(enum) {
    fixed: u32,
    per_alt, // "A"
    per_allele, // "R"
    per_genotype, // "G"
    variable, // "."
};

pub const HeaderFieldType = enum {
    integer,
    float,
    flag,
    character,
    string,
};

/// Iterator over comma-separated ALT alleles.
pub const AltIterator = struct {
    data: []const u8,
    pos: usize,

    pub fn next(self: *AltIterator) ?[]const u8 {
        if (self.pos >= self.data.len) return null;
        const start = self.pos;
        while (self.pos < self.data.len and self.data[self.pos] != ',') {
            self.pos += 1;
        }
        const result = self.data[start..self.pos];
        if (self.pos < self.data.len) self.pos += 1; // skip comma
        return result;
    }
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "NativeReader: read from memory" {
    const allocator = std.testing.allocator;

    const vcf_text =
        "##fileformat=VCFv4.2\n" ++
        "##contig=<ID=chr1,length=249250621>\n" ++
        "##contig=<ID=chr2,length=243199373>\n" ++
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" ++
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSampleA\tSampleB\n" ++
        "chr1\t1001\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/1\t1/2\n" ++
        "chr2\t5000\trs456\tTG\tT\t30\t.\tDP=42\tGT\t0/0\t0/1\n";

    var reader = try NativeReader.fromMemory(allocator, vcf_text);
    defer reader.deinit();

    // Header checks
    try std.testing.expectEqual(@as(usize, 2), reader.header.nSamples());
    try std.testing.expectEqualStrings("SampleA", reader.header.sample_names.items[0]);
    try std.testing.expectEqualStrings("SampleB", reader.header.sample_names.items[1]);
    try std.testing.expectEqual(@as(usize, 2), reader.header.contig_names.items.len);
    try std.testing.expectEqualStrings("chr1", reader.header.contig_names.items[0]);
    try std.testing.expectEqualStrings("chr2", reader.header.contig_names.items[1]);

    // Contig map
    try std.testing.expectEqual(@as(i32, 0), reader.header.contig_map.get("chr1").?);
    try std.testing.expectEqual(@as(i32, 1), reader.header.contig_map.get("chr2").?);

    // INFO field metadata
    const dp_info = reader.header.info_fields.get("DP").?;
    try std.testing.expectEqual(HeaderFieldType.integer, dp_info.field_type);
    try std.testing.expectEqual(@as(u32, 1), dp_info.number.fixed);

    // FORMAT field metadata
    _ = reader.header.format_fields.get("GT").?;

    // Record 1
    const ok1 = try reader.next();
    try std.testing.expect(ok1);
    try std.testing.expectEqualStrings("chr1", reader.chrom());
    try std.testing.expectEqual(@as(u32, 1000), reader.pos());
    try std.testing.expectEqualStrings(".", reader.id());
    try std.testing.expectEqualStrings("A", reader.refAllele());
    try std.testing.expectEqual(@as(?f32, null), reader.qual());
    try std.testing.expectEqualStrings("PASS", reader.filter());
    try std.testing.expectEqualStrings("GT", reader.format());

    // ALT alleles
    var alts1 = reader.altAlleles();
    try std.testing.expectEqualStrings("C", alts1.next().?);
    try std.testing.expectEqualStrings("G", alts1.next().?);
    try std.testing.expectEqual(@as(?[]const u8, null), alts1.next());

    // Record 2
    const ok2 = try reader.next();
    try std.testing.expect(ok2);
    try std.testing.expectEqualStrings("chr2", reader.chrom());
    try std.testing.expectEqual(@as(u32, 4999), reader.pos());
    try std.testing.expectEqualStrings("rs456", reader.id());
    try std.testing.expectEqualStrings("TG", reader.refAllele());
    try std.testing.expectEqual(@as(?f32, 30.0), reader.qual());
    try std.testing.expectEqualStrings("DP=42", reader.info());

    var alts2 = reader.altAlleles();
    try std.testing.expectEqualStrings("T", alts2.next().?);
    try std.testing.expectEqual(@as(?[]const u8, null), alts2.next());

    // EOF
    const ok3 = try reader.next();
    try std.testing.expect(!ok3);
}

test "NativeReader: no samples" {
    const allocator = std.testing.allocator;

    const vcf_text =
        "##fileformat=VCFv4.2\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" ++
        "chr1\t100\t.\tA\tT\t.\t.\t.\n";

    var reader = try NativeReader.fromMemory(allocator, vcf_text);
    defer reader.deinit();

    try std.testing.expectEqual(@as(usize, 0), reader.header.nSamples());

    const ok = try reader.next();
    try std.testing.expect(ok);
    try std.testing.expectEqualStrings("chr1", reader.chrom());
    try std.testing.expectEqual(@as(u32, 99), reader.pos());
    try std.testing.expectEqualStrings("A", reader.refAllele());

    var alts = reader.altAlleles();
    try std.testing.expectEqualStrings("T", alts.next().?);
    try std.testing.expectEqual(@as(?[]const u8, null), alts.next());
}

test "NativeReader: lazy field parsing" {
    const allocator = std.testing.allocator;

    const vcf_text =
        "##fileformat=VCFv4.2\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" ++
        "chr1\t200\t.\tG\tA\t50\tPASS\tDP=10\n";

    var reader = try NativeReader.fromMemory(allocator, vcf_text);
    defer reader.deinit();

    const ok = try reader.next();
    try std.testing.expect(ok);

    // Before any access, fields should not be parsed
    try std.testing.expect(!reader.fields_parsed);

    // Accessing a field triggers parsing
    const c = reader.chrom();
    try std.testing.expectEqualStrings("chr1", c);
    try std.testing.expect(reader.fields_parsed);
}

test "NativeReader: no alt alleles (dot)" {
    const allocator = std.testing.allocator;

    const vcf_text =
        "##fileformat=VCFv4.2\n" ++
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" ++
        "chr1\t100\t.\tA\t.\t.\t.\t.\n";

    var reader = try NativeReader.fromMemory(allocator, vcf_text);
    defer reader.deinit();

    const ok = try reader.next();
    try std.testing.expect(ok);

    var alts = reader.altAlleles();
    try std.testing.expectEqual(@as(?[]const u8, null), alts.next());
}

test "AltIterator" {
    // Multiple alleles
    var it = AltIterator{ .data = "C,G,T", .pos = 0 };
    try std.testing.expectEqualStrings("C", it.next().?);
    try std.testing.expectEqualStrings("G", it.next().?);
    try std.testing.expectEqualStrings("T", it.next().?);
    try std.testing.expectEqual(@as(?[]const u8, null), it.next());

    // Single allele
    var it2 = AltIterator{ .data = "C", .pos = 0 };
    try std.testing.expectEqualStrings("C", it2.next().?);
    try std.testing.expectEqual(@as(?[]const u8, null), it2.next());

    // Empty
    var it3 = AltIterator{ .data = &.{}, .pos = 0 };
    try std.testing.expectEqual(@as(?[]const u8, null), it3.next());
}
