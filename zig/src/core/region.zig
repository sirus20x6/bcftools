const std = @import("std");

pub fn RegionIndex(comptime Payload: type) type {
    return struct {
        const Self = @This();

        pub const Interval = struct {
            beg: u32,
            end: u32,
            insert_order: u32 = 0, // preserves GFF insertion order for stable sorting
            payload: Payload,
        };

        pub const OverlapIterator = struct {
            intervals: []const Interval,
            idx: usize,
            query_beg: u32,
            query_end: u32,

            pub fn next(self: *OverlapIterator) ?*const Interval {
                while (self.idx < self.intervals.len) {
                    const iv = &self.intervals[self.idx];
                    self.idx += 1;
                    if (iv.beg > self.query_end) return null; // past query, done
                    if (iv.end >= self.query_beg) return iv; // overlaps
                }
                return null;
            }

            pub fn reset(self: *OverlapIterator, beg: u32, end: u32) void {
                self.query_beg = beg;
                self.query_end = end;
                // binary search for first interval that could overlap
                self.idx = lowerBound(self.intervals, beg);
            }
        };

        sequences: std.StringHashMap(std.ArrayList(Interval)),
        sorted: bool,
        allocator: std.mem.Allocator,
        next_insert_order: u32 = 0,

        pub fn init(alloc: std.mem.Allocator) Self {
            return .{
                .sequences = std.StringHashMap(std.ArrayList(Interval)).init(alloc),
                .sorted = false,
                .allocator = alloc,
            };
        }

        pub fn deinit(self: *Self) void {
            var it = self.sequences.iterator();
            while (it.next()) |entry| {
                entry.value_ptr.deinit(self.allocator);
            }
            self.sequences.deinit();
        }

        pub fn insert(self: *Self, seq: []const u8, beg: u32, end: u32, payload: Payload) !void {
            const result = try self.sequences.getOrPut(seq);
            if (!result.found_existing) {
                result.value_ptr.* = .empty;
            }
            try result.value_ptr.append(self.allocator, .{ .beg = beg, .end = end, .insert_order = self.next_insert_order, .payload = payload });
            self.next_insert_order += 1;
            self.sorted = false;
        }

        pub fn sort(self: *Self) void {
            var it = self.sequences.iterator();
            while (it.next()) |entry| {
                std.mem.sort(Interval, entry.value_ptr.items, {}, struct {
                    fn lessThan(_: void, a: Interval, b: Interval) bool {
                        if (a.beg != b.beg) return a.beg < b.beg;
                        if (a.end != b.end) return a.end < b.end;
                        return a.insert_order < b.insert_order;
                    }
                }.lessThan);
            }
            self.sorted = true;
        }

        pub fn overlap(self: *Self, seq: []const u8, beg: u32, end: u32) OverlapIterator {
            if (!self.sorted) self.sort();
            const intervals = if (self.sequences.get(seq)) |list| list.items else &[_]Interval{};
            // Start from 0 because intervals sorted by beg don't have monotonic
            // end values — a binary search on end can miss earlier intervals whose
            // beg is small but end extends past the query start.  The iterator's
            // beg > query_end check still provides an early exit.
            return .{
                .intervals = intervals,
                .idx = 0,
                .query_beg = beg,
                .query_end = end,
            };
        }

        pub fn hasSeq(self: *const Self, seq: []const u8) bool {
            return self.sequences.contains(seq);
        }

        fn lowerBound(intervals: []const Interval, target_beg: u32) usize {
            // Find first interval whose end >= target_beg
            // (all intervals before this point end before our query begins)
            var lo: usize = 0;
            var hi: usize = intervals.len;
            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                if (intervals[mid].end < target_beg) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            return lo;
        }
    };
}

test "region index overlap" {
    const alloc = std.testing.allocator;
    var idx = RegionIndex(u32).init(alloc);
    defer idx.deinit();

    try idx.insert("chr1", 100, 200, 1);
    try idx.insert("chr1", 150, 300, 2);
    try idx.insert("chr1", 400, 500, 3);

    var itr = idx.overlap("chr1", 120, 160);
    const first = itr.next();
    try std.testing.expect(first != null);
    try std.testing.expectEqual(@as(u32, 1), first.?.payload);

    const second = itr.next();
    try std.testing.expect(second != null);
    try std.testing.expectEqual(@as(u32, 2), second.?.payload);

    const third = itr.next();
    try std.testing.expect(third == null);
}
