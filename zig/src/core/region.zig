const std = @import("std");

pub fn RegionIndex(comptime Payload: type) type {
    return struct {
        const Self = @This();

        pub const Interval = struct {
            beg: u32,
            end: u32,
            insert_order: u32 = 0, // preserves GFF insertion order for stable sorting
            max_end: u32 = 0, // max(end) for all intervals at this index and before (set after sort)
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
                    // intervals sorted by beg: once beg > query_end, no more can overlap
                    if (iv.beg > self.query_end) return null;
                    if (iv.end >= self.query_beg) return iv; // overlaps
                }
                return null;
            }

            pub fn reset(self: *OverlapIterator, beg: u32, end: u32) void {
                self.query_beg = beg;
                self.query_end = end;
                self.idx = 0;
            }
        };

        sequences: std.StringHashMap(std.ArrayList(Interval)),
        sorted: bool,
        allocator: std.mem.Allocator,
        next_insert_order: u32 = 0,

        // Gap cache: when a query returns no overlaps, we cache the gap
        // boundaries so that subsequent queries within the same gap can
        // skip the binary searches entirely.  Genomic data is sorted, so
        // consecutive intergenic variants benefit heavily from this.
        gap_cache_seq: ?[*]const u8 = null, // pointer identity of cached sequence key
        gap_cache_lo: u32 = 0, // lower bound of gap (exclusive: max_end of last interval before gap)
        gap_cache_hi: u32 = 0, // upper bound of gap (exclusive: beg of first interval after gap)

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
                // Build prefix max-end cache for fast overlap scan-start
                var running_max: u32 = 0;
                for (entry.value_ptr.items) |*iv| {
                    running_max = @max(running_max, iv.end);
                    iv.max_end = running_max;
                }
            }
            self.sorted = true;
        }

        /// Empty iterator singleton (returned for gap cache hits).
        const empty_intervals: []const Interval = &[_]Interval{};

        pub fn overlap(self: *Self, seq: []const u8, beg: u32, end: u32) OverlapIterator {
            if (!self.sorted) self.sort();

            // Look up the interval list for this sequence.
            // We use getKeyPtr to get the stable HashMap key pointer for gap caching.
            const entry = self.sequences.getEntry(seq);
            const intervals = if (entry) |e| e.value_ptr.items else &[_]Interval{};
            const stable_key_ptr: ?[*]const u8 = if (entry) |e| e.key_ptr.*.ptr else null;

            // Gap cache fast path: if query falls entirely within a cached
            // gap on the same sequence, return empty immediately (no binary searches).
            if (self.gap_cache_seq) |cached_ptr| {
                if (stable_key_ptr == cached_ptr and beg > self.gap_cache_lo and end < self.gap_cache_hi) {
                    return .{
                        .intervals = empty_intervals,
                        .idx = 0,
                        .query_beg = beg,
                        .query_end = end,
                    };
                }
            }

            // Binary search: find the first interval where beg > query_end.
            const hi = upperBound(intervals, end);
            // Binary search: find first index where max_end >= query_beg.
            const lo = lowerBoundMaxEnd(intervals[0..hi], beg);

            // Update gap cache when the result is empty (no overlaps found).
            if (lo >= hi) {
                if (stable_key_ptr) |skp| {
                    self.gap_cache_seq = skp;
                    self.gap_cache_lo = if (lo > 0) intervals[lo - 1].max_end else 0;
                    self.gap_cache_hi = if (hi < intervals.len) intervals[hi].beg else std.math.maxInt(u32);
                }
            } else {
                // Overlaps found — invalidate cache for this sequence
                if (self.gap_cache_seq) |cached_ptr| {
                    if (stable_key_ptr == cached_ptr) {
                        self.gap_cache_seq = null;
                    }
                }
            }

            return .{
                .intervals = intervals[lo..hi],
                .idx = 0,
                .query_beg = beg,
                .query_end = end,
            };
        }

        pub fn hasSeq(self: *const Self, seq: []const u8) bool {
            return self.sequences.contains(seq);
        }

        /// Find the first index where max_end >= target_beg.
        /// All intervals before this index have max_end < target_beg,
        /// meaning none of them (or any earlier interval) can overlap.
        fn lowerBoundMaxEnd(intervals: []const Interval, target_beg: u32) usize {
            var lo: usize = 0;
            var hi: usize = intervals.len;
            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                if (intervals[mid].max_end < target_beg) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            return lo;
        }

        /// Find the first index where interval.beg > target.
        /// All intervals before this index have beg <= target, so they
        /// *could* overlap a query ending at `target`.
        fn upperBound(intervals: []const Interval, target: u32) usize {
            var lo: usize = 0;
            var hi: usize = intervals.len;
            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                if (intervals[mid].beg > target) {
                    hi = mid;
                } else {
                    lo = mid + 1;
                }
            }
            return lo; // first index where beg > target
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
