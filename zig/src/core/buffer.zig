const std = @import("std");

pub fn RingBuffer(comptime T: type) type {
    return struct {
        const Self = @This();

        items: []T,
        n: usize,      // number of elements
        f: usize,       // front index
        allocator: std.mem.Allocator,

        pub fn init(alloc: std.mem.Allocator) Self {
            return .{
                .items = &.{},
                .n = 0,
                .f = 0,
                .allocator = alloc,
            };
        }

        pub fn deinit(self: *Self) void {
            if (self.items.len > 0) {
                self.allocator.free(self.items);
            }
        }

        pub fn ensureCapacity(self: *Self, min_cap: usize) !void {
            if (self.items.len >= min_cap) return;
            const new_cap = @max(min_cap, if (self.items.len == 0) 8 else self.items.len * 2);
            const new_items = try self.allocator.alloc(T, new_cap);
            // Copy existing elements
            if (self.n > 0) {
                const end = self.f + self.n;
                if (end <= self.items.len) {
                    @memcpy(new_items[0..self.n], self.items[self.f..end]);
                } else {
                    const first_part = self.items.len - self.f;
                    @memcpy(new_items[0..first_part], self.items[self.f..]);
                    @memcpy(new_items[first_part..self.n], self.items[0..end - self.items.len]);
                }
            }
            if (self.items.len > 0) {
                self.allocator.free(self.items);
            }
            self.items = new_items;
            self.f = 0;
        }

        pub fn append(self: *Self) !usize {
            try self.ensureCapacity(self.n + 1);
            const idx = (self.f + self.n) % self.items.len;
            self.n += 1;
            return idx;
        }

        pub fn shift(self: *Self) ?usize {
            if (self.n == 0) return null;
            const idx = self.f;
            self.f = (self.f + 1) % self.items.len;
            self.n -= 1;
            return idx;
        }

        pub fn last(self: *const Self) ?usize {
            if (self.n == 0) return null;
            return (self.f + self.n - 1) % self.items.len;
        }

        pub fn kth(self: *const Self, k: usize) usize {
            return (self.f + k) % self.items.len;
        }
    };
}

test "ring buffer basic operations" {
    const alloc = std.testing.allocator;
    var rb = RingBuffer(i32).init(alloc);
    defer rb.deinit();

    const idx0 = try rb.append();
    rb.items[idx0] = 10;
    const idx1 = try rb.append();
    rb.items[idx1] = 20;

    try std.testing.expectEqual(@as(usize, 2), rb.n);

    const shifted = rb.shift().?;
    try std.testing.expectEqual(@as(i32, 10), rb.items[shifted]);
    try std.testing.expectEqual(@as(usize, 1), rb.n);
}
