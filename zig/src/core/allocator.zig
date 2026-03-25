const std = @import("std");

/// Arena allocator for per-transcript work.
/// All allocations within a transcript's lifetime use this arena,
/// and are freed in bulk when the transcript is flushed.
pub const TranscriptArena = struct {
    arena: std.heap.ArenaAllocator,

    pub fn init(backing: std.mem.Allocator) TranscriptArena {
        return .{ .arena = std.heap.ArenaAllocator.init(backing) };
    }

    pub fn allocator(self: *TranscriptArena) std.mem.Allocator {
        return self.arena.allocator();
    }

    /// Free all transcript-lifetime allocations at once.
    pub fn reset(self: *TranscriptArena) void {
        self.arena.reset(.free_all);
    }

    pub fn deinit(self: *TranscriptArena) void {
        self.arena.deinit();
    }
};

/// Scratch arena for per-record temporaries.
pub const ScratchArena = struct {
    arena: std.heap.ArenaAllocator,

    pub fn init(backing: std.mem.Allocator) ScratchArena {
        return .{ .arena = std.heap.ArenaAllocator.init(backing) };
    }

    pub fn allocator(self: *ScratchArena) std.mem.Allocator {
        return self.arena.allocator();
    }

    pub fn reset(self: *ScratchArena) void {
        self.arena.reset(.retain_capacity);
    }

    pub fn deinit(self: *ScratchArena) void {
        self.arena.deinit();
    }
};
