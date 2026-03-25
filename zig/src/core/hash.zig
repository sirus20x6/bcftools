const std = @import("std");

/// Position to Vbuf mapping (replaces kh_pos2vbuf)
pub fn PosMap(comptime V: type) type {
    return std.AutoHashMap(u32, V);
}

/// String to integer mapping (replaces khash_str2int)
pub fn StringIntMap() type {
    return std.StringHashMap(u32);
}

/// String to string mapping (replaces khash_str2str)
pub fn StringMap() type {
    return std.StringHashMap([]const u8);
}
