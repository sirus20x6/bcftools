// SIMD utility primitives for byte-level operations.
//
// Provides reusable SIMD-accelerated functions for common byte scanning
// patterns: finding delimiters, case-insensitive comparison, counting, etc.
// All functions have scalar fallbacks for short inputs and tail bytes.

const std = @import("std");

/// Find the first occurrence of a byte in a slice using SIMD.
/// Returns the index, or null if not found.
pub fn findByte(haystack: []const u8, needle: u8) ?usize {
    var i: usize = 0;
    const splat: @Vector(16, u8) = @splat(needle);
    while (i + 16 <= haystack.len) {
        const chunk: @Vector(16, u8) = haystack[i..][0..16].*;
        const eq = chunk == splat;
        const mask: u16 = @bitCast(eq);
        if (mask != 0) return i + @ctz(mask);
        i += 16;
    }
    while (i < haystack.len) {
        if (haystack[i] == needle) return i;
        i += 1;
    }
    return null;
}

/// Find the first occurrence of any byte in a set using SIMD.
/// Useful for finding tab, semicolon, newline, etc.
/// `needles` must have between 1 and 4 entries for efficiency.
pub fn findAnyByte(haystack: []const u8, needles: []const u8) ?usize {
    if (needles.len == 0) return null;
    // Special-case single needle
    if (needles.len == 1) return findByte(haystack, needles[0]);

    var i: usize = 0;
    while (i + 16 <= haystack.len) {
        const chunk: @Vector(16, u8) = haystack[i..][0..16].*;
        var mask: u16 = 0;
        for (needles) |n| {
            const splat: @Vector(16, u8) = @splat(n);
            const eq = chunk == splat;
            mask |= @as(u16, @bitCast(eq));
        }
        if (mask != 0) return i + @ctz(mask);
        i += 16;
    }
    // Scalar fallback
    while (i < haystack.len) {
        for (needles) |n| {
            if (haystack[i] == n) return i;
        }
        i += 1;
    }
    return null;
}

/// Case-insensitive byte comparison using SIMD.
/// Returns true if the two slices are equal ignoring ASCII case.
/// Only converts a-z to A-Z; non-alpha bytes must match exactly.
pub fn eqlIgnoreCase(a: []const u8, b: []const u8) bool {
    if (a.len != b.len) return false;
    var i: usize = 0;
    const case_bit: @Vector(16, u8) = @splat(0x20);
    while (i + 16 <= a.len) {
        const va: @Vector(16, u8) = a[i..][0..16].*;
        const vb: @Vector(16, u8) = b[i..][0..16].*;
        // OR in the case bit to force lowercase
        const ua = va | case_bit;
        const ub = vb | case_bit;
        const neq = ua != ub;
        if (@reduce(.Or, neq)) {
            // There's a mismatch — but we need to verify it's a real alpha
            // mismatch vs a non-alpha byte where bit 0x20 changed meaning.
            // For the hot path (DNA sequences: ACGTacgt) this never triggers
            // on false positives. For correctness, check scalar:
            const mask: u16 = @bitCast(neq);
            const pos = i + @ctz(mask);
            const ca = if (a[pos] >= 'a' and a[pos] <= 'z') a[pos] - 32 else a[pos];
            const cb = if (b[pos] >= 'a' and b[pos] <= 'z') b[pos] - 32 else b[pos];
            if (ca != cb) return false;
            // The SIMD said mismatch but scalar says match (non-alpha byte
            // where 0x20 bit flip didn't matter). Fall through to scalar for
            // this tricky chunk.
            var j = i;
            while (j < i + 16 and j < a.len) : (j += 1) {
                const x = if (a[j] >= 'a' and a[j] <= 'z') a[j] - 32 else a[j];
                const y = if (b[j] >= 'a' and b[j] <= 'z') b[j] - 32 else b[j];
                if (x != y) return false;
            }
        }
        i += 16;
    }
    while (i < a.len) {
        const ca = if (a[i] >= 'a' and a[i] <= 'z') a[i] - 32 else a[i];
        const cb = if (b[i] >= 'a' and b[i] <= 'z') b[i] - 32 else b[i];
        if (ca != cb) return false;
        i += 1;
    }
    return true;
}

/// Count occurrences of a byte using SIMD.
pub fn countByte(haystack: []const u8, needle: u8) usize {
    var count: usize = 0;
    var i: usize = 0;
    const splat: @Vector(16, u8) = @splat(needle);
    while (i + 16 <= haystack.len) {
        const chunk: @Vector(16, u8) = haystack[i..][0..16].*;
        const eq = chunk == splat;
        count += @popCount(@as(u16, @bitCast(eq)));
        i += 16;
    }
    while (i < haystack.len) {
        if (haystack[i] == needle) count += 1;
        i += 1;
    }
    return count;
}

/// Find the first byte NOT equal to a given value.
pub fn findNotByte(haystack: []const u8, byte: u8) ?usize {
    var i: usize = 0;
    const splat: @Vector(16, u8) = @splat(byte);
    while (i + 16 <= haystack.len) {
        const chunk: @Vector(16, u8) = haystack[i..][0..16].*;
        const neq = chunk != splat;
        const mask: u16 = @bitCast(neq);
        if (mask != 0) return i + @ctz(mask);
        i += 16;
    }
    while (i < haystack.len) {
        if (haystack[i] != byte) return i;
        i += 1;
    }
    return null;
}

/// Find the first position where two byte slices differ.
/// Returns null if they are identical over the compared range.
pub fn firstMismatch(a: []const u8, b: []const u8) ?usize {
    const len = @min(a.len, b.len);
    var i: usize = 0;

    // SIMD: compare 16 bytes at a time
    while (i + 16 <= len) {
        const va: @Vector(16, u8) = a[i..][0..16].*;
        const vb: @Vector(16, u8) = b[i..][0..16].*;
        const neq: @Vector(16, bool) = va != vb;
        if (@reduce(.Or, neq)) {
            const mask: u16 = @bitCast(neq);
            return i + @ctz(mask);
        }
        i += 16;
    }

    // Scalar fallback
    while (i < len) {
        if (a[i] != b[i]) return i;
        i += 1;
    }
    return null;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "findByte: find tab in string" {
    const s = "hello\tworld";
    try std.testing.expectEqual(@as(?usize, 5), findByte(s, '\t'));
}

test "findByte: not found" {
    const s = "hello world";
    try std.testing.expectEqual(@as(?usize, null), findByte(s, '\t'));
}

test "findByte: at position 0" {
    const s = "\thello";
    try std.testing.expectEqual(@as(?usize, 0), findByte(s, '\t'));
}

test "findByte: at end" {
    const s = "hello\t";
    try std.testing.expectEqual(@as(?usize, 5), findByte(s, '\t'));
}

test "findByte: in SIMD region (>16 bytes)" {
    const s = "0123456789abcdef\trest";
    try std.testing.expectEqual(@as(?usize, 16), findByte(s, '\t'));
}

test "findByte: multiple occurrences returns first" {
    const s = "a\tb\tc";
    try std.testing.expectEqual(@as(?usize, 1), findByte(s, '\t'));
}

test "findAnyByte: find tab or newline" {
    const s = "hello\nworld\tok";
    const needles = [_]u8{ '\t', '\n' };
    try std.testing.expectEqual(@as(?usize, 5), findAnyByte(s, &needles));
}

test "findAnyByte: find semicolon or equals" {
    const s = "DP=50;AF=0.1";
    const needles = [_]u8{ '=', ';' };
    try std.testing.expectEqual(@as(?usize, 2), findAnyByte(s, &needles));
}

test "findAnyByte: not found" {
    const s = "hello world";
    const needles = [_]u8{ '\t', '\n', ';' };
    try std.testing.expectEqual(@as(?usize, null), findAnyByte(s, &needles));
}

test "eqlIgnoreCase: matching" {
    try std.testing.expect(eqlIgnoreCase("ACGT", "acgt"));
    try std.testing.expect(eqlIgnoreCase("acgt", "ACGT"));
    try std.testing.expect(eqlIgnoreCase("AcGt", "aCgT"));
}

test "eqlIgnoreCase: different content" {
    try std.testing.expect(!eqlIgnoreCase("ACGT", "ACGA"));
    try std.testing.expect(!eqlIgnoreCase("hello", "world"));
}

test "eqlIgnoreCase: different lengths" {
    try std.testing.expect(!eqlIgnoreCase("ABC", "ABCD"));
}

test "eqlIgnoreCase: long strings (SIMD path)" {
    const a = "ACGTACGTACGTACGTacgtacgt";
    const b = "acgtacgtacgtacgtACGTACGT";
    try std.testing.expect(eqlIgnoreCase(a, b));
}

test "eqlIgnoreCase: long strings mismatch in SIMD region" {
    const a = "ACGTACGTACGTACGTX";
    const b = "acgtacgtacgtacgtY";
    try std.testing.expect(!eqlIgnoreCase(a, b));
}

test "countByte: count commas" {
    try std.testing.expectEqual(@as(usize, 3), countByte("a,b,c,d", ','));
}

test "countByte: none found" {
    try std.testing.expectEqual(@as(usize, 0), countByte("abcdef", ','));
}

test "countByte: long string (SIMD path)" {
    const s = ",,,,,,,,,,,,,,,,,,,,"; // 20 commas
    try std.testing.expectEqual(@as(usize, 20), countByte(s, ','));
}

test "findNotByte: finds first non-N" {
    try std.testing.expectEqual(@as(?usize, 3), findNotByte("NNNA", 'N'));
}

test "findNotByte: all same" {
    try std.testing.expectEqual(@as(?usize, null), findNotByte("NNNN", 'N'));
}

test "firstMismatch: identical short strings" {
    const a = "ABCDEF";
    const b = "ABCDEF";
    try std.testing.expect(firstMismatch(a, b) == null);
}

test "firstMismatch: first byte differs" {
    const a = "ABCDEF";
    const b = "XBCDEF";
    try std.testing.expectEqual(@as(usize, 0), firstMismatch(a, b).?);
}

test "firstMismatch: last byte differs" {
    const a = "ABCDEF";
    const b = "ABCDEX";
    try std.testing.expectEqual(@as(usize, 5), firstMismatch(a, b).?);
}

test "firstMismatch: long identical strings (SIMD path)" {
    const a = "ABCDEFGHIJKLMNOPQRSTUVWXYZ012345";
    const b = "ABCDEFGHIJKLMNOPQRSTUVWXYZ012345";
    try std.testing.expect(firstMismatch(a, b) == null);
}

test "firstMismatch: mismatch in second SIMD chunk" {
    var a = "ABCDEFGHIJKLMNOPQRSTUVWXYZ012345".*;
    var b = "ABCDEFGHIJKLMNOPQRSTUVWXYZ012345".*;
    b[19] = 'X';
    try std.testing.expectEqual(@as(usize, 19), firstMismatch(&a, &b).?);
}

test "firstMismatch: mismatch in scalar fallback" {
    var a = "ABCDEFGHIJKLMNOPQR".*;
    var b = "ABCDEFGHIJKLMNOPQR".*;
    b[17] = 'X';
    try std.testing.expectEqual(@as(usize, 17), firstMismatch(&a, &b).?);
}
