// Filter expression engine for bcftools -i/-e expressions.
//
// Supports a subset of the full C filter engine covering ~80% of real-world
// usage: QUAL, INFO field comparisons, TYPE, N_ALT, boolean operators,
// comparison operators, and parenthesized grouping.
//
// Architecture: tokenize -> shunting-yard RPN conversion -> per-record evaluation.

const std = @import("std");
const VcfRecord = @import("../vcf/record.zig").VcfRecord;

/// Token types used in the filter expression.
pub const TokenType = enum {
    // Literals
    int_literal,
    float_literal,
    string_literal,

    // Field references
    qual,
    info_int,
    info_float,
    info_string,
    info_flag,
    n_alt,
    type_field,

    // Comparison operators
    op_gt,
    op_gte,
    op_lt,
    op_lte,
    op_eq,
    op_neq,

    // Boolean operators
    op_and,
    op_or,
    op_not,

    // Grouping
    lparen,
    rparen,
};

/// A single token in the filter expression.
pub const Token = struct {
    tok_type: TokenType,
    int_val: i64 = 0,
    float_val: f64 = 0,
    str_val: ?[]const u8 = null,
    /// For INFO field references: the tag name (e.g. "DP" for INFO/DP).
    field_name: ?[]const u8 = null,
};

/// Operator precedence for the shunting-yard algorithm.
fn precedence(tok: TokenType) u8 {
    return switch (tok) {
        .op_not => 6,
        .op_gt, .op_gte, .op_lt, .op_lte, .op_eq, .op_neq => 4,
        .op_and => 2,
        .op_or => 1,
        else => 0,
    };
}

/// Whether the operator is right-associative.
fn isRightAssoc(tok: TokenType) bool {
    return tok == .op_not;
}

/// Filter expression engine.
///
/// Parse an expression string, convert to RPN via shunting-yard, and
/// evaluate against VcfRecord instances.
pub const Filter = struct {
    allocator: std.mem.Allocator,
    tokens: []Token,
    rpn: []Token,
    /// True when this is an exclude filter (negate the result).
    negate: bool,

    /// Compile a filter expression.
    /// `expr` is the expression string (e.g. "QUAL>20 && INFO/DP>30").
    /// If `exclude` is true, the filter result is negated (for -e semantics).
    pub fn init(allocator: std.mem.Allocator, expr: []const u8, exclude: bool) !Filter {
        var tokens: std.ArrayListUnmanaged(Token) = .empty;
        defer tokens.deinit(allocator);

        try tokenize(allocator, expr, &tokens);

        const owned_tokens = try allocator.dupe(Token, tokens.items);
        errdefer allocator.free(owned_tokens);

        var rpn_list: std.ArrayListUnmanaged(Token) = .empty;
        defer rpn_list.deinit(allocator);

        try shuntingYard(allocator, owned_tokens, &rpn_list);

        const rpn = try allocator.dupe(Token, rpn_list.items);

        return .{
            .allocator = allocator,
            .tokens = owned_tokens,
            .rpn = rpn,
            .negate = exclude,
        };
    }

    pub fn deinit(self: *Filter) void {
        self.allocator.free(self.tokens);
        self.allocator.free(self.rpn);
        self.* = undefined;
    }

    /// Evaluate the filter against a VCF record.
    /// Returns true if the record passes the filter.
    pub fn eval(self: *const Filter, rec: *const VcfRecord) bool {
        var stack: [64]Value = undefined;
        var sp: usize = 0;

        for (self.rpn) |tok| {
            switch (tok.tok_type) {
                // --- Push literals ---
                .int_literal => {
                    stack[sp] = .{ .number = @floatFromInt(tok.int_val) };
                    sp += 1;
                },
                .float_literal => {
                    stack[sp] = .{ .number = tok.float_val };
                    sp += 1;
                },
                .string_literal => {
                    stack[sp] = .{ .string = tok.str_val orelse "" };
                    sp += 1;
                },

                // --- Push field values ---
                .qual => {
                    if (rec.qual) |q| {
                        stack[sp] = .{ .number = @floatCast(q) };
                    } else {
                        stack[sp] = .missing;
                    }
                    sp += 1;
                },
                .info_int => {
                    if (tok.field_name) |name| {
                        if (getInfoValue(rec, name)) |val| {
                            stack[sp] = .{ .number = val };
                        } else {
                            stack[sp] = .missing;
                        }
                    } else {
                        stack[sp] = .missing;
                    }
                    sp += 1;
                },
                .info_float => {
                    if (tok.field_name) |name| {
                        if (getInfoValue(rec, name)) |val| {
                            stack[sp] = .{ .number = val };
                        } else {
                            stack[sp] = .missing;
                        }
                    } else {
                        stack[sp] = .missing;
                    }
                    sp += 1;
                },
                .info_string => {
                    if (tok.field_name) |name| {
                        if (getInfoString(rec, name)) |val| {
                            stack[sp] = .{ .string = val };
                        } else {
                            stack[sp] = .missing;
                        }
                    } else {
                        stack[sp] = .missing;
                    }
                    sp += 1;
                },
                .info_flag => {
                    if (tok.field_name) |name| {
                        stack[sp] = .{ .boolean = hasInfoFlag(rec, name) };
                    } else {
                        stack[sp] = .{ .boolean = false };
                    }
                    sp += 1;
                },
                .n_alt => {
                    const n: u32 = if (rec.alt_alleles.items.len > 0)
                        @intCast(rec.alt_alleles.items.len)
                    else
                        0;
                    stack[sp] = .{ .number = @floatFromInt(n) };
                    sp += 1;
                },
                .type_field => {
                    stack[sp] = .{ .string = getVariantType(rec) };
                    sp += 1;
                },

                // --- Comparison operators ---
                .op_gt, .op_gte, .op_lt, .op_lte, .op_eq, .op_neq => {
                    if (sp < 2) {
                        stack[0] = .{ .boolean = false };
                        sp = 1;
                        continue;
                    }
                    sp -= 2;
                    const lhs = stack[sp];
                    const rhs = stack[sp + 1];
                    stack[sp] = .{ .boolean = evalComparison(tok.tok_type, lhs, rhs) };
                    sp += 1;
                },

                // --- Boolean operators ---
                .op_and => {
                    if (sp < 2) {
                        stack[0] = .{ .boolean = false };
                        sp = 1;
                        continue;
                    }
                    sp -= 2;
                    const lhs = toBool(stack[sp]);
                    const rhs = toBool(stack[sp + 1]);
                    stack[sp] = .{ .boolean = lhs and rhs };
                    sp += 1;
                },
                .op_or => {
                    if (sp < 2) {
                        stack[0] = .{ .boolean = false };
                        sp = 1;
                        continue;
                    }
                    sp -= 2;
                    const lhs = toBool(stack[sp]);
                    const rhs = toBool(stack[sp + 1]);
                    stack[sp] = .{ .boolean = lhs or rhs };
                    sp += 1;
                },
                .op_not => {
                    if (sp < 1) {
                        stack[0] = .{ .boolean = true };
                        sp = 1;
                        continue;
                    }
                    sp -= 1;
                    stack[sp] = .{ .boolean = !toBool(stack[sp]) };
                    sp += 1;
                },

                // Parentheses should not appear in RPN
                .lparen, .rparen => {},
            }
        }

        const result = if (sp > 0) toBool(stack[sp - 1]) else false;
        return if (self.negate) !result else result;
    }
};

// ---------------------------------------------------------------------------
// Evaluation helpers
// ---------------------------------------------------------------------------

const Value = union(enum) {
    number: f64,
    string: []const u8,
    boolean: bool,
    missing,
};

fn toBool(v: Value) bool {
    return switch (v) {
        .number => |n| n != 0,
        .string => |s| s.len > 0,
        .boolean => |b| b,
        .missing => false,
    };
}

fn evalComparison(op: TokenType, lhs: Value, rhs: Value) bool {
    // Missing values: comparisons with missing always fail
    if (lhs == .missing or rhs == .missing) return false;

    // String comparison (for TYPE="snp" etc.)
    if (lhs == .string or rhs == .string) {
        const ls = valueToString(lhs);
        const rs = valueToString(rhs);
        return switch (op) {
            .op_eq => std.mem.eql(u8, ls, rs),
            .op_neq => !std.mem.eql(u8, ls, rs),
            else => false, // <, >, <=, >= not meaningful for strings in this subset
        };
    }

    // Numeric comparison
    const ln = valueToNumber(lhs);
    const rn = valueToNumber(rhs);
    return switch (op) {
        .op_gt => ln > rn,
        .op_gte => ln >= rn,
        .op_lt => ln < rn,
        .op_lte => ln <= rn,
        .op_eq => ln == rn,
        .op_neq => ln != rn,
        else => false,
    };
}

fn valueToNumber(v: Value) f64 {
    return switch (v) {
        .number => |n| n,
        .boolean => |b| if (b) @as(f64, 1) else @as(f64, 0),
        .string => |s| std.fmt.parseFloat(f64, s) catch 0,
        .missing => 0,
    };
}

fn valueToString(v: Value) []const u8 {
    return switch (v) {
        .string => |s| s,
        else => "",
    };
}

// ---------------------------------------------------------------------------
// VCF field accessors
// ---------------------------------------------------------------------------

/// Extract a numeric INFO field value from a VcfRecord.
/// Parses the INFO column (field 7) from the record's stored line.
fn getInfoValue(rec: *const VcfRecord, field: []const u8) ?f64 {
    const info = getInfoColumn(rec) orelse return null;
    return findInfoNumeric(info, field);
}

/// Extract a string INFO field value.
fn getInfoString(rec: *const VcfRecord, field: []const u8) ?[]const u8 {
    const info = getInfoColumn(rec) orelse return null;
    return findInfoString(info, field);
}

/// Check if an INFO flag is present.
fn hasInfoFlag(rec: *const VcfRecord, field: []const u8) bool {
    const info = getInfoColumn(rec) orelse return false;
    return findInfoFlag(info, field);
}

/// Get the INFO column (field 7) from the record's stored line.
fn getInfoColumn(rec: *const VcfRecord) ?[]const u8 {
    const storage = rec._storage orelse return null;
    var col: usize = 0;
    var start: usize = 0;
    for (storage, 0..) |c, i| {
        if (c == '\t') {
            if (col == 7) {
                return storage[start..i];
            }
            col += 1;
            start = i + 1;
        }
    }
    // Last column
    if (col == 7) {
        return storage[start..];
    }
    return null;
}

/// Find a numeric value for `tag` in an INFO string like "DP=50;AF=0.1;MQ=30".
fn findInfoNumeric(info: []const u8, tag: []const u8) ?f64 {
    var pos: usize = 0;
    while (pos < info.len) {
        // Find the start of the next key
        const key_start = pos;
        var key_end = key_start;
        while (key_end < info.len and info[key_end] != '=' and info[key_end] != ';') {
            key_end += 1;
        }
        const key = info[key_start..key_end];

        if (std.mem.eql(u8, key, tag)) {
            if (key_end < info.len and info[key_end] == '=') {
                const val_start = key_end + 1;
                var val_end = val_start;
                while (val_end < info.len and info[val_end] != ';') {
                    val_end += 1;
                }
                return std.fmt.parseFloat(f64, info[val_start..val_end]) catch null;
            }
            return null; // Flag field, no numeric value
        }

        // Skip to next field
        var next = key_end;
        if (next < info.len and info[next] == '=') {
            next += 1;
            while (next < info.len and info[next] != ';') next += 1;
        }
        if (next < info.len and info[next] == ';') next += 1;
        pos = next;
    }
    return null;
}

/// Find a string value for `tag` in INFO.
fn findInfoString(info: []const u8, tag: []const u8) ?[]const u8 {
    var pos: usize = 0;
    while (pos < info.len) {
        const key_start = pos;
        var key_end = key_start;
        while (key_end < info.len and info[key_end] != '=' and info[key_end] != ';') {
            key_end += 1;
        }
        const key = info[key_start..key_end];

        if (std.mem.eql(u8, key, tag)) {
            if (key_end < info.len and info[key_end] == '=') {
                const val_start = key_end + 1;
                var val_end = val_start;
                while (val_end < info.len and info[val_end] != ';') {
                    val_end += 1;
                }
                return info[val_start..val_end];
            }
            return null;
        }

        var next = key_end;
        if (next < info.len and info[next] == '=') {
            next += 1;
            while (next < info.len and info[next] != ';') next += 1;
        }
        if (next < info.len and info[next] == ';') next += 1;
        pos = next;
    }
    return null;
}

/// Check if `tag` is present as a flag in INFO.
fn findInfoFlag(info: []const u8, tag: []const u8) bool {
    var pos: usize = 0;
    while (pos < info.len) {
        const key_start = pos;
        var key_end = key_start;
        while (key_end < info.len and info[key_end] != '=' and info[key_end] != ';') {
            key_end += 1;
        }
        const key = info[key_start..key_end];

        if (std.mem.eql(u8, key, tag)) return true;

        var next = key_end;
        if (next < info.len and info[next] == '=') {
            next += 1;
            while (next < info.len and info[next] != ';') next += 1;
        }
        if (next < info.len and info[next] == ';') next += 1;
        pos = next;
    }
    return false;
}

/// Determine variant type from REF/ALT alleles.
/// Returns "snp", "mnp", "indel", "other", or "ref" (no alts).
pub fn getVariantType(rec: *const VcfRecord) []const u8 {
    if (rec.alt_alleles.items.len == 0) return "ref";

    // Use the first alt allele for type determination (matches bcftools behavior)
    const alt = rec.alt_alleles.items[0];
    const ref_allele = rec.ref_allele;

    if (ref_allele.len == alt.len) {
        if (ref_allele.len == 1) return "snp";
        return "mnp";
    }
    return "indel";
}

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------

fn tokenize(allocator: std.mem.Allocator, expr: []const u8, tokens: *std.ArrayListUnmanaged(Token)) !void {
    var i: usize = 0;
    while (i < expr.len) {
        const c = expr[i];

        // Skip whitespace
        if (c == ' ' or c == '\t' or c == '\n' or c == '\r') {
            i += 1;
            continue;
        }

        // Parentheses
        if (c == '(') {
            try tokens.append(allocator,.{ .tok_type = .lparen });
            i += 1;
            continue;
        }
        if (c == ')') {
            try tokens.append(allocator,.{ .tok_type = .rparen });
            i += 1;
            continue;
        }

        // Boolean operators
        if (c == '&' and i + 1 < expr.len and expr[i + 1] == '&') {
            try tokens.append(allocator,.{ .tok_type = .op_and });
            i += 2;
            continue;
        }
        if (c == '|' and i + 1 < expr.len and expr[i + 1] == '|') {
            try tokens.append(allocator,.{ .tok_type = .op_or });
            i += 2;
            continue;
        }
        if (c == '!') {
            // Could be != or !
            if (i + 1 < expr.len and expr[i + 1] == '=') {
                try tokens.append(allocator,.{ .tok_type = .op_neq });
                i += 2;
                continue;
            }
            try tokens.append(allocator,.{ .tok_type = .op_not });
            i += 1;
            continue;
        }

        // Comparison operators
        if (c == '>') {
            if (i + 1 < expr.len and expr[i + 1] == '=') {
                try tokens.append(allocator,.{ .tok_type = .op_gte });
                i += 2;
            } else {
                try tokens.append(allocator,.{ .tok_type = .op_gt });
                i += 1;
            }
            continue;
        }
        if (c == '<') {
            if (i + 1 < expr.len and expr[i + 1] == '=') {
                try tokens.append(allocator,.{ .tok_type = .op_lte });
                i += 2;
            } else {
                try tokens.append(allocator,.{ .tok_type = .op_lt });
                i += 1;
            }
            continue;
        }
        if (c == '=') {
            // "=" or "=="
            if (i + 1 < expr.len and expr[i + 1] == '=') {
                i += 2;
            } else {
                i += 1;
            }
            try tokens.append(allocator,.{ .tok_type = .op_eq });
            continue;
        }

        // String literals: "..." or '...'
        if (c == '"' or c == '\'') {
            const quote = c;
            i += 1;
            const start = i;
            while (i < expr.len and expr[i] != quote) i += 1;
            const str = expr[start..i];
            if (i < expr.len) i += 1; // skip closing quote
            try tokens.append(allocator,.{ .tok_type = .string_literal, .str_val = str });
            continue;
        }

        // INFO/ prefix
        if (startsWithAt(expr, i, "INFO/")) {
            i += 5; // skip "INFO/"
            const name_start = i;
            while (i < expr.len and isIdentChar(expr[i])) i += 1;
            const name = expr[name_start..i];

            // Peek at the next non-whitespace character to decide the token type.
            // If the next meaningful thing is a comparison operator, we'll determine
            // int vs float vs string during evaluation. Default to info_float for
            // numeric comparisons (it can hold ints too), info_string for ==/!=
            // with a string literal, and info_flag if there's no comparison at all.
            var peek = i;
            while (peek < expr.len and (expr[peek] == ' ' or expr[peek] == '\t')) peek += 1;

            if (peek < expr.len and (expr[peek] == '>' or expr[peek] == '<')) {
                try tokens.append(allocator,.{ .tok_type = .info_float, .field_name = name });
            } else if (peek < expr.len and (expr[peek] == '=' or expr[peek] == '!')) {
                // Could be string or numeric equality — peek at the RHS
                var rhs_peek = peek;
                if (rhs_peek < expr.len and expr[rhs_peek] == '!') rhs_peek += 1; // skip !
                if (rhs_peek < expr.len and expr[rhs_peek] == '=') rhs_peek += 1; // skip =
                if (rhs_peek < expr.len and expr[rhs_peek] == '=') rhs_peek += 1; // skip second =
                while (rhs_peek < expr.len and (expr[rhs_peek] == ' ' or expr[rhs_peek] == '\t')) rhs_peek += 1;
                if (rhs_peek < expr.len and (expr[rhs_peek] == '"' or expr[rhs_peek] == '\'')) {
                    try tokens.append(allocator,.{ .tok_type = .info_string, .field_name = name });
                } else {
                    try tokens.append(allocator,.{ .tok_type = .info_float, .field_name = name });
                }
            } else {
                // No comparison operator follows — treat as flag
                try tokens.append(allocator,.{ .tok_type = .info_flag, .field_name = name });
            }
            continue;
        }

        // QUAL
        if (startsWithAt(expr, i, "QUAL") and
            (i + 4 >= expr.len or !isIdentChar(expr[i + 4])))
        {
            try tokens.append(allocator,.{ .tok_type = .qual });
            i += 4;
            continue;
        }

        // N_ALT
        if (startsWithAt(expr, i, "N_ALT") and
            (i + 5 >= expr.len or !isIdentChar(expr[i + 5])))
        {
            try tokens.append(allocator,.{ .tok_type = .n_alt });
            i += 5;
            continue;
        }

        // TYPE
        if (startsWithAt(expr, i, "TYPE") and
            (i + 4 >= expr.len or !isIdentChar(expr[i + 4])))
        {
            try tokens.append(allocator,.{ .tok_type = .type_field });
            i += 4;
            continue;
        }

        // Numeric literals (int or float)
        if (isDigitOrMinus(c, i, expr, tokens)) {
            const start = i;
            var is_float = false;
            if (c == '-') i += 1;
            while (i < expr.len and (std.ascii.isDigit(expr[i]) or expr[i] == '.')) {
                if (expr[i] == '.') is_float = true;
                i += 1;
            }
            // Scientific notation
            if (i < expr.len and (expr[i] == 'e' or expr[i] == 'E')) {
                is_float = true;
                i += 1;
                if (i < expr.len and (expr[i] == '+' or expr[i] == '-')) i += 1;
                while (i < expr.len and std.ascii.isDigit(expr[i])) i += 1;
            }
            const num_str = expr[start..i];
            if (is_float) {
                const val = std.fmt.parseFloat(f64, num_str) catch 0;
                try tokens.append(allocator,.{ .tok_type = .float_literal, .float_val = val });
            } else {
                const val = std.fmt.parseInt(i64, num_str, 10) catch 0;
                try tokens.append(allocator,.{ .tok_type = .int_literal, .int_val = val });
            }
            continue;
        }

        // Unknown character — skip
        i += 1;
    }
}

/// Check if a minus sign is part of a negative numeric literal (not subtraction).
fn isDigitOrMinus(c: u8, pos: usize, expr: []const u8, tokens: *const std.ArrayListUnmanaged(Token)) bool {
    if (std.ascii.isDigit(c)) return true;
    if (c == '-') {
        // Negative number: must be followed by a digit, and preceded by an operator or start
        if (pos + 1 < expr.len and std.ascii.isDigit(expr[pos + 1])) {
            if (tokens.items.len == 0) return true;
            const prev = tokens.items[tokens.items.len - 1].tok_type;
            return switch (prev) {
                .op_gt, .op_gte, .op_lt, .op_lte, .op_eq, .op_neq, .op_and, .op_or, .op_not, .lparen => true,
                else => false,
            };
        }
    }
    return false;
}

fn isIdentChar(c: u8) bool {
    return std.ascii.isAlphanumeric(c) or c == '_';
}

fn startsWithAt(s: []const u8, offset: usize, prefix: []const u8) bool {
    if (offset + prefix.len > s.len) return false;
    return std.mem.eql(u8, s[offset..][0..prefix.len], prefix);
}

// ---------------------------------------------------------------------------
// Shunting-yard algorithm: infix tokens -> RPN
// ---------------------------------------------------------------------------

fn shuntingYard(allocator: std.mem.Allocator, tokens: []const Token, output: *std.ArrayListUnmanaged(Token)) !void {
    var op_stack: [64]Token = undefined;
    var op_sp: usize = 0;

    for (tokens) |tok| {
        switch (tok.tok_type) {
            // Values and field references go straight to output
            .int_literal, .float_literal, .string_literal, .qual, .info_int, .info_float, .info_string, .info_flag, .n_alt, .type_field => {
                try output.append(allocator,tok);
            },

            // Operators
            .op_gt, .op_gte, .op_lt, .op_lte, .op_eq, .op_neq, .op_and, .op_or, .op_not => {
                while (op_sp > 0) {
                    const top = op_stack[op_sp - 1];
                    if (top.tok_type == .lparen) break;
                    const top_prec = precedence(top.tok_type);
                    const cur_prec = precedence(tok.tok_type);
                    if (top_prec > cur_prec or (top_prec == cur_prec and !isRightAssoc(tok.tok_type))) {
                        try output.append(allocator,top);
                        op_sp -= 1;
                    } else {
                        break;
                    }
                }
                op_stack[op_sp] = tok;
                op_sp += 1;
            },

            .lparen => {
                op_stack[op_sp] = tok;
                op_sp += 1;
            },
            .rparen => {
                while (op_sp > 0 and op_stack[op_sp - 1].tok_type != .lparen) {
                    op_sp -= 1;
                    try output.append(allocator,op_stack[op_sp]);
                }
                if (op_sp > 0 and op_stack[op_sp - 1].tok_type == .lparen) {
                    op_sp -= 1; // discard the left paren
                }
            },
        }
    }

    // Pop remaining operators
    while (op_sp > 0) {
        op_sp -= 1;
        try output.append(allocator,op_stack[op_sp]);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

fn makeTestRecord(allocator: std.mem.Allocator, line: []const u8) !VcfRecord {
    var rec = VcfRecord.init(allocator);
    try rec.parseLine(line);
    return rec;
}

test "QUAL > 20 with QUAL=30 passes" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL>20", false);
    defer f.deinit();

    // chr1 100 . A T 30 PASS DP=50
    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=50");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "QUAL > 20 with QUAL=10 fails" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL>20", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t10\tPASS\tDP=50");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "INFO/DP > 30 with DP=50 passes" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/DP>30", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=50;AF=0.5");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "INFO/DP > 30 with DP=20 fails" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/DP>30", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=20;AF=0.5");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "INFO/AF > 0.01 with AF=0.5 passes" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/AF>0.01", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=50;AF=0.5");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "TYPE=snp with A->T passes" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "TYPE=\"snp\"", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "TYPE=indel with A->AT passes" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "TYPE=\"indel\"", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tAT\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "TYPE=snp with A->AT fails" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "TYPE=\"snp\"", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tAT\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "TYPE=mnp with AT->GC passes" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "TYPE=\"mnp\"", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tAT\tGC\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "N_ALT > 1 with 2 alts passes" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "N_ALT>1", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT,G\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "N_ALT > 1 with 1 alt fails" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "N_ALT>1", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "QUAL>20 && INFO/DP>30 compound passes" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL>20 && INFO/DP>30", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=50");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "QUAL>20 && INFO/DP>30 compound fails when DP low" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL>20 && INFO/DP>30", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=10");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "QUAL>20 || TYPE=indel with OR" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL>20 || TYPE=\"indel\"", false);
    defer f.deinit();

    // QUAL=10 but it is an indel
    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tAT\t10\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "!(TYPE=snp) negation" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "!(TYPE=\"snp\")", false);
    defer f.deinit();

    // This is a SNP -> negated -> should fail
    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "!(TYPE=snp) with indel passes" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "!(TYPE=\"snp\")", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tAT\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "exclude mode negates result" {
    const allocator = std.testing.allocator;
    // With exclude=true, QUAL>20 means exclude records with QUAL>20
    var f = try Filter.init(allocator, "QUAL>20", true);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\t.");
    defer rec.deinit();

    // QUAL=30 > 20 is true, but exclude negates, so record does NOT pass
    try std.testing.expect(!f.eval(&rec));
}

test "parenthesized grouping" {
    const allocator = std.testing.allocator;
    // Without parens: QUAL>20 && (INFO/DP>30 || TYPE="indel")
    var f = try Filter.init(allocator, "QUAL>20 && (INFO/DP>30 || TYPE=\"indel\")", false);
    defer f.deinit();

    // QUAL=25 > 20, DP=10 < 30, but TYPE is indel -> should pass due to grouping
    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tAT\t25\tPASS\tDP=10");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "missing QUAL with comparison fails" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL>20", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t.\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "missing INFO field with comparison fails" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/DP>30", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tAF=0.5");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "QUAL >= 30 boundary" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL>=30", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "QUAL < 30 boundary" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL<30", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "QUAL <= 30 boundary" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL<=30", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\t.");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "INFO flag presence" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/DB", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDB;DP=50");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "INFO flag absence" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/DB", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=50");
    defer rec.deinit();

    try std.testing.expect(!f.eval(&rec));
}

test "INFO string equality" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/SVTYPE=\"DEL\"", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tSVTYPE=DEL;DP=50");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "INFO numeric equality" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/DP=50", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=50");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "INFO numeric inequality" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/DP!=50", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=30");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "complex expression with multiple operators" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "QUAL>20 && INFO/DP>30 && TYPE=\"snp\"", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}

test "INFO field at end of INFO column (no trailing semicolon)" {
    const allocator = std.testing.allocator;
    var f = try Filter.init(allocator, "INFO/MQ>40", false);
    defer f.deinit();

    var rec = try makeTestRecord(allocator, "chr1\t100\t.\tA\tT\t30\tPASS\tDP=50;MQ=60");
    defer rec.deinit();

    try std.testing.expect(f.eval(&rec));
}
