const std = @import("std");
const types = @import("types.zig");
const region = @import("../core/region.zig");

const RegionIndex = region.RegionIndex;
const Biotype = types.Biotype;
const Strand = types.Strand;
const CdsPhase = types.CdsPhase;
const UtrType = types.UtrType;
const Trim = types.Trim;
const CdsEntry = types.CdsEntry;
const Gene = types.Gene;
const Exon = types.Exon;
const Utr = types.Utr;
const Transcript = types.Transcript;

const n_splice_region_intron = types.n_splice_region_intron;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

pub const GffError = error{
    ParseError,
    MissingId,
    MissingParent,
    UnknownTranscript,
    NoTranscriptsFound,
};

// ---------------------------------------------------------------------------
// IdTable — maps string IDs to numeric indices and back
// ---------------------------------------------------------------------------

pub const IdTable = struct {
    str2id: std.StringHashMap(u32),
    id2str: std.ArrayList([]const u8),
    backing_alloc: std.mem.Allocator,

    pub fn init(alloc: std.mem.Allocator) IdTable {
        return .{
            .str2id = std.StringHashMap(u32).init(alloc),
            .id2str = .empty,
            .backing_alloc = alloc,
        };
    }

    pub fn deinit(self: *IdTable) void {
        self.str2id.deinit();
        self.id2str.deinit(self.backing_alloc);
    }

    /// Register a name, returning its numeric ID. If already registered, returns
    /// the existing ID. The name is duped into `arena` on first registration.
    pub fn register(self: *IdTable, arena: std.mem.Allocator, name: []const u8) !u32 {
        if (self.str2id.get(name)) |existing_id| {
            return existing_id;
        }
        const id: u32 = @intCast(self.id2str.items.len);
        const duped = try arena.dupe(u8, name);
        try self.id2str.append(self.backing_alloc, duped);
        try self.str2id.put(duped, id);
        return id;
    }

    pub fn getString(self: *const IdTable, id: u32) []const u8 {
        return self.id2str.items[id];
    }
};

// ---------------------------------------------------------------------------
// Temporary feature record (replaces C ftr_t)
// ---------------------------------------------------------------------------

const Feature = struct {
    ftype: FeatureType,
    beg: u32,
    end: u32,
    trid: u32,
    strand: Strand,
    phase: CdsPhase,
    iseq: u32,
};

const FeatureType = enum {
    cds,
    exon,
    utr3,
    utr5,
};

// ---------------------------------------------------------------------------
// Parsed attributes from column 9
// ---------------------------------------------------------------------------

const Attributes = struct {
    id: ?[]const u8 = null,
    parent: ?[]const u8 = null,
    name: ?[]const u8 = null,
    biotype: ?[]const u8 = null,
    is_gene_from_id: bool = false,
};

// ---------------------------------------------------------------------------
// GffParser
// ---------------------------------------------------------------------------

pub const GffParser = struct {
    allocator: std.mem.Allocator,
    arena: std.heap.ArenaAllocator,

    // Indexes built during parse
    idx_cds: RegionIndex(*CdsEntry),
    idx_utr: RegionIndex(*Utr),
    idx_exon: RegionIndex(*Exon),
    idx_tscript: RegionIndex(*Transcript),

    // Transcript and gene storage
    transcripts: std.AutoHashMap(u32, *Transcript),
    genes: std.AutoHashMap(u32, *Gene),

    // ID tables
    tscript_ids: IdTable,
    gene_ids: IdTable,

    // Sequence names
    seq_names: std.StringHashMap(u32),
    seq_list: std.ArrayList([]const u8) = .empty,

    // Options
    verbosity: i32,
    force: bool,

    pub fn init(allocator: std.mem.Allocator) GffParser {
        return .{
            .allocator = allocator,
            .arena = std.heap.ArenaAllocator.init(allocator),
            .idx_cds = RegionIndex(*CdsEntry).init(allocator),
            .idx_utr = RegionIndex(*Utr).init(allocator),
            .idx_exon = RegionIndex(*Exon).init(allocator),
            .idx_tscript = RegionIndex(*Transcript).init(allocator),
            .transcripts = std.AutoHashMap(u32, *Transcript).init(allocator),
            .genes = std.AutoHashMap(u32, *Gene).init(allocator),
            .tscript_ids = IdTable.init(allocator),
            .gene_ids = IdTable.init(allocator),
            .seq_names = std.StringHashMap(u32).init(allocator),
            .seq_list = .empty,
            .verbosity = 0,
            .force = false,
        };
    }

    pub fn deinit(self: *GffParser) void {
        self.idx_cds.deinit();
        self.idx_utr.deinit();
        self.idx_exon.deinit();
        self.idx_tscript.deinit();
        self.transcripts.deinit();
        self.genes.deinit();
        self.tscript_ids.deinit();
        self.gene_ids.deinit();
        self.seq_names.deinit();
        self.seq_list.deinit(self.allocator);
        self.arena.deinit();
    }

    /// Look up the string for a transcript/gene ID.
    pub fn id2string(self: *const GffParser, id: u32) []const u8 {
        return self.tscript_ids.getString(id);
    }

    /// Check whether a sequence name was seen during parsing.
    pub fn hasSeq(self: *const GffParser, chr: []const u8) bool {
        return self.seq_names.contains(chr);
    }

    /// Return the number of distinct sequences seen.
    pub fn nseq(self: *const GffParser) usize {
        return self.seq_list.items.len;
    }

    // -----------------------------------------------------------------------
    // Main parse entry point
    // -----------------------------------------------------------------------

    /// Parse a GFF3 file. For now uses std.fs; no bgzf support.
    pub fn parse(self: *GffParser, fname: []const u8) !void {
        const file = try std.fs.cwd().openFile(fname, .{});
        defer file.close();
        try self.parseReader(file.reader());
    }

    /// Parse GFF3 data from any reader (useful for testing with in-memory data).
    pub fn parseReader(self: *GffParser, reader: anytype) !void {
        const arena_alloc = self.arena.allocator();

        // Phase 1: collect features (exon/CDS/UTR) and register genes + transcripts
        var features: std.ArrayList(Feature) = .empty;
        defer features.deinit(self.allocator);

        var buf: [8192]u8 = undefined;
        while (true) {
            const line = reader.readUntilDelimiter(&buf, '\n') catch |err| switch (err) {
                error.EndOfStream => break,
                else => return err,
            };
            // Strip trailing \r if present
            const clean = if (line.len > 0 and line[line.len - 1] == '\r')
                line[0 .. line.len - 1]
            else
                line;

            if (clean.len == 0) continue;
            if (clean[0] == '#') continue;

            try self.parseLine(arena_alloc, clean, &features);
        }

        // Phase 2: register features (CDS, exon, UTR) into their indexes
        for (features.items) |*ftr| {
            const tr = self.transcripts.get(ftr.trid) orelse continue;

            // Extend transcript bounds if feature is outside
            if (ftr.beg < tr.beg) tr.beg = ftr.beg;
            if (ftr.end > tr.end) tr.end = ftr.end;

            tr.used = true;
            if (tr.gene) |gene| gene.used = true;

            switch (ftr.ftype) {
                .cds => try self.registerCds(arena_alloc, ftr, tr),
                .exon => try self.registerExon(arena_alloc, ftr, tr),
                .utr3 => try self.registerUtr(arena_alloc, ftr, tr, .prime3),
                .utr5 => try self.registerUtr(arena_alloc, ftr, tr, .prime5),
            }
        }

        // Phase 3: finalize transcripts — sort CDS, set offsets, build idx_tscript and idx_cds
        try self.finalizeTscripts(arena_alloc);
    }

    // -----------------------------------------------------------------------
    // Internal: sequence name registration
    // -----------------------------------------------------------------------

    fn registerSeq(self: *GffParser, arena_alloc: std.mem.Allocator, chr: []const u8) !u32 {
        if (self.seq_names.get(chr)) |id| return id;
        const id: u32 = @intCast(self.seq_list.items.len);
        const duped = try arena_alloc.dupe(u8, chr);
        try self.seq_list.append(self.allocator, duped);
        try self.seq_names.put(duped, id);
        return id;
    }

    // -----------------------------------------------------------------------
    // Internal: parse a single GFF line
    // -----------------------------------------------------------------------

    fn parseLine(
        self: *GffParser,
        arena_alloc: std.mem.Allocator,
        line: []const u8,
        features: *std.ArrayList(Feature),
    ) !void {
        // Split into tab-separated columns. GFF has 9 columns.
        var cols: [9][]const u8 = undefined;
        var col_count: usize = 0;
        var rest = line;
        while (col_count < 9) : (col_count += 1) {
            if (col_count == 8) {
                // Last column gets everything remaining
                cols[col_count] = rest;
                col_count += 1;
                break;
            }
            if (std.mem.indexOfScalar(u8, rest, '\t')) |tab_pos| {
                cols[col_count] = rest[0..tab_pos];
                rest = rest[tab_pos + 1 ..];
            } else {
                cols[col_count] = rest;
                col_count += 1;
                break;
            }
        }
        if (col_count < 9) return; // malformed line, skip silently

        const chr = cols[0];
        const type_str = cols[2];
        const beg_str = cols[3];
        const end_str = cols[4];
        const strand_str = cols[6];
        const phase_str = cols[7];
        const attrs_str = cols[8];

        // Parse begin/end (1-based in file -> 0-based internally)
        const beg = (std.fmt.parseUnsigned(u32, beg_str, 10) catch return) -| 1;
        const end = (std.fmt.parseUnsigned(u32, end_str, 10) catch return) -| 1;

        // Parse strand
        const strand: Strand = if (strand_str.len > 0) switch (strand_str[0]) {
            '+' => .forward,
            '-' => .reverse,
            else => .unknown,
        } else .unknown;

        // Parse phase
        const phase: CdsPhase = if (phase_str.len > 0) switch (phase_str[0]) {
            '0' => .phase0,
            '1' => .phase1,
            '2' => .phase2,
            else => .unknown,
        } else .unknown;

        // Parse attributes (column 9)
        const attrs = parseAttributes(attrs_str);

        // Determine what kind of line this is
        var is_gene_line = std.mem.eql(u8, type_str, "gene");
        // The ID=gene: prefix also marks a gene line
        if (attrs.is_gene_from_id) is_gene_line = true;

        const is_special = isSpecialType(type_str);

        if (is_gene_line or attrs.parent == null) {
            // Gene line (or top-level feature with no parent)
            try self.handleGene(arena_alloc, chr, beg, end, strand, attrs);
            return;
        }

        if (is_special) |ftype| {
            // CDS / exon / UTR line
            const parent = attrs.parent orelse return;
            const trid = try self.tscript_ids.register(arena_alloc, parent);
            const iseq = try self.registerSeq(arena_alloc, chr);

            try features.append(self.allocator, .{
                .ftype = ftype,
                .beg = beg,
                .end = end,
                .trid = trid,
                .strand = strand,
                .phase = phase,
                .iseq = iseq,
            });
            return;
        }

        // Otherwise: transcript / mRNA line
        try self.handleTranscript(arena_alloc, chr, beg, end, strand, type_str, attrs);
    }

    /// Check if the type column is one of the special structural types.
    fn isSpecialType(type_str: []const u8) ?FeatureType {
        if (std.mem.eql(u8, type_str, "CDS")) return .cds;
        if (std.mem.eql(u8, type_str, "exon")) return .exon;
        if (std.mem.eql(u8, type_str, "three_prime_UTR")) return .utr3;
        if (std.mem.eql(u8, type_str, "five_prime_UTR")) return .utr5;
        return null;
    }

    // -----------------------------------------------------------------------
    // Internal: attribute parsing
    // -----------------------------------------------------------------------

    fn parseAttributes(attrs_str: []const u8) Attributes {
        var attrs = Attributes{};
        var remaining = attrs_str;

        while (remaining.len > 0) {
            // Find the end of this key=value pair
            const semi_pos = std.mem.indexOfScalar(u8, remaining, ';') orelse remaining.len;
            const pair = remaining[0..semi_pos];
            remaining = if (semi_pos < remaining.len) remaining[semi_pos + 1 ..] else &[_]u8{};

            if (std.mem.startsWith(u8, pair, "ID=")) {
                var val = pair[3..];
                // Strip optional prefixes: "gene:", "transcript:"
                if (std.mem.startsWith(u8, val, "gene:")) {
                    val = val[5..];
                    attrs.is_gene_from_id = true;
                } else if (std.mem.startsWith(u8, val, "transcript:")) {
                    val = val[11..];
                }
                attrs.id = val;
            } else if (std.mem.startsWith(u8, pair, "Parent=")) {
                var val = pair[7..];
                if (std.mem.startsWith(u8, val, "gene:")) {
                    val = val[5..];
                } else if (std.mem.startsWith(u8, val, "transcript:")) {
                    val = val[11..];
                }
                attrs.parent = val;
            } else if (std.mem.startsWith(u8, pair, "Name=")) {
                attrs.name = pair[5..];
            } else if (std.mem.startsWith(u8, pair, "biotype=")) {
                attrs.biotype = pair[8..];
            } else if (std.mem.startsWith(u8, pair, "gene_biotype=")) {
                attrs.biotype = pair[13..];
            } else if (std.mem.startsWith(u8, pair, "transcript_biotype=")) {
                attrs.biotype = pair[19..];
            }
        }

        return attrs;
    }

    // -----------------------------------------------------------------------
    // Internal: gene handling
    // -----------------------------------------------------------------------

    fn handleGene(
        self: *GffParser,
        arena_alloc: std.mem.Allocator,
        chr: []const u8,
        beg: u32,
        end: u32,
        strand: Strand,
        attrs: Attributes,
    ) !void {
        const id_str = attrs.id orelse return; // no ID, skip
        const gene_id = try self.gene_ids.register(arena_alloc, id_str);

        // Get or create gene
        const gop = try self.genes.getOrPut(gene_id);
        if (!gop.found_existing) {
            const gene = try arena_alloc.create(Gene);
            gene.* = .{
                .name = null,
                .iseq = 0,
                .id = gene_id,
                .beg = beg,
                .end = end,
                .strand = strand,
                .used = false,
            };
            gop.value_ptr.* = gene;
        }

        const gene = gop.value_ptr.*;

        // If gene already has a name, this is a duplicate — skip
        if (gene.name != null) return;

        gene.iseq = try self.registerSeq(arena_alloc, chr);
        gene.beg = beg;
        gene.end = end;
        gene.strand = strand;

        // Set gene name: prefer Name= attribute, fall back to ID
        const name_src = attrs.name orelse id_str;
        const duped = try arena_alloc.dupeZ(u8, name_src);
        gene.name = duped;
    }

    // -----------------------------------------------------------------------
    // Internal: transcript handling
    // -----------------------------------------------------------------------

    fn handleTranscript(
        self: *GffParser,
        arena_alloc: std.mem.Allocator,
        _: []const u8, // chr (unused directly; gene carries iseq)
        beg: u32,
        end: u32,
        strand: Strand,
        type_str: []const u8,
        attrs: Attributes,
    ) !void {
        // Determine biotype: first from biotype= attribute, then from type column
        var biotype: ?Biotype = null;
        if (attrs.biotype) |bt_str| {
            biotype = types.biotype_map.get(bt_str);
        }
        if (biotype == null) {
            biotype = types.biotype_map.get(type_str);
        }
        // Also accept mRNA (case insensitive) as protein_coding
        if (biotype == null) {
            if (std.ascii.eqlIgnoreCase(type_str, "mrna")) {
                biotype = .protein_coding;
            }
        }
        if (biotype == null) return; // unknown biotype, skip

        // A special structural type in column 3 should not be treated as transcript
        if (biotype.?.isSpecial()) return;

        const id_str = attrs.id orelse return;
        const parent_str = attrs.parent orelse return;

        const trid = try self.tscript_ids.register(arena_alloc, id_str);
        const gene_id = try self.gene_ids.register(arena_alloc, parent_str);

        // Ensure gene exists
        const gene_gop = try self.genes.getOrPut(gene_id);
        if (!gene_gop.found_existing) {
            const gene = try arena_alloc.create(Gene);
            gene.* = .{
                .name = null,
                .iseq = 0,
                .id = gene_id,
                .beg = beg,
                .end = end,
                .strand = strand,
                .used = false,
            };
            gene_gop.value_ptr.* = gene;
        }

        const tr = try arena_alloc.create(Transcript);
        tr.* = Transcript.init(arena_alloc);
        tr.id = trid;
        tr.beg = beg;
        tr.end = end;
        tr.strand = strand;
        tr.biotype = biotype.?;
        tr.gene = gene_gop.value_ptr.*;

        try self.transcripts.put(trid, tr);
    }

    // -----------------------------------------------------------------------
    // Internal: register CDS/exon/UTR into indexes
    // -----------------------------------------------------------------------

    fn registerCds(
        _: *GffParser,
        arena_alloc: std.mem.Allocator,
        ftr: *const Feature,
        tr: *Transcript,
    ) !void {
        const cds = try arena_alloc.create(CdsEntry);
        cds.* = .{
            .tr = tr,
            .beg = ftr.beg,
            .pos = 0,
            .len = ftr.end - ftr.beg + 1,
            .icds = 0,
            .phase = ftr.phase,
        };
        try tr.cds.append(tr.allocator, cds);
        // idx_cds insertion happens in finalizeTscripts after sorting
    }

    fn registerExon(
        self: *GffParser,
        arena_alloc: std.mem.Allocator,
        ftr: *const Feature,
        tr: *Transcript,
    ) !void {
        const exon_entry = try arena_alloc.create(Exon);
        exon_entry.* = .{
            .beg = ftr.beg,
            .end = ftr.end,
            .tr = tr,
        };

        const chr = self.seqName(tr);
        // Extend exon region by splice region intron padding
        const padded_beg = ftr.beg -| n_splice_region_intron;
        const padded_end = ftr.end + n_splice_region_intron;
        try self.idx_exon.insert(chr, padded_beg, padded_end, exon_entry);
    }

    fn registerUtr(
        self: *GffParser,
        arena_alloc: std.mem.Allocator,
        ftr: *const Feature,
        tr: *Transcript,
        which: UtrType,
    ) !void {
        const utr = try arena_alloc.create(Utr);
        utr.* = .{
            .which = which,
            .beg = ftr.beg,
            .end = ftr.end,
            .tr = tr,
        };

        const chr = self.seqName(tr);
        try self.idx_utr.insert(chr, ftr.beg, ftr.end, utr);
    }

    /// Get the sequence name string for a transcript (via its gene).
    fn seqName(self: *const GffParser, tr: *const Transcript) []const u8 {
        const gene = tr.gene orelse return "";
        return self.seq_list.items[gene.iseq];
    }

    // -----------------------------------------------------------------------
    // Internal: finalize transcripts
    // -----------------------------------------------------------------------

    fn finalizeTscripts(self: *GffParser, arena_alloc: std.mem.Allocator) !void {
        _ = arena_alloc;
        var it = self.transcripts.iterator();
        while (it.next()) |entry| {
            const tr = entry.value_ptr.*;
            const chr = self.seqName(tr);

            // Register transcript in idx_tscript
            try self.idx_tscript.insert(chr, tr.beg, tr.end, tr);

            const cds_items = tr.cds.items;
            if (cds_items.len == 0) continue;

            // Sort CDS by position
            std.mem.sort(*CdsEntry, cds_items, {}, struct {
                fn lessThan(_: void, a: *CdsEntry, b: *CdsEntry) bool {
                    return if (a.beg != b.beg) a.beg < b.beg else a.len < b.len;
                }
            }.lessThan);

            // Trim 5' end based on phase (forward strand: first CDS; reverse: last CDS)
            if (tr.strand == .forward) {
                self.trim5PrimeFwd(tr);
            } else if (tr.strand == .reverse) {
                self.trim5PrimeRev(tr);
            } else {
                continue; // unknown strand, skip
            }

            // Check total length; trim 3' if not multiple of 3
            var total_len: u32 = 0;
            for (cds_items) |cds| {
                total_len += cds.len;
            }
            if (total_len % 3 != 0) {
                tr.trim = .prime3;
                self.trim3Prime(tr, &total_len);
            }

            // Set icds and pos offsets, then insert into idx_cds
            var pos: u32 = 0;
            for (cds_items, 0..) |cds, i| {
                cds.icds = @intCast(i);
                cds.pos = pos;
                pos += cds.len;

                try self.idx_cds.insert(chr, cds.beg, cds.beg + cds.len -| 1, cds);
            }
        }
    }

    fn trim5PrimeFwd(self: *GffParser, tr: *Transcript) void {
        _ = self;
        const cds_items = tr.cds.items;
        if (cds_items.len == 0) return;
        const first = cds_items[0];
        if (first.phase != .unknown) {
            const phase_val = @intFromEnum(first.phase);
            if (phase_val > 0) tr.trim = .prime5;
            first.beg += phase_val;
            first.len -= phase_val;
            first.phase = .phase0;
        }
    }

    fn trim5PrimeRev(self: *GffParser, tr: *Transcript) void {
        _ = self;
        const cds_items = tr.cds.items;
        if (cds_items.len == 0) return;
        var i: usize = cds_items.len - 1;
        const last = cds_items[i];
        if (last.phase != .unknown) {
            var phase_remaining: u32 = @intFromEnum(last.phase);
            if (phase_remaining > 0) tr.trim = .prime5;
            // Phase can span multiple CDS segments
            while (phase_remaining > 0) {
                const cds = cds_items[i];
                if (phase_remaining >= cds.len) {
                    phase_remaining -= cds.len;
                    cds.phase = .phase0;
                    cds.len = 0;
                    if (i == 0) break;
                    i -= 1;
                } else {
                    cds.len -= phase_remaining;
                    cds.phase = .phase0;
                    phase_remaining = 0;
                }
            }
        }
    }

    fn trim3Prime(self: *GffParser, tr: *Transcript, total_len: *u32) void {
        _ = self;
        const cds_items = tr.cds.items;
        if (tr.strand == .forward) {
            // Trim from the last CDS
            var i: usize = cds_items.len;
            while (i > 0 and total_len.* % 3 != 0) {
                i -= 1;
                const cds = cds_items[i];
                const dlen = if (cds.len >= total_len.* % 3) total_len.* % 3 else cds.len;
                cds.len -= dlen;
                total_len.* -= dlen;
            }
        } else if (tr.strand == .reverse) {
            // Trim from the first CDS
            var i: usize = 0;
            while (i < cds_items.len and total_len.* % 3 != 0) {
                const cds = cds_items[i];
                const dlen = if (cds.len >= total_len.* % 3) total_len.* % 3 else cds.len;
                cds.len -= dlen;
                cds.beg += dlen;
                total_len.* -= dlen;
                i += 1;
            }
        }
    }
};

// ===========================================================================
// Tests
// ===========================================================================

test "parse minimal GFF3" {
    const testing = std.testing;
    const alloc = testing.allocator;

    const gff_data =
        "##gff-version 3\n" ++
        "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;biotype=protein_coding\n" ++
        "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=tx1;Parent=gene1;biotype=protein_coding\n" ++
        "chr1\t.\texon\t1000\t1200\t.\t+\t.\tParent=tx1\n" ++
        "chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tParent=tx1\n" ++
        "chr1\t.\texon\t1500\t2000\t.\t+\t.\tParent=tx1\n" ++
        "chr1\t.\tCDS\t1500\t2000\t.\t+\t0\tParent=tx1\n";

    var parser = GffParser.init(alloc);
    defer parser.deinit();

    var stream = std.io.fixedBufferStream(gff_data);
    try parser.parseReader(stream.reader());

    // 1 gene
    try testing.expectEqual(@as(usize, 1), parser.genes.count());

    // 1 transcript
    try testing.expectEqual(@as(usize, 1), parser.transcripts.count());

    // Verify transcript has 2 CDS entries
    var tr_it = parser.transcripts.iterator();
    const tr_entry = tr_it.next().?;
    const tr = tr_entry.value_ptr.*;
    try testing.expectEqual(@as(usize, 2), tr.cds.items.len);

    // CDS entries should be sorted by position (0-based: 999..1199 and 1499..1999)
    const cds0 = tr.cds.items[0];
    const cds1 = tr.cds.items[1];
    try testing.expect(cds0.beg < cds1.beg);
    try testing.expectEqual(@as(u32, 999), cds0.beg);
    try testing.expectEqual(@as(u32, 201), cds0.len);
    try testing.expectEqual(@as(u32, 1499), cds1.beg);
    try testing.expectEqual(@as(u32, 501), cds1.len);

    // CDS icds should be set
    try testing.expectEqual(@as(u30, 0), cds0.icds);
    try testing.expectEqual(@as(u30, 1), cds1.icds);

    // CDS pos offsets
    try testing.expectEqual(@as(u32, 0), cds0.pos);
    try testing.expectEqual(@as(u32, 201), cds1.pos);

    // idx_cds should have entries
    var overlap = parser.idx_cds.overlap("chr1", 999, 1199);
    const first = overlap.next();
    try testing.expect(first != null);
    try testing.expectEqual(@as(u32, 999), first.?.payload.beg);

    // Second CDS region
    var overlap2 = parser.idx_cds.overlap("chr1", 1499, 1999);
    const second = overlap2.next();
    try testing.expect(second != null);
    try testing.expectEqual(@as(u32, 1499), second.?.payload.beg);

    // idx_tscript should have the transcript
    var ts_overlap = parser.idx_tscript.overlap("chr1", 999, 1999);
    const ts_hit = ts_overlap.next();
    try testing.expect(ts_hit != null);

    // Sequence name should be registered
    try testing.expect(parser.hasSeq("chr1"));
    try testing.expectEqual(@as(usize, 1), parser.nseq());

    // Gene should be linked to transcript
    try testing.expect(tr.gene != null);
    try testing.expectEqual(Strand.forward, tr.strand);
    try testing.expectEqual(Biotype.protein_coding, tr.biotype);
}

test "parse attributes" {
    const attrs = GffParser.parseAttributes("ID=gene:ENSG00000001;Name=TP53;biotype=protein_coding");
    try std.testing.expect(attrs.id != null);
    try std.testing.expectEqualStrings("ENSG00000001", attrs.id.?);
    try std.testing.expect(attrs.is_gene_from_id);
    try std.testing.expectEqualStrings("TP53", attrs.name.?);
    try std.testing.expectEqualStrings("protein_coding", attrs.biotype.?);
}

test "parse attributes with transcript prefix" {
    const attrs = GffParser.parseAttributes("ID=transcript:ENST00000001;Parent=gene:ENSG00000001;biotype=lncRNA");
    try std.testing.expectEqualStrings("ENST00000001", attrs.id.?);
    try std.testing.expect(!attrs.is_gene_from_id);
    try std.testing.expectEqualStrings("ENSG00000001", attrs.parent.?);
    try std.testing.expectEqualStrings("lncRNA", attrs.biotype.?);
}

test "IdTable register and getString" {
    const alloc = std.testing.allocator;
    var arena = std.heap.ArenaAllocator.init(alloc);
    defer arena.deinit();

    var tbl = IdTable.init(alloc);
    defer tbl.deinit();

    const id0 = try tbl.register(arena.allocator(), "gene1");
    const id1 = try tbl.register(arena.allocator(), "gene2");
    const id0_again = try tbl.register(arena.allocator(), "gene1");

    try std.testing.expectEqual(@as(u32, 0), id0);
    try std.testing.expectEqual(@as(u32, 1), id1);
    try std.testing.expectEqual(@as(u32, 0), id0_again);
    try std.testing.expectEqualStrings("gene1", tbl.getString(0));
    try std.testing.expectEqualStrings("gene2", tbl.getString(1));
}

test "reverse strand CDS trimming" {
    const testing = std.testing;
    const alloc = testing.allocator;

    // A reverse-strand transcript with phase=1 on last CDS
    const gff_data =
        "##gff-version 3\n" ++
        "chr1\t.\tgene\t1000\t3000\t.\t-\t.\tID=gene1;biotype=protein_coding\n" ++
        "chr1\t.\tmRNA\t1000\t3000\t.\t-\t.\tID=tx1;Parent=gene1;biotype=protein_coding\n" ++
        "chr1\t.\tCDS\t1000\t1200\t.\t-\t0\tParent=tx1\n" ++
        "chr1\t.\tCDS\t2000\t3000\t.\t-\t1\tParent=tx1\n";

    var parser = GffParser.init(alloc);
    defer parser.deinit();

    var stream = std.io.fixedBufferStream(gff_data);
    try parser.parseReader(stream.reader());

    var tr_it = parser.transcripts.iterator();
    const tr = tr_it.next().?.value_ptr.*;

    // Should have 2 CDS entries
    try testing.expectEqual(@as(usize, 2), tr.cds.items.len);

    // The last CDS (highest position = index 1 after sort) had phase=1,
    // so its length should be reduced by 1
    const last_cds = tr.cds.items[1];
    try testing.expectEqual(@as(u32, 1000), last_cds.len); // was 1001, minus 1 for phase
}

test "gene_biotype and transcript_biotype attributes" {
    const attrs = GffParser.parseAttributes("ID=gene1;gene_biotype=protein_coding");
    try std.testing.expectEqualStrings("protein_coding", attrs.biotype.?);

    const attrs2 = GffParser.parseAttributes("ID=tx1;Parent=gene1;transcript_biotype=lncRNA");
    try std.testing.expectEqualStrings("lncRNA", attrs2.biotype.?);
}

test "skip comment and blank lines" {
    const alloc = std.testing.allocator;

    const gff_data =
        "##gff-version 3\n" ++
        "# this is a comment\n" ++
        "\n" ++
        "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;biotype=protein_coding\n" ++
        "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=tx1;Parent=gene1;biotype=protein_coding\n";

    var parser = GffParser.init(alloc);
    defer parser.deinit();

    var stream = std.io.fixedBufferStream(gff_data);
    try parser.parseReader(stream.reader());

    try std.testing.expectEqual(@as(usize, 1), parser.genes.count());
    try std.testing.expectEqual(@as(usize, 1), parser.transcripts.count());
}
