// Main CSQ processing pipeline — haplotype-aware consequence caller.
//
// This is the Zig port of the core logic from csq.c, wiring together
// the GFF annotation index, haplotype tree, splice analysis, and
// translation modules.
//
// The CsqContext replaces the C `args_t` struct and owns the full
// lifecycle: init -> process (per-record) -> flush -> deinit.

const std = @import("std");
const format = @import("format.zig");

// Re-export format types used in public API
pub const Vcsq = format.Vcsq;
pub const CsqType = format.CsqType;
pub const FormatOptions = format.FormatOptions;
pub const formatVcsq = format.formatVcsq;
pub const formatVcsqList = format.formatVcsqList;
pub const formatAaPrediction = format.formatAaPrediction;

// Import consequence type constants
pub const CSQ_PRINTED_UPSTREAM = format.CSQ_PRINTED_UPSTREAM;
pub const CSQ_SYNONYMOUS_VARIANT = format.CSQ_SYNONYMOUS_VARIANT;
pub const CSQ_MISSENSE_VARIANT = format.CSQ_MISSENSE_VARIANT;
pub const CSQ_STOP_LOST = format.CSQ_STOP_LOST;
pub const CSQ_STOP_GAINED = format.CSQ_STOP_GAINED;
pub const CSQ_INFRAME_DELETION = format.CSQ_INFRAME_DELETION;
pub const CSQ_INFRAME_INSERTION = format.CSQ_INFRAME_INSERTION;
pub const CSQ_FRAMESHIFT_VARIANT = format.CSQ_FRAMESHIFT_VARIANT;
pub const CSQ_SPLICE_ACCEPTOR = format.CSQ_SPLICE_ACCEPTOR;
pub const CSQ_SPLICE_DONOR = format.CSQ_SPLICE_DONOR;
pub const CSQ_START_LOST = format.CSQ_START_LOST;
pub const CSQ_SPLICE_REGION = format.CSQ_SPLICE_REGION;
pub const CSQ_STOP_RETAINED = format.CSQ_STOP_RETAINED;
pub const CSQ_UTR5 = format.CSQ_UTR5;
pub const CSQ_UTR3 = format.CSQ_UTR3;
pub const CSQ_NON_CODING = format.CSQ_NON_CODING;
pub const CSQ_INTRON = format.CSQ_INTRON;
pub const CSQ_INFRAME_ALTERING = format.CSQ_INFRAME_ALTERING;
pub const CSQ_UPSTREAM_STOP = format.CSQ_UPSTREAM_STOP;
pub const CSQ_INCOMPLETE_CDS = format.CSQ_INCOMPLETE_CDS;
pub const CSQ_CODING_SEQUENCE = format.CSQ_CODING_SEQUENCE;
pub const CSQ_ELONGATION = format.CSQ_ELONGATION;
pub const CSQ_TRUNCATION = format.CSQ_TRUNCATION;
pub const CSQ_START_RETAINED = format.CSQ_START_RETAINED;
pub const CSQ_COMPOUND = format.CSQ_COMPOUND;
pub const CSQ_START_STOP = format.CSQ_START_STOP;

/// Padding bases on each side of reference sequences to avoid boundary effects.
pub const N_REF_PAD = 10;

/// Sentinel position meaning "flush everything".
pub const POS_MAX: u32 = std.math.maxInt(u32);

// ---------------------------------------------------------------------------
// Phase handling
// ---------------------------------------------------------------------------

pub const Phase = enum {
    /// Require phased genotypes (error on unphased hets).
    require,
    /// Merge all GTs into a single haplotype.
    merge,
    /// Take GTs as-is regardless of phasing.
    as_is,
    /// Skip unphased heterozygous sites.
    skip,
    /// Create non-reference haplotypes if possible.
    non_ref,
    /// Drop genotypes entirely (e.g. --samples -).
    drop_gt,
};

// ---------------------------------------------------------------------------
// Placeholder VCF record (until full htslib wrapper exists)
// ---------------------------------------------------------------------------

/// Minimal VCF record representation for the CSQ pipeline.
/// This will be replaced by a proper htslib binding wrapper.
pub const VcfRecord = struct {
    pos: u32,
    rid: i32,
    n_allele: u32,
    alleles: []const []const u8,
    rlen: u32,

    pub fn seqname(self: *const VcfRecord) []const u8 {
        // TODO: get from header via htslib binding
        _ = self;
        return "unknown";
    }
};

// ---------------------------------------------------------------------------
// VCF record buffering types (port of vrec_t / vbuf_t / csq_t)
// ---------------------------------------------------------------------------

/// A single VCF record with attached consequences, analogous to vrec_t.
pub const Vrec = struct {
    /// The VCF record. Owned (or swapped) pointer.
    rec: ?*const VcfRecord = null,
    /// Bitmask of sample consequences (first/second haplotype interleaved).
    fmt_bm: ?[]u32 = null,
    /// Number of fmt integers per sample actually used.
    nfmt: u32 = 0,
    /// Consequences attached to this record.
    vcsqs: std.ArrayList(Vcsq) = .empty,

    pub fn deinit(self: *Vrec, allocator: std.mem.Allocator) void {
        self.vcsqs.deinit(allocator);
        if (self.fmt_bm) |bm| {
            allocator.free(bm);
            self.fmt_bm = null;
        }
    }
};

/// Buffer of VCF records at the same position, analogous to vbuf_t.
pub const Vbuf = struct {
    vrecs: std.ArrayList(Vrec) = .empty,
    /// Maximum transcript end overlapping this buffer; controls when we can flush.
    keep_until: u32 = 0,

    pub fn deinit(self: *Vbuf, allocator: std.mem.Allocator) void {
        for (self.vrecs.items) |*vrec| {
            vrec.deinit(allocator);
        }
        self.vrecs.deinit(allocator);
    }

    pub fn pos(self: *const Vbuf) ?u32 {
        if (self.vrecs.items.len > 0) {
            if (self.vrecs.items[0].rec) |rec| return rec.pos;
        }
        return null;
    }
};

/// A consequence tied to a haplotype, analogous to csq_t.
pub const Csq = struct {
    pos: u32 = 0,
    /// Back-pointer to the vrec this consequence is attached to.
    vrec_idx: ?usize = null,
    /// Index of this consequence within the vrec's vcsq list.
    csq_idx: ?usize = null,
    /// The consequence type/annotation data.
    vcsq: Vcsq = .{},
};

// ---------------------------------------------------------------------------
// Ring buffer for Vbuf pointers
// ---------------------------------------------------------------------------

/// Simple ring buffer for managing ordered Vbuf entries.
fn RingBuffer(comptime T: type) type {
    return struct {
        const Self = @This();
        items: []T,
        head: usize = 0,
        len: usize = 0,
        allocator: std.mem.Allocator,

        pub fn init(allocator: std.mem.Allocator, capacity: usize) !Self {
            const items = try allocator.alloc(T, capacity);
            return .{
                .items = items,
                .allocator = allocator,
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.items);
        }

        /// Append an element, growing if necessary. Returns the index.
        pub fn append(self: *Self, value: T) !usize {
            if (self.len == self.items.len) {
                try self.grow();
            }
            const idx = (self.head + self.len) % self.items.len;
            self.items[idx] = value;
            self.len += 1;
            return idx;
        }

        /// Remove and return the first element.
        pub fn shift(self: *Self) ?T {
            if (self.len == 0) return null;
            const val = self.items[self.head];
            self.head = (self.head + 1) % self.items.len;
            self.len -= 1;
            return val;
        }

        /// Peek at the first element without removing.
        pub fn front(self: *const Self) ?T {
            if (self.len == 0) return null;
            return self.items[self.head];
        }

        /// Peek at the last element.
        pub fn last(self: *const Self) ?T {
            if (self.len == 0) return null;
            const idx = (self.head + self.len - 1) % self.items.len;
            return self.items[idx];
        }

        /// Get the k-th element (0-indexed from head).
        pub fn kth(self: *const Self, k: usize) T {
            std.debug.assert(k < self.len);
            return self.items[(self.head + k) % self.items.len];
        }

        fn grow(self: *Self) !void {
            const old_cap = self.items.len;
            const new_cap = if (old_cap == 0) 8 else old_cap * 2;
            const new_items = try self.allocator.alloc(T, new_cap);
            // Copy elements in order
            for (0..self.len) |k| {
                new_items[k] = self.items[(self.head + k) % old_cap];
            }
            self.allocator.free(self.items);
            self.items = new_items;
            self.head = 0;
        }
    };
}

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

pub const Options = struct {
    gff_fname: []const u8,
    fasta_fname: []const u8 = "",
    phase: Phase = .require,
    local_csq: bool = false,
    verbosity: i32 = 1,
    force: bool = false,
    bcsq_tag: []const u8 = "BCSQ",
    ncsq2_max: u32 = 15 * 2,
    gencode_id: i32 = 0,
    brief_predictions: u32 = 0,
};

// ---------------------------------------------------------------------------
// CsqContext — main state container, replaces args_t
// ---------------------------------------------------------------------------

pub const CsqContext = struct {
    allocator: std.mem.Allocator,

    // GFF annotation — TODO: replace with gff_mod.GffParser once available
    // gff: gff_mod.GffParser,

    // Region indexes — TODO: populated from GFF parser
    // idx_cds: *RegionIndex,
    // idx_utr: *RegionIndex,
    // idx_exon: *RegionIndex,
    // idx_tscript: *RegionIndex,

    // Haplotype processing — TODO: wire haplotype.HapContext
    // hap_ctx: haplotype.HapContext,

    // VCF record buffering
    pos2vbuf: std.AutoHashMap(u32, usize), // pos -> ring buffer index (for existence check)
    vcf_rbuf: RingBuffer(*Vbuf),

    // Transcript management — TODO: wire gff_types.Transcript
    // active_transcripts: heap of active transcripts for flushing
    // rm_transcripts: list of transcripts pending cleanup

    // CSQ buffer for non-CDS consequences
    csq_buf: std.ArrayList(Csq),

    // Output
    output: std.ArrayList(u8),

    // Options
    phase: Phase,
    local_csq: bool,
    verbosity: i32,
    force: bool,
    bcsq_tag: []const u8,
    ncsq2_max: u32,
    nfmt_bcsq: u32,
    brief_predictions: u32,

    // State
    current_rid: i32,
    prev_rid: i32,
    prev_pos: i32,

    // Warnings (emit once)
    warned_faidx_fetch_failed: bool,
    warned_ref_allele_mismatch: bool,

    /// Initialize the CSQ context with the given options.
    ///
    /// Corresponds to init_data() in csq.c. The GFF parsing, fasta index
    /// loading, and VCF header setup are deferred to when the htslib
    /// bindings are available.
    pub fn init(allocator: std.mem.Allocator, options: Options) !CsqContext {
        // ncsq2_max -> nfmt_bcsq: see ncsq2_to_nfmt in csq.c
        const nfmt = ncsq2ToNfmt(options.ncsq2_max);

        return CsqContext{
            .allocator = allocator,
            .pos2vbuf = std.AutoHashMap(u32, usize).init(allocator),
            .vcf_rbuf = try RingBuffer(*Vbuf).init(allocator, 64),
            .csq_buf = .empty,
            .output = .empty,
            .phase = options.phase,
            .local_csq = options.local_csq,
            .verbosity = options.verbosity,
            .force = options.force,
            .bcsq_tag = options.bcsq_tag,
            .ncsq2_max = options.ncsq2_max,
            .nfmt_bcsq = nfmt,
            .brief_predictions = options.brief_predictions,
            .current_rid = -1,
            .prev_rid = -1,
            .prev_pos = -1,
            .warned_faidx_fetch_failed = false,
            .warned_ref_allele_mismatch = false,
        };
    }

    /// Release all resources.  Corresponds to destroy_data() in csq.c.
    pub fn deinit(self: *CsqContext) void {
        // Free all vbufs still in the ring buffer
        for (0..self.vcf_rbuf.len) |k| {
            var vbuf = self.vcf_rbuf.kth(k);
            vbuf.deinit(self.allocator);
            self.allocator.destroy(vbuf);
        }
        self.vcf_rbuf.deinit();
        self.pos2vbuf.deinit();
        self.csq_buf.deinit(self.allocator);
        self.output.deinit(self.allocator);
    }

    // -----------------------------------------------------------------
    // vbufPush — add a record to the position-indexed buffer
    // -----------------------------------------------------------------

    /// Buffer a VCF record. Records at the same position share a Vbuf.
    /// Returns a pointer to the Vbuf containing the record.
    ///
    /// Port of vbuf_push() from csq.c (line 2675).
    pub fn vbufPush(self: *CsqContext, rec: *const VcfRecord) !*Vbuf {
        // Check if the last buffered vbuf has the same position
        const last_vbuf: ?*Vbuf = self.vcf_rbuf.last();
        const same_pos = if (last_vbuf) |vb| blk: {
            break :blk if (vb.pos()) |p| p == rec.pos else false;
        } else false;

        var vbuf: *Vbuf = undefined;
        if (same_pos) {
            vbuf = last_vbuf.?;
        } else {
            // Allocate a new Vbuf for this position
            vbuf = try self.allocator.create(Vbuf);
            vbuf.* = .{};
            _ = try self.vcf_rbuf.append(vbuf);
        }

        // Add the record as a new Vrec
        var vrec = Vrec{};
        vrec.rec = rec;
        try vbuf.vrecs.append(self.allocator, vrec);

        // Register in pos2vbuf for O(1) existence check by position
        try self.pos2vbuf.put(rec.pos, 0);

        return vbuf;
    }

    // -----------------------------------------------------------------
    // vbufFlush — flush records up to a position
    // -----------------------------------------------------------------

    /// Flush all buffered VCF records whose keep_until <= pos.
    /// Formats BCSQ strings and writes output.
    ///
    /// Port of vbuf_flush() from csq.c (line 2715).
    pub fn vbufFlush(self: *CsqContext, pos: u32) !void {
        while (self.vcf_rbuf.len > 0) {
            const vbuf = self.vcf_rbuf.front().?;

            // Cannot flush if transcript still active beyond this position
            if (!self.local_csq and vbuf.keep_until > pos) break;

            _ = self.vcf_rbuf.shift();

            // Remove from pos2vbuf
            if (vbuf.pos()) |vpos| {
                _ = self.pos2vbuf.remove(vpos);
            }

            // Format consequences for each record in the vbuf
            for (vbuf.vrecs.items) |*vrec| {
                if (vrec.vcsqs.items.len == 0) {
                    // No consequences — in VCF mode we'd write the record as-is.
                    // TODO: write unmodified record via htslib
                    continue;
                }

                // Format the BCSQ INFO string
                self.output.clearRetainingCapacity();
                try formatVcsqList(
                    vrec.vcsqs.items,
                    .{ .brief_predictions = self.brief_predictions },
                    self.output.writer(self.allocator),
                );

                // TODO: bcf_update_info_string(hdr, rec, bcsq_tag, output)
                // TODO: bcf_update_format_int32(hdr, rec, bcsq_tag, fmt_bm, ...)
                // TODO: bcf_write(out_fh, hdr, rec)
            }

            vbuf.deinit(self.allocator);
            self.allocator.destroy(vbuf);
        }

        // TODO: when active_transcripts heap is empty, clean up rm_transcripts
        // (destroy haplotype trees, free reference sequences, etc.)
        self.csq_buf.clearRetainingCapacity();
    }

    // -----------------------------------------------------------------
    // hapFlush — flush completed transcript haplotypes
    // -----------------------------------------------------------------

    /// Flush haplotypes for transcripts that end at or before `pos`.
    ///
    /// Port of hap_flush() from csq.c (line 2627).
    /// TODO: requires active_transcripts min-heap and hap_finalize.
    pub fn hapFlush(self: *CsqContext, pos: u32) !void {
        _ = self;
        _ = pos;
        // TODO: while active_transcripts heap has transcripts ending <= pos:
        //   1. Pop transcript from heap
        //   2. Call hap_finalize to walk the haplotype tree and emit consequences
        //   3. For VCF output with genotypes, call hap_stage_vcf per sample
        //   4. Mark transcript for deferred cleanup in rm_transcripts
    }

    // -----------------------------------------------------------------
    // process — main per-record dispatch
    // -----------------------------------------------------------------

    /// Process a single VCF record through the consequence pipeline.
    ///
    /// Port of process() from csq.c (line 3602).
    pub fn process(self: *CsqContext, rec: *const VcfRecord) !void {
        // Validate sort order
        if (self.prev_rid == rec.rid and self.prev_pos > @as(i32, @intCast(rec.pos))) {
            return error.UnsortedInput;
        }

        // Chromosome change
        if (self.prev_rid != rec.rid) {
            self.prev_rid = rec.rid;
            self.prev_pos = @intCast(rec.pos);
            // TODO: validate chromosome exists in fasta and GFF
        }

        // Check if this record has callable alt alleles
        var call_csq = true;
        if (rec.n_allele < 2) {
            call_csq = false;
        } else if (rec.n_allele == 2 and rec.alleles.len >= 2) {
            const alt = rec.alleles[1];
            if (alt.len > 0 and (alt[0] == '*')) {
                call_csq = false; // gVCF, not a real alt
            }
        }

        if (!call_csq) {
            // Still buffer the record for pass-through output
            _ = try self.vbufPush(rec);
            if (rec.pos > 0) {
                try self.hapFlush(rec.pos - 1);
                try self.vbufFlush(rec.pos - 1);
            }
            return;
        }

        // Flush on chromosome change
        if (self.current_rid != rec.rid) {
            try self.hapFlush(POS_MAX);
            try self.vbufFlush(POS_MAX);
        }
        self.current_rid = rec.rid;

        const vbuf = try self.vbufPush(rec);

        // Check for symbolic ALTs
        if (rec.alleles.len >= 2 and rec.alleles[1].len > 0 and rec.alleles[1][0] == '<') {
            // TODO: test_symbolic_alt
        } else {
            // Annotation lookup cascade: CDS -> UTR -> splice -> transcript
            var hit: bool = false;
            if (self.local_csq) {
                hit = try self.testCdsLocal(rec);
            } else {
                hit = try self.testCds(rec, vbuf);
            }
            if (!hit) {
                const utr_hit = try self.testUtr(rec);
                hit = hit or utr_hit;
            }
            {
                const splice_hit = try self.testSplice(rec);
                hit = hit or splice_hit;
            }
            if (!hit) {
                _ = try self.testTscript(rec);
            }
        }

        // Flush completed haplotypes and records
        if (rec.pos > 0) {
            try self.hapFlush(rec.pos - 1);
            try self.vbufFlush(rec.pos - 1);
        }

        self.prev_pos = @intCast(rec.pos);
    }

    /// Flush all remaining records. Call at end of input.
    ///
    /// Corresponds to process(args, NULL) in csq.c.
    pub fn flush(self: *CsqContext) !void {
        try self.hapFlush(POS_MAX);
        try self.vbufFlush(POS_MAX);
    }

    // -----------------------------------------------------------------
    // testCds — look up CDS regions overlapping the variant
    // -----------------------------------------------------------------

    /// Check if the variant overlaps coding sequences and build haplotype nodes.
    ///
    /// Port of test_cds() from csq.c (line 3075).
    /// TODO: requires GFF RegionIndex (idx_cds) and haplotype tree operations.
    fn testCds(self: *CsqContext, rec: *const VcfRecord, vbuf: *Vbuf) !bool {
        _ = self;
        _ = rec;
        _ = vbuf;
        // TODO: Implementation steps:
        // 1. regidx_overlap(idx_cds, chr, rec.pos, rec.pos + rec.rlen)
        // 2. For each overlapping CDS:
        //    a. Get transcript, check if coding
        //    b. Set vbuf.keep_until = max(keep_until, tr.end)
        //    c. Initialize transcript aux if needed (fetch ref, build haplotype root)
        //    d. Sanity-check ref allele
        //    e. For phase==drop_gt: single haplotype path
        //    f. Otherwise: iterate samples/haplotypes, extend tree
        return false;
    }

    /// Local (non-haplotype-aware) CDS consequence calling.
    ///
    /// Port of test_cds_local() from csq.c (line 2872).
    fn testCdsLocal(self: *CsqContext, rec: *const VcfRecord) !bool {
        _ = self;
        _ = rec;
        // TODO: Implementation steps:
        // 1. regidx_overlap(idx_cds, chr, rec.pos, rec.pos + rec.rlen)
        // 2. For each overlapping CDS:
        //    a. For each alt allele, call hap_init to get single-variant node
        //    b. If HAP_SSS: stage start/stop/splice consequence
        //    c. Otherwise: translate ref and alt, compare amino acids
        //    d. Build csq with protein/DNA change string, call csqStage
        return false;
    }

    // -----------------------------------------------------------------
    // testUtr — look up UTR regions
    // -----------------------------------------------------------------

    /// Check if the variant overlaps UTR regions.
    ///
    /// Port of test_utr() from csq.c (line 3360).
    fn testUtr(self: *CsqContext, rec: *const VcfRecord) !bool {
        _ = self;
        _ = rec;
        // TODO: Implementation steps:
        // 1. regidx_overlap(idx_utr, chr, rec.pos, rec.pos + rec.rlen)
        // 2. For each overlapping UTR:
        //    a. For each alt allele, check splice consequences
        //    b. If inside/overlap: create UTR5 or UTR3 consequence
        //    c. Call csqStage
        return false;
    }

    // -----------------------------------------------------------------
    // testSplice — check for splice site variants
    // -----------------------------------------------------------------

    /// Check if the variant affects splice sites.
    ///
    /// Port of test_splice() from csq.c (line 3400).
    fn testSplice(self: *CsqContext, rec: *const VcfRecord) !bool {
        _ = self;
        _ = rec;
        // TODO: Implementation steps:
        // 1. regidx_overlap(idx_exon, chr, rec.pos, rec.pos + rec.rlen)
        // 2. For each overlapping exon:
        //    a. Skip non-coding transcripts (ncds == 0)
        //    b. Check region boundaries for acceptor/donor
        //    c. For each alt allele, call splice_csq
        return false;
    }

    // -----------------------------------------------------------------
    // testTscript — check for intronic / non-coding variants
    // -----------------------------------------------------------------

    /// Check if the variant falls within a transcript (intron or non-coding).
    ///
    /// Port of test_tscript() from csq.c (line 3434).
    fn testTscript(self: *CsqContext, rec: *const VcfRecord) !bool {
        _ = self;
        _ = rec;
        // TODO: Implementation steps:
        // 1. regidx_overlap(idx_tscript, chr, rec.pos, rec.pos + rec.rlen)
        // 2. For each overlapping transcript:
        //    a. For each alt allele, call splice_csq to check boundaries
        //    b. If inside/overlap: create INTRON (coding) or NON_CODING consequence
        //    c. Call csqStage
        return false;
    }

    // -----------------------------------------------------------------
    // csqPush — add consequence to vbuf with deduplication
    // -----------------------------------------------------------------

    /// Add a consequence to the appropriate vrec in the vbuf.
    /// Returns true if the consequence was already present (duplicate).
    ///
    /// Port of csq_push() from csq.c (line 1996).
    /// Uses field-level matching for dedup (audit fix #30): checks
    /// transcript, biotype, gene, allele, and vstr to decide merging.
    pub fn csqPush(self: *CsqContext, csq: *Csq, rec: *const VcfRecord) !bool {
        // Look up the vbuf for this position
        _ = self.pos2vbuf.get(csq.pos) orelse {
            return error.VbufNotFound;
        };

        // Find the vbuf at this position by scanning the ring buffer
        var vbuf: ?*Vbuf = null;
        for (0..self.vcf_rbuf.len) |k| {
            const candidate = self.vcf_rbuf.kth(k);
            if (candidate.pos()) |p| {
                if (p == csq.pos) {
                    vbuf = candidate;
                    break;
                }
            }
        }
        const vb = vbuf orelse return error.VbufNotFound;

        // Find the vrec matching this record
        var vrec_idx: ?usize = null;
        for (vb.vrecs.items, 0..) |*vrec, idx| {
            if (vrec.rec == rec) {
                vrec_idx = idx;
                break;
            }
        }
        const vi = vrec_idx orelse return error.VrecNotFound;
        const vrec = &vb.vrecs.items[vi];

        // Apply type masking rules
        var t = csq.vcsq.csq_type;
        if (t & CSQ_INFRAME_INSERTION != 0 and t & CSQ_ELONGATION != 0) t &= ~CSQ_INFRAME_INSERTION;
        if (t & CSQ_INFRAME_DELETION != 0 and t & CSQ_TRUNCATION != 0) t &= ~CSQ_INFRAME_DELETION;
        csq.vcsq.csq_type = t;

        // Splice region / donor+acceptor masking
        if (t & CSQ_SPLICE_REGION != 0 and t & (CSQ_SPLICE_DONOR | CSQ_SPLICE_ACCEPTOR) != 0) {
            csq.vcsq.csq_type &= ~CSQ_SPLICE_REGION;
            t = csq.vcsq.csq_type;
        }

        // Deduplication: scan existing consequences for a match
        for (vrec.vcsqs.items, 0..) |*existing, idx| {
            if (isDuplicate(existing, &csq.vcsq)) {
                // Merge type bits into existing consequence
                existing.csq_type |= t;
                csq.vrec_idx = vi;
                csq.csq_idx = idx;
                return true; // duplicate
            }
        }

        // New consequence — append
        csq.vrec_idx = vi;
        csq.csq_idx = vrec.vcsqs.items.len;
        try vrec.vcsqs.append(self.allocator, csq.vcsq);
        return false;
    }

    /// Check if an existing consequence is a duplicate of a new one.
    fn isDuplicate(existing: *const Vcsq, new: *const Vcsq) bool {
        // Different transcript IDs with printable transcripts -> not duplicate
        if (existing.trid != new.trid and
            (existing.csq_type | new.csq_type) & format.CSQ_PRN_TSCRIPT != 0)
            return false;

        if (existing.biotype != new.biotype) return false;

        // For compound consequences, also check gene, vcf_ial, and vstr
        if (new.csq_type & CSQ_COMPOUND != 0) {
            if (!strEql(existing.gene, new.gene)) return false;
            if (existing.vcf_ial != new.vcf_ial) return false;

            // Both have vstr: must match
            if (existing.vstr != null and new.vstr != null) {
                if (!std.mem.eql(u8, existing.vstr.?, new.vstr.?)) return false;
            } else if (existing.vstr != null or new.vstr != null) {
                // One has vstr, the other doesn't: special START_STOP merging
                if (existing.csq_type & CSQ_START_STOP != 0 and new.csq_type & CSQ_START_STOP != 0) {
                    return true; // will merge
                }
                return false;
            }
            return true;
        }

        // Non-compound: simpler check
        if (existing.csq_type & CSQ_COMPOUND == 0) {
            return true; // same biotype + trid, merge
        }

        // Existing is compound, new is not: check if the merge would be redundant
        return existing.csq_type == (existing.csq_type | new.csq_type);
    }

    fn strEql(a: ?[]const u8, b: ?[]const u8) bool {
        if (a == null and b == null) return true;
        if (a == null or b == null) return false;
        return std.mem.eql(u8, a.?, b.?);
    }

    // -----------------------------------------------------------------
    // csqStage — stage a consequence and handle genotype assignment
    // -----------------------------------------------------------------

    /// Stage a consequence: push it and assign to samples via genotypes.
    ///
    /// Port of csq_stage() from csq.c (line 3290).
    pub fn csqStage(self: *CsqContext, csq: *Csq, rec: *const VcfRecord) !void {
        const is_dup = try self.csqPush(csq, rec);
        if (is_dup and self.phase == .drop_gt) return;

        // TODO: genotype-aware sample assignment
        // 1. Get genotypes from the record
        // 2. For each sample, check if the allele matches csq.vcsq.vcf_ial
        // 3. For tab output: call csq_print_text
        // 4. For VCF output: set bits in vrec.fmt_bm
    }

    // -----------------------------------------------------------------
    // Utility
    // -----------------------------------------------------------------

    /// Convert ncsq2_max to the number of FORMAT integers needed.
    /// Port of ncsq2_to_nfmt() from csq.c.
    fn ncsq2ToNfmt(ncsq2_max: u32) u32 {
        return (ncsq2_max + 30) / 31; // 31 usable bits per int32 (1 reserved for BCF missing)
    }

    pub const ProcessError = error{
        UnsortedInput,
        VbufNotFound,
        VrecNotFound,
        OutOfMemory,
    };
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "vbufPush: two records at same position share vbuf" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "test.gff",
        .phase = .drop_gt,
    });
    defer ctx.deinit();

    const alleles = [_][]const u8{ "A", "T" };
    var rec1 = VcfRecord{
        .pos = 100,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
    };
    var rec2 = VcfRecord{
        .pos = 100,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
    };

    const vbuf1 = try ctx.vbufPush(&rec1);
    const vbuf2 = try ctx.vbufPush(&rec2);

    // Both should return the same vbuf
    try std.testing.expectEqual(vbuf1, vbuf2);
    // The vbuf should contain 2 records
    try std.testing.expectEqual(@as(usize, 2), vbuf1.vrecs.items.len);
    // Ring buffer should have exactly 1 entry
    try std.testing.expectEqual(@as(usize, 1), ctx.vcf_rbuf.len);
}

test "vbufPush: records at different positions get separate vbufs" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "test.gff",
        .phase = .drop_gt,
    });
    defer ctx.deinit();

    const alleles = [_][]const u8{ "A", "T" };
    var rec1 = VcfRecord{
        .pos = 100,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
    };
    var rec2 = VcfRecord{
        .pos = 200,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
    };

    const vbuf1 = try ctx.vbufPush(&rec1);
    const vbuf2 = try ctx.vbufPush(&rec2);

    // Should be different vbufs
    try std.testing.expect(vbuf1 != vbuf2);
    // Each should have 1 record
    try std.testing.expectEqual(@as(usize, 1), vbuf1.vrecs.items.len);
    try std.testing.expectEqual(@as(usize, 1), vbuf2.vrecs.items.len);
    // Ring buffer should have 2 entries
    try std.testing.expectEqual(@as(usize, 2), ctx.vcf_rbuf.len);
}

test "csqPush: dedup pushes same consequence only once" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "test.gff",
        .phase = .drop_gt,
    });
    defer ctx.deinit();

    const alleles = [_][]const u8{ "A", "T" };
    var rec = VcfRecord{
        .pos = 42,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
    };

    _ = try ctx.vbufPush(&rec);

    // Push the same consequence twice
    var csq1 = Csq{
        .pos = 42,
        .vcsq = .{
            .csq_type = CSQ_MISSENSE_VARIANT,
            .trid = 1,
            .biotype = 100,
            .vcf_ial = 1,
            .gene = "BRCA1",
        },
    };
    var csq2 = Csq{
        .pos = 42,
        .vcsq = .{
            .csq_type = CSQ_MISSENSE_VARIANT,
            .trid = 1,
            .biotype = 100,
            .vcf_ial = 1,
            .gene = "BRCA1",
        },
    };

    const dup1 = try ctx.csqPush(&csq1, &rec);
    const dup2 = try ctx.csqPush(&csq2, &rec);

    try std.testing.expect(!dup1); // first push: not a duplicate
    try std.testing.expect(dup2); // second push: is a duplicate

    // Only one consequence stored in the vrec
    const vbuf = ctx.vcf_rbuf.front().?;
    try std.testing.expectEqual(@as(usize, 1), vbuf.vrecs.items[0].vcsqs.items.len);
}

test "vbufFlush: flushes all records" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "test.gff",
        .phase = .drop_gt,
        .local_csq = true, // so keep_until is not checked against active transcripts
    });
    defer ctx.deinit();

    const alleles = [_][]const u8{ "A", "T" };
    var rec1 = VcfRecord{ .pos = 10, .rid = 0, .n_allele = 2, .alleles = &alleles, .rlen = 1 };
    var rec2 = VcfRecord{ .pos = 20, .rid = 0, .n_allele = 2, .alleles = &alleles, .rlen = 1 };
    var rec3 = VcfRecord{ .pos = 30, .rid = 0, .n_allele = 2, .alleles = &alleles, .rlen = 1 };

    _ = try ctx.vbufPush(&rec1);
    _ = try ctx.vbufPush(&rec2);
    _ = try ctx.vbufPush(&rec3);

    try std.testing.expectEqual(@as(usize, 3), ctx.vcf_rbuf.len);

    // Flush everything
    try ctx.vbufFlush(POS_MAX);

    try std.testing.expectEqual(@as(usize, 0), ctx.vcf_rbuf.len);
}

test "ncsq2ToNfmt calculation" {
    try std.testing.expectEqual(@as(u32, 1), CsqContext.ncsq2ToNfmt(1));
    try std.testing.expectEqual(@as(u32, 1), CsqContext.ncsq2ToNfmt(30));
    try std.testing.expectEqual(@as(u32, 2), CsqContext.ncsq2ToNfmt(32));
}
