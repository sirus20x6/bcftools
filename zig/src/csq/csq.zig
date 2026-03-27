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
const splice_mod = @import("splice.zig");
const gff_mod = @import("../gff/gff.zig");
const gff_types = @import("../gff/types.zig");
const region = @import("../core/region.zig");
const types = @import("types.zig");

const haplotype_mod = @import("haplotype.zig");
const translate = @import("translate.zig");

const Splice = splice_mod.Splice;
const SpliceResult = splice_mod.SpliceResult;
const GffParser = gff_mod.GffParser;
const Transcript = gff_types.Transcript;
const CdsEntry = gff_types.CdsEntry;
const Utr = gff_types.Utr;
const Exon = gff_types.Exon;
const Tscript = types.Tscript;
const HapNode = types.HapNode;
const HapNodeType = types.HapNodeType;
const HapInitResult = haplotype_mod.HapInitResult;
const HapContext = haplotype_mod.HapContext;

// ---------------------------------------------------------------------------
// ActiveTranscriptQueue -- min-heap of transcripts sorted by end position
// ---------------------------------------------------------------------------

/// Priority queue of active transcripts, ordered by ascending `.end` position.
/// Used to determine which transcripts have been fully traversed at a given
/// genomic position so their haplotype trees can be finalized.
fn trLessThan(_: void, a: *Transcript, b: *Transcript) std.math.Order {
    return std.math.order(a.end, b.end);
}
const ActiveTranscriptQueue = std.PriorityQueue(*Transcript, void, trLessThan);

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
// Module-level GFF pointer for format callbacks
// ---------------------------------------------------------------------------

/// Module-level GFF parser reference used by the trid/biotype format callbacks.
/// Set by CsqContext.init when the GFF is available.
var g_gff_ptr: ?*GffParser = null;

fn tridToStringCallback(trid: u32) []const u8 {
    if (g_gff_ptr) |gff| {
        return gff.id2string(trid);
    }
    return "";
}

fn biotypeToStringCallback(biotype_id: u32) []const u8 {
    const biotype: gff_types.Biotype = @enumFromInt(biotype_id);
    return biotype.toGffString();
}

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

/// Parsed genotype for a single sample.
pub const Genotype = struct {
    /// Allele indices: -1 for missing, 0 for ref, 1+ for alt.
    alleles: [2]i32,
    /// Whether the genotype separator was '|' (phased).
    phased: bool,
    /// Number of alleles present (1 for haploid, 2 for diploid).
    ploidy: u8,

    /// Parse a single GT string like "0/1", "1|0", "./.", "0", "1", ".".
    pub fn parse(gt_str: []const u8) Genotype {
        if (gt_str.len == 0) return .{ .alleles = .{ -1, -1 }, .phased = false, .ploidy = 0 };

        var result = Genotype{ .alleles = .{ -1, -1 }, .phased = false, .ploidy = 1 };
        var sep_pos: ?usize = null;

        for (gt_str, 0..) |c, idx| {
            if (c == '/' or c == '|') {
                if (c == '|') result.phased = true;
                sep_pos = idx;
                result.ploidy = 2;
                break;
            }
        }

        // Parse first allele
        const first_end = sep_pos orelse gt_str.len;
        const first = gt_str[0..first_end];
        if (first.len == 1 and first[0] == '.') {
            result.alleles[0] = -1;
        } else {
            result.alleles[0] = std.fmt.parseInt(i32, first, 10) catch -1;
        }

        // Parse second allele if present
        if (sep_pos) |sp| {
            const second = gt_str[sp + 1 ..];
            if (second.len == 1 and second[0] == '.') {
                result.alleles[1] = -1;
            } else {
                result.alleles[1] = std.fmt.parseInt(i32, second, 10) catch -1;
            }
        }

        return result;
    }
};

/// Minimal VCF record representation for the CSQ pipeline.
/// This will be replaced by a proper htslib binding wrapper.
pub const VcfRecord = struct {
    pos: u32,
    rid: i32,
    n_allele: u32,
    alleles: []const []const u8,
    rlen: u32,
    /// Chromosome/sequence name for region index lookups.
    chr: []const u8 = "unknown",
    /// Raw VCF line for genotype parsing (optional, used for text VCF input).
    raw_line: ?[]const u8 = null,
    /// Cached parsed genotypes (one per sample, lazily populated).
    gt_cache: ?[]Genotype = null,

    pub fn seqname(self: *const VcfRecord) []const u8 {
        return self.chr;
    }

    /// Return the chromosome name as a null-terminated pointer.
    /// The underlying `chr` slice is expected to originate from a
    /// null-terminated source (e.g. VCF header seqname).  If it does
    /// not end with a sentinel zero we fall back to a comptime default.
    /// Scratch buffer for null-terminated chromosome name.
    var chr_z_buf: [256]u8 = undefined;

    pub fn chrZ(self: *const VcfRecord) [*:0]const u8 {
        if (self.chr.len == 0) return "unknown";
        // Try sentinel-terminated reinterpret first
        if (self.chr.ptr[self.chr.len] == 0) {
            return self.chr.ptr[0..self.chr.len :0];
        }
        // Fall back to copy into scratch buffer
        const len = @min(self.chr.len, chr_z_buf.len - 1);
        @memcpy(chr_z_buf[0..len], self.chr[0..len]);
        chr_z_buf[len] = 0;
        return chr_z_buf[0..len :0];
    }

    /// Parse genotypes from the raw VCF line for all samples.
    /// Returns a slice of Genotype, one per sample, or null if no genotypes.
    /// The GT field must be the first subfield in FORMAT.
    pub fn parseGenotypes(self: *const VcfRecord, allocator: std.mem.Allocator) !?[]Genotype {
        // If already cached, return cached
        if (self.gt_cache) |cached| return cached;

        const line = self.raw_line orelse return null;

        // Find FORMAT column (index 8) and sample columns (index 9+)
        var col: usize = 0;
        var col_start: usize = 0;
        var format_start: usize = 0;
        var format_end: usize = 0;
        var samples_start: usize = 0;

        for (line, 0..) |c, idx| {
            if (c == '\t') {
                if (col == 8) {
                    format_start = col_start;
                    format_end = idx;
                    samples_start = idx + 1;
                }
                col += 1;
                col_start = idx + 1;
            }
        }
        // Handle if FORMAT is the last column found
        if (col == 8) {
            format_start = col_start;
            format_end = line.len;
            return null; // No sample columns
        }
        if (col < 9) return null; // Not enough columns for FORMAT + samples

        // Check that GT is the first FORMAT subfield
        const format_field = line[format_start..format_end];
        const gt_ok = if (format_field.len >= 2)
            (std.mem.startsWith(u8, format_field, "GT\t") or
                std.mem.startsWith(u8, format_field, "GT:") or
                std.mem.eql(u8, format_field, "GT"))
        else
            false;
        if (!gt_ok) return null;

        // Count samples
        var n_samples: usize = 1;
        for (line[samples_start..]) |c| {
            if (c == '\t') n_samples += 1;
        }

        var genotypes = try allocator.alloc(Genotype, n_samples);

        // Parse each sample's GT field (first subfield before ':')
        var smpl_idx: usize = 0;
        var pos: usize = samples_start;
        while (smpl_idx < n_samples) : (smpl_idx += 1) {
            // Find the end of the GT subfield (first ':' or '\t' or end of line)
            var gt_end: usize = pos;
            while (gt_end < line.len and line[gt_end] != ':' and line[gt_end] != '\t' and line[gt_end] != '\n') {
                gt_end += 1;
            }
            genotypes[smpl_idx] = Genotype.parse(line[pos..gt_end]);

            // Advance to the next sample (skip to next '\t')
            var next_pos = gt_end;
            while (next_pos < line.len and line[next_pos] != '\t') {
                next_pos += 1;
            }
            pos = if (next_pos < line.len) next_pos + 1 else next_pos;
        }

        return genotypes;
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
        // Free any owned vstr strings in consequences
        for (self.vcsqs.items) |vcsq| {
            if (vcsq.vstr) |vs| allocator.free(vs);
        }
        self.vcsqs.deinit(allocator);
        if (self.fmt_bm) |bm| {
            allocator.free(bm);
            self.fmt_bm = null;
        }
        // Free the owned VcfRecord (heap-allocated by vbufPush)
        if (self.rec) |rec_ptr| {
            // Free deep-copied slices
            for (rec_ptr.alleles) |a| {
                allocator.free(a);
            }
            if (rec_ptr.alleles.len > 0) {
                allocator.free(rec_ptr.alleles);
            }
            if (rec_ptr.chr.len > 0) {
                allocator.free(rec_ptr.chr);
            }
            if (rec_ptr.raw_line) |rl| {
                allocator.free(rl);
            }
            // rec is *const VcfRecord, need to cast to free
            const mutable: *VcfRecord = @constCast(rec_ptr);
            allocator.destroy(mutable);
            self.rec = null;
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

/// Opaque type-erased pointer to an HtsFaidx instance (from vcf/htslib.zig).
/// Stored as `*anyopaque` so that this module does not require htslib C headers
/// at compile time.  The caller (main.zig) is responsible for opening the faidx
/// and passing it in via Options.fai_ptr.
///
/// To perform a reference fetch through this pointer, use `faidxFetchSeq` below
/// which casts back to the concrete HtsFaidx type via the htslib module.
pub const FaidxPtr = *anyopaque;

/// Function pointer type for fetching a reference sequence region.
/// This abstracts over the concrete htslib faidx so the pipeline can be
/// tested without linking htslib.
///
/// Parameters:
///   ctx      — opaque context (e.g. FaidxPtr)
///   chr      — chromosome name (null-terminated)
///   beg, end — 0-based, inclusive coordinates
///
/// Returns the fetched sequence as an owned slice, or null on failure.
pub const FetchSeqFn = *const fn (ctx: *anyopaque, allocator: std.mem.Allocator, chr: [*:0]const u8, beg: i64, end: i64) ?[]u8;

pub const Options = struct {
    gff_fname: []const u8 = "",
    gff_ptr: ?*GffParser = null, // pre-built GFF parser (tests)
    fasta_fname: []const u8 = "",
    /// Opaque pointer to an opened HtsFaidx.  May be null when running without
    /// htslib (text-only testing path).
    fai_ptr: ?FaidxPtr = null,
    /// Optional function for fetching reference sequences.  When non-null this
    /// is called instead of going through fai_ptr directly, allowing test stubs.
    fetch_seq_fn: ?FetchSeqFn = null,
    phase: Phase = .require,
    local_csq: bool = false,
    verbosity: i32 = 1,
    force: bool = false,
    bcsq_tag: []const u8 = "BCSQ",
    ncsq2_max: u32 = 15 * 2,
    gencode_id: i32 = 0,
    brief_predictions: u32 = 0,
    /// Number of samples in the VCF.
    n_samples: u32 = 0,
    /// Indices of selected samples.  If null, all samples are used.
    sample_indices: ?[]const u32 = null,
};

// ---------------------------------------------------------------------------
// CsqContext — main state container, replaces args_t
// ---------------------------------------------------------------------------

pub const CsqContext = struct {
    allocator: std.mem.Allocator,

    // GFF annotation — provides region indexes for CDS, UTR, exon, transcript lookups
    gff: ?*GffParser,
    owns_gff: bool,

    // Haplotype processing
    hap_ctx: HapContext,

    // FASTA reference access
    fasta_fname: []const u8,
    fai_ptr: ?FaidxPtr,
    fetch_seq_fn: ?FetchSeqFn,

    // VCF record buffering
    pos2vbuf: std.AutoHashMap(u32, usize), // pos -> ring buffer index (for existence check)
    vcf_rbuf: RingBuffer(*Vbuf),

    // Transcript management -- min-heap of active transcripts sorted by end position
    active_transcripts: ActiveTranscriptQueue,
    // Transcripts pending cleanup after vbuf flush (cannot delete immediately because
    // by-position VCF output needs them when flushed by vbufFlush)
    rm_transcripts: std.ArrayList(*Transcript),

    // CSQ buffer for non-CDS consequences
    csq_buf: std.ArrayList(Csq),

    // Output
    output: std.ArrayList(u8),

    // Flushed records — populated by vbufFlush, consumed by the caller
    flushed_records: std.ArrayList(FlushedRecord) = .empty,

    // Sample management
    n_samples: u32,
    /// Indices of selected samples (maps from 0..smpl_n to VCF sample column).
    /// If null, all samples are used (identity mapping).
    sample_indices: ?[]const u32,

    // Cached genotype array (reused across csqStage calls for the same record)
    gt_cache_rec: ?*const VcfRecord,
    gt_cache: ?[]Genotype,

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

    // Warnings (emit once, count for verbosity > 1)
    warned_faidx_fetch_failed: u32,
    warned_ref_allele_mismatch: u32,

    /// Initialize the CSQ context with the given options.
    ///
    /// Corresponds to init_data() in csq.c. The GFF parsing, fasta index
    /// loading, and VCF header setup are deferred to when the htslib
    /// bindings are available.
    pub fn init(allocator: std.mem.Allocator, options: Options) !CsqContext {
        // ncsq2_max -> nfmt_bcsq: see ncsq2_to_nfmt in csq.c
        const nfmt = ncsq2ToNfmt(options.ncsq2_max);

        const gencode = translate.findGeneticCode(options.gencode_id) orelse
            translate.findGeneticCode(0).?;

        // Parse GFF annotation file (skip if empty fname or pre-built gff provided)
        var gff_ptr: ?*GffParser = options.gff_ptr;
        var owns_gff = false;
        if (gff_ptr == null and options.gff_fname.len > 0) {
            const gff_obj = try allocator.create(GffParser);
            gff_obj.* = GffParser.init(allocator);
            gff_obj.parse(options.gff_fname) catch |e| {
                allocator.destroy(gff_obj);
                std.debug.print("Error: failed to parse GFF '{s}': {}\n", .{ options.gff_fname, e });
                return e;
            };
            gff_ptr = gff_obj;
            owns_gff = true;
        }

        // Set module-level GFF pointer for format callbacks
        g_gff_ptr = gff_ptr;

        return CsqContext{
            .allocator = allocator,
            .gff = gff_ptr,
            .owns_gff = owns_gff,
            .hap_ctx = HapContext.init(allocator, gencode),
            .fasta_fname = options.fasta_fname,
            .fai_ptr = options.fai_ptr,
            .fetch_seq_fn = options.fetch_seq_fn,
            .pos2vbuf = std.AutoHashMap(u32, usize).init(allocator),
            .vcf_rbuf = try RingBuffer(*Vbuf).init(allocator, 64),
            .active_transcripts = ActiveTranscriptQueue.init(allocator, {}),
            .rm_transcripts = .empty,
            .csq_buf = .empty,
            .output = .empty,
            .flushed_records = .empty,
            .n_samples = options.n_samples,
            .sample_indices = options.sample_indices,
            .gt_cache_rec = null,
            .gt_cache = null,
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
            .warned_faidx_fetch_failed = 0,
            .warned_ref_allele_mismatch = 0,
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
        // Free any duped bcsq_value strings and fmt_bm in flushed records
        for (self.flushed_records.items) |fr| {
            if (fr.bcsq_value) |bv| self.allocator.free(bv);
            if (fr.fmt_bm) |bm| self.allocator.free(bm);
        }
        self.flushed_records.deinit(self.allocator);
        // Free cached genotypes
        if (self.gt_cache) |gc| self.allocator.free(gc);
        self.gt_cache = null;
        self.gt_cache_rec = null;
        // Clean up haplotype context
        self.hap_ctx.deinit();
        // Clean up any remaining transcripts in the removal list
        for (self.rm_transcripts.items) |tr| {
            self.destroyTranscriptAux(tr);
        }
        self.rm_transcripts.deinit(self.allocator);
        self.active_transcripts.deinit();
        // Free the GFF parser (only if we created it)
        if (self.owns_gff) {
            if (self.gff) |g| {
                g.deinit();
                self.allocator.destroy(g);
            }
        }
    }


    // -----------------------------------------------------------------
    // FASTA reference — fetch, init, splice, sanity check
    // -----------------------------------------------------------------

    /// Fetch a reference sequence region via the configured faidx.
    ///
    /// Returns the sequence as an owned slice, or null if no faidx is
    /// available or the fetch fails.
    fn fetchRefSeq(self: *CsqContext, chr: [*:0]const u8, beg: i64, end: i64) ?[]u8 {
        // Use the function-pointer path (allows test stubs and htslib wrappers)
        if (self.fetch_seq_fn) |fetch_fn| {
            if (self.fai_ptr) |fai| {
                return fetch_fn(fai, self.allocator, chr, beg, end);
            }
        }
        // No faidx available — running without htslib
        return null;
    }

    /// Uppercase a byte slice in-place (ASCII only).
    fn uppercaseInPlace(seq: []u8) void {
        for (seq) |*c| {
            if (c.* >= 'a' and c.* <= 'z') c.* -= 32;
        }
    }

    /// Initialize the reference sequence for a transcript.
    ///
    /// Port of tscript_init_ref() from csq.c (line 2797).
    /// Fetches the genomic region [tr.beg - N_REF_PAD, tr.end + N_REF_PAD]
    /// from the FASTA reference and stores it in tscript.ref_seq.  If the
    /// transcript is close to the start of the chromosome, the left padding
    /// is filled with 'N' characters.
    ///
    /// Returns error.FaidxFetchFailed if the sequence cannot be fetched and
    /// --force is not set; returns error.FaidxSkipped if --force is set and
    /// the fetch fails (caller should skip this transcript).
    pub fn tscriptInitRef(self: *CsqContext, tr: *Transcript, chr: [*:0]const u8) !void {
        const tscript = try self.getOrCreateTscript(tr);

        const pad: u32 = N_REF_PAD;
        const pad_beg: u32 = if (tr.beg >= pad) pad else tr.beg;
        const beg: i64 = @as(i64, @intCast(tr.beg)) - @as(i64, @intCast(pad_beg));
        const end: i64 = @as(i64, @intCast(tr.end)) + @as(i64, @intCast(pad));

        const raw_seq = self.fetchRefSeq(chr, beg, end) orelse {
            // Fetch failed
            if (!self.force) {
                std.log.err("unable to fetch the region of the fasta reference {s}:{d}-{d}", .{
                    std.mem.span(chr), tr.beg + 1, tr.end + 1,
                });
                return error.FaidxFetchFailed;
            }
            if (self.verbosity > 0 and (self.warned_faidx_fetch_failed == 0 or self.verbosity > 1)) {
                std.log.warn("unable to fetch the region of the fasta reference {s}:{d}-{d}", .{
                    std.mem.span(chr), tr.beg + 1, tr.end + 1,
                });
                if (self.verbosity < 2) {
                    std.log.warn("This message is printed only once, the verbosity can be increased with `--verbosity 2`", .{});
                }
            }
            self.warned_faidx_fetch_failed += 1;
            return error.FaidxSkipped;
        };
        defer self.allocator.free(raw_seq);

        const raw_len: u32 = @intCast(raw_seq.len);
        const tr_len: u32 = tr.end - tr.beg + 1;

        // Determine actual padding achieved on each side
        const pad_end: u32 = if (raw_len > tr_len + pad_beg) raw_len - tr_len - pad_beg else 0;

        // If we got full padding on both sides, use the sequence directly
        if (pad_beg == pad and pad_end == pad) {
            tscript.ref_seq = try self.allocator.alloc(u8, raw_seq.len);
            @memcpy(tscript.ref_seq.?, raw_seq);
        } else {
            // Need to pad with N characters to reach N_REF_PAD on each side
            const total_len: usize = tr_len + 2 * pad;
            const ref = try self.allocator.alloc(u8, total_len);

            // Left N-padding
            const left_pad = pad - pad_beg;
            @memset(ref[0..left_pad], 'N');
            var pos: usize = left_pad;

            // Copy the fetched sequence
            @memcpy(ref[pos .. pos + raw_seq.len], raw_seq);
            pos += raw_seq.len;

            // Right N-padding
            const right_pad = pad - pad_end;
            if (pos + right_pad <= ref.len) {
                @memset(ref[pos .. pos + right_pad], 'N');
            }

            tscript.ref_seq = ref;
        }

        // Uppercase the entire reference
        if (tscript.ref_seq) |ref| {
            uppercaseInPlace(ref);
        }
    }

    /// Build the spliced reference sequence from the genomic reference and CDS entries.
    ///
    /// Port of tscript_splice_ref() from csq.c (line 1971).
    /// Concatenates N_REF_PAD bases of upstream context, all CDS segments,
    /// and N_REF_PAD bases of downstream context into tscript.sref.
    pub fn tscriptSpliceRef(self: *CsqContext, tr: *Transcript) !void {
        const tscript = self.getTscript(tr) orelse return error.TscriptNotInitialized;
        const ref = tscript.ref_seq orelse return error.RefNotLoaded;

        if (tr.cds.items.len == 0) return error.NoCdsEntries;

        // Total spliced length = sum of CDS lengths + 2 * N_REF_PAD
        var cds_total_len: u32 = 0;
        for (tr.cds.items) |cds| {
            cds_total_len += cds.len;
        }

        const pad: u32 = N_REF_PAD;
        const total_len: usize = cds_total_len + 2 * pad;
        const sref = try self.allocator.alloc(u8, total_len);

        var pos: usize = 0;

        // Copy N_REF_PAD bases upstream of the first CDS.
        // In the genomic ref, the first CDS starts at offset (cds[0].beg - tr.beg + N_REF_PAD).
        // We want N_REF_PAD bases before that: offset (cds[0].beg - tr.beg).
        const first_cds = tr.cds.items[0];
        const first_offset = first_cds.beg - tr.beg;
        if (first_offset + pad <= ref.len) {
            @memcpy(sref[0..pad], ref[first_offset .. first_offset + pad]);
        }
        pos = pad;

        // Copy each CDS segment from the genomic reference
        for (tr.cds.items) |cds| {
            const cds_offset: usize = pad + cds.beg - tr.beg;
            if (cds_offset + cds.len <= ref.len) {
                @memcpy(sref[pos .. pos + cds.len], ref[cds_offset .. cds_offset + cds.len]);
            }
            pos += cds.len;
        }

        // Copy N_REF_PAD bases downstream of the last CDS
        const last_cds = tr.cds.items[tr.cds.items.len - 1];
        const last_offset: usize = pad + last_cds.beg - tr.beg + last_cds.len;
        if (last_offset + pad <= ref.len) {
            @memcpy(sref[pos .. pos + pad], ref[last_offset .. last_offset + pad]);
        }

        tscript.sref = sref;
        tscript.nsref = @intCast(total_len);
    }

    /// Verify that the VCF REF allele matches the FASTA reference.
    ///
    /// Port of sanity_check_ref() from csq.c (line 2840).
    /// Returns error.RefAlleleMismatch on mismatch when --force is not set.
    /// Returns error.RefMismatchSkipped when --force is set (caller should
    /// skip this variant).  Returns success (void) on match.
    pub fn sanityCheckRef(self: *CsqContext, tr: *Transcript, rec: *const VcfRecord) !void {
        const tscript = self.getTscript(tr) orelse return error.TscriptNotInitialized;
        const ref = tscript.ref_seq orelse return error.RefNotLoaded;

        if (rec.alleles.len == 0) return;
        const vcf_ref = rec.alleles[0];
        if (vcf_ref.len == 0) return;

        // Calculate offset into the padded reference
        var vbeg: usize = 0;
        var rbeg_signed: i64 = @as(i64, @intCast(rec.pos)) - @as(i64, @intCast(tr.beg)) + @as(i64, N_REF_PAD);
        if (rbeg_signed < 0) {
            vbeg = @intCast(-rbeg_signed);
            rbeg_signed = 0;
        }
        const rbeg: usize = @intCast(rbeg_signed);

        if (rbeg >= ref.len or vbeg >= vcf_ref.len) return;

        // Compare character by character
        var i: usize = 0;
        while (rbeg + i < ref.len and vbeg + i < vcf_ref.len) : (i += 1) {
            const rc = std.ascii.toUpper(ref[rbeg + i]);
            const vc = std.ascii.toUpper(vcf_ref[vbeg + i]);
            if (rc != vc) {
                if (!self.force) {
                    std.log.err("the fasta reference does not match the VCF REF allele at {s}:{d} .. fasta={c} vcf={c}", .{
                        rec.seqname(), rec.pos + @as(u32, @intCast(vbeg)) + 1, rc, vc,
                    });
                    return error.RefAlleleMismatch;
                }

                if (self.verbosity > 0 and (self.warned_ref_allele_mismatch == 0 or self.verbosity > 1)) {
                    std.log.warn("the fasta reference does not match the VCF REF allele at {s}:{d} .. fasta={c} vcf={c}", .{
                        rec.seqname(), rec.pos + @as(u32, @intCast(vbeg)) + 1, rc, vc,
                    });
                    if (self.verbosity < 2) {
                        std.log.warn("This message is printed only once, the verbosity can be increased with `--verbosity 2`", .{});
                    }
                }
                self.warned_ref_allele_mismatch += 1;
                return error.RefMismatchSkipped;
            }
        }
    }

    /// Get the Tscript auxiliary data for a transcript, or null if not yet initialized.
    fn getTscript(self: *CsqContext, tr: *Transcript) ?*Tscript {
        _ = self;
        const aux = tr.aux orelse return null;
        return @as(*Tscript, @ptrCast(@alignCast(aux)));
    }

    /// Get or create the Tscript auxiliary data for a transcript.
    fn getOrCreateTscript(self: *CsqContext, tr: *Transcript) !*Tscript {
        if (tr.aux) |aux| {
            return @as(*Tscript, @ptrCast(@alignCast(aux)));
        }
        const tscript = try self.allocator.create(Tscript);
        tscript.* = .{};
        tr.aux = tscript;
        return tscript;
    }

    /// Initialize transcript auxiliary data: fetch reference, build spliced
    /// reference, and create the haplotype tree root node.
    ///
    /// This is the entry point called from testCds/testCdsLocal when a
    /// transcript is first encountered.  Corresponds to the transcript
    /// initialization block in test_cds() / test_cds_local() in csq.c.
    ///
    /// Errors from faidx fetch are propagated; the caller decides whether
    /// to skip the transcript or abort.
    pub fn initTranscriptAux(self: *CsqContext, tr: *Transcript, chr: [*:0]const u8) !void {
        // Already initialized?
        if (self.getTscript(tr)) |ts| {
            if (ts.ref_seq != null) return;
        }

        // 1. Fetch genomic reference around the transcript
        try self.tscriptInitRef(tr, chr);

        // 2. Build the spliced reference from CDS segments
        if (tr.cds.items.len > 0) {
            self.tscriptSpliceRef(tr) catch |err| {
                std.log.warn("failed to build spliced reference for transcript {d}: {}", .{ tr.id, err });
            };
        }

        // 3. Create the haplotype tree root node
        const tscript = self.getTscript(tr).?;
        if (tscript.root == null) {
            const root = try self.allocator.create(HapNode);
            root.* = HapNode.init(.root);
            tscript.root = root;
        }
    }

    /// Free transcript auxiliary data (Tscript and its owned allocations).
    ///
    /// Called when a transcript is removed from the active set (after all
    /// overlapping variants have been processed and flushed).
    pub fn destroyTranscriptAux(self: *CsqContext, tr: *Transcript) void {
        const tscript = self.getTscript(tr) orelse return;

        if (tscript.ref_seq) |ref| {
            self.allocator.free(ref);
            tscript.ref_seq = null;
        }
        if (tscript.sref) |sref| {
            self.allocator.free(sref);
            tscript.sref = null;
        }
        if (tscript.root) |root| {
            root.deinit(self.allocator);
            self.allocator.destroy(root);
            tscript.root = null;
        }
        tscript.hap.deinit(self.allocator);

        self.allocator.destroy(tscript);
        tr.aux = null;
    }

    // -----------------------------------------------------------------
    // vbufPush — add a record to the position-indexed buffer
    // -----------------------------------------------------------------

    /// Buffer a VCF record. Records at the same position share a Vbuf.
    /// Returns a pointer to the Vbuf containing the record.
    ///
    /// Port of vbuf_push() from csq.c (line 2675).
    pub const VbufPushResult = struct {
        vbuf: *Vbuf,
        owned_rec: *const VcfRecord,
    };

    pub fn vbufPush(self: *CsqContext, rec: *const VcfRecord) !VbufPushResult {
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

        // Heap-allocate a copy of the VcfRecord so it outlives the caller's stack.
        // The pipeline buffers records and flushes them later, so the record must
        // remain valid until vbufFlush processes it.
        const owned_rec = try self.allocator.create(VcfRecord);
        owned_rec.* = rec.*;
        // Deep-copy the alleles slice (points into caller's scratch buffer)
        if (rec.alleles.len > 0) {
            const alleles_copy = try self.allocator.alloc([]const u8, rec.alleles.len);
            for (rec.alleles, 0..) |a, ai| {
                alleles_copy[ai] = try self.allocator.dupe(u8, a);
            }
            owned_rec.alleles = alleles_copy;
        }
        // Deep-copy chr slice
        if (rec.chr.len > 0) {
            owned_rec.chr = try self.allocator.dupe(u8, rec.chr);
        }
        // Deep-copy raw_line
        if (rec.raw_line) |rl| {
            owned_rec.raw_line = try self.allocator.dupe(u8, rl);
        }

        // Add the record as a new Vrec
        var vrec = Vrec{};
        vrec.rec = owned_rec;
        try vbuf.vrecs.append(self.allocator, vrec);

        // Register in pos2vbuf for O(1) existence check by position
        try self.pos2vbuf.put(rec.pos, 0);

        return .{ .vbuf = vbuf, .owned_rec = owned_rec };
    }

    // -----------------------------------------------------------------
    // vbufFlush — flush records up to a position
    // -----------------------------------------------------------------

    /// A flushed record with an optional BCSQ annotation string.
    /// Callers (e.g. writeVcfRecord in main.zig) iterate flushed_records
    /// after each vbufFlush call and inject the BCSQ value into the VCF
    /// line before writing.
    pub const FlushedRecord = struct {
        /// Position and rid for looking up the original VCF line.
        pos: u32 = 0,
        rid: i32 = 0,
        /// Formatted BCSQ value, e.g. "missense|GENE|TR|protein_coding|+|5T>5I|100A>G".
        /// Null means no consequences — write the record as-is.
        bcsq_value: ?[]const u8 = null,
        /// Per-sample FORMAT/BCSQ bitmask integers (interleaved first/second haplotype).
        /// Null means no sample annotation.
        fmt_bm: ?[]const u32 = null,
        /// Number of FORMAT integers per sample.
        nfmt: u32 = 0,
    };

    /// Comparison function for sorting consequences in BCSQ output.
    /// Ordering priority (lower = first):
    ///   0: PRINTED_UPSTREAM (@-references)
    ///   1: Compound consequences (CDS-level with vstr)
    ///   2: Non-compound splice-only consequences (splice_donor, splice_acceptor, splice_region)
    ///   3: Non-compound (UTR, intron, non_coding, etc.)
    /// This matches the C code's output ordering.
    fn vcsqCmpLessThan(_: void, a: Vcsq, b: Vcsq) bool {
        return vcsqSortKey(a) < vcsqSortKey(b);
    }

    fn vcsqSortKey(v: Vcsq) u3 {
        // Match C's output order: non-CDS consequences first (pushed during
        // process by testUtr/testSplice/testTscript), then CDS consequences
        // (appended later by hapFlush/transferTreeCsqToVbuf).
        // Preserve insertion order — don't sort. Return 0 for all entries
        // so the stable sort maintains the original push order.
        // The C code outputs consequences in csq_push insertion order.
        _ = v;
        return 0;
    }

    /// Flush all buffered VCF records whose keep_until <= pos.
    /// Formats BCSQ strings and populates flushed_records.
    ///
    /// Port of vbuf_flush() from csq.c (line 2715).
    pub fn vbufFlush(self: *CsqContext, pos: u32) !void {
        // NOTE: Do NOT clear flushed_records here.  process() may call
        // vbufFlush multiple times (e.g. chromosome-change flush followed
        // by position-based flush) and the caller is responsible for
        // consuming and clearing flushed_records between process() calls.

        while (self.vcf_rbuf.len > 0) {
            const vbuf = self.vcf_rbuf.front().?;

            // Cannot flush if there are active transcripts and the vbuf's
            // keep_until extends beyond the current position (C line 2721-2726)
            if (!self.local_csq and self.active_transcripts.count() > 0) {
                if (vbuf.keep_until > pos) break;
            }

            _ = self.vcf_rbuf.shift();

            // Remove from pos2vbuf
            if (vbuf.pos()) |vpos| {
                _ = self.pos2vbuf.remove(vpos);
            }

            // Format consequences for each record in the vbuf
            for (vbuf.vrecs.items) |*vrec| {
                const rec_ptr = vrec.rec orelse continue;

                if (vrec.vcsqs.items.len == 0) {
                    // No consequences — record passes through unmodified
                    try self.flushed_records.append(self.allocator, .{
                        .pos = rec_ptr.pos,
                        .rid = rec_ptr.rid,
                    });
                    continue;
                }

                // Sort consequences: non-compound (UTR, intron, splice-only)
                // before compound (CDS-level). This matches C's output order
                // where non-CDS consequences are pushed first during process(),
                // and CDS consequences are appended later from hapFlush.
                if (vrec.vcsqs.items.len > 1) {
                    std.mem.sort(Vcsq, vrec.vcsqs.items, {}, vcsqCmpLessThan);
                }

                // Format the BCSQ INFO string
                self.output.clearRetainingCapacity();
                try formatVcsqList(
                    vrec.vcsqs.items,
                    .{
                        .brief_predictions = self.brief_predictions,
                        .trid_to_string = &tridToStringCallback,
                        .biotype_to_string = &biotypeToStringCallback,
                    },
                    self.output.writer(self.allocator),
                );

                // Dupe the formatted string so it outlives the output buffer reuse
                const bcsq_str = try self.allocator.dupe(u8, self.output.items);

                try self.flushed_records.append(self.allocator, .{
                    .pos = rec_ptr.pos,
                    .rid = rec_ptr.rid,
                    .bcsq_value = bcsq_str,
                    .fmt_bm = vrec.fmt_bm,
                    .nfmt = vrec.nfmt,
                });
                // Transfer ownership of fmt_bm to the flushed record
                // so that vrec.deinit() won't free it.
                vrec.fmt_bm = null;
            }

            vbuf.deinit(self.allocator);
            self.allocator.destroy(vbuf);
        }

        // When all active transcripts have been flushed, clean up the removal list
        if (self.active_transcripts.count() == 0) {
            for (self.rm_transcripts.items) |tr| {
                self.destroyTranscriptAux(tr);
            }
            self.rm_transcripts.clearRetainingCapacity();
        }
        self.csq_buf.clearRetainingCapacity();
    }

    // -----------------------------------------------------------------
    // hapFlush — flush completed transcript haplotypes
    // -----------------------------------------------------------------

    /// Flush haplotypes for transcripts that end at or before `pos`.
    ///
    /// Port of hap_flush() from csq.c (line 2627).
    /// Pops transcripts from the active heap whose end <= pos, finalizes
    /// their haplotype trees, stages per-sample VCF consequences, and
    /// defers transcript cleanup until after vbufFlush.
    pub fn hapFlush(self: *CsqContext, pos: u32) !void {
        while (self.active_transcripts.count() > 0) {
            const tr = self.active_transcripts.peek().?;
            if (tr.end > pos) break;

            // Pop the transcript with the smallest end position
            _ = self.active_transcripts.remove();

            // Point the haplotype context at this transcript
            self.hap_ctx.tr = tr;

            const taux: *Tscript = @ptrCast(@alignCast(tr.aux orelse {
                // No aux data -- nothing to finalize, just mark for removal
                try self.rm_transcripts.append(self.allocator, tr);
                continue;
            }));

            const root = taux.root orelse {
                try self.rm_transcripts.append(self.allocator, tr);
                continue;
            };

            if (root.children.items.len > 0) {
                // Finalize the haplotype tree: DFS traversal, translation,
                // consequence determination for each leaf path
                haplotype_mod.hapFinalize(&self.hap_ctx) catch |err| {
                    std.log.warn("hapFinalize failed for transcript {d}: {}", .{ tr.id, err });
                };

                // Transfer consequences from haplotype tree leaf nodes into vbuf vrecs.
                // hapFinalize populates csq_list on leaf nodes (types.Csq), but the
                // pipeline vbuf uses format.Vcsq. We walk the tree and push each
                // consequence into the matching vrec by position.
                self.transferTreeCsqToVbuf(root, tr) catch |err| {
                    std.log.warn("transferTreeCsqToVbuf failed for transcript {d}: {}", .{ tr.id, err });
                };

                // Stage per-sample VCF consequences (unless DROP_GT mode)
                if (self.phase != .drop_gt) {
                    const n_smpl = self.n_samples;
                    var i: u32 = 0;
                    while (i < n_smpl) : (i += 1) {
                        const ismpl: i32 = if (self.sample_indices) |idx|
                            @intCast(idx[i])
                        else
                            @intCast(i);

                        // Two haplotypes per sample
                        if (taux.hap.items.len > i * 2) {
                            self.hapStageVcf(ismpl, 0, taux.hap.items[i * 2]);
                        }
                        if (taux.hap.items.len > i * 2 + 1) {
                            self.hapStageVcf(ismpl, 1, taux.hap.items[i * 2 + 1]);
                        }
                    }
                }
            }

            // Mark transcript for deferred cleanup (cannot delete now because
            // vbuf_flush still needs the transcript data for by-position output)
            try self.rm_transcripts.append(self.allocator, tr);
        }
    }

    /// Walk the haplotype tree (DFS) and transfer consequences from leaf nodes
    /// into the vbuf vrecs.
    ///
    /// After hapFinalize, leaf nodes have csq_list populated with types.Csq entries.
    /// These need to be converted to format.Vcsq and pushed into the pipeline's
    /// Vrec.vcsqs so that vbufFlush can format the BCSQ string.
    fn transferTreeCsqToVbuf(self: *CsqContext, root: *HapNode, tr: *Transcript) !void {
        // DFS traversal using an explicit stack
        var stack: std.ArrayList(*HapNode) = .empty;
        defer stack.deinit(self.allocator);
        try stack.append(self.allocator, root);

        while (stack.items.len > 0) {
            const node = stack.pop() orelse break;

            // Push children in reverse order so they pop in left-to-right order
            // (matching the C code's istack-based DFS traversal).
            {
                var ci: usize = node.children.items.len;
                while (ci > 0) {
                    ci -= 1;
                    try stack.append(self.allocator, node.children.items[ci]);
                }
            }

            // Process consequences on this node
            if (node.csq_list.items.len == 0) continue;

            for (node.csq_list.items) |csq_entry| {
                // Convert types.Vcsq -> format.Vcsq
                const vcsq = Vcsq{
                    .csq_type = csq_entry.type_info.csq_type.toInt(),
                    .biotype = csq_entry.type_info.biotype,
                    .strand = if (csq_entry.type_info.strand) .fwd else .rev,
                    .trid = csq_entry.type_info.trid,
                    .vcf_ial = csq_entry.type_info.vcf_ial,
                    .gene = csq_entry.type_info.gene,
                    .ref_pos = csq_entry.ref_pos,
                    .vstr = if (csq_entry.type_info.vstr.items.len > 0)
                        csq_entry.type_info.vstr.items
                    else
                        null,
                };

                // Find the vbuf at this position and push the consequence
                try self.pushCsqToVbufByPos(csq_entry.pos, vcsq, tr);
            }
        }
    }

    /// Push a consequence into the vbuf vrec matching the given position.
    /// Unlike csqPush which requires pointer identity, this matches by position
    /// and allele index to find the correct vrec.
    fn pushCsqToVbufByPos(self: *CsqContext, rec_pos: u32, vcsq: Vcsq, tr: *Transcript) !void {
        _ = tr;

        // Find the vbuf at this position
        var vbuf: ?*Vbuf = null;
        for (0..self.vcf_rbuf.len) |k| {
            const candidate = self.vcf_rbuf.kth(k);
            if (candidate.pos()) |p| {
                if (p == rec_pos) {
                    vbuf = candidate;
                    break;
                }
            }
        }
        const vb = vbuf orelse return; // record may have already been flushed

        // Find the vrec matching this allele.
        // For multi-allelic sites, match by vcf_ial; for biallelic, use the first vrec.
        var target_vrec: ?*Vrec = null;
        for (vb.vrecs.items) |*vrec| {
            const rec_ptr = vrec.rec orelse continue;
            // Match: the vrec's record must have enough alleles for this ial
            if (vcsq.vcf_ial < rec_ptr.n_allele) {
                target_vrec = vrec;
                break;
            }
        }
        if (target_vrec == null and vb.vrecs.items.len > 0) {
            // Fallback: use the first vrec at this position
            target_vrec = &vb.vrecs.items[0];
        }
        const vrec = target_vrec orelse return;

        // Apply type masking rules (same as csqPush)
        var t = vcsq.csq_type;
        if (t & CSQ_INFRAME_INSERTION != 0 and t & CSQ_ELONGATION != 0) t &= ~CSQ_INFRAME_INSERTION;
        if (t & CSQ_INFRAME_DELETION != 0 and t & CSQ_TRUNCATION != 0) t &= ~CSQ_INFRAME_DELETION;
        if (t & CSQ_SPLICE_REGION != 0 and t & (CSQ_SPLICE_DONOR | CSQ_SPLICE_ACCEPTOR) != 0) {
            t &= ~CSQ_SPLICE_REGION;
        }
        // Remove stop_lost & synonymous if stop_retained set
        if (t & CSQ_STOP_RETAINED != 0) t &= ~(CSQ_STOP_LOST | CSQ_SYNONYMOUS_VARIANT);
        // Remove start_lost & synonymous if start_retained set
        if (t & CSQ_START_RETAINED != 0) t &= ~(CSQ_START_LOST | CSQ_SYNONYMOUS_VARIANT);
        var masked_vcsq = vcsq;
        masked_vcsq.csq_type = t;

        // Deduplication: special handling for CSQ_PRINTED_UPSTREAM (C csq.c lines 2018-2033)
        if (t & CSQ_PRINTED_UPSTREAM != 0) {
            for (vrec.vcsqs.items) |*existing| {
                // START_STOP replaces START_STOP
                if (t & CSQ_START_STOP != 0 and existing.csq_type & CSQ_START_STOP != 0) {
                    existing.* = masked_vcsq;
                    return;
                }
                // Only match existing PRINTED_UPSTREAM with same ref_pos
                if (existing.csq_type & CSQ_PRINTED_UPSTREAM == 0) continue;
                if (existing.ref_pos != null and masked_vcsq.ref_pos != null and
                    existing.ref_pos.? == masked_vcsq.ref_pos.?)
                {
                    return; // duplicate
                }
            }
            // Not a duplicate: append
            try vrec.vcsqs.append(self.allocator, masked_vcsq);
            return;
        }

        // Standard deduplication: check if an identical consequence already exists
        for (vrec.vcsqs.items) |*existing| {
            if (isDuplicate(existing, &masked_vcsq)) {
                existing.csq_type |= t;

                // Remove stop_lost & synonymous if stop_retained set (C line 2059-2060)
                if (existing.csq_type & CSQ_STOP_RETAINED != 0)
                    existing.csq_type &= ~(CSQ_STOP_LOST | CSQ_SYNONYMOUS_VARIANT);

                // Remove start_lost & synonymous if start_retained set (C line 2062-2063)
                if (existing.csq_type & CSQ_START_RETAINED != 0)
                    existing.csq_type &= ~(CSQ_START_LOST | CSQ_SYNONYMOUS_VARIANT);

                // Copy vstr from new to existing only for compound merges (C line 2065).
                // Non-compound merges (C line 2082-2085) do NOT copy vstr.
                if (t & CSQ_COMPOUND != 0) {
                    if (existing.vstr == null and masked_vcsq.vstr != null) {
                        existing.vstr = masked_vcsq.vstr;
                    }
                }
                return;
            }
        }

        // Append new consequence
        try vrec.vcsqs.append(self.allocator, masked_vcsq);
    }

    /// Stage VCF consequence bitmask bits for a single sample/haplotype leaf node.
    ///
    /// Port of hap_stage_vcf() from csq.c (line 2597).
    /// For each consequence in the leaf node's csq_list, find the matching
    /// consequence in the pipeline vrec and set the appropriate bitmask bit.
    fn hapStageVcf(self: *CsqContext, ismpl: i32, ihap: u1, node: *HapNode) void {
        if (ismpl < 0) return;
        if (node.csq_list.items.len == 0) return;

        for (node.csq_list.items) |csq| {
            // Find the pipeline vrec at this position
            var target_vrec: ?*Vrec = null;
            for (0..self.vcf_rbuf.len) |k| {
                const candidate = self.vcf_rbuf.kth(k);
                if (candidate.pos()) |p| {
                    if (p == csq.pos) {
                        if (candidate.vrecs.items.len > 0) {
                            target_vrec = &candidate.vrecs.items[0];
                        }
                        break;
                    }
                }
            }
            const vrec = target_vrec orelse continue;

            // Find the matching consequence index in the vrec's vcsqs list
            const csq_type_raw = csq.type_info.csq_type.toInt();
            var csq_idx: ?u32 = null;
            for (vrec.vcsqs.items, 0..) |*vcsq, idx| {
                // Match by type (including upstream_stop) and vstr
                const existing_type = vcsq.csq_type;
                // For PRINTED_UPSTREAM, match by ref_pos
                if (csq_type_raw & CSQ_PRINTED_UPSTREAM != 0) {
                    if (existing_type & CSQ_PRINTED_UPSTREAM != 0) {
                        csq_idx = @intCast(idx);
                        break;
                    }
                    continue;
                }
                // Match type bits (mask out printed_upstream for comparison)
                const type_mask = ~CSQ_PRINTED_UPSTREAM;
                if ((existing_type & type_mask) != (csq_type_raw & type_mask)) continue;
                // Match vstr
                const csq_vstr = if (csq.type_info.vstr.items.len > 0) csq.type_info.vstr.items else "";
                const existing_vstr = vcsq.vstr orelse "";
                if (!std.mem.eql(u8, csq_vstr, existing_vstr)) continue;
                csq_idx = @intCast(idx);
                break;
            }

            const ci = csq_idx orelse continue;
            const icsq2: u32 = ci * 2 + @as(u32, ihap);

            if (icsq2 >= self.ncsq2_max) break;

            const ival: u32 = icsq2 / 30;
            const ibit: u5 = @intCast(icsq2 % 30);
            if (vrec.nfmt < 1 + ival) vrec.nfmt = 1 + ival;

            // Allocate fmt_bm if needed
            if (vrec.fmt_bm == null) {
                const bm_size = self.n_samples * self.nfmt_bcsq;
                if (bm_size > 0) {
                    vrec.fmt_bm = self.allocator.alloc(u32, bm_size) catch continue;
                    @memset(vrec.fmt_bm.?, 0);
                }
            }

            // Set the bit
            if (vrec.fmt_bm) |bm| {
                const sample_u: usize = @intCast(ismpl);
                const offset = sample_u * self.nfmt_bcsq + ival;
                if (offset < bm.len) {
                    bm[offset] |= @as(u32, 1) << ibit;
                }
            }
        }
    }

    /// Convert a doubled consequence index to the (ival, ibit) pair for
    /// indexing into the fmt_bm bitmask array.
    ///
    /// Port of icsq2_to_bit() from csq.c (line 585).
    pub fn icsq2ToBit(icsq2: u32) struct { ival: u32, ibit: u5 } {
        return .{
            .ival = icsq2 / 30,
            .ibit = @intCast(icsq2 % 30),
        };
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

        const push_result = try self.vbufPush(rec);
        const vbuf = push_result.vbuf;
        // Use the heap-owned record for haplotype tree identity checks.
        // The caller's `rec` may be a stack variable reused across calls,
        // so pointer identity wouldn't distinguish different records.
        const owned_rec = push_result.owned_rec;

        // Check for symbolic ALTs
        if (rec.alleles.len >= 2 and rec.alleles[1].len > 0 and rec.alleles[1][0] == '<') {
            try self.testSymbolicAlt(rec);
        } else {
            // Annotation lookup: CDS, UTR, splice are all checked independently;
            // only tscript (intron/non-coding) is skipped if any of the above hit.
            // This matches the C code: hit = test_cds(); hit += test_utr(); hit += test_splice();
            var hit: bool = false;
            if (self.local_csq) {
                hit = try self.testCdsLocal(rec);
            } else {
                hit = try self.testCds(owned_rec, vbuf);
            }
            {
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
    /// For the haplotype-aware path, this queries the CDS region index and
    /// extends the per-transcript haplotype tree. Currently implements the
    /// drop_gt (no-genotype) simplified path; full sample-aware haplotype
    /// extension requires htslib genotype access.
    fn testCds(self: *CsqContext, rec: *const VcfRecord, vbuf: *Vbuf) !bool {
        const gff = self.gff orelse return false;
        const chr = rec.seqname();

        // Note: off-by-one extension of rlen is deliberate to account for insertions
        var itr = gff.idx_cds.overlap(chr, rec.pos, rec.pos + rec.rlen);

        var ret = false;
        while (itr.next()) |interval| {
            const cds: *CdsEntry = interval.payload;
            const tr: *Transcript = cds.tr;
            if (!tr.biotype.isCoding()) continue;

            // Extend the vbuf keep_until to cover the full transcript
            if (vbuf.keep_until < tr.end) vbuf.keep_until = tr.end;
            ret = true;

            // Initialize transcript aux if first time
            if (tr.aux == null) {
                const taux = try self.allocator.create(Tscript);
                taux.* = .{};
                tr.aux = taux;
                // Create haplotype tree root node
                const root = try self.allocator.create(HapNode);
                root.* = HapNode.init(.root);
                taux.root = root;
                const nhap: u32 = if (self.phase == .drop_gt) 1 else 2 * self.n_samples;
                root.nend = nhap;
                // Fetch the FASTA reference for this transcript (if faidx available)
                self.tscriptInitRef(tr, rec.chrZ()) catch |err| switch (err) {
                    error.FaidxSkipped => {
                        // --force: skip this transcript, clean up aux
                        root.deinit(self.allocator);
                        self.allocator.destroy(root);
                        self.allocator.destroy(taux);
                        tr.aux = null;
                        continue;
                    },
                    error.FaidxFetchFailed => {
                        // No fasta: clean up aux and skip (matches C: free(tr->aux); tr->aux=NULL; continue)
                        root.deinit(self.allocator);
                        self.allocator.destroy(root);
                        self.allocator.destroy(taux);
                        tr.aux = null;
                        continue;
                    },
                    else => return err,
                };
                // Build the spliced reference from CDS segments
                self.tscriptSpliceRef(tr) catch {};
                try self.active_transcripts.add(tr);
            }

            const taux_ptr: *Tscript = @ptrCast(@alignCast(tr.aux orelse continue));
            // Verify VCF REF allele matches the FASTA reference
            self.sanityCheckRef(tr, rec) catch |err| switch (err) {
                error.RefMismatchSkipped => continue,
                error.RefAlleleMismatch => return err,
                error.TscriptNotInitialized, error.RefNotLoaded => {},
                else => return err,
            };

            if (self.phase == .drop_gt) {
                // Simplified path: single haplotype, no genotype tracking.
                if (rec.alleles.len < 2) continue;
                const alt = rec.alleles[1];
                if (alt.len > 0 and (alt[0] == '<' or alt[0] == '*')) continue;

                // Get current leaf or root
                const parent: *HapNode = if (taux_ptr.hap.items.len > 0)
                    taux_ptr.hap.items[0]
                else
                    (taux_ptr.root orelse continue);

                var child = try self.allocator.create(HapNode);
                child.* = HapNode.init(.cds);
                const ref_allele = rec.alleles[0];

                const hap_ret = haplotype_mod.hapInit(
                    self.allocator,
                    parent,
                    child,
                    cds,
                    rec.pos,
                    ref_allele,
                    alt,
                    1,
                    taux_ptr,
                ) catch {
                    self.allocator.destroy(child);
                    continue;
                };

                switch (hap_ret.kind) {
                    .overlapping => {
                        child.deinit(self.allocator);
                        self.allocator.destroy(child);
                        continue;
                    },
                    .discarded => {
                        child.deinit(self.allocator);
                        self.allocator.destroy(child);
                        ret = true;
                        continue;
                    },
                    .added => {},
                }

                // Stage the splice consequence (mirrors C csq_stage_splice called
                // from within splice_csq during hap_init).  This creates an early
                // "placeholder" entry that later gets merged with the compound
                // consequence from hapFinalize, producing combined annotations
                // like "start_lost&splice_region".
                // Use the pre-clearing splice_csq from hapInit which preserves
                // synonymous_variant (C stages this via csq_stage_splice before
                // clearing synonymous for the CDS path).
                if (hap_ret.splice_csq.toInt() != 0) {
                    var splice_csq = Csq{
                        .pos = rec.pos,
                        .vcsq = .{
                            .csq_type = hap_ret.splice_csq.toInt(),
                            .biotype = @intFromEnum(tr.biotype),
                            .strand = if (tr.strand == .forward) .fwd else .rev,
                            .trid = tr.id,
                            .vcf_ial = 1,
                            .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                        },
                    };
                    _ = self.csqStage(&splice_csq, rec) catch {};
                }

                // Splice-only (HAP_SSS): stage the splice consequence directly
                if (child.payload == .sss) {
                    var csq = Csq{
                        .pos = rec.pos,
                        .vcsq = .{
                            .csq_type = child.csq.toInt(),
                            .biotype = @intFromEnum(tr.biotype),
                            .strand = if (tr.strand == .forward) .fwd else .rev,
                            .trid = tr.id,
                            .vcf_ial = 1,
                            .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                        },
                    };
                    try self.csqStage(&csq, rec);
                    child.deinit(self.allocator);
                    self.allocator.destroy(child);
                    ret = true;
                    continue;
                }

                // Attach child to parent in the haplotype tree
                parent.nend -= 1;
                try parent.children.append(self.allocator, child);
                if (taux_ptr.hap.items.len == 0) {
                    try taux_ptr.hap.append(self.allocator, child);
                } else {
                    taux_ptr.hap.items[0] = child;
                }
                taux_ptr.hap.items[0].nend = 1;
                continue;
            }

            // ── Genotype-aware per-sample haplotype path ──────────────
            // Port of C lines 3168-3288: iterate samples & haplotypes,
            // extend the per-transcript haplotype tree for each non-ref allele.
            const genotypes = rec.parseGenotypes(self.allocator) catch null;
            if (genotypes == null) continue;
            const gts = genotypes.?;
            defer self.allocator.free(gts);

            if (gts.len == 0) continue;

            // Use maximum ploidy across all samples (C uses ngt from
            // bcf_get_genotypes which pads haploid samples with missing).
            // This ensures diploid samples are fully processed even when
            // some samples are haploid.
            var ngts: u32 = 0;
            for (gts) |g| {
                if (g.ploidy > ngts) ngts = g.ploidy;
            }
            if (ngts != 1 and ngts != 2) {
                // Non-haploid/diploid: skip (warn once)
                if (self.verbosity > 0) {
                    std.log.warn("Skipping site with non-diploid/non-haploid genotypes at {s}:{d}", .{ chr, rec.pos + 1 });
                }
                continue;
            }

            // Ensure the hap array is large enough: 2 * n_samples entries
            const n_smpl = if (self.sample_indices) |si| @as(u32, @intCast(si.len)) else @as(u32, @intCast(gts.len));
            const nhap_needed: usize = 2 * @as(usize, n_smpl);
            while (taux_ptr.hap.items.len < nhap_needed) {
                try taux_ptr.hap.append(self.allocator, taux_ptr.root orelse continue);
            }
            // Note: root.nend is initialized at transcript creation (nhap = 2 * n_samples)
            // and should NOT be re-initialized here, as haplotypes may have already
            // moved from root to child nodes, decrementing root.nend correctly.

            for (0..n_smpl) |ismpl_idx| {
                const ismpl: usize = if (self.sample_indices) |si| @as(usize, si[ismpl_idx]) else ismpl_idx;
                if (ismpl >= gts.len) continue;

                var gt = gts[ismpl];
                if (gt.alleles[0] < 0) continue; // first allele missing

                // Handle unphased heterozygous
                if (ngts > 1 and gt.alleles[1] >= 0 and gt.alleles[0] != gt.alleles[1]) {
                    if (self.phase == .merge) {
                        if (gt.alleles[0] == 0) gt.alleles[0] = gt.alleles[1];
                    }
                    if (!gt.phased) {
                        switch (self.phase) {
                            .require => return error.UnphasedHeterozygous,
                            .skip => continue,
                            .non_ref => {
                                if (gt.alleles[0] == 0) {
                                    gt.alleles[0] = gt.alleles[1];
                                } else if (gt.alleles[1] == 0) {
                                    gt.alleles[1] = gt.alleles[0];
                                }
                            },
                            else => {},
                        }
                    }
                }

                var ihap: u32 = 0;
                while (ihap < ngts) : (ihap += 1) {
                    if (gt.alleles[ihap] <= 0) continue; // missing or ref
                    const ial: u32 = @intCast(gt.alleles[ihap]);
                    if (ial >= rec.n_allele) continue;
                    if (ial >= rec.alleles.len) continue;

                    const alt = rec.alleles[ial];
                    if (alt.len > 0 and (alt[0] == '<' or alt[0] == '*')) continue;

                    const i: usize = 2 * ismpl_idx + ihap;
                    const root_ptr = taux_ptr.root orelse continue;
                    const parent: *HapNode = if (i < taux_ptr.hap.items.len and taux_ptr.hap.items[i] != root_ptr)
                        taux_ptr.hap.items[i]
                    else
                        root_ptr;

                    // Check if this haplotype already seen for another sample at this record
                    const rec_opaque: *const anyopaque = @ptrCast(rec);
                    if (parent.cur_rec != null and parent.cur_rec.? == rec_opaque) {
                        // Look up the cached child for this allele
                        if (ial < parent.cur_child.items.len) {
                            const cached_idx = parent.cur_child.items[ial];
                            if (cached_idx >= 0 and @as(usize, @intCast(cached_idx)) < parent.children.items.len) {
                                taux_ptr.hap.items[i] = parent.children.items[@intCast(cached_idx)];
                                taux_ptr.hap.items[i].nend += 1;
                                parent.nend -|= 1;
                                continue;
                            }
                        }
                    }

                    var child = try self.allocator.create(HapNode);
                    child.* = HapNode.init(.cds);
                    const ref_allele = rec.alleles[0];

                    const hap_ret = haplotype_mod.hapInit(
                        self.allocator,
                        parent,
                        child,
                        cds,
                        rec.pos,
                        ref_allele,
                        alt,
                        ial,
                        taux_ptr,
                    ) catch {
                        self.allocator.destroy(child);
                        continue;
                    };

                    switch (hap_ret.kind) {
                        .overlapping => {
                            child.deinit(self.allocator);
                            self.allocator.destroy(child);
                            continue;
                        },
                        .discarded => {
                            child.deinit(self.allocator);
                            self.allocator.destroy(child);
                            continue;
                        },
                        .added => {},
                    }

                    // Stage the pre-clearing splice consequence (mirrors C csq_stage_splice)
                    if (hap_ret.splice_csq.toInt() != 0) {
                        var splice_csq_entry = Csq{
                            .pos = rec.pos,
                            .vcsq = .{
                                .csq_type = hap_ret.splice_csq.toInt(),
                                .biotype = @intFromEnum(tr.biotype),
                                .strand = if (tr.strand == .forward) .fwd else .rev,
                                .trid = tr.id,
                                .vcf_ial = ial,
                                .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                            },
                        };
                        _ = self.csqStage(&splice_csq_entry, rec) catch {};
                    }

                    // Splice-only (HAP_SSS): stage the splice consequence directly
                    if (child.payload == .sss) {
                        var csq_sss = Csq{
                            .pos = rec.pos,
                            .vcsq = .{
                                .csq_type = child.csq.toInt(),
                                .biotype = @intFromEnum(tr.biotype),
                                .strand = if (tr.strand == .forward) .fwd else .rev,
                                .trid = tr.id,
                                .vcf_ial = ial,
                                .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                            },
                        };
                        try self.csqStage(&csq_sss, rec);
                        child.deinit(self.allocator);
                        self.allocator.destroy(child);
                        continue;
                    }

                    // Initialize cur_child tracking on the parent for this record
                    if (parent.cur_rec == null or parent.cur_rec.? != rec_opaque) {
                        parent.cur_child.clearRetainingCapacity();
                        while (parent.cur_child.items.len < rec.n_allele) {
                            try parent.cur_child.append(self.allocator, -1);
                        }
                        parent.cur_rec = rec_opaque;
                    }

                    // Attach child to parent in the haplotype tree
                    const child_idx: i32 = @intCast(parent.children.items.len);
                    if (ial < parent.cur_child.items.len) {
                        parent.cur_child.items[ial] = child_idx;
                    }
                    try parent.children.append(self.allocator, child);
                    taux_ptr.hap.items[i] = child;
                    taux_ptr.hap.items[i].nend += 1;
                    parent.nend -|= 1;
                    parent.nend -|= 1;
                }
            }
        }
        return ret;
    }

    /// Local (non-haplotype-aware) CDS consequence calling.
    ///
    /// Port of test_cds_local() from csq.c (line 2872).
    /// Queries the CDS region index and, for each overlapping CDS and each alt
    /// allele, determines the coding consequence. This is the simplified path
    /// that does not track per-sample haplotypes.
    fn testCdsLocal(self: *CsqContext, rec: *const VcfRecord) !bool {
        const gff = self.gff orelse return false;
        const chr = rec.seqname();

        var itr = gff.idx_cds.overlap(chr, rec.pos, rec.pos + rec.rlen);

        // Working buffers for translation (reused across iterations)
        var tref_buf: std.ArrayList(u8) = .empty;
        defer tref_buf.deinit(self.allocator);
        var tseq_buf: std.ArrayList(u8) = .empty;
        defer tseq_buf.deinit(self.allocator);
        var tref_stop_buf: std.ArrayList(u8) = .empty;
        defer tref_stop_buf.deinit(self.allocator);
        var tseq_stop_buf: std.ArrayList(u8) = .empty;
        defer tseq_stop_buf.deinit(self.allocator);

        var ret = false;
        while (itr.next()) |interval| {
            const cds: *CdsEntry = interval.payload;
            const tr: *Transcript = cds.tr;
            if (!tr.biotype.isCoding()) continue;
            ret = true;

            // Initialize transcript aux if first time
            if (tr.aux == null) {
                const taux_new = try self.allocator.create(Tscript);
                taux_new.* = .{};
                tr.aux = taux_new;
                // Fetch FASTA reference and build spliced reference
                self.tscriptInitRef(tr, rec.chrZ()) catch |err| switch (err) {
                    error.FaidxSkipped, error.FaidxFetchFailed => {
                        // No fasta: clean up aux (matches C: free(tr->aux); tr->aux=NULL; continue)
                        self.allocator.destroy(taux_new);
                        tr.aux = null;
                        continue;
                    },
                    else => return err,
                };
                self.tscriptSpliceRef(tr) catch {};
                try self.active_transcripts.add(tr);
            }

            const taux: *Tscript = @ptrCast(@alignCast(tr.aux orelse continue));
            // Verify VCF REF allele matches the FASTA reference
            self.sanityCheckRef(tr, rec) catch |err| switch (err) {
                error.RefMismatchSkipped => continue,
                error.RefAlleleMismatch => return err,
                error.TscriptNotInitialized, error.RefNotLoaded => {},
                else => return err,
            };

            // For each alt allele
            var ial: u32 = 1;
            while (ial < rec.n_allele) : (ial += 1) {
                if (ial >= rec.alleles.len) break;
                const alt = rec.alleles[ial];
                if (alt.len > 0 and (alt[0] == '<' or alt[0] == '*')) continue;

                const ref_allele = rec.alleles[0];

                // Use a temporary root to do single-variant hapInit
                var tmp_root = HapNode.init(.root);
                defer tmp_root.deinit(self.allocator);
                var node = HapNode.init(.cds);

                const hap_ret = haplotype_mod.hapInit(
                    self.allocator,
                    &tmp_root,
                    &node,
                    cds,
                    rec.pos,
                    ref_allele,
                    alt,
                    ial,
                    taux,
                ) catch continue;

                if (hap_ret.kind != .added) continue;
                defer {
                    if (node.payload == .cds) {
                        if (node.payload.cds.seq) |seq| self.allocator.free(seq);
                    }
                    if (node.var_str) |vs| self.allocator.free(vs);
                }

                var csq = Csq{
                    .pos = rec.pos,
                    .vcsq = .{
                        .biotype = @intFromEnum(tr.biotype),
                        .strand = if (tr.strand == .forward) .fwd else .rev,
                        .trid = tr.id,
                        .vcf_ial = ial,
                        .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                    },
                };

                var csq_type: u32 = node.csq.toInt();

                // Splice-only node (HAP_SSS): stage directly
                if (node.payload == .sss) {
                    csq.vcsq.csq_type = csq_type;
                    try self.csqStage(&csq, rec);
                    continue;
                }

                // CDS node: translate ref and alt, compare amino acids
                const sref = taux.sref orelse {
                    // No spliced reference available; fall back to coding_sequence
                    csq.vcsq.csq_type = CSQ_CODING_SEQUENCE;
                    try self.csqStage(&csq, rec);
                    continue;
                };
                const sref_len: usize = @intCast(taux.nsref);
                const gencode = translate.findGeneticCode(0) orelse {
                    csq.vcsq.csq_type = CSQ_CODING_SEQUENCE;
                    try self.csqStage(&csq, rec);
                    continue;
                };

                // Translate the alt allele
                const node_seq = if (node.payload == .cds) node.payload.cds.seq else null;
                const alen: usize = if (node_seq) |s| s.len else 0;
                const fill_val: i32 = if (@rem(node.dlen, 3) != 0 and alen > 0) 1 else 0;

                haplotype_mod.cdsTranslate(
                    self.allocator,
                    sref,
                    sref_len,
                    node_seq orelse &[_]u8{},
                    if (sref_len >= 2 * N_REF_PAD) sref_len - 2 * N_REF_PAD + @as(usize, @intCast(@max(node.dlen, 0))) else 0,
                    node.sbeg,
                    node.sbeg,
                    node.sbeg + @as(u32, @intCast(@max(@as(i32, 0), node.rlen))),
                    tr.strand,
                    &tseq_buf,
                    &tseq_stop_buf,
                    fill_val,
                    gencode,
                ) catch {
                    csq.vcsq.csq_type = CSQ_CODING_SEQUENCE;
                    try self.csqStage(&csq, rec);
                    continue;
                };

                // Translate the reference
                {
                    const ref_start = N_REF_PAD + node.sbeg;
                    const ref_len_u: u32 = @intCast(@max(@as(i32, 0), node.rlen));
                    const ref_slice = if (ref_start + ref_len_u <= sref.len)
                        sref[ref_start .. ref_start + ref_len_u]
                    else
                        &[_]u8{};

                    haplotype_mod.cdsTranslate(
                        self.allocator,
                        sref,
                        sref_len,
                        ref_slice,
                        if (sref_len >= 2 * N_REF_PAD) sref_len - 2 * N_REF_PAD else 0,
                        node.sbeg,
                        node.sbeg,
                        node.sbeg + ref_len_u,
                        tr.strand,
                        &tref_buf,
                        &tref_stop_buf,
                        fill_val,
                        gencode,
                    ) catch {
                        csq.vcsq.csq_type = CSQ_CODING_SEQUENCE;
                        try self.csqStage(&csq, rec);
                        continue;
                    };
                }

                // Use hapAddCsq to determine the consequence type
                const csq_result = haplotype_mod.hapAddCsq(
                    tref_buf.items,
                    tref_stop_buf.items,
                    tseq_buf.items,
                    tseq_stop_buf.items,
                    node.dlen,
                    node.dlen != 0,
                    types.CsqType.fromInt(csq_type),
                    false, // not sss
                    false, // not compound (single variant)
                    false, // no upstream stop
                );
                csq_type = csq_result.csq_type.toInt();

                // Stage compound consequences (with variant string)
                if (csq_type & format.CSQ_COMPOUND != 0) {
                    var vstr_buf = std.ArrayList(u8).empty;
                    defer vstr_buf.deinit(self.allocator);

                    const aa_rbeg: usize = if (tr.strand == .forward)
                        node.sbeg / 3 + 1
                    else blk: {
                        const nsref_coding = if (sref_len >= 2 * N_REF_PAD) sref_len - 2 * N_REF_PAD else 0;
                        const rlen_u: usize = @intCast(@max(@as(i32, 0), node.rlen));
                        break :blk (nsref_coding -| node.sbeg -| rlen_u) / 3 + 1;
                    };

                    try vstr_buf.append(self.allocator, '|');
                    {
                        var num_buf: [32]u8 = undefined;
                        const s = std.fmt.bufPrint(&num_buf, "{d}", .{aa_rbeg}) catch unreachable;
                        try vstr_buf.appendSlice(self.allocator, s);
                    }
                    try vstr_buf.appendSlice(self.allocator, tref_buf.items);
                    if (csq_type & CSQ_SYNONYMOUS_VARIANT == 0) {
                        try vstr_buf.append(self.allocator, '>');
                        const aa_sbeg: usize = if (tr.strand == .forward)
                            node.sbeg / 3 + 1
                        else blk: {
                            const nsref_coding_dlen = if (sref_len >= 2 * N_REF_PAD)
                                @as(i64, @intCast(sref_len - 2 * N_REF_PAD)) + node.dlen
                            else
                                @as(i64, node.dlen);
                            break :blk @as(usize, @intCast(@max(nsref_coding_dlen - @as(i64, @intCast(node.sbeg)) - @as(i64, @intCast(alen)), 0))) / 3 + 1;
                        };
                        var num_buf: [32]u8 = undefined;
                        const s = std.fmt.bufPrint(&num_buf, "{d}", .{aa_sbeg}) catch unreachable;
                        try vstr_buf.appendSlice(self.allocator, s);
                        try vstr_buf.appendSlice(self.allocator, tseq_buf.items);
                    }
                    try vstr_buf.append(self.allocator, '|');
                    {
                        var num_buf: [32]u8 = undefined;
                        const s = std.fmt.bufPrint(&num_buf, "{d}", .{rec.pos + 1}) catch unreachable;
                        try vstr_buf.appendSlice(self.allocator, s);
                    }
                    if (node.var_str) |vs| try vstr_buf.appendSlice(self.allocator, vs);

                    csq.vcsq.vstr = try self.allocator.dupe(u8, vstr_buf.items);
                    csq.vcsq.csq_type = csq_type & format.CSQ_COMPOUND;
                    try self.csqStage(&csq, rec);

                    // Track vstr for cleanup
                    if (taux.root == null) {
                        const root = try self.allocator.create(HapNode);
                        root.* = HapNode.init(.root);
                        taux.root = root;
                    }
                }

                // Stage non-compound consequences separately
                if (csq_type & ~format.CSQ_COMPOUND != 0) {
                    csq.vcsq.csq_type = csq_type & ~format.CSQ_COMPOUND;
                    csq.vcsq.vstr = null;
                    try self.csqStage(&csq, rec);
                }
            }
        }
        return ret;
    }

    // -----------------------------------------------------------------
    // testUtr — look up UTR regions
    // -----------------------------------------------------------------

    /// Check if the variant overlaps UTR regions.
    ///
    /// Port of test_utr() from csq.c (line 3360).
    /// Queries the UTR region index. For each overlapping UTR and each alt
    /// allele, runs splice analysis and stages a UTR5 or UTR3 consequence
    /// if the variant falls inside the UTR.
    fn testUtr(self: *CsqContext, rec: *const VcfRecord) !bool {
        const gff = self.gff orelse return false;
        const chr = rec.seqname();

        var itr = gff.idx_utr.overlap(chr, rec.pos, rec.pos + rec.rlen);

        var ret = false;
        while (itr.next()) |interval| {
            const utr: *Utr = interval.payload;
            const tr: *Transcript = utr.tr;

            // For each alt allele
            var ial: u32 = 1;
            while (ial < rec.n_allele) : (ial += 1) {
                if (ial >= rec.alleles.len) break;
                const alt = rec.alleles[ial];
                if (alt.len > 0 and (alt[0] == '<' or alt[0] == '*')) continue;

                const ref_allele = rec.alleles[0];

                // Run splice analysis
                var splice = Splice.init(self.allocator, tr);
                defer splice.deinit();

                splice.reset(
                    @intCast(rec.pos),
                    @intCast(ref_allele.len),
                    @intCast(ial),
                    ref_allele,
                    alt,
                );

                const splice_ret = splice.spliceCsq(utr.beg, utr.end);
                if (splice_ret != .inside and splice_ret != .overlap) {
                    // For insertions at the exact CDS/UTR boundary (e.g. last CDS base
                    // where the insertion position equals the UTR start), spliceCsq
                    // returns .outside but C stages the UTR consequence internally
                    // (csq.c line 1101-1117, check_utr path from hapInit).
                    // This only applies when the variant also overlaps a CDS (i.e.,
                    // the insertion is at a CDS boundary, not intergenic).
                    if (splice_ret == .outside and ref_allele.len < alt.len) {
                        const last_ref_base = rec.pos + @as(u32, @intCast(ref_allele.len)) - 1;
                        // Check: last ref base is adjacent to UTR AND overlaps a CDS
                        var has_adjacent_cds = false;
                        if (last_ref_base + 1 == utr.beg) {
                            var cds_itr = gff.idx_cds.overlap(chr, rec.pos, rec.pos + rec.rlen);
                            if (cds_itr.next() != null) has_adjacent_cds = true;
                        }
                        if (has_adjacent_cds) {
                            // Stage UTR consequence for CDS/UTR boundary insertion
                        } else {
                            continue;
                        }
                    } else {
                        continue;
                    }
                }

                // Determine UTR type
                const utr_csq: CsqType = if (utr.which == .prime5) CSQ_UTR5 else CSQ_UTR3;

                var csq = Csq{
                    .pos = rec.pos,
                    .vcsq = .{
                        .csq_type = utr_csq,
                        .biotype = @intFromEnum(tr.biotype),
                        .strand = if (tr.strand == .forward) .fwd else .rev,
                        .trid = tr.id,
                        .vcf_ial = ial,
                        .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                    },
                };
                try self.csqStage(&csq, rec);
                ret = true;
            }
        }
        return ret;
    }

    // -----------------------------------------------------------------
    // testSplice — check for splice site variants
    // -----------------------------------------------------------------

    /// Check if the variant affects splice sites.
    ///
    /// Port of test_splice() from csq.c (line 3400).
    /// Queries the exon region index. For each overlapping exon in a coding
    /// transcript, runs splice analysis with donor/acceptor checking enabled.
    /// Splice consequences are staged within the splice analysis itself
    /// (via the csq flags on the Splice struct); here we just check if any
    /// consequence was set and return accordingly.
    fn testSplice(self: *CsqContext, rec: *const VcfRecord) !bool {
        const gff = self.gff orelse return false;
        const chr = rec.seqname();

        var itr = gff.idx_exon.overlap(chr, rec.pos, rec.pos + rec.rlen);

        var ret = false;
        while (itr.next()) |interval| {
            const exon: *Exon = interval.payload;
            const tr: *Transcript = exon.tr;

            // Skip non-coding transcripts (no CDS entries)
            if (tr.cds.items.len == 0) continue;

            // Determine whether to check region boundaries.
            // If the exon starts/ends at the transcript boundary, there is no
            // intron on that side, so don't check for splice region there.
            const check_region_beg = tr.beg != exon.beg;
            const check_region_end = tr.end != exon.end;

            // For each alt allele
            var ial: u32 = 1;
            while (ial < rec.n_allele) : (ial += 1) {
                if (ial >= rec.alleles.len) break;
                const alt = rec.alleles[ial];
                if (alt.len > 0 and (alt[0] == '<' or alt[0] == '*')) continue;

                const ref_allele = rec.alleles[0];

                var splice = Splice.init(self.allocator, tr);
                defer splice.deinit();

                splice.reset(
                    @intCast(rec.pos),
                    @intCast(ref_allele.len),
                    @intCast(ial),
                    ref_allele,
                    alt,
                );

                // Enable donor/acceptor checking (the key difference from testUtr)
                splice.flags.check_donor = true;
                splice.flags.check_acceptor = true;
                splice.flags.check_region_beg = check_region_beg;
                splice.flags.check_region_end = check_region_end;

                _ = splice.spliceCsq(exon.beg, exon.end);

                // If any splice consequence was set, stage it
                if (splice.csq.toInt() != 0) {
                    var csq = Csq{
                        .pos = rec.pos,
                        .vcsq = .{
                            .csq_type = splice.csq.toInt(),
                            .biotype = @intFromEnum(tr.biotype),
                            .strand = if (tr.strand == .forward) .fwd else .rev,
                            .trid = tr.id,
                            .vcf_ial = ial,
                            .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                        },
                    };
                    try self.csqStage(&csq, rec);
                    ret = true;
                }
            }
        }
        return ret;
    }

    // -----------------------------------------------------------------
    // testSymbolicAlt — handle <INS:*> and <DEL> symbolic alleles
    // -----------------------------------------------------------------

    /// Handle symbolic ALT alleles like <INS:ME:ALU> and <DEL>.
    ///
    /// Port of test_symbolic_alt() from csq.c (line 3472).
    /// Checks CDS, UTR, exon (splice), and transcript indices.
    fn testSymbolicAlt(self: *CsqContext, rec: *const VcfRecord) !void {
        const gff = self.gff orelse return;
        const chr = rec.seqname();

        if (rec.alleles.len < 2) return;
        const alt = rec.alleles[1];

        // Determine elongation or truncation
        var csq_class: CsqType = 0;
        if (alt.len >= 4 and std.ascii.eqlIgnoreCase(alt[0..4], "<INS")) {
            csq_class = CSQ_ELONGATION;
        } else if (alt.len >= 4 and std.ascii.eqlIgnoreCase(alt[0..4], "<DEL")) {
            csq_class = CSQ_TRUNCATION;
        } else return;

        // Symbolic ALTs use pos+1 as the query position (C: beg = rec->pos + 1)
        const beg = rec.pos + 1;
        const end = beg;

        var hit = false;

        // Check CDS index
        {
            var itr = gff.idx_cds.overlap(chr, beg, end);
            while (itr.next()) |interval| {
                const cds: *CdsEntry = interval.payload;
                const tr: *Transcript = cds.tr;
                const coding_csq: CsqType = if (tr.biotype.isCoding()) CSQ_CODING_SEQUENCE else CSQ_NON_CODING;
                var csq = Csq{
                    .pos = rec.pos,
                    .vcsq = .{
                        .csq_type = coding_csq | csq_class,
                        .biotype = @intFromEnum(tr.biotype),
                        .strand = if (tr.strand == .forward) .fwd else .rev,
                        .trid = tr.id,
                        .vcf_ial = 1,
                        .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                    },
                };
                try self.csqStage(&csq, rec);
                hit = true;
            }
        }

        // Check UTR index
        {
            var itr = gff.idx_utr.overlap(chr, beg, end);
            while (itr.next()) |interval| {
                const utr: *Utr = interval.payload;
                const tr: *Transcript = utr.tr;
                const utr_csq: CsqType = if (utr.which == .prime5) CSQ_UTR5 else CSQ_UTR3;
                var csq = Csq{
                    .pos = rec.pos,
                    .vcsq = .{
                        .csq_type = utr_csq | csq_class,
                        .biotype = @intFromEnum(tr.biotype),
                        .strand = if (tr.strand == .forward) .fwd else .rev,
                        .trid = tr.id,
                        .vcf_ial = 1,
                        .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                    },
                };
                try self.csqStage(&csq, rec);
                hit = true;
            }
        }

        // Check exon index for splice consequences
        {
            var itr = gff.idx_exon.overlap(chr, beg, end);
            while (itr.next()) |interval| {
                const exon: *Exon = interval.payload;
                const tr: *Transcript = exon.tr;
                if (tr.cds.items.len == 0) continue;

                const check_region_beg = tr.beg != exon.beg;
                const check_region_end = tr.end != exon.end;

                var splice = Splice.init(self.allocator, tr);
                defer splice.deinit();

                const ref_allele = rec.alleles[0];
                splice.reset(
                    @intCast(rec.pos),
                    @intCast(ref_allele.len),
                    1,
                    ref_allele,
                    alt,
                );

                splice.flags.check_donor = true;
                splice.flags.check_acceptor = true;
                splice.flags.check_region_beg = check_region_beg;
                splice.flags.check_region_end = check_region_end;
                // Pre-set csq to csq_class so splice adds to it (C: splice.csq = csq_class)
                splice.csq = types.CsqType.fromInt(csq_class);

                _ = splice.spliceCsq(exon.beg, exon.end);

                if (splice.csq.toInt() != 0) {
                    var csq = Csq{
                        .pos = rec.pos,
                        .vcsq = .{
                            .csq_type = splice.csq.toInt(),
                            .biotype = @intFromEnum(tr.biotype),
                            .strand = if (tr.strand == .forward) .fwd else .rev,
                            .trid = tr.id,
                            .vcf_ial = 1,
                            .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                        },
                    };
                    try self.csqStage(&csq, rec);
                    hit = true;
                }
            }
        }

        // Check transcript index if nothing else hit
        if (!hit) {
            var itr = gff.idx_tscript.overlap(chr, beg, end);
            while (itr.next()) |interval| {
                const tr: *Transcript = interval.payload;
                const csq_type: CsqType = if (tr.biotype.isCoding()) CSQ_INTRON else CSQ_NON_CODING;
                var csq = Csq{
                    .pos = rec.pos,
                    .vcsq = .{
                        .csq_type = csq_type | csq_class,
                        .biotype = @intFromEnum(tr.biotype),
                        .strand = if (tr.strand == .forward) .fwd else .rev,
                        .trid = tr.id,
                        .vcf_ial = 1,
                        .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                    },
                };
                try self.csqStage(&csq, rec);
            }
        }
    }

    // -----------------------------------------------------------------
    // testTscript — check for intronic / non-coding variants
    // -----------------------------------------------------------------

    /// Check if the variant falls within a transcript (intron or non-coding).
    ///
    /// Port of test_tscript() from csq.c (line 3434).
    /// Queries the transcript region index. For each overlapping transcript
    /// and each alt allele, runs splice analysis against the full transcript
    /// span. If the variant is inside/overlapping, it is classified as INTRON
    /// (for coding transcripts) or NON_CODING (for non-coding transcripts).
    fn testTscript(self: *CsqContext, rec: *const VcfRecord) !bool {
        const gff = self.gff orelse return false;
        const chr = rec.seqname();

        var itr = gff.idx_tscript.overlap(chr, rec.pos, rec.pos + rec.rlen);

        var ret = false;
        while (itr.next()) |interval| {
            const tr: *Transcript = interval.payload;

            // For each alt allele
            var ial: u32 = 1;
            while (ial < rec.n_allele) : (ial += 1) {
                if (ial >= rec.alleles.len) break;
                const alt = rec.alleles[ial];
                if (alt.len > 0 and (alt[0] == '<' or alt[0] == '*')) continue;

                const ref_allele = rec.alleles[0];

                var splice = Splice.init(self.allocator, tr);
                defer splice.deinit();

                splice.reset(
                    @intCast(rec.pos),
                    @intCast(ref_allele.len),
                    @intCast(ial),
                    ref_allele,
                    alt,
                );

                const splice_ret = splice.spliceCsq(tr.beg, tr.end);
                if (splice_ret != .inside and splice_ret != .overlap) continue;

                // Coding transcript -> INTRON; non-coding -> NON_CODING
                const csq_type: CsqType = if (tr.biotype.isCoding()) CSQ_INTRON else CSQ_NON_CODING;

                var csq = Csq{
                    .pos = rec.pos,
                    .vcsq = .{
                        .csq_type = csq_type,
                        .biotype = @intFromEnum(tr.biotype),
                        .strand = if (tr.strand == .forward) .fwd else .rev,
                        .trid = tr.id,
                        .vcf_ial = ial,
                        .gene = if (tr.gene) |g| @as(?[]const u8, if (g.name) |n| std.mem.span(n) else null) else null,
                    },
                };
                try self.csqStage(&csq, rec);
                ret = true;
            }
        }
        return ret;
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

        // Find the vrec matching this record.
        // Match by position + allele count since the stored record is a
        // heap-allocated copy (not the same pointer as the caller's rec).
        var vrec_idx: ?usize = null;
        for (vb.vrecs.items, 0..) |*vrec, idx| {
            if (vrec.rec) |stored_rec| {
                if (stored_rec.pos == rec.pos and stored_rec.n_allele == rec.n_allele and stored_rec.rid == rec.rid) {
                    vrec_idx = idx;
                    break;
                }
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

        // Remove stop_lost & synonymous if stop_retained set (C line 2059-2060)
        if (t & CSQ_STOP_RETAINED != 0) {
            t &= ~(CSQ_STOP_LOST | CSQ_SYNONYMOUS_VARIANT);
            csq.vcsq.csq_type = t;
        }
        // Remove start_lost & synonymous if start_retained set (C line 2062-2063)
        if (t & CSQ_START_RETAINED != 0) {
            t &= ~(CSQ_START_LOST | CSQ_SYNONYMOUS_VARIANT);
            csq.vcsq.csq_type = t;
        }

        // Deduplication: scan existing consequences for a match
        for (vrec.vcsqs.items, 0..) |*existing, idx| {
            if (isDuplicate(existing, &csq.vcsq)) {
                // Merge type bits into existing consequence
                existing.csq_type |= t;

                // Remove stop_lost & synonymous if stop_retained set (C line 2059-2060)
                if (existing.csq_type & CSQ_STOP_RETAINED != 0)
                    existing.csq_type &= ~(CSQ_STOP_LOST | CSQ_SYNONYMOUS_VARIANT);

                // Remove start_lost & synonymous if start_retained set (C line 2062-2063)
                if (existing.csq_type & CSQ_START_RETAINED != 0)
                    existing.csq_type &= ~(CSQ_START_LOST | CSQ_SYNONYMOUS_VARIANT);

                // Copy vstr from new to existing only for compound merges (C line 2065).
                // Non-compound merges (C line 2082-2085) do NOT copy vstr.
                if (t & CSQ_COMPOUND != 0) {
                    if (existing.vstr == null and csq.vcsq.vstr != null) {
                        existing.vstr = csq.vcsq.vstr;
                        csq.vcsq.vstr = null; // transfer ownership
                    }
                }

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

        // For compound consequences, also check gene, vcf_ial, upstream_stop, and vstr
        if (new.csq_type & CSQ_COMPOUND != 0) {
            if (!strEql(existing.gene, new.gene)) return false;
            if (existing.vcf_ial != new.vcf_ial) return false;
            // Both must or mustn't have upstream_stop (C line 2043)
            if ((existing.csq_type & CSQ_UPSTREAM_STOP) ^ (new.csq_type & CSQ_UPSTREAM_STOP) != 0) return false;

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

        if (self.phase == .drop_gt) {
            // No genotype handling needed in drop_gt mode.
            return;
        }

        // Genotype-aware sample assignment.
        // Port of csq_stage() from csq.c lines 3296-3358.

        // Parse (or reuse cached) genotypes
        var ngt: u32 = 0;
        var gts: ?[]Genotype = null;

        if (self.gt_cache_rec == rec and self.gt_cache != null) {
            gts = self.gt_cache;
        } else {
            // Free previously cached genotypes
            if (self.gt_cache) |gc| self.allocator.free(gc);
            self.gt_cache = null;
            self.gt_cache_rec = null;

            gts = rec.parseGenotypes(self.allocator) catch null;
            if (gts != null) {
                self.gt_cache = gts;
                self.gt_cache_rec = rec;
            }
        }

        if (gts == null or gts.?.len == 0) {
            // No genotypes: output with no sample (tab text mode would print here)
            return;
        }
        const genotypes = gts.?;

        // Use max ploidy across all samples (matches testCds)
        ngt = 0;
        for (genotypes) |g| {
            if (g.ploidy > ngt) ngt = g.ploidy;
        }
        if (ngt == 0 or ngt > 2) return;

        // VCF output: set bits in vrec.fmt_bm for matching samples
        const n_smpl = if (self.sample_indices) |si| @as(u32, @intCast(si.len)) else @as(u32, @intCast(genotypes.len));
        const csq_idx: u32 = if (csq.csq_idx) |ci| @intCast(ci) else return;

        // Find the vrec for this consequence
        const vrec_idx = csq.vrec_idx orelse return;
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
        const vb = vbuf orelse return;
        if (vrec_idx >= vb.vrecs.items.len) return;
        var vrec = &vb.vrecs.items[vrec_idx];

        // Ensure fmt_bm is allocated: n_smpl * nfmt_bcsq u32s
        if (vrec.fmt_bm == null) {
            const bm_size = n_smpl * self.nfmt_bcsq;
            vrec.fmt_bm = try self.allocator.alloc(u32, bm_size);
            @memset(vrec.fmt_bm.?, 0);
        }

        for (0..n_smpl) |ismpl_idx| {
            const ismpl: usize = if (self.sample_indices) |si| @as(usize, si[ismpl_idx]) else ismpl_idx;
            if (ismpl >= genotypes.len) continue;

            const gt = genotypes[ismpl];

            var j: u32 = 0;
            while (j < ngt) : (j += 1) {
                if (gt.alleles[j] <= 0) continue; // missing or ref
                const ial: u32 = @intCast(gt.alleles[j]);
                if (ial != csq.vcsq.vcf_ial) continue;

                // icsq2 = 2 * csq_idx + haplotype (interleave first/second haplotype)
                const icsq2: u32 = 2 * csq_idx + j;
                if (icsq2 >= self.ncsq2_max) {
                    if (self.verbosity > 0) {
                        std.log.warn("Too many consequences at pos {d}, keeping first {d}", .{ csq.pos + 1, icsq2 + 1 });
                    }
                    break;
                }

                // Convert icsq2 to (ival, ibit) pair
                const ival: u32 = icsq2 / 31;
                const ibit: u5 = @intCast(icsq2 % 31);

                if (vrec.nfmt < 1 + ival) vrec.nfmt = 1 + ival;

                const bm_idx = @as(usize, ismpl_idx) * self.nfmt_bcsq + ival;
                if (bm_idx < vrec.fmt_bm.?.len) {
                    vrec.fmt_bm.?[bm_idx] |= @as(u32, 1) << ibit;
                }
            }
        }
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
        UnphasedHeterozygous,
    };
};

// ---------------------------------------------------------------------------
// BCSQ injection into text VCF lines
// ---------------------------------------------------------------------------

/// Inject a BCSQ annotation into a text VCF line by appending `tag=value`
/// to the INFO column (column index 7, 0-based).
///
/// If the INFO field is "." it is replaced entirely; otherwise the tag is
/// appended with a semicolon separator.
///
/// Returns a newly-allocated line (without trailing newline).
pub fn injectBcsq(
    allocator: std.mem.Allocator,
    original_line: []const u8,
    bcsq_value: []const u8,
    bcsq_tag: []const u8,
) ![]u8 {
    // Strip trailing newline/CR
    var line = original_line;
    if (line.len > 0 and line[line.len - 1] == '\n') line = line[0 .. line.len - 1];
    if (line.len > 0 and line[line.len - 1] == '\r') line = line[0 .. line.len - 1];

    // Find tab-delimited column boundaries.
    // We need to locate column 7 (INFO).
    var col_starts: [9]usize = undefined; // cols 0..8
    var col_ends: [9]usize = undefined;
    var col: usize = 0;
    var start: usize = 0;
    for (line, 0..) |c, i| {
        if (c == '\t') {
            if (col < 9) {
                col_starts[col] = start;
                col_ends[col] = i;
            }
            col += 1;
            start = i + 1;
            if (col >= 9) break;
        }
    }
    // Handle last/remaining field
    if (col < 9) {
        col_starts[col] = start;
        col_ends[col] = line.len;
        col += 1;
    }

    if (col < 8) return error.TooFewColumns;

    const info_start = col_starts[7];
    const info_end = col_ends[7];
    const info_field = line[info_start..info_end];

    // Calculate result size
    const is_dot = std.mem.eql(u8, info_field, ".");
    const inject_len = bcsq_tag.len + 1 + bcsq_value.len; // "TAG=VALUE"
    const separator: usize = if (is_dot) 0 else 1; // ";" before tag

    const new_info_len = if (is_dot)
        inject_len
    else
        info_field.len + separator + inject_len;

    const result_len = line.len - info_field.len + new_info_len;
    const result = try allocator.alloc(u8, result_len);

    // Copy: [before INFO] [new INFO] [after INFO]
    var pos: usize = 0;
    // Everything up to (but not including) INFO content
    @memcpy(result[pos .. pos + info_start], line[0..info_start]);
    pos += info_start;

    if (is_dot) {
        // Replace "." with "TAG=VALUE"
        @memcpy(result[pos .. pos + bcsq_tag.len], bcsq_tag);
        pos += bcsq_tag.len;
        result[pos] = '=';
        pos += 1;
        @memcpy(result[pos .. pos + bcsq_value.len], bcsq_value);
        pos += bcsq_value.len;
    } else {
        // Keep existing INFO, append ";TAG=VALUE"
        @memcpy(result[pos .. pos + info_field.len], info_field);
        pos += info_field.len;
        result[pos] = ';';
        pos += 1;
        @memcpy(result[pos .. pos + bcsq_tag.len], bcsq_tag);
        pos += bcsq_tag.len;
        result[pos] = '=';
        pos += 1;
        @memcpy(result[pos .. pos + bcsq_value.len], bcsq_value);
        pos += bcsq_value.len;
    }

    // Copy everything after INFO (including trailing columns)
    const after = line[info_end..];
    @memcpy(result[pos .. pos + after.len], after);
    pos += after.len;

    std.debug.assert(pos == result_len);
    return result;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "vbufPush: two records at same position share vbuf" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
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

    const r1 = try ctx.vbufPush(&rec1);
    const r2 = try ctx.vbufPush(&rec2);

    // Both should return the same vbuf
    try std.testing.expectEqual(r1.vbuf, r2.vbuf);
    // The vbuf should contain 2 records
    try std.testing.expectEqual(@as(usize, 2), r1.vbuf.vrecs.items.len);
    // Ring buffer should have exactly 1 entry
    try std.testing.expectEqual(@as(usize, 1), ctx.vcf_rbuf.len);
}

test "vbufPush: records at different positions get separate vbufs" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
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

    const r1 = try ctx.vbufPush(&rec1);
    const r2 = try ctx.vbufPush(&rec2);

    // Should be different vbufs
    try std.testing.expect(r1.vbuf != r2.vbuf);
    // Each should have 1 record
    try std.testing.expectEqual(@as(usize, 1), r1.vbuf.vrecs.items.len);
    try std.testing.expectEqual(@as(usize, 1), r2.vbuf.vrecs.items.len);
    // Ring buffer should have 2 entries
    try std.testing.expectEqual(@as(usize, 2), ctx.vcf_rbuf.len);
}

test "csqPush: dedup pushes same consequence only once" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
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
        .gff_fname = "",
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

// ---------------------------------------------------------------------------
// Integration tests — GFF + CSQ pipeline wiring
// ---------------------------------------------------------------------------

/// Helper: create a minimal GffParser with a single coding transcript on chr1
/// spanning [100, 900], a CDS at [200, 400], an exon at [200, 400], a UTR5
/// at [100, 199], and the transcript itself at [100, 900].
const TestGff = struct {
    gff: *GffParser,
    tr: *Transcript,
};

fn makeTestGff(allocator: std.mem.Allocator) !TestGff {
    var gff = try allocator.create(GffParser);
    gff.* = GffParser.init(allocator);
    const arena = gff.arena.allocator();

    // Create gene
    const gene = try arena.create(gff_types.Gene);
    gene.* = .{
        .name = null,
        .iseq = 0,
        .id = 0,
        .beg = 100,
        .end = 900,
        .strand = .forward,
        .used = true,
    };

    // Create transcript
    const tr = try arena.create(Transcript);
    tr.* = Transcript.init(arena);
    tr.id = 0;
    tr.beg = 100;
    tr.end = 900;
    tr.strand = .forward;
    tr.biotype = .protein_coding;
    tr.gene = gene;

    // Create CDS entry
    const cds = try arena.create(CdsEntry);
    cds.* = .{
        .tr = tr,
        .beg = 200,
        .pos = 0,
        .len = 201,
        .icds = 0,
        .phase = .phase0,
    };
    try tr.cds.append(arena, cds);

    // Insert into region indexes
    try gff.idx_cds.insert("chr1", 200, 400, cds);

    const utr = try arena.create(Utr);
    utr.* = .{
        .which = .prime5,
        .beg = 100,
        .end = 199,
        .tr = tr,
    };
    try gff.idx_utr.insert("chr1", 100, 199, utr);

    const exon = try arena.create(Exon);
    exon.* = .{
        .beg = 200,
        .end = 400,
        .tr = tr,
    };
    try gff.idx_exon.insert("chr1", 200, 400, exon);

    try gff.idx_tscript.insert("chr1", 100, 900, tr);

    return .{ .gff = gff, .tr = tr };
}

fn destroyTestGff(allocator: std.mem.Allocator, tgff: TestGff) void {
    cleanupTranscriptAux(allocator, tgff.tr);
    tgff.gff.deinit();
    allocator.destroy(tgff.gff);
}

/// Clean up transcript aux data (Tscript + HapNode root + hap array)
/// that was allocated by testCds/testCdsLocal during testing.
fn cleanupTranscriptAux(allocator: std.mem.Allocator, tr: *Transcript) void {
    if (tr.aux) |aux_raw| {
        const taux: *Tscript = @ptrCast(@alignCast(aux_raw));
        if (taux.root) |root| {
            // Free children recursively (shallow: only direct children in test scenarios)
            for (root.children.items) |child_node| {
                if (child_node.payload == .cds) {
                    if (child_node.payload.cds.seq) |seq| allocator.free(seq);
                }
                if (child_node.var_str) |vs| allocator.free(vs);
                child_node.deinit(allocator);
                allocator.destroy(child_node);
            }
            root.deinit(allocator);
            allocator.destroy(root);
        }
        taux.hap.deinit(allocator);
        allocator.destroy(taux);
        tr.aux = null;
    }
}

test "testCds: variant overlapping CDS is detected" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGff(allocator);
    defer destroyTestGff(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
        .force = true,
        .verbosity = 0,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    const alleles = [_][]const u8{ "A", "T" };
    var rec = VcfRecord{
        .pos = 250,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    const vbuf = (try ctx.vbufPush(&rec)).vbuf;
    const hit = try ctx.testCds(&rec, vbuf);

    try std.testing.expect(hit);
    // keep_until should be extended to transcript end
    try std.testing.expectEqual(@as(u32, 900), vbuf.keep_until);
}

test "testCds: variant outside CDS is not detected" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGff(allocator);
    defer destroyTestGff(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    const alleles = [_][]const u8{ "A", "T" };
    // Position 500 is outside the CDS [200, 400]
    var rec = VcfRecord{
        .pos = 500,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    const vbuf = (try ctx.vbufPush(&rec)).vbuf;
    const hit = try ctx.testCds(&rec, vbuf);

    try std.testing.expect(!hit);
}

test "testCdsLocal: variant overlapping CDS is detected (no fasta)" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGff(allocator);
    defer destroyTestGff(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
        .local_csq = true,
        .force = true,
        .verbosity = 0,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    const alleles = [_][]const u8{ "A", "T" };
    var rec = VcfRecord{
        .pos = 300,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    _ = try ctx.vbufPush(&rec);
    const hit = try ctx.testCdsLocal(&rec);

    // CDS overlap is detected even without fasta (ret=true set before init)
    try std.testing.expect(hit);

    // Without fasta, transcript aux cannot be initialized, so no consequence
    // is staged. The CDS hit still prevents fallthrough to UTR/intron.
}

test "testUtr: variant in UTR5 region is detected" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGff(allocator);
    defer destroyTestGff(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    const alleles = [_][]const u8{ "A", "T" };
    // Position 150 is inside UTR5 [100, 199]
    var rec = VcfRecord{
        .pos = 150,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    _ = try ctx.vbufPush(&rec);
    const hit = try ctx.testUtr(&rec);

    try std.testing.expect(hit);

    // Check that UTR5 consequence was staged
    const vbuf = ctx.vcf_rbuf.front().?;
    const vrec = &vbuf.vrecs.items[0];
    try std.testing.expect(vrec.vcsqs.items.len > 0);
    try std.testing.expect(vrec.vcsqs.items[0].csq_type & CSQ_UTR5 != 0);
}

test "testTscript: intronic variant in coding transcript gets INTRON" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGff(allocator);
    defer destroyTestGff(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    const alleles = [_][]const u8{ "A", "T" };
    // Position 500 is inside transcript [100, 900] but outside CDS [200, 400]
    var rec = VcfRecord{
        .pos = 500,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    _ = try ctx.vbufPush(&rec);
    const hit = try ctx.testTscript(&rec);

    try std.testing.expect(hit);

    // Check that INTRON consequence was staged (protein_coding is coding)
    const vbuf = ctx.vcf_rbuf.front().?;
    const vrec = &vbuf.vrecs.items[0];
    try std.testing.expect(vrec.vcsqs.items.len > 0);
    try std.testing.expect(vrec.vcsqs.items[0].csq_type & CSQ_INTRON != 0);
}

test "testTscript: variant in non-coding transcript gets NON_CODING" {
    const allocator = std.testing.allocator;

    // Build a custom GFF with a non-coding transcript
    var gff = try allocator.create(GffParser);
    defer {
        gff.deinit();
        allocator.destroy(gff);
    }
    gff.* = GffParser.init(allocator);
    const arena = gff.arena.allocator();

    const gene = try arena.create(gff_types.Gene);
    gene.* = .{
        .name = null,
        .iseq = 0,
        .id = 0,
        .beg = 100,
        .end = 900,
        .strand = .forward,
        .used = true,
    };

    const tr = try arena.create(Transcript);
    tr.* = Transcript.init(arena);
    tr.id = 0;
    tr.beg = 100;
    tr.end = 900;
    tr.strand = .forward;
    tr.biotype = .lncRNA; // non-coding
    tr.gene = gene;

    try gff.idx_tscript.insert("chr1", 100, 900, tr);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
    });
    defer ctx.deinit();
    ctx.gff = gff;

    const alleles = [_][]const u8{ "A", "T" };
    var rec = VcfRecord{
        .pos = 500,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    _ = try ctx.vbufPush(&rec);
    const hit = try ctx.testTscript(&rec);

    try std.testing.expect(hit);

    const vbuf = ctx.vcf_rbuf.front().?;
    const vrec = &vbuf.vrecs.items[0];
    try std.testing.expect(vrec.vcsqs.items.len > 0);
    try std.testing.expect(vrec.vcsqs.items[0].csq_type & CSQ_NON_CODING != 0);
}

test "testSplice: variant near exon boundary sets splice consequence" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGff(allocator);
    defer destroyTestGff(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    const alleles = [_][]const u8{ "A", "T" };
    // Position 399 is within the last 3bp of exon [200, 400], which triggers
    // splice_region via checkExonEnd. The exon end != transcript end, so
    // check_region_end is set.
    var rec = VcfRecord{
        .pos = 399,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    _ = try ctx.vbufPush(&rec);
    const hit = try ctx.testSplice(&rec);

    try std.testing.expect(hit);

    // Check that a splice consequence was staged
    const vbuf = ctx.vcf_rbuf.front().?;
    const vrec = &vbuf.vrecs.items[0];
    try std.testing.expect(vrec.vcsqs.items.len > 0);
    // Should have splice_region set (within last 3bp of exon, check_region_end enabled)
    try std.testing.expect(vrec.vcsqs.items[0].csq_type & CSQ_SPLICE_REGION != 0);
}

test "testSplice: variant far from exon boundary has no splice consequence" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGff(allocator);
    defer destroyTestGff(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    const alleles = [_][]const u8{ "A", "T" };
    // Position 500 is well outside the exon [200, 400] and beyond the
    // splice region window (8bp)
    var rec = VcfRecord{
        .pos = 500,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    _ = try ctx.vbufPush(&rec);
    const hit = try ctx.testSplice(&rec);

    // No exon overlaps position 500, so no splice hit
    try std.testing.expect(!hit);
}

test "process: full cascade CDS -> UTR -> splice -> tscript" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGff(allocator);
    defer destroyTestGff(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
        .local_csq = true,
        .force = true,
        .verbosity = 0,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    // Variant in CDS region
    const alleles_cds = [_][]const u8{ "A", "T" };
    var rec_cds = VcfRecord{
        .pos = 300,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles_cds,
        .rlen = 1,
        .chr = "chr1",
    };
    try ctx.process(&rec_cds);

    // Variant in intron (outside CDS, inside transcript)
    const alleles_intr = [_][]const u8{ "G", "C" };
    var rec_intron = VcfRecord{
        .pos = 500,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles_intr,
        .rlen = 1,
        .chr = "chr1",
    };
    try ctx.process(&rec_intron);

    // Flush everything
    try ctx.flush();

    // Verify that both records were processed (ring buffer is empty after flush)
    try std.testing.expectEqual(@as(usize, 0), ctx.vcf_rbuf.len);
}

/// Helper: create a test GFF with pre-populated reference sequences for
/// CDS translation tests. The transcript spans chr1:[100,900], forward strand,
/// with a single CDS at [200,400] (len=201). The reference sequence starts
/// with ATG (Met) at position 200-202, followed by GAA (Glu) at 203-205, etc.
fn makeTestGffWithRef(allocator: std.mem.Allocator) !TestGff {
    const tgff = try makeTestGff(allocator);
    const tr = tgff.tr;

    // Build a mock reference sequence covering [tr.beg-10, tr.end+10] = [90, 910]
    // ref_seq length = (tr.end - tr.beg + 1) + 2*N_REF_PAD = 801 + 20 = 821
    const ref_len: usize = @as(usize, tr.end - tr.beg + 1) + 2 * N_REF_PAD;
    const ref_seq = try allocator.alloc(u8, ref_len);
    @memset(ref_seq, 'A'); // default all A

    // Set up the CDS region: starts at ref_seq[N_REF_PAD + (cds.beg - tr.beg)]
    // = ref_seq[10 + 100] = ref_seq[110]
    // CDS positions [200, 400] => ref_seq[110..311]
    //
    // First codon (pos 200-202): ATG = Met (start codon)
    ref_seq[110] = 'A';
    ref_seq[111] = 'T';
    ref_seq[112] = 'G';
    // Second codon (pos 203-205): GAA = Glu
    ref_seq[113] = 'G';
    ref_seq[114] = 'A';
    ref_seq[115] = 'A';
    // Third codon (pos 206-208): TGC = Cys
    ref_seq[116] = 'T';
    ref_seq[117] = 'G';
    ref_seq[118] = 'C';
    // Fill rest of CDS with AAA (Lys) codons
    var i: usize = 119;
    while (i < 311) : (i += 3) {
        ref_seq[i] = 'A';
        if (i + 1 < 311) ref_seq[i + 1] = 'A';
        if (i + 2 < 311) ref_seq[i + 2] = 'A';
    }

    // Build the spliced reference: N_REF_PAD + CDS + N_REF_PAD
    const cds_len: usize = 201;
    const sref_len: usize = 2 * N_REF_PAD + cds_len;
    const sref = try allocator.alloc(u8, sref_len);
    // Left padding: from ref_seq at cds_offset - N_REF_PAD = 100
    @memcpy(sref[0..N_REF_PAD], ref_seq[100 .. 100 + N_REF_PAD]);
    // CDS region
    @memcpy(sref[N_REF_PAD .. N_REF_PAD + cds_len], ref_seq[110 .. 110 + cds_len]);
    // Right padding: from ref_seq after CDS
    @memcpy(sref[N_REF_PAD + cds_len .. sref_len], ref_seq[110 + cds_len .. 110 + cds_len + N_REF_PAD]);

    // Create the Tscript (aux data) with the reference
    const taux = try allocator.create(Tscript);
    taux.* = .{
        .ref_seq = ref_seq,
        .sref = sref,
        .nsref = @intCast(sref_len),
    };

    // Create haplotype tree root node (needed for testCds)
    const root = try allocator.create(HapNode);
    root.* = HapNode.init(.root);
    root.nend = 1; // single haplotype for DROP_GT
    taux.root = root;

    tr.aux = taux;

    return tgff;
}

fn destroyTestGffWithRef(allocator: std.mem.Allocator, tgff: TestGff) void {
    // Clean up the ref_seq and sref we allocated
    if (tgff.tr.aux) |aux_raw| {
        const taux: *Tscript = @ptrCast(@alignCast(aux_raw));
        if (taux.ref_seq) |rs| allocator.free(rs);
        if (taux.sref) |sr| allocator.free(sr);
        // Clean up any children created by hapInit
        if (taux.root) |root| {
            for (root.children.items) |child_node| {
                if (child_node.payload == .cds) {
                    if (child_node.payload.cds.seq) |seq| allocator.free(seq);
                }
                if (child_node.var_str) |vs| allocator.free(vs);
                child_node.deinit(allocator);
                allocator.destroy(child_node);
            }
            root.deinit(allocator);
            allocator.destroy(root);
        }
        taux.hap.deinit(allocator);
        allocator.destroy(taux);
        tgff.tr.aux = null;
    }
    tgff.gff.deinit();
    allocator.destroy(tgff.gff);
}

test "testCdsLocal with ref: missense consequence for T>A at second codon position" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGffWithRef(allocator);
    defer destroyTestGffWithRef(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
        .local_csq = true,
        .force = true,
        .verbosity = 0,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    // Variant at position 201 (second base of first codon ATG):
    // REF=T, ALT=A => codon changes ATG(Met) -> AAG(Lys) = missense
    const alleles = [_][]const u8{ "T", "A" };
    var rec = VcfRecord{
        .pos = 201,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    _ = try ctx.vbufPush(&rec);
    const hit = try ctx.testCdsLocal(&rec);

    try std.testing.expect(hit);

    // Check that a consequence was staged
    const vbuf = ctx.vcf_rbuf.front().?;
    const vrec = &vbuf.vrecs.items[0];
    try std.testing.expect(vrec.vcsqs.items.len > 0);

    // The consequence should contain missense_variant (bit 2)
    const csq_type = vrec.vcsqs.items[0].csq_type;
    try std.testing.expect(csq_type & CSQ_MISSENSE_VARIANT != 0);
}

test "testCds DROP_GT: haplotype tree node created for CDS variant" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGffWithRef(allocator);
    defer destroyTestGffWithRef(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
        .force = true,
        .verbosity = 0,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    // Variant at position 204 (second base of second codon GAA):
    // REF=A, ALT=T => codon changes GAA(Glu) -> GTA(Val) = missense
    const alleles = [_][]const u8{ "A", "T" };
    var rec = VcfRecord{
        .pos = 204,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
    };

    const vbuf = (try ctx.vbufPush(&rec)).vbuf;
    const hit = try ctx.testCds(&rec, vbuf);

    try std.testing.expect(hit);
    try std.testing.expectEqual(@as(u32, 900), vbuf.keep_until);

    // The transcript should have aux data with a haplotype tree
    const taux: *Tscript = @ptrCast(@alignCast(tgff.tr.aux orelse unreachable));
    try std.testing.expect(taux.root != null);

    // The root should have a child (the variant node)
    const root = taux.root.?;
    try std.testing.expect(root.children.items.len > 0);

    // The child should be a CDS node with the variant applied
    const child = root.children.items[0];
    try std.testing.expect(child.payload == .cds);
    try std.testing.expectEqual(@as(u32, 204), child.rbeg);

    // The haplotype leaf should point to the child
    try std.testing.expect(taux.hap.items.len > 0);
    try std.testing.expectEqual(child, taux.hap.items[0]);
    try std.testing.expectEqual(@as(u32, 1), child.nend);
}

// ---------------------------------------------------------------------------
// injectBcsq tests
// ---------------------------------------------------------------------------

test "injectBcsq: replace dot INFO with BCSQ" {
    const allocator = std.testing.allocator;
    const line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1";
    const result = try injectBcsq(allocator, line, "missense|GENE|TR|protein_coding|+", "BCSQ");
    defer allocator.free(result);
    try std.testing.expectEqualStrings(
        "chr1\t100\t.\tA\tT\t.\tPASS\tBCSQ=missense|GENE|TR|protein_coding|+\tGT\t0/1",
        result,
    );
}

test "injectBcsq: append to existing INFO" {
    const allocator = std.testing.allocator;
    const line = "chr1\t100\t.\tA\tT\t.\tPASS\tDP=30;AF=0.5\tGT\t0/1";
    const result = try injectBcsq(allocator, line, "intron|G||lncRNA", "BCSQ");
    defer allocator.free(result);
    try std.testing.expectEqualStrings(
        "chr1\t100\t.\tA\tT\t.\tPASS\tDP=30;AF=0.5;BCSQ=intron|G||lncRNA\tGT\t0/1",
        result,
    );
}

test "injectBcsq: custom tag name" {
    const allocator = std.testing.allocator;
    const line = "chr1\t100\t.\tA\tT\t.\tPASS\t.";
    const result = try injectBcsq(allocator, line, "synonymous|X||", "MY_CSQ");
    defer allocator.free(result);
    try std.testing.expectEqualStrings(
        "chr1\t100\t.\tA\tT\t.\tPASS\tMY_CSQ=synonymous|X||",
        result,
    );
}

test "injectBcsq: line with trailing newline" {
    const allocator = std.testing.allocator;
    const line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\n";
    const result = try injectBcsq(allocator, line, "stop_gained|G|T|pc", "BCSQ");
    defer allocator.free(result);
    try std.testing.expectEqualStrings(
        "chr1\t100\t.\tA\tT\t.\tPASS\tBCSQ=stop_gained|G|T|pc",
        result,
    );
}

test "injectBcsq: minimal line without FORMAT/samples" {
    const allocator = std.testing.allocator;
    const line = "chr1\t100\t.\tA\tT\t.\tPASS\t.";
    const result = try injectBcsq(allocator, line, "missense|G||pc|+", "BCSQ");
    defer allocator.free(result);
    try std.testing.expectEqualStrings(
        "chr1\t100\t.\tA\tT\t.\tPASS\tBCSQ=missense|G||pc|+",
        result,
    );
}

// ---------------------------------------------------------------------------
// Genotype parsing tests
// ---------------------------------------------------------------------------

test "Genotype.parse: 0/1 -> alleles=[0,1], phased=false, ploidy=2" {
    const gt = Genotype.parse("0/1");
    try std.testing.expectEqual(@as(i32, 0), gt.alleles[0]);
    try std.testing.expectEqual(@as(i32, 1), gt.alleles[1]);
    try std.testing.expect(!gt.phased);
    try std.testing.expectEqual(@as(u8, 2), gt.ploidy);
}

test "Genotype.parse: 1|0 -> alleles=[1,0], phased=true, ploidy=2" {
    const gt = Genotype.parse("1|0");
    try std.testing.expectEqual(@as(i32, 1), gt.alleles[0]);
    try std.testing.expectEqual(@as(i32, 0), gt.alleles[1]);
    try std.testing.expect(gt.phased);
    try std.testing.expectEqual(@as(u8, 2), gt.ploidy);
}

test "Genotype.parse: ./. -> missing alleles" {
    const gt = Genotype.parse("./.");
    try std.testing.expectEqual(@as(i32, -1), gt.alleles[0]);
    try std.testing.expectEqual(@as(i32, -1), gt.alleles[1]);
    try std.testing.expectEqual(@as(u8, 2), gt.ploidy);
}

test "Genotype.parse: 0 -> haploid" {
    const gt = Genotype.parse("0");
    try std.testing.expectEqual(@as(i32, 0), gt.alleles[0]);
    try std.testing.expectEqual(@as(i32, -1), gt.alleles[1]);
    try std.testing.expectEqual(@as(u8, 1), gt.ploidy);
}

test "Genotype.parse: 1 -> haploid alt" {
    const gt = Genotype.parse("1");
    try std.testing.expectEqual(@as(i32, 1), gt.alleles[0]);
    try std.testing.expectEqual(@as(i32, -1), gt.alleles[1]);
    try std.testing.expectEqual(@as(u8, 1), gt.ploidy);
}

test "Genotype.parse: . -> missing haploid" {
    const gt = Genotype.parse(".");
    try std.testing.expectEqual(@as(i32, -1), gt.alleles[0]);
    try std.testing.expectEqual(@as(u8, 1), gt.ploidy);
}

test "Genotype.parse: 0/2 -> multi-allelic" {
    const gt = Genotype.parse("0/2");
    try std.testing.expectEqual(@as(i32, 0), gt.alleles[0]);
    try std.testing.expectEqual(@as(i32, 2), gt.alleles[1]);
    try std.testing.expect(!gt.phased);
}

test "VcfRecord.parseGenotypes: two samples from raw VCF line" {
    const allocator = std.testing.allocator;
    const alleles = [_][]const u8{ "A", "T" };
    const rec = VcfRecord{
        .pos = 100,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
        .raw_line = "chr1\t101\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t1|0",
    };
    const gts = try rec.parseGenotypes(allocator);
    try std.testing.expect(gts != null);
    defer allocator.free(gts.?);

    try std.testing.expectEqual(@as(usize, 2), gts.?.len);

    // Sample 0: 0/1
    try std.testing.expectEqual(@as(i32, 0), gts.?[0].alleles[0]);
    try std.testing.expectEqual(@as(i32, 1), gts.?[0].alleles[1]);
    try std.testing.expect(!gts.?[0].phased);

    // Sample 1: 1|0
    try std.testing.expectEqual(@as(i32, 1), gts.?[1].alleles[0]);
    try std.testing.expectEqual(@as(i32, 0), gts.?[1].alleles[1]);
    try std.testing.expect(gts.?[1].phased);
}

test "VcfRecord.parseGenotypes: no raw_line returns null" {
    const allocator = std.testing.allocator;
    const alleles = [_][]const u8{ "A", "T" };
    const rec = VcfRecord{
        .pos = 100,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
    };
    const gts = try rec.parseGenotypes(allocator);
    try std.testing.expect(gts == null);
}

test "testCds genotype-aware: two samples heterozygous get tree nodes" {
    const allocator = std.testing.allocator;

    const tgff = try makeTestGffWithRef(allocator);
    defer destroyTestGffWithRef(allocator, tgff);

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .as_is, // allow unphased
        .force = true,
        .verbosity = 0,
        .n_samples = 2,
    });
    defer ctx.deinit();
    ctx.gff = tgff.gff;

    // Two-sample VCF line with GT: sample 0 is 0/1, sample 1 is 1|0
    const alleles = [_][]const u8{ "A", "T" };
    var rec = VcfRecord{
        .pos = 204,
        .rid = 0,
        .n_allele = 2,
        .alleles = &alleles,
        .rlen = 1,
        .chr = "chr1",
        .raw_line = "chr1\t205\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t1|0",
    };

    // Pre-set root nend to 4 (2 samples * 2 haplotypes) since we're not in drop_gt mode
    const taux: *Tscript = @ptrCast(@alignCast(tgff.tr.aux orelse unreachable));
    const root = taux.root.?;
    root.nend = 4;

    const vbuf = (try ctx.vbufPush(&rec)).vbuf;
    const hit = try ctx.testCds(&rec, vbuf);

    try std.testing.expect(hit);

    // The root should have at least one child node (both samples share same allele 1)
    try std.testing.expect(root.children.items.len > 0);

    // The child should be a CDS node
    const child = root.children.items[0];
    try std.testing.expect(child.payload == .cds);
    try std.testing.expectEqual(@as(u32, 204), child.rbeg);

    // At least some haplotypes should point to the child node
    // Sample 0 hap[1] (second allele=1) and sample 1 hap[2] (first allele=1)
    try std.testing.expect(taux.hap.items.len >= 4);

    // Both haplotypes carrying alt allele should point to the same child
    // (sharing optimization: second sample reuses the node from the first)
    var alt_hap_count: u32 = 0;
    for (taux.hap.items) |h| {
        if (h == child) alt_hap_count += 1;
    }
    // At least 2 haplotypes should point to the child (one from each sample)
    try std.testing.expect(alt_hap_count >= 2);
}

// ---------------------------------------------------------------------------
// hapFlush tests
// ---------------------------------------------------------------------------

test "hapFlush: transcript ending before pos is flushed from active set" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
        .local_csq = true,
    });
    defer ctx.deinit();

    // Create a transcript ending at position 500
    var tr = Transcript.init(allocator);
    defer tr.deinit();
    tr.id = 1;
    tr.beg = 100;
    tr.end = 500;
    tr.strand = .forward;
    tr.biotype = .protein_coding;

    // Add to active transcripts
    try ctx.active_transcripts.add(&tr);
    try std.testing.expectEqual(@as(usize, 1), ctx.active_transcripts.count());

    // Flush with pos=600 (past the transcript end)
    try ctx.hapFlush(600);

    // Transcript should have been removed from active set
    try std.testing.expectEqual(@as(usize, 0), ctx.active_transcripts.count());

    // And added to the removal list
    try std.testing.expectEqual(@as(usize, 1), ctx.rm_transcripts.items.len);

    // Clear rm_transcripts without destroying aux (there is none)
    ctx.rm_transcripts.clearRetainingCapacity();
}

test "hapFlush: transcript ending after pos is NOT flushed" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
        .local_csq = true,
    });
    defer ctx.deinit();

    // Create a transcript ending at position 500
    var tr = Transcript.init(allocator);
    defer tr.deinit();
    tr.id = 1;
    tr.beg = 100;
    tr.end = 500;
    tr.strand = .forward;
    tr.biotype = .protein_coding;

    try ctx.active_transcripts.add(&tr);

    // Flush with pos=300 (before the transcript end)
    try ctx.hapFlush(300);

    // Transcript should still be in the active set
    try std.testing.expectEqual(@as(usize, 1), ctx.active_transcripts.count());
    // Nothing in removal list
    try std.testing.expectEqual(@as(usize, 0), ctx.rm_transcripts.items.len);
}

test "hapFlush: two transcripts flush in order of end position" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
        .local_csq = true,
    });
    defer ctx.deinit();

    // Transcript A ends at 400, transcript B ends at 600
    var tr_a = Transcript.init(allocator);
    defer tr_a.deinit();
    tr_a.id = 1;
    tr_a.beg = 100;
    tr_a.end = 400;
    tr_a.strand = .forward;
    tr_a.biotype = .protein_coding;

    var tr_b = Transcript.init(allocator);
    defer tr_b.deinit();
    tr_b.id = 2;
    tr_b.beg = 200;
    tr_b.end = 600;
    tr_b.strand = .forward;
    tr_b.biotype = .protein_coding;

    // Add in reverse order to test heap ordering
    try ctx.active_transcripts.add(&tr_b);
    try ctx.active_transcripts.add(&tr_a);
    try std.testing.expectEqual(@as(usize, 2), ctx.active_transcripts.count());

    // Flush at pos=500: only tr_a (end=400) should be flushed
    try ctx.hapFlush(500);

    try std.testing.expectEqual(@as(usize, 1), ctx.active_transcripts.count());
    try std.testing.expectEqual(@as(usize, 1), ctx.rm_transcripts.items.len);
    try std.testing.expectEqual(@as(u32, 1), ctx.rm_transcripts.items[0].id); // tr_a

    // Flush at pos=700: tr_b (end=600) should be flushed
    try ctx.hapFlush(700);

    try std.testing.expectEqual(@as(usize, 0), ctx.active_transcripts.count());
    try std.testing.expectEqual(@as(usize, 2), ctx.rm_transcripts.items.len);
    try std.testing.expectEqual(@as(u32, 2), ctx.rm_transcripts.items[1].id); // tr_b

    // Clear rm_transcripts without destroying aux
    ctx.rm_transcripts.clearRetainingCapacity();
}

test "hapFlush: POS_MAX drains all active transcripts" {
    const allocator = std.testing.allocator;

    var ctx = try CsqContext.init(allocator, .{
        .gff_fname = "",
        .phase = .drop_gt,
        .local_csq = true,
    });
    defer ctx.deinit();

    var tr1 = Transcript.init(allocator);
    defer tr1.deinit();
    tr1.id = 1;
    tr1.beg = 0;
    tr1.end = 100;
    tr1.biotype = .protein_coding;

    var tr2 = Transcript.init(allocator);
    defer tr2.deinit();
    tr2.id = 2;
    tr2.beg = 50;
    tr2.end = std.math.maxInt(u32) - 1;
    tr2.biotype = .protein_coding;

    try ctx.active_transcripts.add(&tr1);
    try ctx.active_transcripts.add(&tr2);

    try ctx.hapFlush(POS_MAX);

    try std.testing.expectEqual(@as(usize, 0), ctx.active_transcripts.count());
    try std.testing.expectEqual(@as(usize, 2), ctx.rm_transcripts.items.len);

    ctx.rm_transcripts.clearRetainingCapacity();
}

test "icsq2ToBit: basic calculations" {
    const r0 = CsqContext.icsq2ToBit(0);
    try std.testing.expectEqual(@as(u32, 0), r0.ival);
    try std.testing.expectEqual(@as(u5, 0), r0.ibit);

    const r29 = CsqContext.icsq2ToBit(29);
    try std.testing.expectEqual(@as(u32, 0), r29.ival);
    try std.testing.expectEqual(@as(u5, 29), r29.ibit);

    const r30 = CsqContext.icsq2ToBit(30);
    try std.testing.expectEqual(@as(u32, 1), r30.ival);
    try std.testing.expectEqual(@as(u5, 0), r30.ibit);

    const r61 = CsqContext.icsq2ToBit(61);
    try std.testing.expectEqual(@as(u32, 2), r61.ival);
    try std.testing.expectEqual(@as(u5, 1), r61.ibit);
}
