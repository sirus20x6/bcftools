// haplotype.zig — Haplotype tree construction and finalization logic.
// Ported from bcftools csq.c (lines 1596-2670).
//
// This module manages per-transcript haplotype trees and their finalization
// into consequence annotations (protein-level changes such as missense,
// synonymous, stop_gained, frameshift, etc.).

const std = @import("std");
const Allocator = std.mem.Allocator;
const ArrayList = std.ArrayList;

const types = @import("types.zig");
const gff_types = @import("../gff/types.zig");
const translate = @import("translate.zig");

const CsqType = types.CsqType;
const HapNode = types.HapNode;
const HapNodeType = types.HapNodeType;
const Hstack = types.Hstack;
const Csq = types.Csq;
const Vcsq = types.Vcsq;
const n_ref_pad = types.n_ref_pad;

// ---------------------------------------------------------------------------
// HapInitResult — return value from hapInit
// ---------------------------------------------------------------------------

pub const HapInitResult = enum {
    /// Variant was added to the haplotype tree.
    added,
    /// Variant overlaps a previous variant on this haplotype.
    overlapping,
    /// Variant was silently discarded (intronic, alt=ref, etc.).
    discarded,
};

// ---------------------------------------------------------------------------
// HapContext — traversal context for a transcript's haplotype tree
// ---------------------------------------------------------------------------

/// Context for haplotype tree operations on a single transcript.
/// Holds working buffers for DFS traversal, translation, and consequence
/// determination.  Mirrors the C `hap_t` struct.
pub const HapContext = struct {
    allocator: Allocator,
    stack: ArrayList(Hstack),
    tr: ?*gff_types.Transcript,
    /// Spliced haplotype sequence (reference strand).
    sseq: ArrayList(u8),
    /// Translated haplotype protein sequence.
    tseq: ArrayList(u8),
    /// Translated reference protein sequence.
    tref: ArrayList(u8),
    /// Stop/start codon annotations for tseq.
    tseq_stop: ArrayList(u8),
    /// Stop/start codon annotations for tref.
    tref_stop: ArrayList(u8),
    /// Splice-offset for the first node on the stack.
    sbeg: u32,
    /// Whether an upstream premature stop was detected.
    upstream_stop: bool,
    /// Genetic code table to use for translation.
    gencode: *const translate.GeneticCode,

    pub fn init(allocator: Allocator, gencode: *const translate.GeneticCode) HapContext {
        return .{
            .allocator = allocator,
            .stack = .empty,
            .tr = null,
            .sseq = .empty,
            .tseq = .empty,
            .tref = .empty,
            .tseq_stop = .empty,
            .tref_stop = .empty,
            .sbeg = 0,
            .upstream_stop = false,
            .gencode = gencode,
        };
    }

    pub fn deinit(self: *HapContext) void {
        self.stack.deinit(self.allocator);
        self.sseq.deinit(self.allocator);
        self.tseq.deinit(self.allocator);
        self.tref.deinit(self.allocator);
        self.tseq_stop.deinit(self.allocator);
        self.tref_stop.deinit(self.allocator);
    }
};

// ---------------------------------------------------------------------------
// cdsTranslate — translate a (possibly indel-modified) CDS to protein
// ---------------------------------------------------------------------------

/// Translate a spliced CDS sequence to amino acids.
///
/// This is the core translation function that handles both forward and reverse
/// strand, with left/right padding from the reference to complete partial
/// codons at the edges of the variant region.
///
/// Parameters:
///   allocator  — allocator for the output ArrayLists
///   sref       — full spliced reference sequence (with n_ref_pad padding on each end)
///   sref_len   — length of sref
///   seq        — the query (possibly variant-containing) sequence fragment
///   seq_total  — total length of the complete spliced query transcript
///                (needed for reverse strand codon boundary calculation)
///   seq_beg    — offset of `seq` within the spliced query transcript (0-based)
///   ref_beg    — offset of `seq` within sref (excluding padding), 0-based
///   ref_end    — one past the last base of the reference region corresponding to seq
///   strand     — coding strand (.forward or .reverse)
///   result     — output: translated amino acid sequence
///   result_stop — output: stop/start annotation per codon position
///   fill       — if non-zero, extend translation to end (fwd) or start (rev) of transcript
///   gencode    — genetic code table for codon-to-AA lookup
pub fn cdsTranslate(
    allocator: Allocator,
    sref: []const u8,
    sref_len: usize,
    seq: []const u8,
    seq_total: usize,
    seq_beg: u32,
    ref_beg: u32,
    ref_end: u32,
    strand: gff_types.Strand,
    result: *ArrayList(u8),
    result_stop: *ArrayList(u8),
    fill: i32,
    gencode: *const translate.GeneticCode,
) !void {
    result.clearRetainingCapacity();
    result_stop.clearRetainingCapacity();

    if (seq.len == 0) {
        try result.append(allocator, '?');
        try result_stop.append(allocator, '?');
        return;
    }

    if (strand == .forward) {
        try cdsTranslateFwd(allocator, sref, sref_len, seq, seq_beg, ref_beg, ref_end, result, result_stop, fill, gencode);
    } else if (strand == .reverse) {
        try cdsTranslateRev(allocator, sref, seq, seq_total, seq_beg, ref_beg, ref_end, result, result_stop, fill, gencode);
    } else {
        return error.InvalidStrand;
    }
}

/// Forward-strand translation.  Mirrors the `strand==STRAND_FWD` branch
/// of the C cds_translate().
fn cdsTranslateFwd(
    allocator: Allocator,
    sref: []const u8,
    sref_len: usize,
    seq: []const u8,
    seq_beg: u32,
    ref_beg: u32,
    ref_end: u32,
    result: *ArrayList(u8),
    result_stop: *ArrayList(u8),
    fill: i32,
    gencode: *const translate.GeneticCode,
) !void {
    var tmp: [3]u8 = undefined;

    // left padding — number of bases to borrow from reference before seq
    const npad: u32 = seq_beg % 3;
    std.debug.assert(npad <= ref_beg);

    var i: u32 = 0;
    while (i < npad) : (i += 1) {
        tmp[i] = sref[ref_beg + i - npad + n_ref_pad];
    }
    while (i < 3 and (i - npad) < seq.len) : (i += 1) {
        tmp[i] = seq[i - npad];
    }
    const remaining: usize = if (seq.len + npad >= i) seq.len - (i - npad) else 0;

    // Position within seq for continuing after first codon
    const seq_pos: usize = if (i >= npad) i - npad else 0;

    // How many leftover bases from seq go into the trailing partial codon
    var trailing: u32 = 0;

    if (i == 3) {
        // Translate the first (left-padded) codon
        try appendCodonFwd(allocator, result, result_stop, &tmp, gencode);

        // Translate full codons from seq
        const full_codon_bytes = remaining - (remaining % 3);
        var offset: usize = seq_pos;
        while (offset + 3 <= seq_pos + full_codon_bytes) {
            const codon: *const [3]u8 = seq[offset..][0..3];
            try appendCodonFwd(allocator, result, result_stop, codon, gencode);
            offset += 3;
        }

        // Collect trailing partial codon bases from seq
        trailing = @intCast(remaining % 3);
        for (0..trailing) |t| {
            tmp[t] = seq[offset + t];
        }
    } else {
        // seq was shorter than 3 - npad bases; everything is in tmp[0..i]
        trailing = i;
    }

    // right padding — fill trailing partial codon from reference after seq
    var ref_pos: usize = ref_end + n_ref_pad;
    if (trailing > 0) {
        var t = trailing;
        while (t < 3) : (t += 1) {
            if (ref_pos < sref.len) {
                tmp[t] = sref[ref_pos];
                ref_pos += 1;
            }
        }
        try appendCodonFwd(allocator, result, result_stop, &tmp, gencode);
    }

    // fill extension (for frameshift: continue translating reference codons)
    if (fill != 0) {
        const end_pos: usize = if (sref_len >= n_ref_pad) sref_len - n_ref_pad else sref_len;
        while (ref_pos + 3 <= end_pos) {
            const codon: *const [3]u8 = sref[ref_pos..][0..3];
            try appendCodonFwd(allocator, result, result_stop, codon, gencode);
            ref_pos += 3;
        }
    }
}

/// Reverse-strand translation.  Mirrors the `strand==STRAND_REV` branch
/// of the C cds_translate().
fn cdsTranslateRev(
    allocator: Allocator,
    sref: []const u8,
    seq: []const u8,
    seq_total: usize,
    seq_beg: u32,
    ref_beg: u32,
    ref_end: u32,
    result: *ArrayList(u8),
    result_stop: *ArrayList(u8),
    fill: i32,
    gencode: *const translate.GeneticCode,
) !void {
    var tmp: [3]u8 = undefined;

    // right padding — number of bases to borrow from reference after seq
    const tail = seq_beg + seq.len;
    const npad_raw: usize = if (seq_total >= tail) (seq_total - tail) % 3 else 0;
    const npad: u32 = @intCast(npad_raw);

    // Set up the initial partial codon from reference padding + end of seq
    var seq_end: usize = seq.len; // pointer moving backwards through seq

    var i: i32 = 2; // index into tmp, filling right-to-left

    if (npad == 2) {
        tmp[1] = sref[ref_end + n_ref_pad];
        tmp[2] = sref[ref_end + n_ref_pad + 1];
        i = 0;
    } else if (npad == 1) {
        tmp[2] = sref[ref_end + n_ref_pad];
        i = 1;
    } else {
        i = 2;
    }

    // Fill the rest of the first codon from the end of seq (right-to-left)
    while (i >= 0 and seq_end > 0) {
        seq_end -= 1;
        tmp[@intCast(i)] = seq[seq_end];
        i -= 1;
    }

    // How many trailing (leftward) bases remain for a partial codon
    var trailing: i32 = -1; // -1 means "no trailing partial codon"

    if (i == -1) {
        // Full first codon — translate it (reverse complement)
        try appendCodonRev(allocator, result, result_stop, &tmp, gencode);

        // Translate full codons moving leftward through seq
        while (seq_end >= 3) {
            seq_end -= 3;
            const codon: *const [3]u8 = seq[seq_end..][0..3];
            try appendCodonRev(allocator, result, result_stop, codon, gencode);
        }

        // Collect any remaining partial codon at the left of seq
        if (seq_end == 1) {
            tmp[2] = seq[0];
            trailing = 1;
        } else if (seq_end == 2) {
            tmp[1] = seq[0];
            tmp[2] = seq[1];
            trailing = 0;
        } else {
            trailing = -1;
        }
    } else {
        // Didn't fill a complete first codon — everything is partial
        trailing = i;
    }

    // left padding — fill remaining partial codon from reference before seq
    var ref_left: usize = n_ref_pad + ref_beg;
    if (trailing >= 0) {
        var t = trailing;
        while (t >= 0 and ref_left > 0) {
            ref_left -= 1;
            tmp[@intCast(t)] = sref[ref_left];
            t -= 1;
        }
        try appendCodonRev(allocator, result, result_stop, &tmp, gencode);
    }

    // fill extension (for frameshift: continue translating reference codons leftward)
    if (fill != 0) {
        while (ref_left >= n_ref_pad + 3) {
            ref_left -= 3;
            const codon: *const [3]u8 = sref[ref_left..][0..3];
            try appendCodonRev(allocator, result, result_stop, codon, gencode);
        }
    }
}

/// Append the forward-strand translation of a codon to the result buffers.
fn appendCodonFwd(
    allocator: Allocator,
    result: *ArrayList(u8),
    result_stop: *ArrayList(u8),
    codon: *const [3]u8,
    gencode: *const translate.GeneticCode,
) !void {
    const aa = translate.dna2aa(gencode, codon) orelse '?';
    const stop = translate.dna2stop(gencode, codon) orelse '?';
    try result.append(allocator, aa);
    try result_stop.append(allocator, stop);
}

/// Append the reverse-complement translation of a codon to the result buffers.
fn appendCodonRev(
    allocator: Allocator,
    result: *ArrayList(u8),
    result_stop: *ArrayList(u8),
    codon: *const [3]u8,
    gencode: *const translate.GeneticCode,
) !void {
    const aa = translate.cdna2aa(gencode, codon) orelse '?';
    const stop = translate.cdna2stop(gencode, codon) orelse '?';
    try result.append(allocator, aa);
    try result_stop.append(allocator, stop);
}

// ---------------------------------------------------------------------------
// hapAddCsq — determine consequence type from translated ref vs alt protein
// ---------------------------------------------------------------------------

/// Result of consequence comparison between reference and alternate protein.
pub const CsqResult = struct {
    csq_type: CsqType,
    upstream_stop: bool,
};

/// Determine the consequence type by comparing translated reference and
/// alternate protein sequences.
///
/// This implements the core logic of the C `hap_add_csq()` function:
///   - Compares ref vs alt amino acids character by character
///   - Detects synonymous, missense, stop_gained, stop_lost
///   - Handles frameshifts and inframe indels based on dlen
///   - Checks for upstream premature stops
///
/// Parameters:
///   tref       — translated reference protein
///   tref_stop  — stop/start annotations for reference
///   tseq       — translated alternate protein
///   tseq_stop  — stop/start annotations for alternate
///   dlen       — net length difference (alt_len - ref_len in bases)
///   indel      — true if any variant in the group is an indel
///   node_csq   — combined CSQ_COMPOUND flags from constituent nodes
///   is_sss     — true if the first node in the group is a splice-only node
///   is_compound — true if multiple variants are grouped together
///   has_upstream_stop — whether a premature stop was already detected upstream
pub fn hapAddCsq(
    tref: []const u8,
    tref_stop: []const u8,
    tseq: []const u8,
    tseq_stop: []const u8,
    dlen: i32,
    indel: bool,
    node_csq: CsqType,
    is_sss: bool,
    is_compound: bool,
    has_upstream_stop: bool,
) CsqResult {
    var csq = node_csq;
    var rm_csq = CsqType{};
    var upstream_stop = has_upstream_stop;

    if (dlen == 0 and indel) {
        csq.inframe_altering = true;
    }

    if (!is_sss) {
        // Truncate at first stop codon in reference
        var ref_len = tref.len;
        for (tref_stop, 0..) |ch, idx| {
            if (ch == '*') {
                ref_len = idx + 1;
                break;
            }
        }

        // Truncate at first stop codon in alternate
        var seq_len = tseq.len;
        for (tseq_stop, 0..) |ch, idx| {
            if (ch == '*') {
                seq_len = idx + 1;
                upstream_stop = true;
                break;
            }
        }

        const tref_trunc = tref[0..ref_len];
        const tseq_trunc = tseq[0..seq_len];
        const tref_stop_trunc = tref_stop[0..ref_len];
        const tseq_stop_trunc = tseq_stop[0..seq_len];

        // Check stop_lost: does the reference stop codon survive?
        if (csq.stop_lost) {
            if (ref_len > 0 and seq_len > 0 and
                tref_stop_trunc[ref_len - 1] == '*' and
                tref_stop_trunc[ref_len - 1] == tseq_stop_trunc[seq_len - 1])
            {
                rm_csq.stop_lost = true;
                csq.stop_retained = true;
            } else if (ref_len > 0 and tref_stop_trunc[ref_len - 1] != '*') {
                // Incomplete CDS (3' end)
                if (seq_len > 0 and tseq_stop_trunc[seq_len - 1] == '*') {
                    rm_csq.stop_gained = true;
                    csq.stop_retained = true;
                } else {
                    csq.incomplete_cds = true;
                }
            }
        }

        // Check start_lost: does the reference start codon survive?
        if (csq.start_lost) {
            if (ref_len > 0 and seq_len > 0 and
                tref_stop_trunc[ref_len - 1] == 'M' and
                tref_stop_trunc[ref_len - 1] == tseq_stop_trunc[seq_len - 1])
            {
                rm_csq.start_lost = true;
                csq.start_retained = true;
            }
        }

        if (dlen != 0) {
            // Indel consequence classification
            if (@mod(dlen, 3) != 0) {
                csq.frameshift_variant = true;
            } else if (dlen < 0) {
                csq.inframe_deletion = true;
            } else {
                csq.inframe_insertion = true;
            }
            // Check for stop gained by indel
            if (ref_len > 0 and seq_len > 0 and
                tref_stop_trunc[ref_len - 1] != '*' and
                tseq_stop_trunc[seq_len - 1] == '*')
            {
                csq.stop_gained = true;
            }
        } else {
            // Substitution: compare amino acids one by one
            var aa_change = false;
            const cmp_len = @min(tref_trunc.len, tseq_trunc.len);
            for (0..cmp_len) |idx| {
                if (tref_trunc[idx] == tseq_trunc[idx]) continue;
                aa_change = true;
                if (tref_stop_trunc[idx] == '*') {
                    csq.stop_lost = true;
                } else if (tseq_stop_trunc[idx] == '*') {
                    csq.stop_gained = true;
                } else {
                    csq.missense_variant = true;
                }
            }
            if (!aa_change) {
                csq.synonymous_variant = true;
            }
        }
    }

    // Compound inframe variants that actually introduce a premature stop
    // are really frameshifts + stop_gained
    if (is_compound and
        (csq.inframe_deletion or csq.inframe_insertion or csq.inframe_altering))
    {
        if (tseq_stop.len > 0 and tseq_stop[tseq_stop.len - 1] == '*') {
            rm_csq.inframe_deletion = true;
            rm_csq.inframe_insertion = true;
            rm_csq.inframe_altering = true;
            csq.frameshift_variant = true;
            csq.stop_gained = true;
        }
    }

    // Frameshift + start_lost: demote to splice-only
    if (csq.frameshift_variant and csq.start_lost) {
        rm_csq.frameshift_variant = true;
    }

    if (has_upstream_stop) csq.upstream_stop = true;

    // Apply the removal mask
    const csq_raw = csq.toInt() & ~rm_csq.toInt();
    return .{
        .csq_type = CsqType.fromInt(csq_raw),
        .upstream_stop = upstream_stop,
    };
}

// ---------------------------------------------------------------------------
// hapInit — initialize a haplotype tree node for a variant (skeleton)
// ---------------------------------------------------------------------------

/// Initialize a haplotype node for a variant overlapping a CDS exon.
///
/// Ported from the C `hap_init()` function (csq.c lines 1596-1743).
///
/// Steps:
///   1. Run splice analysis to check donor/acceptor/region/start/stop.
///   2. If no coding impact (splice-only), create HAP_SSS node and return .added.
///   3. If the variant overlaps coding sequence, build the spliced CDS sequence
///      incorporating the variant and create a HAP_CDS node.
///   4. If overlapping variants are detected, return .overlapping.
pub fn hapInit(
    allocator: Allocator,
    parent: *HapNode,
    child: *HapNode,
    cds: *const gff_types.CdsEntry,
    rec_pos: u32,
    ref_allele: []const u8,
    alt_allele: []const u8,
    ial: u32,
    tscript_aux: *const types.Tscript,
) !HapInitResult {
    const tr = cds.tr;
    child.icds = cds.icds;
    child.vcf_ial = @intCast(ial);

    // ── Step 1: splice analysis ────────────────────────────────────

    const splice_mod = @import("splice.zig");
    var splice = splice_mod.Splice.init(allocator, tr);
    defer splice.deinit();

    splice.reset(
        @intCast(rec_pos),
        @intCast(ref_allele.len),
        @intCast(ial),
        ref_allele,
        alt_allele,
    );
    splice.vcf.alen = @intCast(alt_allele.len);
    splice.flags = .{
        .check_acceptor = true,
        .check_donor = true,
        .set_refalt = true,
        .check_utr = true,
        .check_start = false,
        .check_stop = false,
        .check_region_beg = cds.icds != 0,
        .check_region_end = cds.icds != tr.cds.items.len - 1,
    };

    // Check start codon: first exon on the coding strand
    if (tr.trim != .prime5) {
        if (tr.strand == .forward and cds.icds == 0)
            splice.flags.check_start = true;
        if (tr.strand == .reverse and cds.icds == tr.cds.items.len - 1)
            splice.flags.check_start = true;
    }
    // Check stop codon: last exon on the coding strand
    if (tr.trim != .prime3) {
        if (tr.strand == .forward and cds.icds == tr.cds.items.len - 1)
            splice.flags.check_stop = true;
        if (tr.strand == .reverse and cds.icds == 0)
            splice.flags.check_stop = true;
    }

    // Verify start codon is actually M before checking for start_lost
    if (splice.flags.check_start) {
        if (tscript_aux.ref_seq) |ref_seq| {
            const translate_mod = @import("translate.zig");
            if (tr.strand == .forward) {
                const off = n_ref_pad + cds.beg -| tr.beg;
                if (off + 3 <= ref_seq.len) {
                    const stop_ch = translate_mod.dna2stop(translate_mod.findGeneticCode(0) orelse unreachable, ref_seq[off..][0..3]);
                    if (stop_ch == null or stop_ch.? != 'M')
                        splice.flags.check_start = false;
                }
            } else if (tr.strand == .reverse) {
                const off = n_ref_pad + cds.beg -| tr.beg + cds.len -| 3;
                if (off + 3 <= ref_seq.len) {
                    const stop_ch = translate_mod.cdna2stop(translate_mod.findGeneticCode(0) orelse unreachable, ref_seq[off..][0..3]);
                    if (stop_ch == null or stop_ch.? != 'M')
                        splice.flags.check_start = false;
                }
            }
        }
    }

    // Run splice consequence analysis
    const ret = splice.spliceCsq(cds.beg, cds.beg + cds.len - 1);

    // ── Step 2: handle non-coding results ──────────────────────────

    if (ret == .var_ref) return .discarded; // not a variant

    if (ret == .outside or ret == .overlap) {
        if (splice.csq.toInt() == 0) return .discarded; // fully intronic

        // Splice region/acceptor/donor: create HAP_SSS node
        child.payload = .{ .sss = {} };
        child.sbeg = 0;
        child.rbeg = rec_pos;
        child.rlen = 0;
        child.dlen = 0;

        // Build "ref>alt" string
        const var_str = try allocator.alloc(u8, ref_allele.len + 1 + alt_allele.len);
        @memcpy(var_str[0..ref_allele.len], ref_allele);
        var_str[ref_allele.len] = '>';
        @memcpy(var_str[ref_allele.len + 1 ..], alt_allele);
        child.var_str = var_str;

        child.csq = splice.csq;
        return .added;
    }

    // Clear synonymous if set by splice (will be re-evaluated after translation)
    if (splice.csq.synonymous_variant)
        splice.csq.synonymous_variant = false;

    // ── Step 3: build the spliced CDS sequence ─────────────────────

    // Handle variant overlapping exon boundary: trim to exon
    var dbeg: u32 = 0;
    if (splice.ref_beg < cds.beg) {
        dbeg = cds.beg - splice.ref_beg;
        splice.ref_beg = cds.beg;
    }

    // Parent must not be HAP_SSS for CDS sequence building
    std.debug.assert(parent.payload != .sss);

    var seq_buf: ArrayList(u8) = .empty;

    const ref_seq = tscript_aux.ref_seq orelse return .discarded;

    if (parent.payload == .cds) {
        const parent_icds = parent.icds;

        if (parent_icds != cds.icds) {
            // Variant is on a new exon: finish the previous exon
            const prev_cds = tr.cds.items[parent_icds];
            const prev_exon_end = prev_cds.beg + prev_cds.len;
            const parent_var_end = parent.rbeg + @as(u32, @intCast(@max(@as(i32, 0), parent.rlen)));
            if (prev_exon_end > parent_var_end) {
                const len = prev_exon_end - parent_var_end;
                const src_off = n_ref_pad + parent_var_end -| tr.beg;
                if (src_off + len <= ref_seq.len)
                    try seq_buf.appendSlice(allocator, ref_seq[src_off .. src_off + len]);
            }

            // Append any skipped non-variant exons
            var i: u32 = parent_icds + 1;
            while (i < cds.icds) : (i += 1) {
                const skip_cds = tr.cds.items[i];
                const src_off = n_ref_pad + skip_cds.beg -| tr.beg;
                if (src_off + skip_cds.len <= ref_seq.len)
                    try seq_buf.appendSlice(allocator, ref_seq[src_off .. src_off + skip_cds.len]);
            }
        }

        if (parent_icds == child.icds) {
            // Same exon: append reference gap between parent variant end and this variant
            const parent_var_end = parent.rbeg + @as(u32, @intCast(@max(@as(i32, 0), parent.rlen)));
            if (splice.ref_beg < parent_var_end) {
                // Overlapping variants
                seq_buf.deinit(allocator);
                return .overlapping;
            }
            const gap = splice.ref_beg - parent_var_end;
            if (gap > 0) {
                const src_off = n_ref_pad + parent_var_end -| tr.beg;
                if (src_off + gap <= ref_seq.len)
                    try seq_buf.appendSlice(allocator, ref_seq[src_off .. src_off + gap]);
            }
        } else {
            // Different exon: reference from start of new exon to variant
            const gap = splice.ref_beg - cds.beg;
            if (gap > 0) {
                const src_off = n_ref_pad + cds.beg -| tr.beg;
                if (src_off + gap <= ref_seq.len)
                    try seq_buf.appendSlice(allocator, ref_seq[src_off .. src_off + gap]);
            }
        }
    }

    // Append the alternate allele (trimmed by dbeg for exon-boundary overlap)
    if (splice.kalt.items.len > dbeg)
        try seq_buf.appendSlice(allocator, splice.kalt.items[dbeg..]);

    // Populate the child node as HAP_CDS
    const owned_seq = try seq_buf.toOwnedSlice(allocator);
    child.payload = .{ .cds = .{ .seq = owned_seq } };
    child.sbeg = cds.pos + (splice.ref_beg - cds.beg);
    child.rbeg = splice.ref_beg;
    child.rlen = @intCast(splice.kref.items.len);
    child.prev = parent;
    child.csq = splice.csq;

    // Set dlen and build "ref>alt" string
    child.dlen = @as(i32, @intCast(alt_allele.len)) - @as(i32, @intCast(ref_allele.len));
    const var_str = try allocator.alloc(u8, ref_allele.len + 1 + alt_allele.len);
    @memcpy(var_str[0..ref_allele.len], ref_allele);
    var_str[ref_allele.len] = '>';
    @memcpy(var_str[ref_allele.len + 1 ..], alt_allele);
    child.var_str = var_str;

    // If the whole CDS is modified/deleted, demote to HAP_SSS
    if (child.rbeg + @as(u32, @intCast(@max(@as(i32, 0), child.rlen))) > cds.beg + cds.len) {
        child.payload = .{ .sss = {} };
        if (child.csq.toInt() == 0) child.csq.coding_sequence = true;
    }

    return .added;
}

// ---------------------------------------------------------------------------
// hapFinalize — DFS traversal of haplotype tree (skeleton)
// ---------------------------------------------------------------------------

/// Finalize all haplotypes for a transcript by performing a DFS traversal
/// of the haplotype tree.
///
/// For each leaf node reached during traversal:
///   1. The spliced alt sequence is reconstructed from the stack
///   2. It is broken into independent parts by codon boundaries
///   3. Each part is translated (both ref and alt)
///   4. hapAddCsq determines the consequence
///
/// This is a skeleton. The full implementation requires:
///   - Building the spliced reference (tscript_splice_ref)
///   - The explicit DFS stack traversal matching the C code's break-by-codon logic
///   - Calling cdsTranslate for each independent part
///   - Calling hapAddCsq and pushing consequences to the output buffer
///
/// TODO:
///   - Implement tscript_splice_ref (build padded spliced reference)
///   - Implement the full DFS with codon-boundary partitioning
///   - Forward strand: walk stack indices 1..istack, group by codon boundaries
///   - Reverse strand: walk stack indices istack..1, group by codon boundaries
///   - For each group: translate alt and ref, call hapAddCsq
pub fn hapFinalize(ctx: *HapContext) !void {
    _ = ctx;
    // TODO: Full DFS implementation
    //
    // Pseudocode from the C version:
    //
    //   1. Ensure sref is built (tscript_splice_ref)
    //   2. Push root onto stack at index 0
    //   3. While stack is not empty:
    //      a. Advance to next unvisited child
    //      b. If no more children, pop (istack--)
    //      c. Otherwise push child, append its seq to sseq
    //      d. If child is a leaf (nend > 0):
    //         - For forward strand: walk i=1..istack, group variants by codon boundary
    //         - For reverse strand: walk i=istack..1, group variants by codon boundary
    //         - For each group:
    //           * Extract alt sequence from sseq
    //           * Extract ref sequence from sref
    //           * Call cdsTranslate for both
    //           * Call hapAddCsq
}

// ---------------------------------------------------------------------------
// hapFlush — flush completed transcripts (skeleton)
// ---------------------------------------------------------------------------

/// Flush completed transcripts from the active-transcript heap.
///
/// TODO:
///   - Pop transcripts whose end <= pos from the heap
///   - Call hapFinalize for each
///   - Stage consequences for VCF or text output
///   - Mark transcripts for deferred deletion
pub fn hapFlush(ctx: *HapContext, pos: u32) !void {
    _ = ctx;
    _ = pos;
    // TODO: implement heap-based flushing
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

pub const HapError = error{
    InvalidStrand,
};

// ===========================================================================
// Tests
// ===========================================================================

test "cdsTranslate forward strand — ATGCCCAGATAA translates to MPR*" {
    const allocator = std.testing.allocator;
    const gencode = translate.findGeneticCode(0) orelse unreachable;

    // Build a padded reference: 10 N's + "ATGCCCAGATAA" + 10 N's
    const coding = "ATGCCCAGATAA";
    const pad = "N" ** n_ref_pad;
    const sref = pad ++ coding ++ pad;

    var result: ArrayList(u8) = .empty;
    defer result.deinit(allocator);
    var result_stop: ArrayList(u8) = .empty;
    defer result_stop.deinit(allocator);

    try cdsTranslate(
        allocator,
        sref,
        sref.len,
        coding, // seq = the coding region itself
        coding.len, // seq_total = total transcript length
        0, // seq_beg = start of transcript
        0, // ref_beg = start of sref (excluding pad)
        @intCast(coding.len), // ref_end = end of coding region
        .forward,
        &result,
        &result_stop,
        0, // no fill
        gencode,
    );

    try std.testing.expectEqualStrings("MPR*", result.items);
}

test "cdsTranslate forward strand — partial codon left padding" {
    const allocator = std.testing.allocator;
    const gencode = translate.findGeneticCode(0) orelse unreachable;

    // Full transcript: ATGCCCAGATAA (4 codons)
    // We test translating just "CCCAGATAA" starting at seq_beg=3 (after ATG).
    // Left pad: seq_beg%3 = 0, so no padding needed — codons are CCC AGA TAA
    const coding = "ATGCCCAGATAA";
    const pad = "N" ** n_ref_pad;
    const sref = pad ++ coding ++ pad;
    const seq = "CCCAGATAA"; // positions 3..12 in the transcript

    var result: ArrayList(u8) = .empty;
    defer result.deinit(allocator);
    var result_stop: ArrayList(u8) = .empty;
    defer result_stop.deinit(allocator);

    try cdsTranslate(
        allocator,
        sref,
        sref.len,
        seq,
        coding.len,
        3, // seq_beg
        3, // ref_beg
        12, // ref_end
        .forward,
        &result,
        &result_stop,
        0,
        gencode,
    );

    try std.testing.expectEqualStrings("PR*", result.items);
}

test "cdsTranslate reverse strand — TTATCTGGGCAT translates to MPR*" {
    const allocator = std.testing.allocator;
    const gencode = translate.findGeneticCode(0) orelse unreachable;

    // On the reverse strand, the sequence TTATCTGGGCAT is read as
    // reverse complement: ATG CCC AGA TAA -> M P R *
    //
    // The reverse strand translation reads codons right-to-left from
    // the sequence, taking the reverse complement of each triplet.
    const coding = "TTATCTGGGCAT";
    const pad = "N" ** n_ref_pad;
    const sref = pad ++ coding ++ pad;

    var result: ArrayList(u8) = .empty;
    defer result.deinit(allocator);
    var result_stop: ArrayList(u8) = .empty;
    defer result_stop.deinit(allocator);

    try cdsTranslate(
        allocator,
        sref,
        sref.len,
        coding,
        coding.len,
        0, // seq_beg
        0, // ref_beg
        @intCast(coding.len), // ref_end
        .reverse,
        &result,
        &result_stop,
        0,
        gencode,
    );

    try std.testing.expectEqualStrings("MPR*", result.items);
}

test "cdsTranslate empty sequence returns ?" {
    const allocator = std.testing.allocator;
    const gencode = translate.findGeneticCode(0) orelse unreachable;

    const pad = "N" ** n_ref_pad;
    const sref = pad ++ pad;

    var result: ArrayList(u8) = .empty;
    defer result.deinit(allocator);
    var result_stop: ArrayList(u8) = .empty;
    defer result_stop.deinit(allocator);

    try cdsTranslate(
        allocator,
        sref,
        sref.len,
        "", // empty seq
        0, // seq_total
        0, // seq_beg
        0, // ref_beg
        0, // ref_end
        .forward,
        &result,
        &result_stop,
        0, // no fill
        gencode,
    );

    try std.testing.expectEqualStrings("?", result.items);
    try std.testing.expectEqualStrings("?", result_stop.items);
}

test "hapAddCsq — missense: ref MPR vs alt MHR" {
    const tref = "MPR";
    const tref_stop = "M--";
    const tseq = "MHR";
    const tseq_stop = "M--";

    const res = hapAddCsq(
        tref,
        tref_stop,
        tseq,
        tseq_stop,
        0, // dlen=0: substitution
        false, // not an indel
        CsqType{}, // no pre-existing node csq
        false, // not splice-only
        false, // not compound
        false, // no upstream stop
    );

    try std.testing.expect(res.csq_type.missense_variant);
    try std.testing.expect(!res.csq_type.synonymous_variant);
    try std.testing.expect(!res.csq_type.stop_gained);
    try std.testing.expect(!res.csq_type.stop_lost);
}

test "hapAddCsq — synonymous: ref MPR vs alt MPR" {
    const tref = "MPR";
    const tref_stop = "M--";
    const tseq = "MPR";
    const tseq_stop = "M--";

    const res = hapAddCsq(
        tref,
        tref_stop,
        tseq,
        tseq_stop,
        0,
        false,
        CsqType{},
        false,
        false,
        false,
    );

    try std.testing.expect(res.csq_type.synonymous_variant);
    try std.testing.expect(!res.csq_type.missense_variant);
}

test "hapAddCsq — stop_gained: ref MPR vs alt M*R" {
    const tref = "MPR";
    const tref_stop = "M--";
    const tseq = "M*R";
    const tseq_stop = "M*-";

    const res = hapAddCsq(
        tref,
        tref_stop,
        tseq,
        tseq_stop,
        0,
        false,
        CsqType{},
        false,
        false,
        false,
    );

    try std.testing.expect(res.csq_type.stop_gained);
    try std.testing.expect(!res.csq_type.synonymous_variant);
}

test "hapAddCsq — stop_lost: ref MPR* vs alt MPRW" {
    const tref = "MPR*";
    const tref_stop = "M--*";
    const tseq = "MPRW";
    const tseq_stop = "M---";

    const res = hapAddCsq(
        tref,
        tref_stop,
        tseq,
        tseq_stop,
        0,
        false,
        CsqType{},
        false,
        false,
        false,
    );

    try std.testing.expect(res.csq_type.stop_lost);
    try std.testing.expect(!res.csq_type.synonymous_variant);
}

test "hapAddCsq — frameshift: dlen not divisible by 3" {
    const tref = "MPR";
    const tref_stop = "M--";
    const tseq = "MI";
    const tseq_stop = "M-";

    const res = hapAddCsq(
        tref,
        tref_stop,
        tseq,
        tseq_stop,
        1, // insertion of 1 base
        true,
        CsqType{},
        false,
        false,
        false,
    );

    try std.testing.expect(res.csq_type.frameshift_variant);
    try std.testing.expect(!res.csq_type.inframe_insertion);
}

test "hapAddCsq — inframe deletion: dlen=-3" {
    const tref = "MPRS";
    const tref_stop = "M---";
    const tseq = "MPS";
    const tseq_stop = "M--";

    const res = hapAddCsq(
        tref,
        tref_stop,
        tseq,
        tseq_stop,
        -3,
        true,
        CsqType{},
        false,
        false,
        false,
    );

    try std.testing.expect(res.csq_type.inframe_deletion);
    try std.testing.expect(!res.csq_type.frameshift_variant);
}

test "hapAddCsq — inframe insertion: dlen=3" {
    const tref = "MPR";
    const tref_stop = "M--";
    const tseq = "MPAR";
    const tseq_stop = "M---";

    const res = hapAddCsq(
        tref,
        tref_stop,
        tseq,
        tseq_stop,
        3,
        true,
        CsqType{},
        false,
        false,
        false,
    );

    try std.testing.expect(res.csq_type.inframe_insertion);
    try std.testing.expect(!res.csq_type.frameshift_variant);
}

test "hapAddCsq — upstream_stop propagation" {
    const tref = "MPR";
    const tref_stop = "M--";
    const tseq = "MPR";
    const tseq_stop = "M--";

    const res = hapAddCsq(
        tref,
        tref_stop,
        tseq,
        tseq_stop,
        0,
        false,
        CsqType{},
        false,
        false,
        true, // has_upstream_stop = true
    );

    try std.testing.expect(res.csq_type.upstream_stop);
    try std.testing.expect(res.csq_type.synonymous_variant);
}

test "hapAddCsq — splice-only node sets no coding consequence" {
    const tref = "";
    const tref_stop = "";
    const tseq = "";
    const tseq_stop = "";

    var node_csq = CsqType{};
    node_csq.splice_donor = true;

    const res = hapAddCsq(
        tref,
        tref_stop,
        tseq,
        tseq_stop,
        0,
        false,
        node_csq,
        true, // is_sss = true (splice-only)
        false,
        false,
    );

    // SSS nodes skip the protein comparison entirely
    try std.testing.expect(res.csq_type.splice_donor);
    try std.testing.expect(!res.csq_type.synonymous_variant);
    try std.testing.expect(!res.csq_type.missense_variant);
}

test "HapContext init and deinit" {
    const allocator = std.testing.allocator;
    const gencode = translate.findGeneticCode(1) orelse unreachable;

    var ctx = HapContext.init(allocator, gencode);
    defer ctx.deinit();

    try std.testing.expect(ctx.tr == null);
    try std.testing.expect(!ctx.upstream_stop);
    try std.testing.expectEqual(@as(usize, 0), ctx.sseq.items.len);
}
