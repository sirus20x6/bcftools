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
const splice_mod = @import("splice.zig");

const CsqType = types.CsqType;
const HapNode = types.HapNode;
const HapNodeType = types.HapNodeType;
const Hstack = types.Hstack;
const Csq = types.Csq;
const Vcsq = types.Vcsq;
const n_ref_pad = types.n_ref_pad;

/// Extract the Tscript auxiliary data from a Transcript's opaque `aux` pointer.
/// Returns null if aux is not set.
fn getTscriptAux(tr: *const gff_types.Transcript) ?*types.Tscript {
    return @ptrCast(@alignCast(tr.aux orelse return null));
}

// ---------------------------------------------------------------------------
// HapInitResult — return value from hapInit
// ---------------------------------------------------------------------------

pub const HapInitResultKind = enum {
    /// Variant was added to the haplotype tree.
    added,
    /// Variant overlaps a previous variant on this haplotype.
    overlapping,
    /// Variant was silently discarded (intronic, alt=ref, etc.).
    discarded,
};

pub const HapInitResult = struct {
    kind: HapInitResultKind,
    /// The full splice consequence from splice_csq, BEFORE synonymous is cleared
    /// for the CDS path.  In the C code, splice_csq_del stages this via
    /// csq_stage_splice before returning SPLICE_INSIDE, so the caller needs the
    /// pre-clearing value for its own staging.
    splice_csq: CsqType = .{},
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
    /// When true, the ibeg node should be treated as SSS (splice-only):
    /// skip the variant string and OR in the node's splice consequence bits.
    /// This happens when frameshift+start_lost demotes the node.
    demote_to_sss: bool = false,
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
    // are really frameshifts + stop_gained.
    // C code (line 2280): uses the truncated tseq_stop length (after
    // truncation at the first stop codon in the !is_sss block above).
    // We check `upstream_stop` which is set to true when truncation happened.
    if (is_compound and
        (csq.inframe_deletion or csq.inframe_insertion or csq.inframe_altering))
    {
        // After truncation, the last character of tseq_stop is '*' if
        // a premature stop was found.  We detect this via upstream_stop
        // which was set during truncation, or by checking the original
        // tseq_stop for any stop codon.
        var trunc_last_is_stop = false;
        for (tseq_stop) |ch| {
            if (ch == '*') {
                trunc_last_is_stop = true;
                break;
            }
        }
        if (trunc_last_is_stop) {
            rm_csq.inframe_deletion = true;
            rm_csq.inframe_insertion = true;
            rm_csq.inframe_altering = true;
            csq.frameshift_variant = true;
            csq.stop_gained = true;
        }
    }

    // Frameshift + start_lost: demote to splice-only (HAP_SSS)
    var demote_to_sss = false;
    if (csq.frameshift_variant and csq.start_lost) {
        rm_csq.frameshift_variant = true;
        demote_to_sss = true;
    }

    if (has_upstream_stop) csq.upstream_stop = true;

    // Apply the removal mask
    const csq_raw = csq.toInt() & ~rm_csq.toInt();
    return .{
        .csq_type = CsqType.fromInt(csq_raw),
        .upstream_stop = upstream_stop,
        .demote_to_sss = demote_to_sss,
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
            if (tr.strand == .forward) {
                const off = n_ref_pad + cds.beg -| tr.beg;
                if (off + 3 <= ref_seq.len) {
                    const stop_ch = translate.dna2stop(translate.findGeneticCode(0) orelse unreachable, ref_seq[off..][0..3]);
                    if (stop_ch == null or stop_ch.? != 'M')
                        splice.flags.check_start = false;
                }
            } else if (tr.strand == .reverse) {
                const off = n_ref_pad + cds.beg -| tr.beg + cds.len -| 3;
                if (off + 3 <= ref_seq.len) {
                    const stop_ch = translate.cdna2stop(translate.findGeneticCode(0) orelse unreachable, ref_seq[off..][0..3]);
                    if (stop_ch == null or stop_ch.? != 'M')
                        splice.flags.check_start = false;
                }
            }
        }
    }

    // Set transcript reference for shifted_del_synonymous and build_hap
    splice.tr_ref = tscript_aux.ref_seq;

    // Run splice consequence analysis
    const ret = splice.spliceCsq(cds.beg, cds.beg + cds.len - 1);

    // ── Step 2: handle non-coding results ──────────────────────────

    if (ret == .var_ref) return .{ .kind = .discarded }; // not a variant

    if (ret == .outside or ret == .overlap) {
        if (splice.csq.toInt() == 0) return .{ .kind = .discarded }; // fully intronic

        // Splice region/acceptor/donor: create HAP_SSS node
        child.payload = .{ .sss = {} };
        child.sbeg = 0;
        child.rbeg = rec_pos;
        child.rec_pos = rec_pos;
        child.rlen = 0;
        child.dlen = 0;

        // Build "ref>alt" string
        const var_str = try allocator.alloc(u8, ref_allele.len + 1 + alt_allele.len);
        @memcpy(var_str[0..ref_allele.len], ref_allele);
        var_str[ref_allele.len] = '>';
        @memcpy(var_str[ref_allele.len + 1 ..], alt_allele);
        child.var_str = var_str;

        child.csq = splice.csq;
        return .{ .kind = .added, .splice_csq = splice.csq };
    }

    // Save the full splice CSQ before clearing synonymous.
    // In the C code, splice_csq_del stages this via csq_stage_splice before
    // returning SPLICE_INSIDE, so the caller needs the pre-clearing value.
    const full_splice_csq = splice.csq;

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

    const ref_seq = tscript_aux.ref_seq orelse return .{ .kind = .discarded };

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
                return .{ .kind = .overlapping };
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
    child.rec_pos = rec_pos; // original VCF position for consequence output
    // C: splice.kref.l -= dbeg; child->rlen = splice.kref.l;
    // Subtract dbeg from kref length when variant overlaps exon boundary
    child.rlen = @intCast(splice.kref.items.len -| dbeg);
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

    return .{ .kind = .added, .splice_csq = full_splice_csq };
}

// ---------------------------------------------------------------------------
// hapFinalize — DFS traversal of haplotype tree (skeleton)
// ---------------------------------------------------------------------------

/// Finalize all haplotypes for a transcript by performing a DFS traversal
/// of the haplotype tree.
///
/// Ported from the C `hap_finalize()` function (csq.c lines 2374-2551).
///
/// For each leaf node reached during traversal:
///   1. The spliced alt sequence is reconstructed from the stack
///   2. It is broken into independent parts by codon boundaries:
///      - Forward strand: left-to-right, break when dlen%3==0 AND
///        consecutive variants are in different codons
///      - Reverse strand: right-to-left with the same logic
///   3. Each part is translated (both ref and alt) via cdsTranslate
///   4. hapAddCsq determines the consequence type
pub fn hapFinalize(ctx: *HapContext) !void {
    const tr_opaque = ctx.tr orelse return;
    // The GFF transcript is stored as an opaque pointer in HapContext.
    // We need the Tscript (aux data) which holds the ref/sref and root.
    // By convention, the transcript's .aux field points to the Tscript.
    const tr_ptr: *const gff_types.Transcript = @ptrCast(@alignCast(tr_opaque));
    const tscript_aux: *types.Tscript = getTscriptAux(tr_ptr) orelse return;
    const allocator = ctx.allocator;

    // Build spliced reference if not done yet
    if (tscript_aux.sref == null)
        try tscriptSpliceRef(allocator, tscript_aux, tr_ptr);

    const sref = tscript_aux.sref orelse return;
    const sref_len: usize = @intCast(tscript_aux.nsref);

    // Initialize traversal stack with root
    ctx.stack.clearRetainingCapacity();
    try ctx.stack.append(allocator, .{
        .node = tscript_aux.root,
        .ichild = -1,
        .slen = 0,
        .dlen = 0,
    });

    ctx.sseq.clearRetainingCapacity();

    var istack: usize = 0;

    while (true) {
        if (istack >= ctx.stack.items.len) break;

        const node = ctx.stack.items[istack].node orelse break;

        // Find next non-null child
        var found_child = false;
        {
            var ichild = ctx.stack.items[istack].ichild + 1;
            while (ichild < @as(i32, @intCast(node.children.items.len))) : (ichild += 1) {
                ctx.stack.items[istack].ichild = ichild;
                found_child = true;
                break;
            }
            if (!found_child)
                ctx.stack.items[istack].ichild = @intCast(node.children.items.len);
        }

        if (!found_child) {
            if (istack == 0) break;
            istack -= 1;
            continue;
        }

        const child_idx: usize = @intCast(ctx.stack.items[istack].ichild);
        const child_node = node.children.items[child_idx];

        istack += 1;

        // Ensure stack capacity
        while (ctx.stack.items.len <= istack)
            try ctx.stack.append(allocator, .{});

        const parent_slen = ctx.stack.items[istack - 1].slen;
        const parent_dlen = ctx.stack.items[istack - 1].dlen;

        ctx.stack.items[istack] = .{
            .node = child_node,
            .ichild = -1,
            .slen = 0, // will be set below
            .dlen = parent_dlen + child_node.dlen,
        };

        // Build spliced sequence up to this point
        ctx.sseq.shrinkRetainingCapacity(parent_slen);
        if (child_node.payload == .cds) {
            if (child_node.payload.cds.seq) |seq|
                try ctx.sseq.appendSlice(allocator, seq);
        }
        ctx.stack.items[istack].slen = ctx.sseq.items.len;

        if (child_node.nend == 0) continue; // not a leaf

        // ── Leaf node: break into independent parts and translate ────

        const total_dlen = ctx.stack.items[istack].dlen;
        const sref_coding_len: i64 = @as(i64, @intCast(sref_len)) - 2 * @as(i64, n_ref_pad);
        const seq_m: usize = @intCast(@max(0, sref_coding_len + total_dlen));
        ctx.upstream_stop = false;

        // Set sbeg from the first real node (index 1)
        if (istack >= 1 and ctx.stack.items.len > 1) {
            if (ctx.stack.items[1].node) |n1|
                ctx.sbeg = n1.sbeg;
        }

        // If the leaf node is SSS-only (splice consequence, no CDS overlap),
        // handle it separately: no translation/vstr needed, just record the
        // splice consequence directly.  This mirrors the C code at csq.c:2296.
        if (child_node.payload == .sss) {
            // SSS-only leaf: record splice consequence without variant string
            const csq_entry_sss: types.Csq = .{
                .pos = child_node.rec_pos,
                .type_info = .{
                    .csq_type = child_node.csq,
                    .trid = tr_ptr.id,
                    .vcf_ial = @intCast(child_node.vcf_ial),
                    .gene = if (tr_ptr.gene) |g| blk: {
                        break :blk if (g.name) |n| std.mem.span(n) else null;
                    } else null,
                    .strand = tr_ptr.strand == .forward,
                    .biotype = @intFromEnum(tr_ptr.biotype),
                },
            };
            try child_node.csq_list.append(allocator, csq_entry_sss);
            istack -= 1;
            continue;
        }
        const stack = ctx.stack.items;

        if (tr_ptr.strand == .forward) {
            var i: usize = 0;
            var ibeg_s: i64 = -1;
            var dlen_acc: i32 = 0;
            var indel_flag = false;

            while (true) {
                i += 1;
                if (i > istack) break;

                std.debug.assert(stack[i].node.?.payload != .sss);

                dlen_acc += stack[i].node.?.dlen;
                if (stack[i].node.?.dlen != 0) indel_flag = true;

                // Decide whether to flush this portion
                if (i < istack) {
                    if (@rem(dlen_acc, 3) != 0) {
                        if (ibeg_s == -1) ibeg_s = @intCast(i);
                        continue;
                    }
                    // Same codon check (forward strand)
                    var icur = ctx.sbeg + (stack[i].slen -| rlenPlusDlen(stack[i].node.?));
                    const inext = ctx.sbeg + (stack[i + 1].slen -| rlenPlusDlen(stack[i + 1].node.?));
                    if (stack[i].node.?.dlen > 0)
                        icur += @as(usize, @intCast(stack[i].node.?.dlen))
                    else if (stack[i].node.?.dlen < 0)
                        icur += 1;
                    if (icur / 3 == inext / 3) {
                        if (ibeg_s == -1) ibeg_s = @intCast(i);
                        continue;
                    }
                }
                if (ibeg_s < 0) ibeg_s = @intCast(i);

                const ibeg_u: usize = @intCast(ibeg_s);
                const ioff = stack[ibeg_u].slen -| rlenPlusDlen(stack[ibeg_u].node.?);
                const icur_pos: u32 = @intCast(ctx.sbeg + ioff);
                const rbeg_val = stack[ibeg_u].node.?.sbeg;
                const rend_val = stack[i].node.?.sbeg + @as(u32, @intCast(@max(@as(i32, 0), stack[i].node.?.rlen)));
                const fill: i32 = @rem(dlen_acc, 3);

                // Translate alt
                if (ctx.sseq.items.len > 0) {
                    const alt_end = stack[i].slen;
                    const alt_seq = if (alt_end > ioff) ctx.sseq.items[ioff..alt_end] else &[_]u8{};
                    try cdsTranslate(allocator, sref, sref_len, alt_seq, seq_m, icur_pos, rbeg_val, rend_val, .forward, &ctx.tseq, &ctx.tseq_stop, fill, ctx.gencode);
                } else {
                    try cdsTranslate(allocator, sref, sref_len, &[_]u8{}, seq_m, icur_pos, rbeg_val, rend_val, .forward, &ctx.tseq, &ctx.tseq_stop, 0, ctx.gencode);
                }

                // Translate ref
                {
                    const ref_start = n_ref_pad + rbeg_val;
                    const ref_len = rend_val - rbeg_val;
                    const ref_slice = if (ref_start + ref_len <= sref.len) sref[ref_start .. ref_start + ref_len] else &[_]u8{};
                    try cdsTranslate(allocator, sref, sref_len, ref_slice, sref_len - 2 * n_ref_pad, rbeg_val, rbeg_val, rend_val, .forward, &ctx.tref, &ctx.tref_stop, fill, ctx.gencode);
                }

                // Determine consequence
                const csq_result = hapAddCsq(
                    ctx.tref.items,
                    ctx.tref_stop.items,
                    ctx.tseq.items,
                    ctx.tseq_stop.items,
                    dlen_acc,
                    indel_flag,
                    accumulateCompound(stack, ibeg_u, i),
                    stack[ibeg_u].node.?.payload == .sss,
                    ibeg_u != i,
                    ctx.upstream_stop,
                );
                ctx.upstream_stop = csq_result.upstream_stop;

                // Record consequence on the leaf node
                var merged_csq = csq_result.csq_type;

                // SSS demotion: OR in the ibeg node's own splice bits, skip vstr
                if (csq_result.demote_to_sss or stack[ibeg_u].node.?.payload == .sss) {
                    const node_splice_bits = stack[ibeg_u].node.?.csq.toInt();
                    merged_csq = CsqType.fromInt(merged_csq.toInt() | node_splice_bits);
                    var csq_entry: types.Csq = .{
                        .pos = stack[ibeg_u].node.?.rec_pos,
                        .type_info = .{
                            .csq_type = merged_csq,
                            .trid = tr_ptr.id,
                            .vcf_ial = @intCast(child_node.vcf_ial),
                            .gene = if (tr_ptr.gene) |g| blk: {
                                break :blk if (g.name) |n| std.mem.span(n) else null;
                            } else null,
                            .strand = tr_ptr.strand == .forward,
                            .biotype = @intFromEnum(tr_ptr.biotype),
                        },
                    };
                    _ = &csq_entry;
                    try child_node.csq_list.append(allocator, csq_entry);

                    ibeg_s = -1;
                    dlen_acc = 0;
                    indel_flag = false;
                    continue;
                }

                // Truncate tref/tseq at first stop codon for buildVstr (C lines 2203-2221).
                // hapAddCsq already determined consequences using its own truncation;
                // this truncation ensures buildVstr uses the short AA strings.
                for (ctx.tref_stop.items, 0..) |ch, tidx| {
                    if (ch == '*') {
                        ctx.tref.shrinkRetainingCapacity(tidx + 1);
                        ctx.tref_stop.shrinkRetainingCapacity(tidx + 1);
                        break;
                    }
                }
                for (ctx.tseq_stop.items, 0..) |ch, tidx| {
                    if (ch == '*') {
                        ctx.tseq.shrinkRetainingCapacity(tidx + 1);
                        ctx.tseq_stop.shrinkRetainingCapacity(tidx + 1);
                        break;
                    }
                }

                var csq_entry: types.Csq = .{
                    .pos = stack[if (tr_ptr.strand == .forward) ibeg_u else i].node.?.rec_pos,
                    .type_info = .{
                        .csq_type = csq_result.csq_type,
                        .trid = tr_ptr.id,
                        .vcf_ial = @intCast(child_node.vcf_ial),
                        .gene = if (tr_ptr.gene) |g| blk: {
                            break :blk if (g.name) |n| std.mem.span(n) else null;
                        } else null,
                        .strand = tr_ptr.strand == .forward,
                        .biotype = @intFromEnum(tr_ptr.biotype),
                    },
                };
                // Build variant string in vstr
                try buildVstr(allocator, &csq_entry.type_info.vstr, stack, ibeg_u, i, ctx, sref_len, seq_m, tr_ptr, csq_result.csq_type);
                try child_node.csq_list.append(allocator, csq_entry);

                // For compound variants (ibeg != iend), create CSQ_PRINTED_UPSTREAM
                // entries at non-ref-node positions (C csq.c lines 2334-2370).
                // ref_node = ibeg for forward strand.
                if (ibeg_u != i) {
                    const ref_node_pos_1based: u32 = stack[ibeg_u].node.?.rbeg + 1;
                    var j: usize = ibeg_u;
                    while (j <= i) : (j += 1) {
                        if (j == ibeg_u) continue; // skip ref_node
                        const node_j = stack[j].node.?;
                        // Create CSQ_PRINTED_UPSTREAM entry
                        var upstream_csq = types.CsqType{};
                        upstream_csq.printed_upstream = true;
                        // Also include the node's own CSQ bits
                        const node_csq_raw = node_j.csq.toInt();
                        upstream_csq = types.CsqType.fromInt(upstream_csq.toInt() | node_csq_raw);
                        const upstream_entry: types.Csq = .{
                            .pos = node_j.rec_pos,
                            .ref_pos = ref_node_pos_1based,
                            .type_info = .{
                                .csq_type = upstream_csq,
                                .trid = tr_ptr.id,
                                .gene = if (tr_ptr.gene) |g| blk: {
                                    break :blk if (g.name) |n| std.mem.span(n) else null;
                                } else null,
                                .strand = tr_ptr.strand == .forward,
                                .biotype = @intFromEnum(tr_ptr.biotype),
                            },
                        };
                        try child_node.csq_list.append(allocator, upstream_entry);
                    }
                }

                ibeg_s = -1;
                dlen_acc = 0;
                indel_flag = false;
            }
        } else if (tr_ptr.strand == .reverse) {
            var i: usize = istack + 1;
            var ibeg_s: i64 = -1;
            var dlen_acc: i32 = 0;
            var indel_flag = false;

            while (i > 1) {
                i -= 1;

                std.debug.assert(stack[i].node.?.payload != .sss);

                dlen_acc += stack[i].node.?.dlen;
                if (stack[i].node.?.dlen != 0) indel_flag = true;

                if (i > 1) {
                    if (@rem(dlen_acc, 3) != 0) {
                        if (ibeg_s == -1) ibeg_s = @intCast(i);
                        continue;
                    }
                    // Same codon check (reverse strand)
                    var icur_rev: i64 = @as(i64, @intCast(seq_m)) - 1 - @as(i64, @intCast(ctx.sbeg + (stack[i].slen -| rlenPlusDlen(stack[i].node.?))));
                    var inext_rev: i64 = @as(i64, @intCast(seq_m)) - 1 - @as(i64, @intCast(ctx.sbeg + (stack[i - 1].slen -| rlenPlusDlen(stack[i - 1].node.?))));
                    if (stack[i].node.?.dlen > 0) icur_rev += stack[i].node.?.dlen - 1 else if (stack[i].node.?.dlen < 0) icur_rev -= stack[i].node.?.dlen;
                    if (stack[i - 1].node.?.dlen > 0) inext_rev -= stack[i - 1].node.?.dlen;
                    if (icur_rev >= 0 and inext_rev >= 0) {
                        if (@as(u64, @intCast(icur_rev)) / 3 == @as(u64, @intCast(inext_rev)) / 3) {
                            if (ibeg_s == -1) ibeg_s = @intCast(i);
                            continue;
                        }
                    }
                }
                if (ibeg_s < 0) ibeg_s = @intCast(i);

                const ibeg_u: usize = @intCast(ibeg_s);
                const ioff = stack[i].slen -| rlenPlusDlen(stack[i].node.?);
                const icur_pos: u32 = @intCast(ctx.sbeg + ioff);
                const rbeg_val = stack[i].node.?.sbeg;
                const rend_val = stack[ibeg_u].node.?.sbeg + @as(u32, @intCast(@max(@as(i32, 0), stack[ibeg_u].node.?.rlen)));
                const fill: i32 = @rem(dlen_acc, 3);

                // Translate alt
                if (ctx.sseq.items.len > 0) {
                    const alt_end = stack[ibeg_u].slen;
                    const alt_seq = if (alt_end > ioff) ctx.sseq.items[ioff..alt_end] else &[_]u8{};
                    try cdsTranslate(allocator, sref, sref_len, alt_seq, seq_m, icur_pos, rbeg_val, rend_val, .reverse, &ctx.tseq, &ctx.tseq_stop, fill, ctx.gencode);
                } else {
                    try cdsTranslate(allocator, sref, sref_len, &[_]u8{}, seq_m, icur_pos, rbeg_val, rend_val, .reverse, &ctx.tseq, &ctx.tseq_stop, 0, ctx.gencode);
                }

                // Translate ref
                {
                    const ref_start = n_ref_pad + rbeg_val;
                    const ref_len = rend_val - rbeg_val;
                    const ref_slice = if (ref_start + ref_len <= sref.len) sref[ref_start .. ref_start + ref_len] else &[_]u8{};
                    try cdsTranslate(allocator, sref, sref_len, ref_slice, sref_len - 2 * n_ref_pad, rbeg_val, rbeg_val, rend_val, .reverse, &ctx.tref, &ctx.tref_stop, fill, ctx.gencode);
                }

                // Determine consequence
                const csq_result = hapAddCsq(
                    ctx.tref.items,
                    ctx.tref_stop.items,
                    ctx.tseq.items,
                    ctx.tseq_stop.items,
                    dlen_acc,
                    indel_flag,
                    accumulateCompound(stack, i, ibeg_u),
                    stack[i].node.?.payload == .sss,
                    i != ibeg_u,
                    ctx.upstream_stop,
                );
                ctx.upstream_stop = csq_result.upstream_stop;

                // Record consequence on the leaf node
                var merged_csq_rev = csq_result.csq_type;

                // SSS demotion (reverse strand): ibeg for reverse is 'i'
                if (csq_result.demote_to_sss or stack[i].node.?.payload == .sss) {
                    const node_splice_bits = stack[i].node.?.csq.toInt();
                    merged_csq_rev = CsqType.fromInt(merged_csq_rev.toInt() | node_splice_bits);
                    var csq_entry_sss: types.Csq = .{
                        .pos = stack[i].node.?.rec_pos,
                        .type_info = .{
                            .csq_type = merged_csq_rev,
                            .trid = tr_ptr.id,
                            .vcf_ial = @intCast(child_node.vcf_ial),
                            .gene = if (tr_ptr.gene) |g| blk: {
                                break :blk if (g.name) |n| std.mem.span(n) else null;
                            } else null,
                            .strand = tr_ptr.strand == .forward,
                            .biotype = @intFromEnum(tr_ptr.biotype),
                        },
                    };
                    _ = &csq_entry_sss;
                    try child_node.csq_list.append(allocator, csq_entry_sss);

                    ibeg_s = -1;
                    dlen_acc = 0;
                    indel_flag = false;
                    continue;
                }

                // Truncate tref/tseq at first stop for buildVstr (C lines 2203-2221)
                for (ctx.tref_stop.items, 0..) |ch, tidx| {
                    if (ch == '*') {
                        ctx.tref.shrinkRetainingCapacity(tidx + 1);
                        ctx.tref_stop.shrinkRetainingCapacity(tidx + 1);
                        break;
                    }
                }
                for (ctx.tseq_stop.items, 0..) |ch, tidx| {
                    if (ch == '*') {
                        ctx.tseq.shrinkRetainingCapacity(tidx + 1);
                        ctx.tseq_stop.shrinkRetainingCapacity(tidx + 1);
                        break;
                    }
                }

                var csq_entry: types.Csq = .{
                    .pos = stack[ibeg_u].node.?.rec_pos,
                    .type_info = .{
                        .csq_type = csq_result.csq_type,
                        .trid = tr_ptr.id,
                        .vcf_ial = @intCast(child_node.vcf_ial),
                        .gene = if (tr_ptr.gene) |g| blk: {
                            break :blk if (g.name) |n| std.mem.span(n) else null;
                        } else null,
                        .strand = tr_ptr.strand == .forward,
                        .biotype = @intFromEnum(tr_ptr.biotype),
                    },
                };
                try buildVstr(allocator, &csq_entry.type_info.vstr, stack, i, ibeg_u, ctx, sref_len, seq_m, tr_ptr, csq_result.csq_type);
                try child_node.csq_list.append(allocator, csq_entry);

                // For compound variants (i != ibeg_u), create CSQ_PRINTED_UPSTREAM
                // entries at non-ref-node positions (C csq.c lines 2334-2370).
                // ref_node = iend = ibeg_u for reverse strand.
                if (i != ibeg_u) {
                    const ref_node_pos_1based: u32 = stack[ibeg_u].node.?.rbeg + 1;
                    var j: usize = i;
                    while (j <= ibeg_u) : (j += 1) {
                        if (j == ibeg_u) continue; // skip ref_node
                        const node_j = stack[j].node.?;
                        var upstream_csq = types.CsqType{};
                        upstream_csq.printed_upstream = true;
                        const node_csq_raw = node_j.csq.toInt();
                        upstream_csq = types.CsqType.fromInt(upstream_csq.toInt() | node_csq_raw);
                        const upstream_entry: types.Csq = .{
                            .pos = node_j.rec_pos,
                            .ref_pos = ref_node_pos_1based,
                            .type_info = .{
                                .csq_type = upstream_csq,
                                .trid = tr_ptr.id,
                                .gene = if (tr_ptr.gene) |g| blk: {
                                    break :blk if (g.name) |n| std.mem.span(n) else null;
                                } else null,
                                .strand = tr_ptr.strand == .forward,
                                .biotype = @intFromEnum(tr_ptr.biotype),
                            },
                        };
                        try child_node.csq_list.append(allocator, upstream_entry);
                    }
                }

                ibeg_s = -1;
                dlen_acc = 0;
                indel_flag = false;
            }
        }
    }
}

/// Build the spliced reference for a transcript by concatenating CDS exons.
/// Corresponds to C function tscript_splice_ref().
fn tscriptSpliceRef(allocator: std.mem.Allocator, tscript_aux: *types.Tscript, tr: *const gff_types.Transcript) !void {
    const ref_seq = tscript_aux.ref_seq orelse return error.InvalidStrand;

    var total_len: usize = 0;
    for (tr.cds.items) |cds| total_len += cds.len;

    const sref_len = total_len + 2 * n_ref_pad;
    const sref_buf = try allocator.alloc(u8, sref_len);
    var pos: usize = 0;

    // Left padding
    const first_cds = tr.cds.items[0];
    const pad_start = first_cds.beg -| tr.beg;
    if (pad_start + n_ref_pad <= ref_seq.len)
        @memcpy(sref_buf[0..n_ref_pad], ref_seq[pad_start .. pad_start + n_ref_pad]);
    pos += n_ref_pad;

    // Copy each exon
    for (tr.cds.items) |cds| {
        const src_off = n_ref_pad + cds.beg -| tr.beg;
        if (src_off + cds.len <= ref_seq.len)
            @memcpy(sref_buf[pos .. pos + cds.len], ref_seq[src_off .. src_off + cds.len]);
        pos += cds.len;
    }

    // Right padding
    const last_cds = tr.cds.items[tr.cds.items.len - 1];
    const right_start = n_ref_pad + last_cds.beg -| tr.beg + last_cds.len;
    if (right_start + n_ref_pad <= ref_seq.len)
        @memcpy(sref_buf[pos .. pos + n_ref_pad], ref_seq[right_start .. right_start + n_ref_pad]);
    pos += n_ref_pad;

    tscript_aux.sref = sref_buf[0..pos];
    tscript_aux.nsref = @intCast(pos);
}

/// Helper: compute rlen + dlen for a node (used in soff calculation).
/// This is the signed sum rlen + dlen, clamped to 0 if negative.
/// For deletions (dlen < 0), the result is the alt allele length.
fn rlenPlusDlen(node: *const HapNode) usize {
    const sum: i64 = @as(i64, node.rlen) + @as(i64, node.dlen);
    return if (sum > 0) @intCast(sum) else 0;
}

/// The compound-consequence bitmask, matching the C CSQ_COMPOUND definition.
/// Duplicated here because CsqType.compound_mask is private to types.zig.
const compound_mask: u32 = (1 << 1) | // synonymous_variant
    (1 << 2) | // missense_variant
    (1 << 3) | // stop_lost
    (1 << 4) | // stop_gained
    (1 << 5) | // inframe_deletion
    (1 << 6) | // inframe_insertion
    (1 << 7) | // frameshift_variant
    (1 << 10) | // start_lost
    (1 << 12) | // stop_retained
    (1 << 18) | // inframe_altering
    (1 << 19) | // upstream_stop
    (1 << 20) | // incomplete_cds
    (1 << 22) | // elongation
    (1 << 23) | // truncation
    (1 << 24); // start_retained

/// Accumulate compound consequence flags from stack entries [ibeg..iend].
fn accumulateCompound(stack: []const Hstack, ibeg: usize, iend: usize) CsqType {
    var raw: u32 = 0;
    var i = ibeg;
    while (i <= iend) : (i += 1) {
        raw |= stack[i].node.?.csq.toInt() & compound_mask;
    }
    return CsqType.fromInt(raw);
}

/// Build the variant string (vstr) for a consequence: |aa_pos ref_aa>alt_aa|dna_variants
fn buildVstr(
    allocator: Allocator,
    vstr: *ArrayList(u8),
    stack: []const Hstack,
    ibeg: usize,
    iend: usize,
    ctx: *const HapContext,
    sref_len: usize,
    seq_m: usize,
    tr: *const gff_types.Transcript,
    csq_type: CsqType,
) !void {
    const rbeg_val = stack[ibeg].node.?.sbeg;
    const rend_val = stack[iend].node.?.sbeg + @as(u32, @intCast(@max(@as(i32, 0), stack[iend].node.?.rlen)));
    _ = rend_val;

    // node2soff(i) = stack[i].slen - (stack[i].node.rlen + stack[i].node.dlen)
    // node2sbeg(i) = ctx.sbeg + node2soff(i)
    // node2send(i) = ctx.sbeg + stack[i].slen
    const aa_rbeg: usize = if (tr.strand == .forward)
        rbeg_val / 3 + 1
    else
        (sref_len -| 2 * n_ref_pad -| (stack[iend].node.?.sbeg + @as(u32, @intCast(@max(@as(i32, 0), stack[iend].node.?.rlen))))) / 3 + 1;

    // aa_sbeg uses the spliced (alt) position, which accounts for indels
    const soff_ibeg = stack[ibeg].slen -| rlenPlusDlen(stack[ibeg].node.?);
    const send_iend = stack[iend].slen;
    const aa_sbeg: usize = if (tr.strand == .forward)
        (ctx.sbeg + soff_ibeg) / 3 + 1
    else
        (seq_m -| (ctx.sbeg + send_iend)) / 3 + 1;

    try vstr.append(allocator, '|');
    {
        var buf: [32]u8 = undefined;
        const s = std.fmt.bufPrint(&buf, "{d}", .{aa_rbeg}) catch return;
        try vstr.appendSlice(allocator, s);
    }
    try vstr.appendSlice(allocator, ctx.tref.items);
    // For synonymous variants, omit the ">alt_aa" part (C: csq.c line 2285)
    if (!csq_type.synonymous_variant) {
        try vstr.append(allocator, '>');
        {
            var buf: [32]u8 = undefined;
            const s = std.fmt.bufPrint(&buf, "{d}", .{aa_sbeg}) catch return;
            try vstr.appendSlice(allocator, s);
        }
        try vstr.appendSlice(allocator, ctx.tseq.items);
    }
    try vstr.append(allocator, '|');

    // DNA variant string: position + var for each node
    var first = true;
    var i = ibeg;
    while (i <= iend) : (i += 1) {
        if (!first) try vstr.append(allocator, '+');
        first = false;
        const n = stack[i].node.?;
        var buf: [32]u8 = undefined;
        const s = std.fmt.bufPrint(&buf, "{d}", .{n.rec_pos + 1}) catch return;
        try vstr.appendSlice(allocator, s);
        if (n.var_str) |vs| try vstr.appendSlice(allocator, vs);
    }
}

// ---------------------------------------------------------------------------
// hapFlush — flush completed transcripts
// ---------------------------------------------------------------------------
// Note: transcript flushing is implemented in CsqContext.flushTranscripts()
// (csq.zig) which owns the active-transcript priority queue and calls
// hapFinalize for each completed transcript.  This stub is retained for
// reference but is not called by the pipeline.
// ---------------------------------------------------------------------------

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

test "hapInit returns discarded for ref==alt" {
    const allocator = std.testing.allocator;

    // Set up a minimal transcript with one CDS exon
    var tr = gff_types.Transcript.init(allocator);
    defer tr.deinit();
    tr.id = 1;
    tr.beg = 100;
    tr.end = 108;
    tr.strand = .forward;

    // Build a reference: pad + ATGAAATTT (MKF) + pad
    const ref_str = "NNNNNNNNNN" ++ "ATGAAATTT" ++ "NNNNNNNNNN";
    var ref_buf: [ref_str.len]u8 = ref_str.*;

    var tscript_data: types.Tscript = .{
        .ref_seq = &ref_buf,
    };

    var cds_entry = gff_types.CdsEntry{
        .tr = &tr,
        .beg = 100,
        .pos = 0,
        .len = 9,
        .icds = 0,
        .phase = .phase0,
    };
    try tr.cds.append(allocator, &cds_entry);

    var root = HapNode.init(.root);
    var child = HapNode.init(.cds);

    // ref == alt: should be discarded
    const result = try hapInit(
        allocator,
        &root,
        &child,
        &cds_entry,
        103,
        "A",
        "A",
        1,
        &tscript_data,
    );

    try std.testing.expectEqual(HapInitResultKind.discarded, result.kind);
}

test "hapInit creates CDS node for coding SNP" {
    const allocator = std.testing.allocator;

    var tr = gff_types.Transcript.init(allocator);
    defer tr.deinit();
    tr.id = 1;
    tr.beg = 100;
    tr.end = 108;
    tr.strand = .forward;

    const ref_str = "NNNNNNNNNN" ++ "ATGAAATTT" ++ "NNNNNNNNNN";
    var ref_buf: [ref_str.len]u8 = ref_str.*;

    var tscript_data: types.Tscript = .{
        .ref_seq = &ref_buf,
    };

    var cds_entry = gff_types.CdsEntry{
        .tr = &tr,
        .beg = 100,
        .pos = 0,
        .len = 9,
        .icds = 0,
        .phase = .phase0,
    };
    try tr.cds.append(allocator, &cds_entry);

    var root = HapNode.init(.root);
    var child = HapNode.init(.cds);
    defer {
        // Free allocations made by hapInit
        if (child.var_str) |vs| allocator.free(vs);
        if (child.payload == .cds) {
            if (child.payload.cds.seq) |seq| allocator.free(seq);
        }
    }

    // SNP inside the exon at position 103: A>T
    const result = try hapInit(
        allocator,
        &root,
        &child,
        &cds_entry,
        103,
        "A",
        "T",
        1,
        &tscript_data,
    );

    try std.testing.expectEqual(HapInitResultKind.added, result.kind);
    try std.testing.expect(child.payload == .cds);
    try std.testing.expectEqual(@as(i32, 0), child.dlen); // SNP: no length change
    try std.testing.expect(child.var_str != null);
    try std.testing.expectEqualStrings("A>T", child.var_str.?);
}

test "hapFinalize produces consequence for simple 2-node tree" {
    const allocator = std.testing.allocator;
    const gencode = translate.findGeneticCode(0) orelse unreachable;

    // Build a minimal transcript: one exon, ATGAAATTT (M K F)
    var tr = gff_types.Transcript.init(allocator);
    defer tr.deinit();
    tr.id = 1;
    tr.beg = 100;
    tr.end = 108;
    tr.strand = .forward;
    tr.biotype = .protein_coding;

    const ref_str = "NNNNNNNNNN" ++ "ATGAAATTT" ++ "NNNNNNNNNN";
    var ref_buf: [ref_str.len]u8 = ref_str.*;

    // Create root node
    var root = HapNode.init(.root);
    defer root.deinit(allocator);

    var tscript_data: types.Tscript = .{
        .ref_seq = &ref_buf,
        .root = &root,
    };
    defer if (tscript_data.sref) |s| allocator.free(s);
    tr.aux = @ptrCast(&tscript_data);

    var cds_entry = gff_types.CdsEntry{
        .tr = &tr,
        .beg = 100,
        .pos = 0,
        .len = 9,
        .icds = 0,
        .phase = .phase0,
    };
    try tr.cds.append(allocator, &cds_entry);

    // CDS child: SNP at spliced pos 3, changing codon AAA -> TAA (K -> stop)
    var child_node = HapNode.init(.cds);
    child_node.payload = .{ .cds = .{ .seq = @constCast("T") } };
    child_node.var_str = "A>T";
    child_node.dlen = 0;
    child_node.rbeg = 3;
    child_node.rlen = 1;
    child_node.sbeg = 3;
    child_node.icds = 0;
    child_node.vcf_ial = 1;
    child_node.nend = 1; // leaf
    child_node.prev = &root;

    try root.children.append(allocator, &child_node);

    // Create HapContext and run finalize
    var ctx = HapContext.init(allocator, gencode);
    defer ctx.deinit();
    ctx.tr = @ptrCast(&tr);

    try hapFinalize(&ctx);

    // The child node should have at least one consequence
    try std.testing.expect(child_node.csq_list.items.len > 0);

    // Clean up csq_list
    for (child_node.csq_list.items) |*c| c.deinit(allocator);
    child_node.csq_list.deinit(allocator);
}

test "hapFinalize merges compound variants in same codon" {
    const allocator = std.testing.allocator;
    const gencode = translate.findGeneticCode(0) orelse unreachable;

    // Transcript with one exon: ATGAAATTTCCC (M K F P)
    var tr = gff_types.Transcript.init(allocator);
    defer tr.deinit();
    tr.id = 1;
    tr.beg = 100;
    tr.end = 111;
    tr.strand = .forward;
    tr.biotype = .protein_coding;

    const ref_str = "NNNNNNNNNN" ++ "ATGAAATTTCCC" ++ "NNNNNNNNNN";
    var ref_buf: [ref_str.len]u8 = ref_str.*;

    var root = HapNode.init(.root);
    defer root.deinit(allocator);

    var tscript_data: types.Tscript = .{
        .ref_seq = &ref_buf,
        .root = &root,
    };
    defer if (tscript_data.sref) |s| allocator.free(s);
    tr.aux = @ptrCast(&tscript_data);

    var cds_entry = gff_types.CdsEntry{
        .tr = &tr,
        .beg = 100,
        .pos = 0,
        .len = 12,
        .icds = 0,
        .phase = .phase0,
    };
    try tr.cds.append(allocator, &cds_entry);

    // Two SNPs in same codon (AAA at positions 3,4): pos 3 A>T, pos 4 A>T
    var child1 = HapNode.init(.cds);
    child1.payload = .{ .cds = .{ .seq = @constCast("T") } };
    child1.var_str = "A>T";
    child1.dlen = 0;
    child1.rbeg = 3;
    child1.rlen = 1;
    child1.sbeg = 3;
    child1.icds = 0;
    child1.vcf_ial = 1;
    child1.nend = 0; // not a leaf
    child1.prev = &root;
    try root.children.append(allocator, &child1);

    var child2 = HapNode.init(.cds);
    child2.payload = .{ .cds = .{ .seq = @constCast("T") } };
    child2.var_str = "A>T";
    child2.dlen = 0;
    child2.rbeg = 4;
    child2.rlen = 1;
    child2.sbeg = 4;
    child2.icds = 0;
    child2.vcf_ial = 1;
    child2.nend = 1; // leaf
    child2.prev = &child1;
    try child1.children.append(allocator, &child2);

    var ctx = HapContext.init(allocator, gencode);
    defer ctx.deinit();
    ctx.tr = @ptrCast(&tr);

    try hapFinalize(&ctx);

    // child2 (the leaf) should have consequence(s) that merge both variants
    try std.testing.expect(child2.csq_list.items.len > 0);

    // The consequence vstr should contain '+' indicating merged DNA variant entries
    var found_compound = false;
    for (child2.csq_list.items) |c| {
        for (c.type_info.vstr.items) |ch| {
            if (ch == '+') {
                found_compound = true;
                break;
            }
        }
        if (found_compound) break;
    }
    try std.testing.expect(found_compound);

    // Clean up: HapNode.deinit only frees its own arrays, not children's
    for (child2.csq_list.items) |*c| c.deinit(allocator);
    child2.csq_list.deinit(allocator);
    child1.children.deinit(allocator);
    child1.csq_list.deinit(allocator);
    child1.cur_child.deinit(allocator);
}
