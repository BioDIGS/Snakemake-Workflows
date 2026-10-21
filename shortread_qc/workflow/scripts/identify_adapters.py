#!/usr/bin/env python3
"""
identify_adapters.py — figure out what adapters (and which instrument) a
paired-end Illumina-style FASTQ pair actually carries, before trimming.

Three independent lines of evidence:

  1. Instrument / chemistry, read from the FASTQ headers (decides whether
     poly-G tails from 2-colour chemistry are expected).
  2. Known-adapter signatures: how many reads contain the first bases of
     TruSeq, Nextera/Tn5, small-RNA and MGI/DNBSEQ adapters, and where.
  3. De novo adapter inference from read-pair overlap: for pairs whose insert
     is shorter than the read, R1 and reverse-complement(R2) overlap exactly,
     and whatever follows the insert in each read IS the adapter. A per-position
     consensus of those tails gives the adapter sequence without assuming a kit.

It also reports whether reads look already trimmed (variable read lengths),
poly-G tail frequency, adapter-dimer frequency and an insert-size summary,
and prints the fastp options that match what it found.

Usage:
  identify_adapters.py R1.fastq.gz R2.fastq.gz [--pairs 200000] [--skip 10000]
                       [--tsv out.tsv] [--txt out.txt] [--sample NAME]
"""
import argparse
import gzip
import re
import statistics
import sys
from collections import Counter

# ---------------------------------------------------------------------------
# Reference adapter sequences (read-through orientation, as seen in the read)
# ---------------------------------------------------------------------------
KNOWN_ADAPTERS = {
    "Illumina TruSeq R1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "Illumina TruSeq R2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "Nextera / Tn5 (R1 & R2)": "CTGTCTCTTATACACATCT",
    "Illumina small RNA 3'": "TGGAATTCTCGGGTGCCAAGG",
    "MGI/DNBSEQ R1": "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
    "MGI/DNBSEQ R2": "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG",
}
# Short signatures used for counting reads that contain each adapter family.
SIGNATURES = {
    "TruSeq (AGATCGGAAGAGC)": "AGATCGGAAGAGC",
    "Nextera/Tn5 (CTGTCTCTTATACACATCT)": "CTGTCTCTTATACACATCT",
    "small RNA (TGGAATTCTCGG)": "TGGAATTCTCGG",
    "MGI R1 (AAGTCGGAGGCC)": "AAGTCGGAGGCC",
    "MGI R2 (AAGTCGGATCGT)": "AAGTCGGATCGT",
}

# Illumina instrument-ID prefixes (heuristic; first field of the read header).
# Longest prefixes first. two_colour=True means poly-G tails are expected.
INSTRUMENTS = [
    ("HWI-", "HiSeq 2000/2500 (older ID format)", False),
    ("LH", "NovaSeq X / X Plus", True),
    ("VH", "NextSeq 1000/2000", True),
    ("VL", "NextSeq 1000/2000", True),
    ("NB", "NextSeq 500/550", True),
    ("NS", "NextSeq 500/550", True),
    ("MN", "MiniSeq", True),
    ("FS", "iSeq 100 (one-channel chemistry)", False),
    ("SN", "HiSeq 2000/2500", False),
    ("A", "NovaSeq 6000", True),
    ("M", "MiSeq", False),
    ("D", "HiSeq 2500", False),
    ("C", "HiSeq 1500/2500", False),
    ("J", "HiSeq 3000", False),
    ("K", "HiSeq 4000", False),
    ("E", "HiSeq X", False),
]
MGI_HEADER = re.compile(r"^[A-Z]\d{9,}L\d")  # e.g. V350012345L1C001R0010000001
SRA_HEADER = re.compile(r"^[SED]RR\d+")

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def open_fq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_pairs(p1, p2, n, skip):
    """Yield (header1, seq1, seq2) for up to n pairs after skipping `skip`."""
    with open_fq(p1) as f1, open_fq(p2) as f2:
        i = 0
        while True:
            h1 = f1.readline()
            if not h1:
                return
            s1 = f1.readline().rstrip()
            f1.readline()
            f1.readline()
            h2 = f2.readline()
            s2 = f2.readline().rstrip()
            f2.readline()
            f2.readline()
            if not h2:
                sys.exit("ERROR: R2 has fewer records than R1 — files are not a matched pair.")
            i += 1
            if i <= skip:
                continue
            if i > skip + n:
                return
            yield h1.rstrip(), s1.upper(), s2.upper()


def guess_instrument(header):
    first = header[1:].split()[0]
    inst_id = first.split(":")[0]
    if SRA_HEADER.match(first):
        return inst_id, "unknown (SRA-style header; instrument ID stripped)", None
    if MGI_HEADER.match(first):
        return inst_id, "MGI/DNBSEQ (BGI) — not Illumina; MGI adapters apply", False
    for prefix, name, two_colour in INSTRUMENTS:
        if inst_id.startswith(prefix):
            return inst_id, name, two_colour
    return inst_id, "unknown", None


def overlap_insert(r1, r2, seed_len=16):
    """If the pair overlaps, return the insert length L (bases of genomic
    fragment); R1[:L] == revcomp(R2)[-L:]. Returns None if no overlap found."""
    seed = r1[:seed_len]
    if len(seed) < seed_len or "N" in seed:
        return None
    rc2 = revcomp(r2)
    p = rc2.find(seed)
    if p < 0:
        return None
    L = len(rc2) - p
    a, b = r1[:L], rc2[p:]
    n = min(len(a), len(b))
    mism = sum(1 for x, y in zip(a[:n], b[:n]) if x != y and x != "N" and y != "N")
    if mism > max(2, int(0.1 * n)):
        return None
    return L


def consensus(tails, max_len=40, min_support=20, min_frac=0.6):
    """Per-position majority consensus of adapter tails."""
    out = []
    support = []
    for i in range(max_len):
        col = Counter(t[i] for t in tails if len(t) > i and t[i] != "N")
        tot = sum(col.values())
        if tot < min_support:
            break
        base, cnt = col.most_common(1)[0]
        if cnt / tot < min_frac:
            break
        out.append(base)
        support.append(cnt / tot)
    seq = "".join(out)
    # Past the end of the adapter, 2-colour instruments read "no signal" = G.
    # Strip that poly-G run so it is not mistaken for adapter sequence.
    stripped = re.sub(r"G{4,}$", "", seq)
    return stripped, support[:len(stripped)]


def best_known_match(seq):
    """Identity of consensus vs each known adapter over the shared prefix."""
    if len(seq) < 8:
        return None, 0.0, 0
    best = (None, 0.0, 0)
    for name, ref in KNOWN_ADAPTERS.items():
        n = min(len(seq), len(ref))
        ident = sum(1 for a, b in zip(seq[:n], ref[:n]) if a == b) / n
        if ident > best[1] or (ident == best[1] and n > best[2]):
            best = (name, ident, n)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("r1")
    ap.add_argument("r2")
    ap.add_argument("--pairs", type=int, default=200000, help="read pairs to inspect (default 200000)")
    ap.add_argument("--skip", type=int, default=10000,
                    help="skip this many leading pairs (first tiles are often atypical; default 10000)")
    ap.add_argument("--sample", default=None, help="sample name for the report")
    ap.add_argument("--tsv", default=None, help="write key<TAB>value summary here")
    ap.add_argument("--txt", default=None, help="write the human-readable report here (also printed)")
    args = ap.parse_args()

    n = 0
    first_header = None
    len1, len2 = Counter(), Counter()
    sig_hits = {k: [0, 0] for k in SIGNATURES}          # reads with signature (R1, R2)
    sig_pos = {k: [] for k in SIGNATURES}
    dimers = 0
    polyg1 = polyg2 = 0
    overlaps = 0
    readthrough = 0
    inserts = []
    tails1, tails2 = [], []

    for h1, s1, s2 in read_pairs(args.r1, args.r2, args.pairs, args.skip):
        if first_header is None:
            first_header = h1
        n += 1
        len1[len(s1)] += 1
        len2[len(s2)] += 1
        if s1.endswith("G" * 10):
            polyg1 += 1
        if s2.endswith("G" * 10):
            polyg2 += 1
        dimer_here = False
        for k, sig in SIGNATURES.items():
            for j, s in enumerate((s1, s2)):
                p = s.find(sig)
                if p >= 0:
                    sig_hits[k][j] += 1
                    sig_pos[k].append(p)
                    if p < 10:
                        dimer_here = True
        dimers += dimer_here

        L = overlap_insert(s1, s2)
        if L is not None:
            overlaps += 1
            inserts.append(L)
            if L < len(s1) or L < len(s2):
                readthrough += 1
                tails1.append(s1[L:])
                tails2.append(s2[L:])

    if n == 0:
        sys.exit(f"ERROR: no read pairs read (file has <= {args.skip} pairs? try --skip 0)")

    inst_id, inst_name, two_colour = guess_instrument(first_header)
    max1, max2 = max(len1), max(len2)
    frac_full1 = len1[max1] / n
    frac_full2 = len2[max2] / n
    already_trimmed = frac_full1 < 0.95 or frac_full2 < 0.95

    cons1, sup1 = consensus(tails1)
    cons2, sup2 = consensus(tails2)
    m1 = best_known_match(cons1)
    m2 = best_known_match(cons2)

    # Which signature family dominates?
    fam = max(SIGNATURES, key=lambda k: sum(sig_hits[k]))
    fam_reads = sum(sig_hits[fam])

    def pct(x, d=n):
        return 100.0 * x / d if d else 0.0

    sample = args.sample or args.r1.split("/")[-1]
    lines = []
    w = lines.append
    w(f"Adapter identification report — {sample}")
    w("=" * 72)
    w(f"Read pairs inspected        : {n:,} (after skipping {args.skip:,})")
    w(f"First header                : {first_header[:70]}")
    w(f"Instrument ID / guess       : {inst_id}  ->  {inst_name}")
    if two_colour is True:
        w("Chemistry                   : 2-colour — poly-G tails expected; keep poly-G trimming ON")
    elif two_colour is False:
        w("Chemistry                   : not 2-colour Illumina — poly-G trimming optional")
    else:
        w("Chemistry                   : unknown — check run metadata with the provider")
    w("")
    w("Read lengths")
    w(f"  R1 max {max1} bp, {pct(len1[max1]):.1f}% of reads at max;  "
      f"R2 max {max2} bp, {pct(len2[max2]):.1f}% at max")
    if already_trimmed:
        w("  NOTE: read lengths vary — these reads look ALREADY trimmed (e.g. by the")
        w("        provider / BaseSpace / an earlier step). Remaining adapter is residual.")
    w("")
    w("Known-adapter signatures (reads containing the signature)")
    for k in SIGNATURES:
        r1h, r2h = sig_hits[k]
        if r1h or r2h:
            med = statistics.median(sig_pos[k]) if sig_pos[k] else float("nan")
            f1 = f"{pct(r1h):6.2f}%" if pct(r1h) >= 0.01 or not r1h else " <0.01%"
            f2 = f"{pct(r2h):6.2f}%" if pct(r2h) >= 0.01 or not r2h else " <0.01%"
            w(f"  {k:<36} R1 {f1}  R2 {f2}  median pos {med:.0f}")
        else:
            w(f"  {k:<36} not seen")
    w(f"  Adapter dimers / near-empty inserts (adapter starts < 10 bp): {pct(dimers):.3f}% of pairs")
    w("")
    w("Poly-G tails (>=10 G at 3' end)")
    w(f"  R1 {pct(polyg1):.2f}%   R2 {pct(polyg2):.2f}%")
    w("")
    w("Overlap-based (kit-agnostic) inference")
    w(f"  Pairs with detectable R1/R2 overlap : {pct(overlaps):.1f}%")
    w(f"  Pairs reading through into adapter  : {pct(readthrough):.1f}%  (insert < read length)")
    if inserts:
        q = statistics.quantiles(inserts, n=4) if len(inserts) > 3 else [inserts[0]] * 3
        w(f"  Insert size among overlapping pairs : median {statistics.median(inserts):.0f} bp "
          f"(IQR {q[0]:.0f}-{q[2]:.0f}); longer inserts do not overlap and are not counted")
    for lab, cons, sup, m in (("R1", cons1, sup1, m1), ("R2", cons2, sup2, m2)):
        if cons:
            w(f"  {lab} adapter consensus ({len(cons)} bp, min support {min(sup):.2f}): {cons}")
            if m[0] and m[1] >= 0.9:
                w(f"     -> matches {m[0]} ({100*m[1]:.0f}% identity over {m[2]} bp)")
            else:
                w(f"     -> no confident match to a known kit (best: {m[0]}, {100*m[1]:.0f}%) "
                  f"— use --adapter_sequence with this consensus")
        else:
            w(f"  {lab} adapter consensus: not enough read-through pairs to call")
    w("")
    w("Suggested fastp adapter options")
    def pick(cons, m):
        # Prefer the full known sequence when the consensus confidently matches it.
        if m[0] and m[1] >= 0.9:
            return KNOWN_ADAPTERS[m[0]]
        return cons

    if cons1 and cons2 and len(cons1) >= 12 and len(cons2) >= 12:
        w(f"  --detect_adapter_for_pe --adapter_sequence {pick(cons1, m1)} "
          f"--adapter_sequence_r2 {pick(cons2, m2)}")
        w("  (explicit sequences make the run independent of fastp's auto-detection;")
        w("   pin them in config.yaml under fastp: adapter_r1 / adapter_r2)")
    else:
        w("  --detect_adapter_for_pe   (too little read-through to pin an explicit sequence)")
    if two_colour:
        w("  --trim_poly_g             (auto-enabled by fastp for NovaSeq/NextSeq IDs; set explicitly)")
    report = "\n".join(lines)
    print(report)
    if args.txt:
        with open(args.txt, "w") as fh:
            fh.write(report + "\n")

    if args.tsv:
        kv = [
            ("sample", sample),
            ("pairs_inspected", n),
            ("instrument_id", inst_id),
            ("instrument_guess", inst_name),
            ("two_colour", two_colour),
            ("r1_max_len", max1),
            ("r2_max_len", max2),
            ("pct_r1_full_length", round(pct(len1[max1]), 2)),
            ("looks_already_trimmed", already_trimmed),
            ("dominant_adapter_family", fam if fam_reads else "none"),
            ("pct_reads_with_dominant_signature", round(pct(fam_reads, 2 * n), 3)),
            ("pct_pairs_readthrough", round(pct(readthrough), 3)),
            ("pct_adapter_dimer", round(pct(dimers), 4)),
            ("pct_r1_polyG", round(pct(polyg1), 3)),
            ("pct_r2_polyG", round(pct(polyg2), 3)),
            ("median_insert_overlapping", statistics.median(inserts) if inserts else "NA"),
            ("adapter_r1_consensus", cons1 or "NA"),
            ("adapter_r1_match", f"{m1[0]} ({100*m1[1]:.0f}%)" if cons1 else "NA"),
            ("adapter_r2_consensus", cons2 or "NA"),
            ("adapter_r2_match", f"{m2[0]} ({100*m2[1]:.0f}%)" if cons2 else "NA"),
        ]
        with open(args.tsv, "w") as fh:
            fh.write("key\tvalue\n")
            for k, v in kv:
                fh.write(f"{k}\t{v}\n")


if __name__ == "__main__":
    main()
