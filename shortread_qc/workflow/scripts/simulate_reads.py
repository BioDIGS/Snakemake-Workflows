#!/usr/bin/env python3
"""
simulate_reads.py — simulate Illumina-style paired-end reads from a FASTA,
including the things that make QC hard: adapter read-through when the insert
is shorter than the read, poly-G fill after the adapter (2-colour chemistry),
and a small substitution error rate.

Used by the workflow's optional validation step to build a human spike-in from
the same reference the host filter uses, so leakage can be measured directly.
Every simulated read carries a tag in the flowcell field of its name
(INSTRUMENT:RUN:<TAG>:LANE:TILE:X:Y) that survives trimming and host removal.

Reads the FASTA through a .fai index (byte offsets), so a 3 Gb reference is
not loaded into memory. If the .fai does not exist it is built once (a single
streaming pass) at --fai.

Usage:
  simulate_reads.py --fasta hs1.fa --pairs 50000 --out-r1 s_R1.fq.gz --out-r2 s_R2.fq.gz
                    [--fai hs1.fa.fai] [--tag HSSPIKE] [--read-len 150]
                    [--insert-mean 280 --insert-sd 90] [--adapter-r1 SEQ --adapter-r2 SEQ]
                    [--instrument A00000] [--seed 1]
"""
import argparse
import bisect
import gzip
import os
import random

TRUSEQ_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRUSEQ_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def build_fai(fasta, fai):
    """Write a samtools-compatible .fai (name, length, offset, linebases, linewidth)."""
    entries = []
    with open(fasta, "rb") as fh:
        name = None
        length = offset = lb = lw = 0
        pos = 0
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    entries.append((name, length, offset, lb, lw))
                name = line[1:].split()[0].decode()
                length = 0
                offset = pos + len(line)
                lb = lw = 0
            else:
                if lb == 0:
                    lb = len(line.rstrip(b"\r\n"))
                    lw = len(line)
                length += len(line.rstrip(b"\r\n"))
            pos += len(line)
        if name is not None:
            entries.append((name, length, offset, lb, lw))
    with open(fai, "w") as out:
        for e in entries:
            out.write("\t".join(map(str, e)) + "\n")


def load_fai(fai):
    recs = []
    with open(fai) as fh:
        for line in fh:
            name, length, offset, lb, lw = line.split("\t")[:5]
            recs.append((name, int(length), int(offset), int(lb), int(lw)))
    return recs


def fetch(fh, rec, start, length):
    """Fetch [start, start+length) from a contig using .fai geometry."""
    _, clen, off, lb, lw = rec
    end = min(start + length, clen)
    first_line, first_col = divmod(start, lb)
    byte_start = off + first_line * lw + first_col
    last_line, last_col = divmod(end, lb)
    byte_end = off + last_line * lw + last_col
    fh.seek(byte_start)
    raw = fh.read(byte_end - byte_start)
    return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def mutate(seq, rate, rng):
    if rate <= 0:
        return seq
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def make_read(fragment_strand, adapter, read_len, rng, err, lowq=0.03, lowq_3p=0.0):
    """Read = fragment, then adapter, then poly-G fill (2-colour 'no signal')."""
    seq = fragment_strand[:read_len]
    if len(seq) < read_len:
        seq = seq + adapter
        seq = seq[:read_len]
        if len(seq) < read_len:
            seq = seq + "G" * (read_len - len(seq))
    seq = mutate(seq, err, rng)
    ins = min(len(fragment_strand), read_len)
    qual = []
    for i in range(read_len):
        if i >= ins + len(adapter):           # poly-G fill region: low quality
            qual.append("#")
        elif rng.random() < lowq + lowq_3p * (i / read_len) ** 2:
            qual.append(",")                  # Q11 (NovaSeq binned); lowq_3p adds 3' decay
        elif rng.random() < 0.05:
            qual.append(":")                  # Q25
        else:
            qual.append("F")                  # Q37
    return seq, "".join(qual)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fai", default=None, help="index path (default: FASTA.fai; built if missing)")
    ap.add_argument("--pairs", type=int, default=50000)
    ap.add_argument("--out-r1", required=True)
    ap.add_argument("--out-r2", required=True)
    ap.add_argument("--tag", default="HSSPIKE", help="flowcell-field tag marking simulated reads")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--insert-mean", type=float, default=280)
    ap.add_argument("--insert-sd", type=float, default=90)
    ap.add_argument("--min-insert", type=int, default=35)
    ap.add_argument("--adapter-r1", default=TRUSEQ_R1)
    ap.add_argument("--adapter-r2", default=TRUSEQ_R2)
    ap.add_argument("--error-rate", type=float, default=0.002)
    ap.add_argument("--instrument", default="A00000", help="instrument ID (A... = NovaSeq 6000)")
    ap.add_argument("--max-n-frac", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--flowcell", default=None, help="flowcell field (default: --tag)")
    ap.add_argument("--run", type=int, default=1)
    ap.add_argument("--lane", type=int, default=1)
    ap.add_argument("--index", default="ACGTACGT+TGCATGCA")
    ap.add_argument("--lowq-rate", type=float, default=0.03, help="baseline fraction of Q11 bases")
    ap.add_argument("--lowq-3prime", type=float, default=0.0, help="extra Q11 fraction reached at the 3' end")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    fai = args.fai or args.fasta + ".fai"
    if not os.path.exists(fai):
        build_fai(args.fasta, fai)
    recs = [r for r in load_fai(fai) if r[1] >= 1000]
    if not recs:
        raise SystemExit("ERROR: no contigs >= 1 kb in the reference")
    cum = []
    tot = 0
    for r in recs:
        tot += r[1]
        cum.append(tot)

    n_written = 0
    with open(args.fasta, "rb") as fa, gzip.open(args.out_r1, "wt") as o1, gzip.open(args.out_r2, "wt") as o2:
        attempts = 0
        while n_written < args.pairs:
            attempts += 1
            if attempts > args.pairs * 50:
                raise SystemExit("ERROR: too many rejected fragments (reference mostly N?)")
            ins = max(args.min_insert, int(rng.gauss(args.insert_mean, args.insert_sd)))
            k = bisect.bisect_left(cum, rng.random() * tot)
            rec = recs[k]
            if rec[1] <= ins:
                continue
            start = rng.randrange(0, rec[1] - ins)
            frag = fetch(fa, rec, start, ins)
            if len(frag) < ins or frag.count("N") > args.max_n_frac * ins:
                continue
            if rng.random() < 0.5:
                frag = revcomp(frag)
            s1, q1 = make_read(frag, args.adapter_r1, args.read_len, rng, args.error_rate,
                               args.lowq_rate, args.lowq_3prime)
            s2, q2 = make_read(revcomp(frag), args.adapter_r2, args.read_len, rng, args.error_rate,
                               args.lowq_rate, args.lowq_3prime * 1.5)
            n_written += 1
            fc = args.flowcell or args.tag
            name = f"{args.instrument}:{args.run}:{fc}:{args.lane}:1101:{n_written}:{ins}"
            o1.write(f"@{name} 1:N:0:{args.index}\n{s1}\n+\n{q1}\n")
            o2.write(f"@{name} 2:N:0:{args.index}\n{s2}\n+\n{q2}\n")
    print(f"wrote {n_written} pairs (tag {args.tag}) to {args.out_r1} / {args.out_r2}")


if __name__ == "__main__":
    main()
