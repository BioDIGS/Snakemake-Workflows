#!/usr/bin/env python3
"""
rename_samples.py — give samples.tsv real sample IDs from your mapping CSV, and
write a metadata table (site, sample ID, recorded instrument) for the sample sheet.

Mapping CSV columns used (same file as populate_data_truth.sh):
  raw_id, site_id, sample_id, sequencing_instrument

Matching (exact only; your old prefix fallback is reported, never applied, because
"A01_10" starts with "A01_1"):
  1. raw_id == the provider stem, e.g. "A01_10" or "HHGCVDRX3-1-AACATCGCGC-ACAAGTGGAC"
  2. raw_id == stem + "_S<n>"          (e.g. "A01_10_S2")
  3. raw_id == stem + "_S<n>_L<lane>"  (any of the sample's lanes)

New name: the plain sample_id, e.g. A01_10 (like the long-read config). If the
same sample_id turns up in more than one delivery, ONLY those samples keep their
delivery suffix (A01_10__20230925, A01_10__20250219) so they stay separate, and a
NOTE lists them. --keep-delivery-suffix keeps the suffix on every sample.

Checks: unmapped samples, raw_ids matching more than one row, two samples mapping
to the same new name, and recorded instrument vs the instrument in each file's
first read header (catches mapping mistakes).

Usage:
  rename_samples.py --samples config/samples.tsv --mapping mapping.csv \\
      --out-samples config/samples.renamed.tsv --out-metadata config/sample_metadata.tsv \\
      [--keep-unmapped]
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from identify_adapters import guess_instrument  # noqa: E402

FILE_PAT = re.compile(r"^(?P<stem>.+?)(?P<s>_S\d+)?(?P<lane>_L\d{3})?[._](?:R)?[12](?:_\d{3})?\.f(?:ast)?q(?:\.gz)?$")
INSTRUMENT_WORDS = {"novaseq 6000": "NovaSeq 6000", "novaseq x": "NovaSeq X"}


def first_header(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return fh.readline().strip()


def norm_inst(s):
    s = (s or "").strip().lower()
    for k, v in INSTRUMENT_WORDS.items():
        if s.startswith(k):
            return v
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--mapping", required=True)
    ap.add_argument("--out-samples", required=True)
    ap.add_argument("--out-metadata", required=True)
    ap.add_argument("--keep-delivery-suffix", action="store_true",
                    help="name every sample <sample_id>__<delivery> instead of plain <sample_id>")
    ap.add_argument("--strip-prefix", default="",
                    help="remove this prefix from sample_id for the new name, e.g. Molecular_")
    ap.add_argument("--instrument-mismatch", choices=["error", "warn"], default="error",
                    help="mapping instrument disagrees with the read header: stop (error) or "
                         "just report (warn). The header's instrument is always recorded.")
    ap.add_argument("--drop-columns", default="filename,size_GB,md5sum,reads,Note",
                    help="mapping columns NOT copied into the metadata (per-file columns); "
                         "all other columns are copied")
    ap.add_argument("--keep-unmapped", action="store_true",
                    help="keep unmapped samples under their current names instead of stopping")
    a = ap.parse_args()

    with open(a.mapping, newline="") as fh:
        rows = list(csv.DictReader(fh))
    need = {"raw_id", "site_id", "sample_id"}
    if not rows or not need <= set(rows[0]):
        sys.exit(f"ERROR: mapping CSV needs columns {sorted(need)}; has {list(rows[0]) if rows else []}")
    by_raw = defaultdict(list)
    for r in rows:
        by_raw[r["raw_id"].strip()].append(r)
    dup_raw = {k for k, v in by_raw.items() if len(v) > 1 and k}
    # Fallback index: '-' and '_' treated as the same (the master sheet writes
    # HHGCVDRX3_1_<i7>-<i5>_S1 where the FASTQ is HHGCVDRX3-1-<i7>-<i5>_S1_L001)
    sep = lambda x: re.sub(r"[-_]", "_", x)
    by_raw_sep = defaultdict(list)
    for k, v in by_raw.items():
        by_raw_sep[sep(k)].append(k)

    with open(a.samples) as fh:
        samples = [r for r in csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")),
                                             delimiter="\t")]

    out, meta, problems, unmapped, near, notes, warns = [], [], [], [], [], [], []
    drop = {c.strip() for c in a.drop_columns.split(",") if c.strip()}
    n_sep = 0      # matched only after treating '-' and '_' as equal
    n_dup_ok = 0   # raw_id on several rows that agree (e.g. one row per R1/R2 file)
    for s in samples:
        old = s["sample"].strip()
        base, _, suffix = old.partition("__")
        r1_files = [p.strip() for p in s["r1"].split(",") if p.strip()]
        cands = [base]
        for f in r1_files:
            m = FILE_PAT.match(os.path.basename(f))
            if m:
                st = m.group("stem")
                if m.group("s"):
                    cands.append(st + m.group("s"))
                    if m.group("lane"):
                        cands.append(st + m.group("s") + m.group("lane"))
        hit = None
        cands = list(dict.fromkeys(cands))
        # pass 1: exact raw_id; pass 2: same candidates with '-' and '_' treated as equal
        c = next((x for x in cands if x in by_raw), None)
        if c is None:
            loose = {by_raw_sep[sep(x)][0] for x in cands if len(by_raw_sep.get(sep(x), [])) == 1}
            if len(loose) == 1:
                c = loose.pop()
                n_sep += 1
            elif len(loose) > 1:
                problems.append(f"{old}: several raw_ids match when ignoring '-' vs '_': {sorted(loose)}")
        if c is not None:
            if c in dup_raw:
                ids = {(r["sample_id"].strip(), r["site_id"].strip()) for r in by_raw[c]}
                if len(ids) == 1:
                    n_dup_ok += 1
                else:
                    problems.append(f"{old}: raw_id '{c}' is in {len(by_raw[c])} mapping rows that DISAGREE: "
                                    + "; ".join(f"{x} (site {y})" for x, y in sorted(ids)))
            hit = by_raw[c][0]
        if hit is None:
            pref = [k for k in by_raw if k and (base.startswith(k) or k.startswith(base))]
            if pref:
                near.append(f"{old}: no exact raw_id; prefix-only candidates (NOT applied): {', '.join(sorted(pref)[:5])}")
            unmapped.append(old)
            if a.keep_unmapped:
                out.append({**s, "sample": old})
                meta.append({"sample": old, "sample_id": "", "site_id": "", "raw_id": "",
                             "sequencing_instrument_recorded": "", "instrument_in_header": "",
                             "delivery_suffix": suffix, "provider_name": base})
            continue

        sid = hit["sample_id"].strip()
        if a.strip_prefix and sid.startswith(a.strip_prefix):
            sid = sid[len(a.strip_prefix):]
        new = re.sub(r"[^A-Za-z0-9._-]", "_", sid)
        if a.keep_delivery_suffix and suffix:
            new = f"{new}__{suffix}"
        _, hdr_model, _ = guess_instrument(first_header(r1_files[0]))
        rec = norm_inst(hit.get("sequencing_instrument", ""))
        hdr = "NovaSeq X" if hdr_model.startswith("NovaSeq X") else hdr_model
        if rec and hdr != "unknown" and rec != hdr and rec.lower() not in hdr.lower():
            (warns if a.instrument_mismatch == "warn" else problems).append(
                f"{old} -> {new}: mapping says '{hit.get('sequencing_instrument')}', "
                            f"read header says '{hdr_model}' — wrong mapping row?")
        out.append({**s, "sample": new})
        row_meta = {"sample": new, "sample_id": hit["sample_id"].strip(), "site_id": hit["site_id"].strip(),
                    "raw_id": hit["raw_id"].strip(),
                    "sequencing_instrument_recorded": hit.get("sequencing_instrument", "").strip(),
                    "instrument_in_header": hdr_model, "delivery_suffix": suffix, "provider_name": base}
        for k, v in hit.items():
            if k and k not in drop and k not in row_meta and k != "sequencing_instrument":
                row_meta[k] = (v or "").strip()
        meta.append(row_meta)

    # Same sample_id in more than one delivery: suffix only those, so they stay separate
    counts = Counter(r["sample"] for r in out)
    for r_out, r_meta in zip(out, meta):
        if counts[r_out["sample"]] > 1 and r_meta["delivery_suffix"] and \
                not r_out["sample"].endswith("__" + r_meta["delivery_suffix"]):
            print(f"NOTE  {r_out['sample']} is in more than one delivery; naming this one "
                  f"{r_out['sample']}__{r_meta['delivery_suffix']}", file=sys.stderr)
            r_out["sample"] = r_meta["sample"] = f"{r_out['sample']}__{r_meta['delivery_suffix']}"
    dupe_new = [n for n, c in Counter(r["sample"] for r in out).items() if c > 1]
    for n in dupe_new:
        problems.append(f"two or more samples would be named '{n}' — check the mapping")

    print(f"{len(samples)} samples: {len(out) - (len(unmapped) if a.keep_unmapped else 0)} renamed, "
          f"{len(unmapped)} unmapped", file=sys.stderr)
    if n_dup_ok:
        print(f"NOTE  {n_dup_ok} samples have their raw_id on several agreeing rows (e.g. one row per "
              f"R1/R2 file); used one row each", file=sys.stderr)
    if n_sep:
        print(f"NOTE  {n_sep} samples matched only after treating '-' and '_' as equal "
              f"(e.g. HHGCVDRX3_1_... in the sheet vs HHGCVDRX3-1-... in the FASTQ names)", file=sys.stderr)
    for line in near + notes:
        print("NOTE  " + line, file=sys.stderr)
    if warns:
        print(f"WARN  instrument in mapping disagrees with read headers for {len(warns)} sample(s); "
              f"continuing (--instrument-mismatch warn). Header instrument is recorded in the "
              f"metadata column 'instrument_in_header'. First few:", file=sys.stderr)
        for w_ in warns[:5]:
            print("WARN    " + w_, file=sys.stderr)
    for u in unmapped:
        print(f"UNMAPPED  {u}", file=sys.stderr)
    for p in problems:
        print("PROBLEM  " + p, file=sys.stderr)
    if problems or (unmapped and not a.keep_unmapped):
        sys.exit("Stopped: fix the problems above (or pass --keep-unmapped for unmapped samples). "
                 "Nothing was written.")

    with open(a.out_samples, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "r1", "r2"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    with open(a.out_metadata, "w", newline="") as fh:
        cols = ["sample", "sample_id", "site_id", "raw_id", "sequencing_instrument_recorded",
                "instrument_in_header", "delivery_suffix", "provider_name"]
        for r in meta:
            for k in r:
                if k not in cols:
                    cols.append(k)
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(meta)
    print(f"wrote {a.out_samples} and {a.out_metadata}", file=sys.stderr)


if __name__ == "__main__":
    main()