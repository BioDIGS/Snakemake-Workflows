#!/usr/bin/env python3
"""
build_sample_sheet.py — one row per short-read sample with everything you'd want
as a record of how it was sequenced and cleaned (and the input for a later
drift analysis):

  * where it was sequenced — instrument, model guess, run number, flowcell(s),
    lane(s) with read shares, index (i7/i5), share of reads on the main
    flowcell/lane, dates parsed from the provider's run-folder name or path;
  * how it was sequenced — read lengths, Q30, per-cycle quality drop at the 3'
    end (R1 and R2), N rate, poly-G, adapter read-through, insert size,
    duplication, GC;
  * what cleaning did to it — fastp losses, % human, final pairs/bases;
  * your own metadata (site, contributor, extraction / library-prep batch,
    dates ...) merged in from an optional TSV/CSV keyed by `sample`.

Writes TSV and CSV with the same content.

Delivery: the provider's dated delivery folder (e.g. 2025.02.19.GRCF.NovoSeq) and
the FASTQ folder inside it are recorded as delivery_folder / delivery_date /
fastq_folder. A delivery date is when data were released, not the run date, but
it is a natural batch label.

Dates: Illumina FASTQ headers carry no date. The run date comes from the
standard run-folder name YYMMDD_<instrument>_<run>_<A|B><flowcell> if it appears
in a file path, else from any YYYYMMDD / YYYY-MM-DD in the path (marked as a hint),
else is left blank; a `seq_date` column in your metadata always wins. File
modification times are recorded but are NOT sequencing dates (copying resets them).
Within one instrument, the run number increases with time, so `run_number`
orders runs even without dates.
"""
import argparse
import csv
import datetime as dt
import json
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from identify_adapters import guess_instrument  # noqa: E402

RUN_FOLDER = re.compile(r"(?<!\d)(\d{6})_([A-Za-z0-9-]+)_(\d{1,5})_([AB]?)([A-Z0-9]{5,12})(?![A-Za-z0-9])")
PATH_DATE = re.compile(r"(?<!\d)(20\d{2})[-_.]?(0[1-9]|1[0-2])[-_.]?(0[1-9]|[12]\d|3[01])(?!\d)")
# Provider delivery folders such as "2025.02.19.GRCF.NovoSeq" or "2023-09-25_delivery"
DELIVERY = re.compile(r"^(20\d{2})[._-](0[1-9]|1[0-2])[._-](0[1-9]|[12]\d|3[01])(?:[._-].*)?$")


def delivery_info(paths):
    """(delivery_folder, delivery_date, fastq_folder) for the first R1 file.
    delivery_date is when the provider released the data, NOT the run date."""
    if not paths:
        return "", "", ""
    parts = os.path.abspath(paths[0]).split(os.sep)
    fastq_folder = parts[-2] if len(parts) > 1 else ""
    for part in reversed(parts[:-1]):
        m = DELIVERY.match(part)
        if m:
            return part, f"{m.group(1)}-{m.group(2)}-{m.group(3)}", fastq_folder
    return "", "", fastq_folder


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_table_any(path):
    """TSV or CSV, decided by extension then by sniffing."""
    with open(path, newline="") as fh:
        head = fh.read(4096)
        fh.seek(0)
        delim = "," if path.endswith(".csv") or (head.count(",") > head.count("\t")) else "\t"
        rows = list(csv.DictReader(fh, delimiter=delim))
    if rows and "sample" not in rows[0]:
        sys.exit(f"ERROR: metadata file {path} needs a 'sample' column")
    return rows


def parse_dates(paths):
    """(run_folder_date, run_folder_name, path_date_hint) from file paths."""
    for p in paths:
        for part in p.split(os.sep):
            m = RUN_FOLDER.search(part)
            if m:
                try:
                    d = dt.datetime.strptime(m.group(1), "%y%m%d").date()
                    return d.isoformat(), m.group(0), ""
                except ValueError:
                    pass
    for p in paths:
        m = PATH_DATE.search(p)
        if m:
            try:
                d = dt.date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
                return "", "", d.isoformat()
            except ValueError:
                pass
    return "", "", ""


def summarize_headers(rows):
    total = sum(int(r["reads"]) for r in rows)
    if total == 0:
        return {"header_reads_scanned": 0}

    def share(key_fn):
        c = {}
        for r in rows:
            k = key_fn(r)
            c[k] = c.get(k, 0) + int(r["reads"])
        return sorted(c.items(), key=lambda kv: -kv[1])

    def fmt(items, min_frac=0.001):
        return ";".join(f"{k}({100*v/total:.1f}%)" for k, v in items if v / total >= min_frac)

    inst = share(lambda r: r["instrument"])
    runs = share(lambda r: f"{r['instrument']}_{r['run']}")
    fcs = share(lambda r: r["flowcell"])
    lanes = share(lambda r: f"{r['flowcell']}:{r['lane']}")
    idx = share(lambda r: r["index"])
    flt = share(lambda r: r["filter"])

    primary_inst = inst[0][0]
    model, two_colour = "NA", "NA"
    if primary_inst != "NON_ILLUMINA_HEADER":
        _, model, tc = guess_instrument(f"@{primary_inst}:1:FC:1:1:1:1")
        two_colour = tc
    i7 = i5 = "NA"
    if idx[0][0] != "NA":
        parts = idx[0][0].split("+")
        i7 = parts[0]
        i5 = parts[1] if len(parts) > 1 else ""
    n_idx_major = sum(1 for _, v in idx if v / total >= 0.01)
    y = dict(flt).get("Y", 0)
    prim_run = runs[0][0].rsplit("_", 1)
    return {
        "header_format": "non-Illumina/SRA" if primary_inst == "NON_ILLUMINA_HEADER" else "Illumina",
        "header_reads_scanned": total,
        "instrument": primary_inst,
        "instrument_model_guess": model,
        "two_colour_chemistry": two_colour,
        "run_number": prim_run[1] if len(prim_run) > 1 else "NA",
        "primary_flowcell": fcs[0][0],
        "primary_flowcell_lane": lanes[0][0],
        "frac_reads_primary_flowcell": round(fcs[0][1] / total, 4),
        "frac_reads_primary_lane": round(lanes[0][1] / total, 4),
        "n_instruments": len(inst),
        "n_runs": len(runs),
        "n_flowcells": len(fcs),
        "n_lanes": len(lanes),
        "runs": fmt(runs),
        "flowcells": fmt(fcs),
        "lanes": fmt(lanes),
        "index_i7": i7,
        "index_i5": i5,
        "frac_reads_main_index": round(idx[0][1] / total, 4),
        "n_indexes_ge1pct": n_idx_major,
        "pct_reads_filter_flag_Y": round(100 * y / total, 4),
    }


def curve_stats(fp, key):
    """Per-cycle mean quality and N% from fastp's before-filtering curves."""
    sect = fp.get(key, {}) or {}
    q = (sect.get("quality_curves") or {}).get("mean") or []
    n = (sect.get("content_curves") or {}).get("N") or []
    out = {}
    if len(q) >= 40:
        head = statistics.mean(q[5:35])           # skip the first few cycles
        tail = statistics.mean(q[-20:])
        out["mean_q_cycles_6_35"] = round(head, 2)
        out["mean_q_last20"] = round(tail, 2)
        out["q_drop_3prime"] = round(head - tail, 2)
    if n:
        out["pct_N"] = round(100 * statistics.mean(n), 4)
    return out


def fastp_extras(path):
    fp = json.load(open(path))
    bf = fp.get("summary", {}).get("before_filtering", {})
    af = fp.get("summary", {}).get("after_filtering", {})
    r1 = curve_stats(fp, "read1_before_filtering")
    r2 = curve_stats(fp, "read2_before_filtering")
    return {
        # GC after trimming: raw GC is inflated by adapter and poly-G tails
        "gc_content_trimmed": af.get("gc_content", "NA"),
        "raw_gc_content": bf.get("gc_content", "NA"),
        "raw_q20_rate": bf.get("q20_rate", "NA"),
        "r1_q_drop_3prime": r1.get("q_drop_3prime", "NA"),
        "r2_q_drop_3prime": r2.get("q_drop_3prime", "NA"),
        "r1_mean_q_last20": r1.get("mean_q_last20", "NA"),
        "r2_mean_q_last20": r2.get("mean_q_last20", "NA"),
        "r1_pct_N_raw": r1.get("pct_N", "NA"),
        "r2_pct_N_raw": r2.get("pct_N", "NA"),
    }


def seqkit_lengths(path):
    rows = read_tsv(path)
    out = {}
    for i, r in enumerate(rows[:2], start=1):
        out[f"r{i}_max_len_raw"] = r.get("max_len", "NA")
        out[f"r{i}_avg_len_raw"] = r.get("avg_len", "NA")
    return out


QC_KEEP = [
    "raw_pairs", "raw_bases", "raw_q30_rate", "raw_looks_already_trimmed",
    "adapter_family_raw", "fastp_r1_adapter", "fastp_r2_adapter",
    "pct_pairs_readthrough_raw", "pct_reads_adapter_signature_raw",
    "pct_polyG_r1_raw", "fastp_duplication_rate", "fastp_insert_peak",
    "pct_pairs_removed_fastp", "pct_pairs_human", "final_pairs", "final_bases",
    "final_q30_pct", "pct_pairs_retained_overall", "counts_consistent",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, help="workflow samples.tsv (sample, r1, r2)")
    ap.add_argument("--qc-summary", required=True)
    ap.add_argument("--headers-dir", required=True)
    ap.add_argument("--fastp-dir", required=True)
    ap.add_argument("--seqkit-dir", required=True)
    ap.add_argument("--metadata", default="", help="optional TSV/CSV keyed by 'sample'")
    ap.add_argument("--exclude", default="", help="comma-separated samples to leave out")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-csv", required=True)
    a = ap.parse_args()

    exclude = {s for s in a.exclude.split(",") if s}
    samples = [r for r in read_tsv(a.samples) if r["sample"].strip() and r["sample"] not in exclude
               and not r["sample"].startswith("#")]
    qc = {r["sample"]: r for r in read_tsv(a.qc_summary)}
    meta = {}
    meta_cols = []
    if a.metadata:
        mrows = read_table_any(a.metadata)
        meta_cols = [c for c in mrows[0].keys() if c != "sample"] if mrows else []
        meta = {r["sample"].strip(): r for r in mrows}
        unknown = sorted(set(meta) - {s["sample"] for s in samples})
        missing = sorted({s["sample"] for s in samples} - set(meta))
        if unknown:
            print(f"WARNING: metadata rows with no matching sample: {', '.join(unknown)}", file=sys.stderr)
        if missing:
            print(f"WARNING: samples with no metadata row: {', '.join(missing)}", file=sys.stderr)

    out_rows = []
    for s in samples:
        name = s["sample"].strip()
        r1_files = [p.strip() for p in s["r1"].split(",") if p.strip()]
        r2_files = [p.strip() for p in s["r2"].split(",") if p.strip()]
        row = {"sample": name}
        for c in meta_cols:
            row[c] = (meta.get(name, {}).get(c) or "").strip()
        row.update(summarize_headers(read_tsv(os.path.join(a.headers_dir, f"{name}.header_counts.tsv"))))

        rf_date, rf_name, hint = parse_dates(r1_files)
        dfolders = sorted({delivery_info([f])[0] for f in r1_files})
        dfolder, ddate, fqfolder = delivery_info(r1_files)
        row["delivery_folder"] = ";".join(d for d in dfolders if d) if len(dfolders) > 1 else dfolder
        row["delivery_date"] = ddate if len(dfolders) == 1 else "MIXED"
        row["fastq_folder"] = fqfolder
        row["run_folder"] = rf_name
        row["run_folder_date"] = rf_date
        row["path_date_hint"] = hint
        try:
            row["raw_file_mtime"] = dt.date.fromtimestamp(os.path.getmtime(r1_files[0])).isoformat()
        except OSError:
            row["raw_file_mtime"] = ""
        if "seq_date" not in meta_cols:
            row["seq_date"] = rf_date
            row["seq_date_source"] = "run folder" if rf_date else ""
        else:
            row["seq_date_source"] = "metadata" if row.get("seq_date") else ""
        if row.get("n_runs", 1) > 1 and row["seq_date_source"] == "run folder":
            row["seq_date_source"] = (f"run folder — but reads come from {row['n_runs']} runs "
                                      f"(see 'runs'); this date covers only some of them")

        row.update(seqkit_lengths(os.path.join(a.seqkit_dir, f"{name}.raw.tsv")))
        row.update(fastp_extras(os.path.join(a.fastp_dir, f"{name}.fastp.json")))
        q = qc.get(name, {})
        for k in QC_KEEP:
            row[k] = q.get(k, "NA")
        row["n_raw_files"] = len(r1_files) + len(r2_files)
        row["raw_r1_files"] = ",".join(r1_files)
        row["raw_r2_files"] = ",".join(r2_files)
        out_rows.append(row)

    cols = []
    for r in out_rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    for path, delim in ((a.out_tsv, "\t"), (a.out_csv, ",")):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter=delim, extrasaction="ignore", restval="")
            w.writeheader()
            w.writerows(out_rows)
    print(f"wrote {len(out_rows)} samples x {len(cols)} columns -> {a.out_tsv}")


if __name__ == "__main__":
    main()
