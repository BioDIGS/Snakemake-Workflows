# BioDIGS short-read QC

This is the Snakemake workflow we use to clean the Illumina short reads from the BioDIGS soil metagenomics project. It trims adapters and poly-G tails, filters out low-quality reads, removes anything that looks human, and documetns what was done.

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [1. What's in the Repository](#1-whats-in-the-repository)
- [2. Requirements](#2-requirements)
- [3. Configuring the Run](#3-configuring-the-run)
- [4. Running the Workflow](#4-running-the-workflow)
- [5. Outputs](#5-outputs)
- [Citing and License](#citing-and-license)

---

## Overview

Every sample goes through these steps:

```
raw paired-end reads
  │
  ├─► adapter check (raw)       which adapter family is actually in the reads?
  │     └─► compared against the adapter set assigned to that sample's delivery
  │
  ├─► fastp                     adapters, poly-G tails, 3' quality trimming,
  │                             quality/N filters, both mates ≥ 50 bp
  │
  ├─► Bowtie2 vs T2T-CHM13      --very-sensitive --reorder
  │     └─► a pair is kept only if BOTH mates fail to align
  │
  └─► clean reads
        ├─► adapter check (clean), seqkit stats, FastQC
        ├─► md5 checksums
        └─► per-sample QC table, sample sheet, MultiQC report
```

A few things are worth knowing before you use it:

**Human removal is strict on purpose.** If either mate of a pair aligns to the human reference, the whole pair is removed. These reads end up public, so we would rather lose a few soil reads than leak human ones. If you want to see what was taken out, set `keep_human_reads: true` and the removed pairs are written to their own folder.

**Adapters are assigned per delivery, then checked.** Our samples came back from the sequencing provider in batches prepared with different library kits (Nextera in 2023, TruSeq in 2025). Rather than relying only on fastp's auto-detection, each delivery gets an explicit adapter set in the config. Before fastp runs, the workflow looks at the raw reads and stops any sample whose adapters don't match the set it was given.

**Everything that affects the output is recorded.** Tool versions, the exact fastp arguments used for each sample, the Bowtie2 index files and their sizes, copies of the config and sample list, and md5 checksums for both raw and clean reads all end up in `results/provenance/`.

---

## Quick Start

1. Clone this repository.
2. Get a Bowtie2 index for T2T-CHM13v2.0 (UCSC `hs1`). See [Requirements](#2-requirements).
3. Copy the example config and edit the lines marked `EDIT`:
   ```bash
   cp config/config.example.yaml config/config.yaml
   ```
4. Make your sample list:
   ```bash
   python workflow/scripts/make_samples_tsv.py /path/to/raw > config/samples.tsv
   # or
   edit the example sample tsv
   ```
5. Build the conda environments (one time only):
   ```bash
   bash run_devnode.sh --envs-only
   ```
6. Check every sample's adapters before committing to the full run.:
   ```bash
   snakemake adapter_audit
   ```
7. Run it:
   ```bash
   snakemake
   # cores and resources can be edited in /profiles/default/config.yaml     
   ```

---

## 1. What's in the Repository

```
.
├── README.md
├── config/
│   ├── config.example.yaml        copy to config.yaml and edit
│   └── samples.tsv.example        format for the sample list
├── profiles/
│   └── default/
│       └── config.yaml            picked up automatically by `snakemake`
└── workflow/
    ├── Snakefile
    ├── multiqc_config.yaml
    ├── envs/                      one conda env per tool
    └── scripts/                   helpers called by the rules
```

- **`workflow/Snakefile`**: all the rules. Snakemake finds it automatically when you run from the repository root.
- **`config/config.example.yaml`**: every setting, commented, with the values we used for BioDIGS. The workflow reads `config/config.yaml`, so copy this first.
- **`profiles/default/config.yaml`**: turns on conda, keeps going past failed samples, and sets conservative default cores and memory. Raise those for a bigger machine.
- **`workflow/envs/`**: conda environments for fastp, Bowtie2/samtools, seqkit, FastQC, MultiQC and Python. Snakemake builds and activates them for you.
- **`workflow/scripts/`**: adapter identification, the adapter-assignment check, per-sample summaries, the sample sheet builder, and the spike-in simulation and validation.

---

## 2. Requirements

**Software**

- Linux
- conda or mamba
- Snakemake v8 or newer. The profile uses `software-deployment-method`, which is v8 syntax. The version we ran with is recorded in `results/provenance/versions.tsv`.

**Input data**

- Paired-end Illumina reads, gzipped FASTQ.
- Samples split across several lanes are fine. List all the lane files and they are concatenated before processing.

**Human reference (Bowtie2 index)**

We filter against T2T-CHM13v2.0, distributed by UCSC as `hs1`. To build the index yourself:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
gunzip hs1.fa.gz
bowtie2-build --threads 16 hs1.fa hs1
```

This takes a few hours. In the config, `bowtie2_index` is the index *prefix* (e.g. `/refs/hs1/hs1`), not the folder. Keep `hs1.fa` around if you want to run the optional spike-in validation.

**Resources**

- Bowtie2 needs about 12 GB of memory per job (`resources.bowtie2_mem_mb` in the config).
- Disk: clean reads come out a little smaller than the raw reads. While a sample is running, there is also a temporary BAM file for it under `results/host/`.

---

## 3. Configuring the Run

You edit two files: `config/config.yaml` (settings) and `config/samples.tsv` (which reads to process).

### config/config.yaml

- **`samples`, `outdir`**: where your sample list is and where results go.
- **`host`**: the Bowtie2 index and the arguments used to align against it. `keep_human_reads` writes the removed pairs to `results/host/human/`. It is off by default. If you turn it on, treat that output as human data.
- **`fastp.adapter_sets`**: the adapters for each delivery (see below).
- **The rest of `fastp`**: trimming and filtering thresholds. One that deserves a closer look is `cut_right`. It slides a 4 bp window from the 5' end and cuts the read at the *first* window with mean quality below 15, dropping everything after it. That is stricter than trimming only the 3' tail, so check the read-length plots in MultiQC if you change it.
- **`threads`, `resources`**: per-tool threads and memory.
- **`sample_sheet.metadata`**: optional; see below.
- **`validation`**: optional human spike-in test; see [Running](#4-running-the-workflow).

### config/samples.tsv

Tab-separated, one row per sample, with columns `sample`, `r1` and `r2`:

```
sample	r1	r2
A01_10	/data/2023.09.25_delivery/A01_10_R1.fastq.gz	/data/2023.09.25_delivery/A01_10_R2.fastq.gz
B12_03	/data/2025.02.19_delivery/B12_03_L001_R1.fastq.gz,/data/2025.02.19_delivery/B12_03_L002_R1.fastq.gz	/data/2025.02.19_delivery/B12_03_L001_R2.fastq.gz,/data/2025.02.19_delivery/B12_03_L002_R2.fastq.gz
```

- Sample names can only use letters, numbers, `.`, `_` and `-`.
- For multi-lane samples, separate the files with commas, keeping R1 and R2 in the same order.
- `workflow/scripts/make_samples_tsv.py` builds this file from a folder of raw reads, which is a lot less error-prone than writing it by hand.

### Adapter sets

Each entry in `fastp.adapter_sets` says which adapters a group of samples should have:

```yaml
adapter_sets:
  - label: nextera_2023
    match_path: "/2023\\.09\\.25[^/]*/"
    r1: CTGTCTCTTATACACATCT
    r2: CTGTCTCTTATACACATCT
```

- **`match_path`**: a regular expression checked against each sample's first R1 path. We match on the delivery folder, since every delivery used one kit.
- **`match`** (optional): a regular expression matched against the sample *name* instead.
- **The first set that matches a sample wins.**

If any sample matches no set, the run refuses to start and lists those samples. That's intentional: a sample with no adapter set is almost always a missing config line, not something you want trimmed on auto-detect alone. If you really do want auto-detection for unmatched samples, set `allow_unmatched_adapter_sets: true`.

Replace our two sets with your own. Running `snakemake adapter_audit` will tell you quickly if you got one wrong.

### Sample metadata (optional)

`results/sample_sheet.tsv` is built mostly from the reads themselves: instrument, run, flowcell, lane and index from the read headers, plus quality, adapter, insert-size and duplication metrics, and what cleaning removed. If you also want field or lab information in that table (site, collection date, extraction batch), point `sample_sheet.metadata` to a TSV with a `sample` column, and the extra columns are added by sample name. Leave it as `""` if you don't have one.

---

## 4. Running the Workflow

Run everything from the repository root.

### Check the adapters first

```bash
snakemake adapter_audit
```

This only reads a couple hundred thousand pairs per sample, so it finishes in minutes. It writes `results/adapters/adapter_audit.tsv` and prints how many samples fall into each status. A sample flagged `LOW_EVIDENCE` had too little adapter read-through to tell either way, which is common for long-insert libraries. Those samples are allowed through with a flag. A real mismatch stops that sample before fastp.

### Do a dry run

```bash
snakemake -n
```

This lists every job Snakemake plans to run without running anything. It's the quickest way to catch a bad path in the config.

### Run on a normal machine

```bash
snakemake --cores 16 --resources mem_mb=32000
```

- **`--cores`**: the total number of threads Snakemake may use at once.
- **`--resources mem_mb`**: the memory budget Snakemake uses to decide how many jobs run in parallel. It doesn't enforce a limit, so set it below the RAM you actually have free.

The default profile already sets high values for both, so plain `snakemake` may not work unless you adjust them.

### Stopping and resuming

If a run stops partway, just run the same command again. Snakemake picks up where it left off and redoes any job that didn't finish. Samples that fail don't stop the others (`keep-going` is on in the profile). Their errors show up in `logs/` and at the end of the run log.

### Optional: human spike-in validation

To check how well the human filter works on your data, set `validation.enabled: true` and point `validation.human_fasta` at the reference FASTA. The workflow then:

1. simulates human read pairs, with adapter read-through and poly-G tails,
2. mixes them into the first `base_pairs` pairs of one real soil sample,
3. runs that mix through the exact same pipeline as an extra sample called `SPIKECHECK`,
4. reports how many human reads leaked through and how many soil reads were lost, in `results/validation/spike_validation.txt`.

`SPIKECHECK` is kept out of the QC table, the sample sheet and the MultiQC report.

One caveat: the simulated reads come from the same reference used for filtering, so this is a best case. Real human contamination from other people carries variation the reference doesn't, so expect slightly more leakage in practice.

---

## 5. Outputs

Everything goes under `results/` (or whatever `outdir` is set to):

```
results/
├── clean/                  clean reads + md5 checksums      ← the main output
├── qc_summary.tsv          one row per sample: read accounting at each step
├── sample_sheet.tsv/.csv   one row per sample: how it was sequenced and what cleaning did
├── multiqc/                MultiQC report
├── adapters/               adapter evidence (raw, clean) and the pre-fastp check
├── fastp/                  fastp JSON + HTML per sample
├── host/                   Bowtie2 logs and human-removal counts per sample
│   └── human/              removed human pairs (only if keep_human_reads: true)
├── fastqc/                 FastQC on raw and clean reads
├── seqkit/                 read statistics before and after
├── runinfo/                instrument/run/flowcell/lane/index counts from read headers
├── summary/                per-sample pieces of qc_summary.tsv
├── provenance/             versions, parameters, checksums
└── validation/             spike-in results (only if enabled)
```

### The files you'll actually look at

- **`clean/{sample}_R1.fastq.gz`, `clean/{sample}_R2.fastq.gz`**: the cleaned read pairs.
- **`clean/{sample}.md5`**: checksums with bare file names, so the file can go straight into an SRA submission.
- **`qc_summary.tsv`**: the read accounting for each sample:
  - `raw_pairs`
  - `pct_pairs_removed_fastp`
  - `pct_pairs_human`
  - `final_pairs`
  - `pct_pairs_retained_overall`
  - adapter signal before and after (`pct_reads_adapter_signature_raw`, `pct_reads_adapter_signature_clean`)
  - `pct_pairs_readthrough_clean`
  - `counts_consistent`, a sanity check that the numbers add up between steps. If a sample says `False` here, look at it before using it.
- **`multiqc/multiqc_report.html`**: everything in one place. In the General Statistics table, `{sample}_raw` rows are FastQC on the raw reads, while `{sample}` rows show FastQC on the clean reads next to the fastp and Bowtie2 numbers.
- **`sample_sheet.tsv`**: useful for spotting batch effects, for example whether a QC metric tracks a flowcell or delivery rather than anything about the soil.

### Provenance

`results/provenance/` is there so a result can be traced back to exactly how it was made:

- **`versions.tsv`**: Snakemake and every tool version.
- **`run_info.tsv`**: run date, host reference, the Bowtie2 index files and their sizes, the Bowtie2 arguments, the shared fastp arguments, and each adapter set.
- **`fastp_adapters_per_sample.tsv`**: the adapter set and full fastp command line used for every sample.
- **`config.yaml`, `samples.tsv`**: copies of what the run used.
- **`raw_md5.tsv`**: checksums of the raw input files.

---

## Citing and License

If you use this workflow, please cite this repository. 
