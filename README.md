# vibrio-genome-qc

A Nextflow (DSL2) pipeline that reassembles and quality-checks a
cohort of *Vibrio* bacterial genome assemblies from paired-end
Illumina short reads:

1. resolves BioSample accessions to SRA runs,
2. downloads and merges the raw reads,
3. re-assembles from *Vibrio*-classified reads,
4. cleans the reassembly with GUNC (per-contig chimera / cross-clade
   contamination gate),
5. compares the original vs. cleaned assembly with CheckM2 / BUSCO /
   SeqKit stats,
6. emits a per-isolate summary and a cohort roll-up.

## Pipeline steps

```
metadata.csv
    │
    ├─► RESOLVE_SRA ─► DOWNLOAD_READS ─► MERGE_RUNS ─► FASTP
    │                                                    │
    │                                                    ▼
    │                                    KRAKEN2 ─► EXTRACT_VIBRIO_READS
    │                                                    │
    │                                                    ▼
    │                                                UNICYCLER
    │                                                    │
    │                          ┌────────── STATS_PRE ────┤
    │                          │           CHECKM2_PRE   │
    │                          │                         ▼
    │                          │                 COV_LENGTH_FILTER
    │                          │                         │
    │                          │                         ▼
    │                          │                       GUNC
    │                          │                         │
    │                          │                         ▼
    │                          │                 APPLY_FCS_POLICY
    │                          │                         │
    │                          │           ┌──── STATS_POST ────┤
    │                          │           │     CHECKM2_POST   │
    │                          │           │     BUSCO_POST     │
    ▼                          ▼           ▼                    ▼
 FETCH_ORIGINAL_ASSEMBLY   STATS_ORIGINAL   CHECKM2_DELTA
      │                    CHECKM2_ORIGINAL   │
      │                    BUSCO_ORIGINAL     ▼
      └──────────────────────► ► ► ► ► SAMPLE_SUMMARY ─► COHORT_SUMMARY
```

## Prerequisites

* Nextflow ≥ 23.10 (DSL2).
* Conda / Mamba. (Docker/Singularity profiles exist in `nextflow.config`
  but no process yet declares a `container` image, so conda is
  currently the only supported execution method.)
* **Hardware recommendations**:
  * **RAM**: ≥ 32 GB recommended (Kraken2 with `--memory-mapping` and Unicycler assembly peaks require ~20–28 GB).
  * **CPU**: ≥ 8 cores recommended.
  * **Disk**: ≥ 60 GB free space for reference databases (~16 GB Kraken2, ~13 GB GUNC, ~4 GB CheckM2/BUSCO) plus temporary FASTQ/assembly work dirs.
* External reference databases (configured in `conf/local.config` or passed via CLI):

### Database Setup & Downloads

#### 1. Kraken2 Database (16 GB Standard Capped)
Instead of building or downloading the ~100+ GB full standard database, use the pre-built **Standard-16** (16 GB hash table limit) provided by the Ben Langmead laboratory. It includes standard Archaea, Bacteria, Viral, plasmid, Human, and UniVec Core sequences:

```bash
mkdir -p databases/kraken2 && cd databases/kraken2
# Download pre-built standard 16 GB database (tar archive ~13-15 GB compressed)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240904.tar.gz
# or fetch the latest standard-16 from: https://benlangmead.github.io/aws-indexes/k2

tar -xzvf k2_standard_16gb_20240904.tar.gz
# Point `kraken2_db` to this directory (containing hash.k2d, opts.k2d, taxo.k2d)
```

#### 2. GUNC Reference Database (~13 GB)
GUNC uses the proGenomes 2.1 Diamond database for chimerism and cross-clade contamination screening:

```bash
mkdir -p databases/gunc && cd databases/gunc
gunc download_db .
# Point `gunc_db` to: <path>/gunc_db_progenomes2.1.dmnd
```

#### 3. BUSCO Lineage Dataset
Download the `vibrionales_odb10` lineage dataset (or allow offline usage):

```bash
mkdir -p databases/busco && cd databases/busco
busco --download vibrionales_odb10 --download_path .
# Point `busco_download_path` to this directory (parent directory containing `lineages/`)
```

#### 4. CheckM2 Database (Optional / Recommended)
CheckM2 requires the `uniref100.KO.1.dmnd` diamond model database:

```bash
mkdir -p databases/checkm2 && cd databases/checkm2
checkm2 database --download --path .
# Point `checkm2_db` to: <path>/CheckM2_database/uniref100.KO.1.dmnd
```

---

## Install

```bash
git clone https://github.com/tetyanaking/vibrio-genome-qc.git
cd vibrio-genome-qc

# Copy the example local config and edit paths for your host:
cp conf/local.config.example conf/local.config
$EDITOR conf/local.config

# Copy the example samplesheet, or point --input at your own:
cp data/metadata_example.csv data/metadata.csv
$EDITOR data/metadata.csv
```

Then run with the conda profile (creates one env per process at first run):

```bash
nextflow run . -profile conda,local
```

## Input manifest

`data/metadata.csv` (or `--input <path>`) — one row per isolate. This
file is host-specific and gitignored; copy the shipped
`data/metadata_example.csv` to `data/metadata.csv` (see Install above)
or pass `--input` to point at your own manifest:

```csv
sample_id,assembly_name,assembly_accession,biosample
Genome1,PDT000000001.1,GCA_000000001.1,SAMN00000001
Genome2,PDT000000002.1,GCA_000000002.1,SAMN00000002
```

Only `sample_id`, `assembly_accession`, and `biosample` are consumed;
`assembly_name` is an optional label kept for reference.

## Running

```bash
# Default run (uses data/metadata.csv, results/ output dir)
nextflow run main.nf -profile local,conda

# Override the input manifest or output directory
nextflow run main.nf -profile local,conda \
    --input  /path/to/other.csv \
    --outdir /path/to/results

# Resume a partially completed run
nextflow run main.nf -profile local,conda -resume
```

Key parameters (all overridable on the command line):

| param                 | default                            | meaning                                  |
| --------------------- | ---------------------------------- | ---------------------------------------- |
| `input` / `samplesheet` | `data/metadata.csv`              | isolate input manifest                   |
| `outdir`              | `results/`                         | output directory                         |
| `vibrio_taxid`        | `662`                              | genus *Vibrio* NCBI taxid                |
| `busco_lineage`       | `vibrionales_odb10`                | BUSCO lineage                            |
| `busco_download_path` | `null` (must be set)               | parent dir of `lineages/…`               |
| `kraken2_db`          | `null` (must be set)               | Kraken2 database dir                     |
| `checkm2_db`          | `null`                             | CheckM2 diamond DB (optional)            |
| `gunc_db`             | `null` (must be set)               | GUNC progenomes DB file                  |

## Outputs (`results/`)

```
results/
├── metadata/                # RESOLVE_SRA: BioSample -> SRA run table per isolate
├── raw_reads/                # DOWNLOAD_READS: downloaded FASTQs per SRA run
├── merged_reads/             # MERGE_RUNS: per-isolate FASTQs (runs concatenated)
├── original_assemblies/      # FETCH_ORIGINAL_ASSEMBLY: original GenBank assembly
├── fastp/                    # FASTP: trimmed reads + HTML/JSON QC report
├── kraken2/                  # KRAKEN2: per-read classification + report
├── vibrio_reads/             # EXTRACT_VIBRIO_READS: Vibrio-classified read pairs
├── unicycler/                # UNICYCLER: reassembled scaffolds + assembly graph/log
├── coverage_filter/          # COV_LENGTH_FILTER: coverage/length-filtered assembly + coverage table
├── gunc/                     # GUNC: per-sample maxCSS TSV + contig assignments
├── fcs_cleaned/              # APPLY_FCS_POLICY: cleaned FASTA + review + decisions per sample
├── seqkit/                   # SEQKIT_STATS: per-assembly stats (original/pre/post)
├── checkm2/                  # CHECKM2: quality reports (original/pre/post)
├── checkm2_delta/            # CHECKM2_DELTA: pre -> post completeness/contamination delta
├── busco/                    # BUSCO: short summaries (original/post)
├── multiqc/                  # MULTIQC: aggregated HTML report + data
└── summary/
    ├── <sample>.summary.tsv  # one row per isolate
    └── cohort_summary.tsv    # concatenated cohort table
```

`cohort_summary.tsv` is the paper-ready artifact — one row per isolate
with `n_contigs`, `total_len`, `N50`, completeness pre/post/delta,
contamination pre/post/delta, BUSCO %C, and the count of contigs the
policy resolved to EXCLUDE / REVIEW.

## Validation Suites

### 1. Kraken2 Classification Validation for *Vibrio*
Before deploying the pipeline on isolate reads, the Kraken2 read classification and extraction step was validated against ground-truth in silico simulated datasets using a Standard-16 GB database (referred to as `k2_standard_16_GB_20260626` in the scripts below — substitute the path to whichever Standard-16 build you downloaded, e.g. the `k2_standard_16gb_20240904` release from the Prerequisites section):

* **Positive *Vibrio* references (6 isolates)**: *V. cholerae*, *V. parahaemolyticus*, *V. vulnificus*, *V. alginolyticus*, *V. anguillarum*, *V. harveyi*.
* **Near-neighbour specificity controls (5 isolates)**: *Photobacterium profundum*, *Aliivibrio salmonicida*, *Aeromonas hydrophila*, *Grimontia hollisae*, *Shewanella oneidensis*.
* **Distant outgroups (3 isolates)**: *Escherichia coli*, *Enterococcus faecalis*, *Staphylococcus aureus*.
* **Mock community evaluation**: 60% *Vibrio* / 20% near-neighbour / 20% distant reads with confidence thresholds swept from `0.00` to `0.20` to verify sensitivity, specificity, and abundance recovery.

**Result**: `--confidence 0.00` (Kraken2's default) was the clear winner and
is what the pipeline uses. Specificity against near-neighbour/distant reads
was already ~99-100% at `0.00` (see
[`metrics/summary.tsv`](validation_kraken2/metrics/summary.tsv)), so raising
the threshold bought essentially nothing on specificity while collapsing
sensitivity for true *Vibrio* reads — from ~97-99.97% at `0.00` down to
~64-69% at `0.10` and under 1% at `0.20`. The pipeline's `KRAKEN2` module
([modules/local/kraken2.nf](modules/local/kraken2.nf)) therefore doesn't
set `--confidence` at all and relies on this default; it isn't exposed as
a `params.*` option.

### 2. Adapting Validation for Other Bacteria
The validation suite in `validation_kraken2/` is modular and can be reused to benchmark any other bacterial genus:

1. **Update Reference Accessions**:
   Edit `validation_kraken2/refs/accessions.tsv` with your target accessions (`category=positive`), closely related phylogenetic neighbors (`category=near`), and outgroup controls (`category=distant`):
   ```tsv
   accession	category	species	taxid	label
   GCF_000006945.2	positive	Salmonella_enterica	28901	S_enterica_LT2
   GCF_000005845.2	near	Escherichia_coli	562	E_coli_K12_MG1655
   GCF_000013425.1	distant	Staphylococcus_aureus	1280	S_aureus_NCTC_8325
   ```
2. **Run Validation Suite** (scripts are numbered and run in order; each is idempotent and skips work already done):
   ```bash
   export KRAKEN2_DB=/path/to/your/kraken2_db
   export TAXDUMP_DIR=/path/to/taxdump   # ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

   cd validation_kraken2/scripts
   ./00_download_refs.sh
   ./10_simulate_reads.sh
   ./20_make_mock.sh
   ./40_confidence_sweep.sh          # sweeps confidence via 30_run_kraken2.sh
   python3 50_evaluate.py
   python3 60_plot.py --confidence 0.00 0.05 0.10 0.15 0.20
   ```
3. **Inspect Metrics**:
   Check `validation_kraken2/metrics/summary.tsv` and `validation_kraken2/plots/sens_vs_spec.png` to confirm classification accuracy and select the optimal confidence threshold.
4. **Configure Pipeline**:
   Update `params.vibrio_taxid` (or pass `--vibrio_taxid <NEW_TAXID>`) and `params.busco_lineage` in `nextflow.config` to match the new organism.
   If your sweep favors a `--confidence` other than Kraken2's default
   (`0.00`, what this pipeline uses for *Vibrio* — see the Result above),
   there's no `params.*` for it: add `--confidence <value>` directly to
   the `kraken2` command in
   [modules/local/kraken2.nf](modules/local/kraken2.nf).

### 3. Pipeline End-to-End Validation
The complete end-to-end workflow was evaluated against a benchmark set of **14 reference isolates** with known ground truth ([reference_accessions.tsv](validation_pipeline/reference_accessions.tsv)):

* **6 True Positives (*Vibrio* targets)**: *V. cholerae*, *V. parahaemolyticus*, *V. vulnificus*, *V. alginolyticus*, *V. anguillarum*, *V. harveyi*.
* **5 Near-Neighbour Controls**: *Photobacterium profundum*, *Aliivibrio salmonicida*, *Aeromonas hydrophila*, *Grimontia hollisae*, *Shewanella oneidensis*.
* **3 Distant Outgroups**: *Escherichia coli*, *Enterococcus faecalis*, *Staphylococcus aureus*.

#### Validation Logic & Results
The evaluation script ([validation_pipeline/evaluate_pipeline.py](validation_pipeline/evaluate_pipeline.py)) checks whether:
1. Kraken2 identifies the expected top genus for pure isolates.
2. Only true *Vibrio* reads pass the taxonomic extraction filter.
3. CheckM2 post-cleaning contamination satisfies MIMAG high-quality caps ($<5\%$).
4. GUNC chimerism gate passes (`pass.GUNC == True`).

**Summary Confusion Matrix** ([validation_pipeline/pipeline_validation_summary.tsv](validation_pipeline/pipeline_validation_summary.tsv)):

| Sample Category | Total ($N$) | Verdict Correct | Verdict Wrong | Genus Correct | Genus Wrong |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`target` (*Vibrio*)** | 6 | **6 (100%)** | 0 | **6 (100%)** | 0 |
| **`non_target_near`** | 5 | **5 (100%)** | 0 | **5 (100%)** | 0 |
| **`non_target_distant`** | 3 | **3 (100%)** | 0 | **3 (100%)** | 0 |

#### Running Pipeline Validation
To run the automated evaluation on completed pipeline results:
```bash
python3 validation_pipeline/evaluate_pipeline.py \
    --pipeline-outdir results \
    --truth validation_pipeline/reference_accessions.tsv \
    --out-report validation_pipeline/pipeline_validation_report.tsv \
    --out-summary validation_pipeline/pipeline_validation_summary.tsv
```

`validation_pipeline/local_samplesheet.csv` (built by `build_local_samplesheet.sh`
from your local `sim/pure/` reads) is a generated, host-specific file and is
gitignored — regenerate it locally rather than expecting it in a fresh clone.

## Quality Standards & Policies

The decontamination and filtering logic in this pipeline adheres to established **NCBI Assembly Quality** and **MIMAG (Bowers et al. 2017)** standards:

1. **NCBI-compliant Contig Length Threshold**: Contigs $< 500\text{ bp}$ are removed ([modules/local/cov_length_filter.nf](modules/local/cov_length_filter.nf)), aligning with NCBI GenBank submission requirements.
2. **Coverage Outlier Filtering**: Contigs with coverage $< 0.25\times$ or $> 3.0\times$ the assembly median are flagged and removed to eliminate spurious or low-confidence fragments.
3. **NCBI FCS-GX / GUNC-style Gating Action Tiers**:
   * **`EXCLUDE`** ($\ge 30\%$ chimeric/foreign signal): Contigs are automatically purged from the final assembly.
   * **`REVIEW`** ($10\% - 30\%$ chimeric signal or genome-level maxCSS failure): Retained but flagged in `<sample>.review.tsv` for human inspection.
   * **`KEEP`** ($< 10\%$ signal): Passed through to the clean assembly.
4. **MIMAG Benchmark Compliance**: Outputs are evaluated against high-quality draft criteria ($\ge 90\%$ CheckM2 / BUSCO completeness, $< 5\%$ contamination).

## GUNC gating

GUNC (Orakov et al. 2021) runs on the coverage/length-filtered assembly
and its per-contig `contamination_portion` plus genome-level `pass.GUNC`
verdict drive the cleaning policy. Thresholds (overridable on the
`apply_fcs_policy.py` command line):

* `--gunc-exclude-threshold` `0.30` — per-contig
  `contamination_portion ≥ 0.30` → contig is dropped from the cleaned
  FASTA (auto EXCLUDE).
* `--gunc-review-threshold` `0.10` — per-contig
  `contamination_portion` in `[0.10, 0.30)` → contig kept but written
  to `<sample>.review.tsv` (manual REVIEW).
* Genome-level `pass.GUNC == False` at the maxCSS level → an extra
  `__genome__` row in `<sample>.review.tsv` so the cohort summary
  surfaces the sample even when no single contig crossed the per-contig
  threshold.

GUNC is the sole per-contig contamination decision maker in this
pipeline.

## Troubleshooting

* **`Missing param 'kraken2_db' / 'gunc_db' / 'busco_download_path'`** — fill
  these in `conf/local.config` (copied from `conf/local.config.example`).
* **BUSCO fails with "Cannot find lineage"** — `busco_download_path`
  must point at the parent that contains `lineages/vibrionales_odb10/`,
  not the lineage directory itself.
* **GUNC OOM** — the process reserves 20 GB / 8 CPU / 4h. Adjust
  `withName: GUNC { … }` in `conf/base.config` if your host has less
  (or needs more — GUNC's Diamond search can exceed 20 GB on some hosts).
* **Kraken2 report empty** — verify `--vibrio_taxid` matches the taxid
  your DB indexes at (default `662` for genus *Vibrio*).
* **Rerunning after fixing config** — use `-resume`; cached tasks are
  keyed on process signature, so config edits that don't change command
  hashes will replay from cache.

## Repo layout

```
vibrio-genome-qc/
├── main.nf                    # workflow wiring
├── nextflow.config            # params, env, profiles
├── conf/
│   ├── base.config            # per-process resources
│   └── local.config.example   # copy to conf/local.config (gitignored)
├── envs/                      # per-tool conda environment YAMLs
├── bin/
│   └── apply_fcs_policy.py    # Biopython policy applier (GUNC-driven)
├── data/
│   └── metadata_example.csv   # copy to data/metadata.csv (gitignored)
├── modules/local/             # one .nf per process
├── assets/                    # MultiQC config, placeholder files
├── validation_kraken2/        # Kraken2 classification validation suite
└── validation_pipeline/       # end-to-end pipeline validation suite
```

## Development

I designed the pipeline architecture, selected and integrated the bioinformatics tools, 
defined the GUNC-based cleaning policy and thresholds, established MIMAG-based quality targets, 
and designed the validation strategy, including the Kraken2 classification sweep 
and a 14-isolate end-to-end benchmark.

Implementation, testing, debugging, and documentation were developed using AI-assisted pair 
programming with Claude (Anthropic).

## Citations

This pipeline wraps a number of tools written by other people —
fastp, Kraken2, Unicycler, GUNC, CheckM2, BUSCO, SeqKit, MultiQC, and
more. Full citations for every tool are in [CITATIONS.md](CITATIONS.md);
please cite them alongside this pipeline if you use it.

It was also applied to clean whole-genome sequencing data for the study: 
King T, Pedrueza M, Rahman M, Oh B, LaMontagne MG. Resolution of MALDI-TOF MS compared to 
whole genome sequencing for the identification of Vibrio parahaemolyticus strains 
isolated from oysters. Arch Microbiol. 2026 Aug 20;208(11):591. 
[https://doi.org/10.1007/s00203-026-05142-8](https://doi.org/10.1007/s00203-026-05142-8)

## License

MIT — see [LICENSE](LICENSE).
