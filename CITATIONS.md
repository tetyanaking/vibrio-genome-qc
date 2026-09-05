# Citations

This pipeline stands on the work of many tool authors. If you use
`vibrio-genome-qc`, please also cite the underlying tools it wraps.

## Workflow manager

* **Nextflow**
  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C.
  > Nextflow enables reproducible computational workflows.
  > *Nat Biotechnol.* 2017;35:316-319. doi: 10.1038/nbt.3820

## Pipeline tools

* **fastp** (read QC/trimming)
  > Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ
  > preprocessor. *Bioinformatics.* 2018;34(17):i884-i890.
  > doi: 10.1093/bioinformatics/bty560

* **Kraken2** (read classification)
  > Wood DE, Lu J, Langmead B. Improved metagenomic analysis with
  > Kraken 2. *Genome Biol.* 2019;20:257.
  > doi: 10.1186/s13059-019-1891-0

* **KrakenTools** (Kraken2 output extraction)
  > Lu J, Rincon N, Wood DE, Breitwieser FP, Pockrandt C, Langmead B,
  > Salzberg SL, Steinegger M. Metagenome analysis using the Kraken
  > software suite. *Nat Protoc.* 2022;17:2815-2839.
  > doi: 10.1038/s41596-022-00738-y
  > (Author Correction: doi: 10.1038/s41596-024-01064-1)

* **Unicycler** (assembly)
  > Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: Resolving
  > bacterial genome assemblies from short and long sequencing reads.
  > *PLoS Comput Biol.* 2017;13(6):e1005595.
  > doi: 10.1371/journal.pcbi.1005595

* **SAMtools** / **minimap2** (read mapping, used for the
  coverage/length filter)
  > Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools
  > and BCFtools. *GigaScience.* 2021;10(2):giab008.
  > doi: 10.1093/gigascience/giab008
  >
  > Li H. Minimap2: pairwise alignment for nucleotide sequences.
  > *Bioinformatics.* 2018;34(18):3094-3100.
  > doi: 10.1093/bioinformatics/bty191

* **GUNC** (chimerism / cross-clade contamination gating)
  > Orakov A, Fullam A, Coelho LP, Khedkar S, Szklarczyk D, Mende DR,
  > Schmidt TSB, Bork P. GUNC: detection of chimerism and contamination
  > in prokaryotic genomes. *Genome Biol.* 2021;22:178.
  > doi: 10.1186/s13059-021-02393-0

* **CheckM2** (genome completeness/contamination)
  > Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid,
  > scalable and accurate tool for assessing microbial genome quality
  > using machine learning. *Nat Methods.* 2023;20:1203-1212.
  > doi: 10.1038/s41592-023-01940-w
  > (Author Correction: doi: 10.1038/s41592-024-02248-z)

* **BUSCO** (gene-completeness benchmarking)
  > Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. BUSCO
  > Update: Novel and Streamlined Workflows along with Broader and
  > Deeper Phylogenetic Coverage for Scoring of Eukaryotic,
  > Prokaryotic, and Viral Genomes. *Mol Biol Evol.* 2021;38(10):4647-4654.
  > doi: 10.1093/molbev/msab199

* **SeqKit** (assembly statistics)
  > Shen W, Le S, Li Y, Hu F. SeqKit: A Cross-Platform and Ultrafast
  > Toolkit for FASTA/Q File Manipulation. *PLoS ONE.* 2016;11(10):e0163962.
  > doi: 10.1371/journal.pone.0163962

* **MultiQC** (aggregated report)
  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize
  > analysis results for multiple tools and samples in a single
  > report. *Bioinformatics.* 2016;32(19):3047-3048.
  > doi: 10.1093/bioinformatics/btw354

* **Biopython** (used by `apply_fcs_policy.py`)
  > Cock PJA, Antao T, Chang JT, et al. Biopython: freely available
  > Python tools for computational molecular biology and
  > bioinformatics. *Bioinformatics.* 2009;25(11):1422-1423.
  > doi: 10.1093/bioinformatics/btp163

* **NCBI Datasets** and **SRA Toolkit** (accession resolution and
  read/assembly download) — National Center for Biotechnology
  Information (NCBI), National Library of Medicine (US).
  https://www.ncbi.nlm.nih.gov/datasets/ ·
  https://github.com/ncbi/sra-tools

## Validation-suite-only tools

* **InSilicoSeq** (simulated reads for the Kraken2 validation suite)
  > Gourlé H, Karlsson-Lindsjö O, Hayer J, Bongcam-Rudloff E.
  > Simulating Illumina metagenomic data with InSilicoSeq.
  > *Bioinformatics.* 2019;35(3):521-522.
  > doi: 10.1093/bioinformatics/bty630

* **seqtk** (read sub-sampling for the mock community) — Heng Li.
  https://github.com/lh3/seqtk

## Standards referenced by the cleaning policy

* **MIMAG** (minimum information standard the quality gates target)
  > Bowers RM, Kyrpides NC, Stepanauskas R, et al. Minimum information
  > about a single amplified genome (MISAG) and a metagenome-assembled
  > genome (MIMAG) of bacteria and archaea. *Nat Biotechnol.*
  > 2017;35:725-731. doi: 10.1038/nbt.3893
  > (Corrigendum: doi: 10.1038/nbt0218-196a)
