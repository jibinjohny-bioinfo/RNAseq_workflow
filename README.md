# RNA-seq Genome-Guided Assembly Pipeline (Nextflow)

This Nextflow DSL2 pipeline performs RNA-seq quality control, genome alignment, genome-guided de novo transcriptome assembly with Trinity, coding sequence prediction with TransDecoder, completeness assessment with BUSCO, protein similarity search with BLASTx, and optional read counting with featureCounts.

## Features
- FastQC per-sample QC and MultiQC aggregation
- HISAT2 genome index and paired-end alignment, sorted/indexed BAMs
- Trinity genome-guided assembly across all BAMs
- TransDecoder ORF prediction
- BUSCO transcriptome completeness
- BLASTx annotation against a reference protein FASTA
- Optional featureCounts quantification (requires GTF/GFF)

## Requirements
- Nextflow 22.10+ (`nextflow -version`)
- Either Conda/Mamba (recommended) or Docker/Singularity for software provisioning

## Input samplesheet
A CSV file with a header and three columns: `sample,fq1,fq2`.
- **sample**: short sample identifier (no spaces)
- **fq1**: path to read 1 FASTQ (gzip accepted)
- **fq2**: path to read 2 FASTQ (gzip accepted)

Example:
```csv
sample,fq1,fq2
F1,data/IT_F1_1.fq.gz,data/IT_F1_2.fq.gz
F2,data/IT_F2_1.fq.gz,data/IT_F2_2.fq.gz
F3,data/IT_F3_1.fq.gz,data/IT_F3_2.fq.gz
```

## Quick start
1. Prepare input files under `data/`:
   - Genome FASTA: `data/genome.fasta`
   - Reference proteins FASTA (for BLASTx DB): `data/ref_proteins.fasta`
   - Samplesheet CSV: `data/samplesheet.csv`

2. Run with Conda (creates per-process environments):

   ```bash
   nextflow run main.nf -profile standard,conda \
     --genome data/genome.fasta \
     --samplesheet data/samplesheet.csv \
     --ref_proteins data/ref_proteins.fasta \
     --busco_db insecta_odb10 \
     --outdir results \
     --trinity_prefix gg_trinity \
     --max_intron 10000
   ```

3. Optional: provide a gene annotation for featureCounts

   ```bash
   nextflow run main.nf -profile standard,conda \
     --genome data/genome.fasta \
     --samplesheet data/samplesheet.csv \
     --ref_proteins data/ref_proteins.fasta \
     --busco_db insecta_odb10 \
     --featurecounts_gtf data/annotation.gtf
   ```

## Outputs
- `results/qc/fastqc/` — FastQC reports per sample
- `results/qc/` — MultiQC report
- `results/genome_index/` — HISAT2 genome index
- `results/*.bam` and `*.bam.bai` — Sorted and indexed BAMs per sample
- `results/gg_trinity_assembly/Trinity-GG.fasta` — Assembled transcripts
- `results/gg_trinity_assembly/Trinity-GG.fasta.transdecoder.*` — TransDecoder outputs
- `results/busco/` — BUSCO results (`busco_out/`)
- `results/blastdb/` — BLAST protein database
- `results/blastx.outfmt6` — BLASTx tabular hits (outfmt 6)
- `results/counts.txt` — featureCounts matrix (if `--featurecounts_gtf` is provided)

## Parameters (selected)
- `--genome` (string): Path to genome FASTA
- `--samplesheet` (string): CSV with `sample,fq1,fq2`
- `--ref_proteins` (string): Path to protein FASTA for BLASTx DB
- `--busco_db` (string): BUSCO lineage name (e.g. `insecta_odb10`) or a local lineage path
- `--outdir` (string): Output directory (default: `results`)
- `--trinity_prefix` (string): Prefix for assembly output directory name
- `--max_intron` (int): Trinity `--genome_guided_max_intron` (default: 10000)
- `--featurecounts_gtf` (string|null): GTF/GFF path for featureCounts; if not set, counting is skipped
- `--max_cpus` (int): Max CPUs per process (default: 16)
- `--max_memory` (string): Max memory per process (default: `100.GB`)

## Notes
- BUSCO may download lineage data on first run; set `--busco_db` to a local path to avoid downloads.
- featureCounts requires a valid annotation (GTF/GFF). Counting against a FASTA is not supported.
- For large datasets, use a batch system or cloud executor via Nextflow profiles.

## Software provisioning tips
- Prefer Mamba for faster solves: run with `-with-mpi -with-conda` and have `mamba` on PATH. With Nextflow 24.04+, set `NXF_CONDA_CHANNELS` to `conda-forge,bioconda,defaults` for robust environment resolution.
- To use containers instead of Conda, append `-profile docker` and set images per process with `process.container`. You can also define a global image if suitable.

## Reproducibility
All software is provisioned per process using Conda when run with `-profile conda`. To use containers instead, supply `-profile docker` and set appropriate `process.container` images per process.