# RNA-seq Genome-Guided Assembly Pipeline (Nextflow)

This Nextflow DSL2 pipeline performs RNA-seq quality control, genome alignment, genome-guided de novo transcriptome assembly with Trinity, coding sequence prediction with TransDecoder, completeness assessment with BUSCO, protein similarity search with BLASTx, sample-level quantification with Salmon, and differential expression with DESeq2 (via tximport). Optional featureCounts is provided if a gene annotation is available.

## Features
- FastQC per-sample QC and MultiQC aggregation
- HISAT2 genome index and paired-end alignment, sorted/indexed BAMs
- Trinity genome-guided assembly across all BAMs
- TransDecoder ORF prediction
- BUSCO transcriptome completeness
- BLASTx annotation against a reference protein FASTA
- Salmon quantification against assembled transcripts
- DESeq2 differential expression using tximport
- Optional featureCounts quantification (requires GTF/GFF)

## Requirements
- Nextflow 22.10+ (`nextflow -version`)
- Either Conda/Mamba (recommended) or Docker/Singularity

## Input samplesheet
CSV with header and four columns: `sample,condition,fq1,fq2`.
- **sample**: short sample identifier (no spaces)
- **condition**: group label (e.g., `M` for male, `F` for female)
- **fq1**: path to read 1 FASTQ (gzip accepted)
- **fq2**: path to read 2 FASTQ (gzip accepted)

Example:
```csv
sample,condition,fq1,fq2
F1,F,data/IT_F1_1.fq.gz,data/IT_F1_2.fq.gz
F2,F,data/IT_F2_1.fq.gz,data/IT_F2_2.fq.gz
F3,F,data/IT_F3_1.fq.gz,data/IT_F3_2.fq.gz
M1,M,data/IT_M1_1.fq.gz,data/IT_M1_2.fq.gz
M2,M,data/IT_M2_1.fq.gz,data/IT_M2_2.fq.gz
M3,M,data/IT_M3_1.fq.gz,data/IT_M3_2.fq.gz
```

## Quick start
Run with Conda:
```bash
nextflow run main.nf -profile standard,conda \
  --genome data/genome.fasta \
  --samplesheet data/samplesheet.csv \
  --ref_proteins data/ref_proteins.fasta \
  --busco_db insecta_odb10 \
  --outdir results \
  --trinity_prefix gg_trinity \
  --max_intron 10000 \
  --de_case_label M --de_control_label F
```

Optional: featureCounts if you have a gene annotation:
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
- `results/salmon_index/` — Salmon index of assembled transcripts
- `results/salmon_quant/` — Per-sample `quant.sf`
- `results/de/` — Differential expression: `deseq2_results.tsv`, `normalized_counts.tsv`, `plots/`
- `results/counts.txt` — featureCounts matrix (if `--featurecounts_gtf` is provided)

## Parameters (selected)
- `--genome` (string): Path to genome FASTA
- `--samplesheet` (string): CSV with `sample,condition,fq1,fq2`
- `--ref_proteins` (string): Path to protein FASTA for BLASTx DB
- `--busco_db` (string): BUSCO lineage name (e.g. `insecta_odb10`) or local lineage path
- `--outdir` (string): Output directory (default: `results`)
- `--trinity_prefix` (string): Prefix for assembly output directory name
- `--max_intron` (int): Trinity `--genome_guided_max_intron` (default: 10000)
- `--de_case_label` (string): Case condition label (default: `M`)
- `--de_control_label` (string): Control condition label (default: `F`)
- `--featurecounts_gtf` (string|null): GTF/GFF path for featureCounts
- `--max_cpus` (int): Max CPUs per process (default: 16)
- `--max_memory` (string): Max memory per process (default: `100.GB`)

## Notes
- DE is performed at transcript-level using Salmon + tximport + DESeq2; adjust `--de_case_label` and `--de_control_label` to match your groups.
- BUSCO may download lineage data on first run; set `--busco_db` to a local path to avoid downloads.
- featureCounts requires a valid annotation (GTF/GFF). Counting against a FASTA is not supported.

