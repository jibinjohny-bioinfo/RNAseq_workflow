# Y-chromosome ONT variant calling pipeline (Nextflow DSL2)

This pipeline filters ONT reads, extracts Y-chromosome reads, aligns them to the reference, calls structural variants (Sniffles2) and SNP/indels (Medaka), and merges SNP VCFs across samples.

## Requirements
- Nextflow 21.10+ (DSL2 enabled)
- Tools in PATH:
  - filtlong
  - minimap2
  - samtools
  - sniffles (Sniffles2)
  - medaka
  - bcftools
  - bgzip, tabix (htslib)

Optionally, configure Conda or containers via a `nextflow.config`.

## Inputs
- FASTQ files: set with `--input_fastq` (glob).
- Reference FASTA: `--reference_fasta` (must be indexed for samtools/minimap2 as needed).
- Y contig name: `--contig_y` (default `chrY`, set to `Y` if reference uses non-`chr` naming).

Directory layout example:
```
ref/
  T2T-CHM13v2.0.fa
  T2T-CHM13v2.0.fa.fai  # recommended
  T2T-CHM13v2.0.mmi     # optional prebuilt minimap2 index for speed
data/
  sample1.fastq.gz
  sample2.fastq.gz
```

## Run
```
nextflow run main.nf \
  --input_fastq "/absolute/path/data/*.{fastq,fq,fastq.gz,fq.gz}" \
  --reference_fasta "/absolute/path/ref/T2T-CHM13v2.0.fa" \
  --contig_y chrY \
  --threads 16 \
  --outdir "/absolute/path/results"
```

## Outputs
- `${outdir}/filtered`: filtered FASTQs
- `${outdir}/y_extract`: Y-only FASTQs
- `${outdir}/aligned`: BAMs aligned to reference (Y reads)
- `${outdir}/variants`: per-sample SV and SNP/indel VCFs (bgzipped + indexed) and `merged_snps.vcf.gz`

## Notes
- If your reference uses `Y` instead of `chrY`, pass `--contig_y Y`.
- Consider adding MAPQ filtering in `ExtractYReads` if needed (e.g., `samtools view -q 10`).
- For GPU Medaka, add `--gpu` to `medaka_variant` and run on a capable executor.
- You can add a `nextflow.config` with conda or container specs for reproducibility.

