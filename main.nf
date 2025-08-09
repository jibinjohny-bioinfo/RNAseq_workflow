nextflow.enable.dsl=2

params {
    input_fastq      = "./data/*.{fastq,fq,fastq.gz,fq.gz}"
    reference_fasta  = "./ref/T2T-CHM13v2.0.fa"
    annotation_gtf   = "./ref/T2T-CHM13v2.0_chrY.gtf"   // Reserved for future use
    contig_y         = "chrY"                             // Set to "Y" if your reference uses non-chr naming
    threads          = 8
    outdir           = "./results"
}

// Utility: derive a clean sample id from filename, handling .fastq/.fq and optional .gz
def inferSampleId(File f) {
    def n = f.getName()
    n = n.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
    return n
}

process FilterReads {
    tag { sample_id }
    publishDir "${params.outdir}/filtered", mode: 'copy'

    cpus params.threads

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.fastq"), emit: filtered

    script:
    """
    filtlong --min_length 5000 --keep_percent 90 ${fastq_file} > ${sample_id}.filtered.fastq
    """
}

process ExtractYReads {
    tag { sample_id }
    publishDir "${params.outdir}/y_extract", mode: 'copy'

    cpus params.threads

    input:
    tuple val(sample_id), path(filtered_fastq_file), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.chrY.fastq"), emit: y_fastq

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${ref_fasta} ${filtered_fastq_file} \
      | samtools sort -@ ${task.cpus} -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    samtools view -@ ${task.cpus} -b ${sample_id}.bam ${params.contig_y} > ${sample_id}.chrY.bam
    samtools fastq ${sample_id}.chrY.bam > ${sample_id}.chrY.fastq
    """
}

process AlignYReads {
    tag { sample_id }
    publishDir "${params.outdir}/aligned", mode: 'copy'

    cpus params.threads

    input:
    tuple val(sample_id), path(yfastq), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.chrY.aligned.bam"), path("${sample_id}.chrY.aligned.bam.bai"), emit: aligned

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${ref_fasta} ${yfastq} \
      | samtools sort -@ ${task.cpus} -o ${sample_id}.chrY.aligned.bam
    samtools index ${sample_id}.chrY.aligned.bam
    """
}

process CallVariants {
    tag { sample_id }
    publishDir "${params.outdir}/variants", mode: 'copy'

    cpus params.threads

    input:
    tuple val(sample_id), path(bam), path(bai), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.sv.vcf.gz"), path("${sample_id}.sv.vcf.gz.tbi"), emit: sv_vcfs
    tuple val(sample_id), path("${sample_id}.snps.vcf.gz"), path("${sample_id}.snps.vcf.gz.tbi"), emit: snp_vcfs

    script:
    """
    # Structural variants with Sniffles2
    sniffles --input ${bam} --vcf ${sample_id}.sv.vcf --reference ${ref_fasta} --threads ${task.cpus}
    bgzip -f ${sample_id}.sv.vcf
    tabix -f -p vcf ${sample_id}.sv.vcf.gz

    # SNP/indel calling with medaka directly from BAM
    medaka_variant -i ${bam} -f ${ref_fasta} -o ${sample_id}.snps.vcf --threads ${task.cpus}
    bgzip -f ${sample_id}.snps.vcf
    tabix -f -p vcf ${sample_id}.snps.vcf.gz
    """
}

process MergeVariants {
    tag { "merge_snps" }
    publishDir "${params.outdir}/variants", mode: 'copy'

    cpus 2

    input:
    val vcf_files

    output:
    path "merged_snps.vcf.gz"
    path "merged_snps.vcf.gz.tbi"

    script:
    """
    bcftools merge -m none -O z -o merged_snps.vcf.gz ${vcf_files.join(' ')}
    tabix -f -p vcf merged_snps.vcf.gz
    """
}

workflow {
    // Build sample tuples and constant channels
    Channel
        .fromPath(params.input_fastq)
        .map { File fq -> tuple(inferSampleId(fq), fq) }
        .set { fastqs }

    def ref_fasta_ch = Channel.fromPath(params.reference_fasta)

    // 1) Filter reads
    filtered = FilterReads(fastqs).filtered

    // 2) Extract chrY reads (pair each sample with the reference)
    def extract_in = filtered.combine(ref_fasta_ch).map { tup, ref -> tuple(tup[0], tup[1], ref) }
    y_fastq = ExtractYReads(extract_in).y_fastq

    // 3) Align chrY reads
    def align_in = y_fastq.combine(ref_fasta_ch).map { tup, ref -> tuple(tup[0], tup[1], ref) }
    aligned = AlignYReads(align_in).aligned

    // 4) Call variants (SV + SNP/indel)
    def call_in = aligned.combine(ref_fasta_ch).map { tup, ref -> tuple(tup[0], tup[1], tup[2], ref) }
    variant_outputs = CallVariants(call_in)
    sv_vcfs   = variant_outputs.sv_vcfs.map { sid, sv, s_tbi -> sv }
    snp_vcfs  = variant_outputs.snp_vcfs.map { sid, snp, s_tbi -> snp }

    // 5) Merge SNP VCFs across samples
    snp_vcfs
        .collect()
        | MergeVariants
}