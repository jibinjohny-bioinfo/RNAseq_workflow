nextflow.enable.dsl=2

// -----------------------------
// Parameters
// -----------------------------
params.outdir          = params.outdir ?: 'results'
params.genome          = params.genome ?: 'data/genome.fasta'
params.samplesheet     = params.samplesheet ?: 'data/samplesheet.csv'  // CSV with: sample,fq1,fq2
params.trinity_prefix  = params.trinity_prefix ?: 'gg_trinity'
params.max_intron      = params.max_intron ?: 10000
params.ref_proteins    = params.ref_proteins ?: 'data/ref_proteins.fasta'
params.busco_db        = params.busco_db ?: 'insecta_odb10'             // lineage name or path
params.featurecounts_gtf = params.featurecounts_gtf ?: null              // optional: GTF/GFF

// Default resources
params.max_cpus   = params.max_cpus ?: 16
params.max_memory = params.max_memory ?: '100.GB'

// -----------------------------
// Helper functions
// -----------------------------
// Parse samplesheet CSV (header: sample,fq1,fq2)
// Returns tuples: [ sample_id, path(fq1), path(fq2) ]
def parseSamplesheet(file) {
    Channel
        .fromPath(file)
        .splitCsv(header:true)
        .map { row ->
            def sample = row.sample as String
            def fq1 = file(row.fq1 as String)
            def fq2 = file(row.fq2 as String)
            if (!fq1.exists()) error "Missing FASTQ file: ${fq1}"
            if (!fq2.exists()) error "Missing FASTQ file: ${fq2}"
            tuple(sample, fq1, fq2)
        }
}

// -----------------------------
// Channels
// -----------------------------
Channel
    .fromPath(params.genome)
    .set { ch_genome }

Channel
    .fromPath(params.ref_proteins)
    .set { ch_ref_proteins }

parseSamplesheet(params.samplesheet)
    .set { ch_reads }

// -----------------------------
// Processes
// -----------------------------

process FASTQC {
    tag "${sample}"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy', overwrite: true

    cpus 2
    memory '2.GB'
    time '2h'
    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(sample), path(fq1), path(fq2)

    output:
    path "${sample}_fastqc", emit: qcdir

    script:
    """
    mkdir -p ${sample}_fastqc
    fastqc -t ${task.cpus} -o ${sample}_fastqc ${fq1} ${fq2}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    cpus 2
    memory '2.GB'
    time '2h'
    conda 'bioconda::multiqc=1.19'

    input:
    path qcdirs

    output:
    path 'multiqc_report.html'
    path 'multiqc_data'

    script:
    def qc_space = qcdirs.collect { it.getName() }.join(' ')
    """
    multiqc -o ./ ${qc_space}
    """
}

process HISAT2_BUILD_INDEX {
    tag 'hisat2_index'
    publishDir "${params.outdir}/genome_index", mode: 'copy', overwrite: true

    cpus params.max_cpus
    memory '8.GB'
    time '6h'
    conda 'bioconda::hisat2=2.2.1'

    input:
    path genome_fa

    output:
    path 'hisat2_index', emit: index_dir

    script:
    """
    mkdir -p hisat2_index
    hisat2-build -p ${task.cpus} ${genome_fa} hisat2_index/genome
    """
}

process HISAT2_ALIGN_SORT_INDEX {
    tag "${sample}"
    publishDir params.outdir, mode: 'copy', overwrite: true

    cpus params.max_cpus
    memory '16.GB'
    time '24h'
    conda 'bioconda::hisat2=2.2.1,bioconda::samtools=1.20'

    input:
    tuple val(sample), path(fq1), path(fq2), path(index_dir)

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), path("${sample}_hisat2.log")

    script:
    """
    hisat2 -p ${task.cpus} -x ${index_dir}/genome -1 ${fq1} -2 ${fq2} 2> ${sample}_hisat2.log \
      | samtools sort -@ ${task.cpus} -o ${sample}.bam
    samtools index -@ ${task.cpus} ${sample}.bam
    """
}

process TRINITY_GENOME_GUIDED {
    tag 'trinity_gg'
    publishDir "${params.outdir}/${params.trinity_prefix}_assembly", mode: 'copy', overwrite: true

    cpus params.max_cpus
    memory params.max_memory
    time '72h'
    conda 'bioconda::trinity=2.15.1,bioconda::samtools=1.20'

    input:
    // List of BAMs staged into work dir
    path bam_list

    output:
    path 'Trinity-GG.fasta'
    path 'trinity_out_dir'

    script:
    def bam_csv = bam_list.collect { it.getName() }.join(',')
    """
    Trinity --genome_guided_bam ${bam_csv} \
            --genome_guided_max_intron ${params.max_intron} \
            --CPU ${task.cpus} --max_memory ${task.memory.toGiga()}G \
            --output trinity_out_dir
    cp trinity_out_dir/Trinity-GG.fasta ./
    """
}

process TRANSDECODER {
    tag 'transdecoder'
    publishDir "${params.outdir}/${params.trinity_prefix}_assembly", mode: 'copy', overwrite: true

    cpus 8
    memory '16.GB'
    time '24h'
    conda 'bioconda::transdecoder=5.7.1'

    input:
    path trinity_fa

    output:
    path 'Trinity-GG.fasta.transdecoder_dir'
    path 'Trinity-GG.fasta.transdecoder.pep'
    path 'Trinity-GG.fasta.transdecoder.bed'

    script:
    """
    TransDecoder.LongOrfs -t ${trinity_fa}
    TransDecoder.Predict -t ${trinity_fa}
    """
}

process BUSCO_TRANSCRIPTOME {
    tag 'busco'
    publishDir "${params.outdir}/busco", mode: 'copy', overwrite: true

    cpus params.max_cpus
    memory '32.GB'
    time '48h'
    conda 'bioconda::busco=5.7.1'

    input:
    path trinity_fa

    output:
    path 'busco_out'

    script:
    """
    busco -i ${trinity_fa} -l ${params.busco_db} -o busco_out -m transcriptome -c ${task.cpus} --download_path .
    """
}

process MAKEBLASTDB {
    tag 'makeblastdb'
    publishDir "${params.outdir}/blastdb", mode: 'copy', overwrite: true

    cpus 2
    memory '2.GB'
    time '6h'
    conda 'bioconda::blast=2.14.1'

    input:
    path ref_prot

    output:
    path 'prot_db'

    script:
    """
    mkdir -p prot_db
    makeblastdb -in ${ref_prot} -dbtype prot -out prot_db/ref_protein_db
    """
}

process BLASTX_ANNOTATION {
    tag 'blastx'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    cpus params.max_cpus
    memory '16.GB'
    time '48h'
    conda 'bioconda::blast=2.14.1'

    input:
    path trinity_fa
    path dbdir

    output:
    path 'blastx.outfmt6'

    script:
    """
    blastx -query ${trinity_fa} -db ${dbdir}/ref_protein_db \
           -out blastx.outfmt6 -evalue 1e-5 -num_threads ${task.cpus} -max_target_seqs 1 -outfmt 6
    """
}

process FEATURECOUNTS {
    tag 'featurecounts'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    when:
    params.featurecounts_gtf != null && params.featurecounts_gtf as String

    cpus params.max_cpus
    memory '16.GB'
    time '24h'
    conda 'bioconda::subread=2.0.6'

    input:
    path gtf
    path bam_files

    output:
    path 'counts.txt'

    script:
    def bam_space = bam_files.collect { it.getName() }.join(' ')
    """
    featureCounts -T ${task.cpus} -a ${gtf} -o counts.txt ${bam_space}
    """
}

// -----------------------------
// Workflow
// -----------------------------
workflow {
    main:
        // QC
        FASTQC(ch_reads)
        MULTIQC(FASTQC.out.qcdir.collect())

        // Genome index and alignments
        HISAT2_BUILD_INDEX(ch_genome)
        def reads_with_index = ch_reads.cross(HISAT2_BUILD_INDEX.out.index_dir)
        HISAT2_ALIGN_SORT_INDEX(reads_with_index)

        // Collect BAMs for downstream
        def bam_list = HISAT2_ALIGN_SORT_INDEX.out.map { it[1] }.collect()

        // Trinity + downstream
        TRINITY_GENOME_GUIDED(bam_list)
        TRANSDECODER(TRINITY_GENOME_GUIDED.out.filter { it.name == 'Trinity-GG.fasta' })
        BUSCO_TRANSCRIPTOME(TRINITY_GENOME_GUIDED.out.filter { it.name == 'Trinity-GG.fasta' })

        // BLASTx annotation
        MAKEBLASTDB(ch_ref_proteins)
        BLASTX_ANNOTATION(TRINITY_GENOME_GUIDED.out.filter { it.name == 'Trinity-GG.fasta' }, MAKEBLASTDB.out)

        // Optional featureCounts
        def ch_gtf = (params.featurecounts_gtf && (params.featurecounts_gtf as String)) ? Channel.fromPath(params.featurecounts_gtf) : Channel.empty()
        FEATURECOUNTS(ch_gtf, bam_list)
}