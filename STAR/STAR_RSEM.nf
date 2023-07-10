/*
 * pipeline input parameters
 */
params.reads = "$projectDir/*_{1,2}.fastq.gz"
params.genome = "$projectDir/GRCh38.primary_assembly.genome.fa"
params.gtf = "$projectDir/gencode.v43.primary_assembly.annotation.gtf"
params.outdir = "results"

log.info """\
    STAR > RSEM  Pipeline
    ===================================
    genome       : ${params.genome}
    gtf          : ${params.gtf}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()


process INDEX {
    tag "STAR index build"
    
    input:
    path genome
    path gtf
    
    output:
    path 'genome_dir'

    script:
    """
    STAR --runMode genomeGenerate --runThreadN $task.cpus --genomeDir genome_dir --genomeFastaFiles $genome --sjdbGTFfile $gtf
    """
}


process ALIGNMENT {
    tag "Alignment of $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)
    val genome_dir
    path rsem_index
    
    output:
    tuple path("${sample_id}.Aligned.toTranscriptome.out.bam"), path("${sample_id}.genes.results")
    
    script:
    """
    STAR --runThreadN $task.cpus --genomeDir $genome_dir --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix ${sample_id}. --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts --readFilesCommand gunzip -c
    
    rsem-calculate-expression -p $task.cpus --paired-end --alignments ${sample_id}.Aligned.toTranscriptome.out.bam rsem_index/$rsem_index ${sample_id}

    """
}

process RSEM_INDEX {
    tag "RSEM index build"
    
    input:
    path genome
    path gtf

    output:
    path 'rsem_index'

    script:
    """
    mkdir rsem_index
    rsem-prepare-reference --gtf $gtf $genome rsem_index/rsem_index
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -t $task.cpus -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    
    index_ch = INDEX(params.genome, params.gtf)
    rsem_index_ch = RSEM_INDEX(params.genome, params.gtf)
    alignment_ch = ALIGNMENT(read_pairs_ch, index_ch, rsem_index_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(alignment_ch.mix(fastqc_ch).collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}

