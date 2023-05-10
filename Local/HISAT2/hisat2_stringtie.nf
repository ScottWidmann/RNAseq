/*
 * pipeline input parameters
 */
params.reads = "$projectDir/*_R{1,2}.fastq.gz"
params.index = "$projectDir/grch38_snp_tran/genome_snp_tran"
params.reference = "$projectDir/Homo_sapiens.GRCh38.109.gtf"
params.outdir = "results"
log.info """\
    Hisat2 > Stringtie Pipeline
    ===================================
    index	 : ${params.index}
    refernce     : ${params.reference}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()


process ALIGNMENT {
    tag "Hisat2 on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)
    val index
    

    output:
    path "${sample_id}.bam"

    script:
    """
    hisat2 -q -t --met-file ${sample_id}.txt -p $task.cpus -x $index -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam
    samtools sort -@ $task.cpus ${sample_id}.sam -o ${sample_id}.bam
    """
}

process QUANTIFY {
    tag "Stringtie on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path bam
    path reference
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}.gtf"
    path "${sample_id}_gene_abund.tab"

    script:
    """
    stringtie -A ${sample_id}_gene_abund.tab -p $task.cpus -e -G $reference -o ${sample_id}.gtf ${sample_id}.bam
    """
}

process MKDIR_MV {
    publishDir params.outdir, mode:'move'
    
    input:
    tuple val(sample_id), path(reads)
    path gtf
    path abund

    output:
    path "${sample_id}"
    path "${sample_id}/${sample_id}.gtf"
    path "${sample_id}/${sample_id}_gene_abund.tab"

    script:
    """
    mkdir ${sample_id}
    mv ${sample_id}.gtf ${sample_id}
    mv ${sample_id}_gene_abund.tab ${sample_id}
    """
}

process PREP {
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)
    path sample
    path gtf
    path abund

    output:
    path "gene_count_matrix.csv"

    script:
    """
    
    python3 $projectDir/prepDE.py3 -i $projectDir/results
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

    alignment_ch = ALIGNMENT(read_pairs_ch, params.index)
    quantify_ch = QUANTIFY(alignment_ch, params.reference, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(alignment_ch.mix(fastqc_ch).collect())
    mkdir_mv_ch = MKDIR_MV(read_pairs_ch, quantify_ch)
    PREP(read_pairs_ch, mkdir_mv_ch)	
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Uh Oh" )
}

