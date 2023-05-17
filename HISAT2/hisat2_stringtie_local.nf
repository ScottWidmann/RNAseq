/*
 * pipeline input parameters
 */
params.reads = "$projectDir/*_{1,2}.fastq.gz"
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
    publishDir params.outdir, mode:'copy'
    tag "Hisat2 on $sample_id"
    

    input:
    tuple val(sample_id), path(reads)
    val index
    path reference

    output:
    tuple path("${sample_id}/${sample_id}.bam"), path("${sample_id}/${sample_id}.gtf"), path("${sample_id}/${sample_id}_gene_abund.tab"), path("${sample_id}/${sample_id}.txt")

    script:
    """
    hisat2 -q -t --met-file ${sample_id}.txt -p $task.cpus -x $index -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam
    samtools sort -@ $task.cpus ${sample_id}.sam -o ${sample_id}.bam
    rm *.sam
    stringtie -A ${sample_id}_gene_abund.tab -p $task.cpus -e -G $reference -o ${sample_id}.gtf ${sample_id}.bam
    mkdir ${sample_id}
    mv ${sample_id}.gtf ${sample_id}
    mv ${sample_id}_gene_abund.tab ${sample_id}
    mv ${sample_id}.txt ${sample_id}
    mv ${sample_id}.bam ${sample_id}
    
    """
}

process PREP {
    publishDir params.outdir, mode:'move'

    input:
    val results

    output:
    tuple path("gene_count_matrix.csv"), path("transcript_count_matrix.csv")

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

    alignment_ch = ALIGNMENT(read_pairs_ch, params.index, params.reference)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(fastqc_ch.collect())
    PREP(alignment_ch.collect())	
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Uh Oh" )
}

