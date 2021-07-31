#!/usr/bin/env nextflow

/*
PIPELINE INSTRUCTIONS HERE

This pipeline has been constructured for illumina sequencing reads
Sample Folder structure <path for reads here>/<sample_id>/<contains fastq.gz read1 and read2>

nextflow run pipeline.nf --reads <s3 path to reads folder> --outdir <s3 path to reports folder>
*/

// PIPELINE PARAMETERS HERE

// Input Files
params.reads = "$baseDir/data/*/*_{R1,R2}_*.fastq.gz"

// Report Directory
params.outdir = 'reports'

println """\
         ===================================
         C F D N A - N F   P I P E L I N E    
         ===================================
         PIPELINE PARAMETERS
         -----------------------------------
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

reads = Channel.fromFilePairs(params.reads)
reads.into { fastqc_reads; reads_for_alignment }


process fastqc_run {
    publishDir "$params.outdir/fastqc/$sample_id/", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'
    tag "$sample_id - FastQC"
    machineType 'e2-highcpu-8'
    cpus 4

    input:
    tuple val(sample_id), file(reads_file) from fastqc_reads

    output:
    file '*.zip' into multiqc_ch

    script:
    """
    fastqc $reads_file -o . --threads ${task.cpus}
    """
}

process multiqc {
    publishDir "$params.outdir/multiqc/", mode: 'copy'
    container 'ewels/multiqc:latest'

    input:
    path '*' from multiqc_ch.collect()

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}
