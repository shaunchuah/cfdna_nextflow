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

// Reference Genomes
params.human_genome_index = "s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/"
human_genome_index_ch = Channel.value(file("${params.human_genome_index}"))

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
    cpus 1

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
    cpus 1

    input:
    path '*' from multiqc_ch.collect()

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

process bwa_mem {
    publishDir "$params.outdir/bwa/"
    container 'biocontainers/bwa:v0.7.17_cv1'
    cpus 4

    input:
    tuple val(sample_id), file(reads_file) from reads_for_alignment
    file reference_genome_index from human_genome_index_ch

    output:
    tuple val(sample_id), file('*.sam') into aligned_ch, stats_ch

    script:
    """
    bwa mem -t ${task.cpus} ${reference_genome_index}/genome.fa ${reads_file[0]} ${reads_file[1]} > ${sample_id}.sam 
    """
}

process samtools_flagstat {
    publishDir "$params.outdir/samtools_flagstat/"
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus 1
    tag "$sample_id"

    input:
    tuple val(sample_id), file(sam_file) from stats_ch

    output:
    path "${sample_id}_flagstat.txt"

    script:
    """
    samtools flagstat $sam_file > ${sample_id}_flagstat.txt
    """
}

process get_unmapped_reads {
    publishDir "$params.outdir/temp/"
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus 1
    tag "$sample_id"

    input:
    tuple val(sample_id), file(sam_file) from aligned_ch

    output:
    tuple val(sample_id), path("${sample_id}_unmapped_reads.fastq") into metaphlan_ch

    script:
    """
    samtools view -bf 4 $sam_file | samtools fastq - > ${sample_id}_unmapped_reads.fastq
    """
}

process run_metaphlan {
    publishDir "$params.outdir/metaphlan/"
    container 'biobakery/metaphlan:3.0.7'
    cpus 2
    tag "$sample_id - metaphlan"
    
    input:
    tuple val(sample_id), file(unmapped_fastq) from metaphlan_ch

    output:
    path "${sample_id}_metaphlan_profile.txt"

    script:
    """
    metaphlan $unmapped_fastq --index latest --nproc ${task.cpus} --input_type fastq -o ${sample_id}_metaphlan_profile.txt
    """ 
}