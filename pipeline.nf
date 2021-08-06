#!/usr/bin/env nextflow

/*
PIPELINE INSTRUCTIONS HERE

This pipeline has been constructured for illumina sequencing reads
Sample Folder structure <path for reads here>/<sample_id>/<contains fastq.gz read1 and read2>

nextflow run pipeline.nf --reads <s3/az/gc path to reads folder> --outdir <s3/az/gc path to reports folder>
alternatively open up the config for the profiles and you can run
nextflow run pipeline.nf -resume -profile az
*/

// PIPELINE PARAMETERS HERE

// Input Files
params.reads = "$baseDir/data/*/*_{R1,R2}_*.fastq.gz"

// Report Directory
params.outdir = 'reports'

// Reference Genomes
params.bowtie2_reference_index = "$baseDir/reference_db/bowtie2/bt2_index.tar.gz"
bowtie2_db_ch = Channel.value(file("${params.bowtie2_reference_index}"))
/* 
Top tip: For azure storage - it does not support folders.
You have to put the file directly in the root of the container otherwise azure throws errors.
So I have zipped up the bowtie2 reference index and placed it at az://<container>/<tar.gz bowtie2 index>
*/

// Kraken2 Database
// By default, the full database is linked here which requires 64GB of RAM
// URL for the 8GB DB = "s3://genome-idx/kraken/k2_standard_8gb_20210517.tar.gz"
params.kraken2_db = "s3://genome-idx/kraken/k2_standard_20210517.tar.gz"
kraken2_db_ch = Channel.value(file("${params.kraken2_db}"))

println """\
         ===================================
         C F D N A - N F   P I P E L I N E    
         ===================================
         PIPELINE PARAMETERS
         -----------------------------------
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         -----------------------------------
         PIPELINE REFERENCE DATABASES
         -----------------------------------
         bowtie2 db   : ${params.bowtie2_reference_index}
         kraken2 db   : ${params.kraken2_db}
         """
         .stripIndent()

reads = Channel.fromFilePairs(params.reads)
reads.into { fastqc_reads; reads_for_alignment }


process fastqc_run {
    publishDir "$params.outdir/fastqc/$sample_id/", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'
    tag "$sample_id - FastQC"
    cpus 16

    input:
    tuple val(sample_id), file(reads_file) from fastqc_reads

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc $reads_file -o . --threads ${task.cpus}
    """
}

/*
process multiqc {
    publishDir "$params.outdir/multiqc/", mode: 'copy'
    container 'ewels/multiqc:v1.11'

    input:
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"

    script:
    """
    multiqc .
    """
}
*/

/*
process bwa_mem {
    publishDir "$params.outdir/bwa/"
    container 'biocontainers/bwa:v0.7.17_cv1'
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), file(reads_file) from reads_for_alignment
    file reference_genome from human_genome_reference_ch

    output:
    tuple val(sample_id), file('*.sam') into aligned_ch, stats_ch

    script:
    """
    bwa index ${reference_genome}/genome.fa
    bwa mem -t ${task.cpus} ${reference_genome}/genome.fa ${reads_file[0]} ${reads_file[1]} > ${sample_id}.sam 
    """
}
*/

process bowtie2 {
    container 'biocontainers/bowtie2:v2.4.1_cv1'
    cpus 16

    input:
    tuple val(sample_id), file(reads_file) from reads_for_alignment
    file db from bowtie2_db_ch

    output:
    tuple val(sample_id), file('*.sam') into aligned_ch, stats_ch

    script:
    """
    tar -xvf $db
    bowtie2 -t -p ${task.cpus} -x bowtie2/GRCh38_bowtie2 -1 ${reads_file[0]} -2 ${reads_file[1]} -S ${sample_id}.sam 
    """
}

process samtools_flagstat {
    publishDir "$params.outdir/samtools_flagstat/"
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus 16
    tag "$sample_id"

    input:
    tuple val(sample_id), file(sam_file) from stats_ch

    output:
    path "${sample_id}_flagstat.txt"

    script:
    """
    samtools flagstat -@ ${task.cpus} $sam_file > ${sample_id}_flagstat.txt
    """
}

process get_unmapped_reads {
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus 16
    tag "$sample_id"

    input:
    tuple val(sample_id), file(sam_file) from aligned_ch

    output:
    tuple val(sample_id), path("${sample_id}_unmapped_reads.fastq") into metaphlan_ch, kraken2_ch

    script:
    """
    samtools view -@ ${task.cpus} -bf 4 $sam_file | samtools fastq -@ ${task.cpus} - > ${sample_id}_unmapped_reads.fastq
    """
}

process run_metaphlan {
    publishDir "$params.outdir/metaphlan/"
    container 'biobakery/metaphlan:3.0.7'
    cpus 16
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

process run_kraken2 {
    publishDir "$params.outdir/kraken2/"
    container 'staphb/kraken2:2.1.2-no-db'
    cpus 16
    tag "$sample_id - kraken2"

    input:
    tuple val(sample_id), file(unmapped_fastq) from kraken2_ch
    file kraken2_db from kraken2_db_ch

    output:
    file '*_kraken2report.txt'
    file '*_kraken2output.txt'

    script:
    """
    tar -xvf $kraken2_db
    kraken2 \
        --db . \
        --report ${sample_id}_kraken2report.txt \
        --output ${sample_id}_kraken2output.txt \
        --use-names \
        --threads ${task.cpus} \
        $unmapped_fastq
    """
}
