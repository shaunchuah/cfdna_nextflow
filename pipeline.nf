#!/usr/bin/env nextflow

/*
==================
PIPELINE STRUCTURE
==================
Illumina reads
    |- Fastqc --> MultiQC
    |- Kraken2 --> Kraken_biom
    |- Bowtie2 against GRCh38 no-alt https://benlangmead.github.io/aws-indexes/bowtie --> samtools ChrM, samtools flagstat

==================
PIPELINE INSTRUCTIONS
==================
This pipeline has been constructured for illumina sequencing reads
Sample Folder structure <path for reads here>/<sample_id>/<contains fastq.gz read1 and read2>

nextflow run pipeline.nf --reads <s3/az/gc path to reads folder> --outdir <s3/az/gc path to reports folder>
alternatively open up the config for the profiles and you can run
nextflow run pipeline.nf -resume -profile az

Azure VM Reference:
Standard_E8a_v4 8cpus 64gb ram
Standard_D16_v3 16cpus 64gb ram

Top tip: For azure storage - it does not support folders.
You have to put the file directly in the root of the container otherwise azure throws errors.
So I have zipped up the bowtie2 reference index and placed it at az://<container>/<tar.gz bowtie2 index>
*/

// PIPELINE PARAMETERS HERE (OVERRIDE WITH CONFIG)
// Input Files
params.reads = "$baseDir/data/*/*_{R1,R2}_*.fastq.gz"
// Report Directory
params.outdir = 'reports'
// Reference Genomes
// params.bowtie2_reference_index = "$baseDir/reference_db/bowtie2/bt2_index.tar.gz"
// CPU configuration
params.cpus = 4
// Kraken2 Database
params.kraken2_db = "$baseDir/reference_db/k2_standard_16gb_20201202.tar.gz"
// Reference Genomes
params.bowtie2_reference_index = "$baseDir/reference_db/bowtie2/bt2_index.tar.gz"
bowtie2_db_ch = Channel.value(file("${params.bowtie2_reference_index}"))


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
         kraken2 db   : ${params.kraken2_db}
         bowtie2 db   : ${params.bowtie2_reference_index}
         """
         .stripIndent()

// SET UP INPUT CHANNELS
// bowtie2_db_ch = Channel.value(file("${params.bowtie2_reference_index}"))
kraken2_db_ch = Channel.value(file("${params.kraken2_db}"))
reads = Channel.fromFilePairs(params.reads)

reads.into { fastqc_reads; reads_for_direct_kraken2; reads_for_alignment}

process fastqc_run {
    publishDir "$params.outdir/fastqc/$sample_id/", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'
    tag "$sample_id - FastQC"
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), file(reads_file) from fastqc_reads

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc $reads_file -o . --threads ${task.cpus}
    """
}


process run_kraken2_direct {
    publishDir "$params.outdir/kraken2/report/", mode: 'copy', pattern: "*_kraken2report.txt"
    publishDir "$params.outdir/kraken2/output/", pattern: "*_kraken2output.txt"
    container 'staphb/kraken2:2.1.2-no-db'
    cpus "$params.cpus".toInteger()
    tag "$sample_id - kraken2_direct"

    input:
    tuple val(sample_id), file(reads_file) from reads_for_direct_kraken2
    file kraken2_db from kraken2_db_ch

    output:
    file '*_kraken2report.txt' into kraken_biom_ch
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
        ${reads_file[0]}
    """
}

process kraken_biom {
    publishDir "$params.outdir/kraken2/biom/", mode: 'copy'
    container 'shaunchuah/kraken_biom'

    input:
    file(kraken2_report_files) from kraken_biom_ch.collect()

    output:
    file "collated_kraken2.biom"

    script:
    """
    kraken-biom ${kraken2_report_files} -o collated_kraken2.biom
    """
}


process bowtie2 {
    container 'biocontainers/bowtie2:v2.4.1_cv1'
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), file(reads_file) from reads_for_alignment
    file db from bowtie2_db_ch

    output:
    tuple val(sample_id), file('*.sam') into stats_ch, chrM_ch

    script:
    """
    tar -xvf $db
    bowtie2 -t -p ${task.cpus} -x GRCh38_noalt_as -1 ${reads_file[0]} -2 ${reads_file[1]} -S ${sample_id}.sam 
    """
}

process samtools_chrM {
    publishDir "$params.outdir/samtools_chrM/", mode: 'copy', pattern: "*_chrM.txt"
    publishDir "$params.outdir/samtools_idxstats/", mode: 'copy', pattern: "*_idxstats.tsv"
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus "$params.cpus".toInteger()
    tag "$sample_id"

    input:
    tuple val(sample_id), file(sam_file) from chrM_ch

    output:
    path "${sample_id}_chrM.txt"
    path "${sample_id}_idxstats.tsv"

    script:
    """
    samtools view -@ ${task.cpus} -S -b $sam_file | samtools sort - > temp.bam
    samtools index -@ ${task.cpus} temp.bam

    echo "${sample_id} - Total Reads vs Mitochondrial Reads" > ${sample_id}_chrM.txt
    echo "Total Reads:" >> ${sample_id}_chrM.txt
    samtools view -@ ${task.cpus} -c temp.bam >> ${sample_id}_chrM.txt

    echo "Mitochondrial Reads:" >> ${sample_id}_chrM.txt
    samtools view -@ ${task.cpus} -c temp.bam chrM >> ${sample_id}_chrM.txt

    samtools idxstats -@ ${task.cpus} temp.bam > ${sample_id}_idxstats.tsv
    """
}

process samtools_flagstat {
    publishDir "$params.outdir/samtools_flagstat/"
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus "$params.cpus".toInteger()
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

/*
process get_mapped_reads {
    publishDir "$params.outdir/mapped_bam_igv/"
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus "$params.cpus".toInteger()
    tag "$sample_id"

    input:
    tuple val(sample_id), file(sam_file) from igv_ch

    output:
    path "${sample_id}_mapped_reads.bam"
    path "${sample_id}_mapped_reads.bam.bai"

    script:
    """
    samtools view -@ ${task.cpus} -bF 4 $sam_file | samtools sort - > ${sample_id}_mapped_reads.bam
    samtools index -@ ${task.cpus} ${sample_id}_mapped_reads.bam
    """
}


process get_unmapped_reads {
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus "$params.cpus".toInteger()
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
*/
