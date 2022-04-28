#!/usr/bin/env nextflow

/*
==================
PIPELINE STRUCTURE
==================
Illumina reads
    |- Fastp --> kraken2 --> bracken --> krakentools to filter host reads --> kraken_biom

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
//
params.cpus = 16
params.kraken2_db = "$baseDir/reference_db/k2_standard_16gb_20201202.tar.gz"
params.bowtie2_reference_index = "$baseDir/reference_db/bowtie2/bt2_index.tar.gz"
params.bowtie2_mito_index = "$baseDir/reference_db/bowtie2/bt2_index.tar.gz"

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
         bt2 mito db  : ${params.bowtie2_mito_index}
         """
         .stripIndent()

// SET UP INPUT CHANNELS
bowtie2_db_ch = Channel.value(file("${params.bowtie2_reference_index}"))
bowtie2_mitodb_ch = Channel.value(file("${params.bowtie2_mito_index}"))
kraken2_db_ch = Channel.value(file("${params.kraken2_db}"))

reads = Channel.fromFilePairs(params.reads)
reads.into { fastp_reads; second_line }

/*
==================
FASTP
==================
*/
process fastp {
    publishDir "$params.outdir/fastp/json/", mode: 'copy', pattern: '*_fastp.json'
    publishDir "$params.outdir/fastp/html/", mode: 'copy', pattern: '*_fastp.html'
    container 'biocontainers/fastp:v0.20.1_cv1'
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), file(reads_file) from fastp_reads

    output:
    tuple val(sample_id), file('*.fq.gz') into kraken2_direct_ch
    file '*fastp.{json,html}'

    script:
    """
    fastp \
    -i ${reads_file[0]} \
    -I ${reads_file[1]} \
    -o ${sample_id}.fastp.R1.fq.gz \
    -O ${sample_id}.fastp.R2.fq.gz \
    -j ${sample_id}.fastp.json \
    -h ${sample_id}.fastp.html \
    -w ${task.cpus}
    """
}

/*
==================
RAW READ FILES AGAINST KRAKEN2
==================
*/
process kraken2_bracken_direct {
    publishDir "$params.outdir/kraken2/report/", mode: 'copy', pattern: '*_kraken2.report'
    publishDir "$params.outdir/kraken2/bracken/", mode: 'copy', pattern: '*_bracken.tsv'
    publishDir "$params.outdir/kraken2/filtered_bracken/", mode: 'copy', pattern: '*.filtered.bracken'
    container 'shaunchuah/kraken_bracken'
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), file(reads_file) from kraken2_direct_ch
    file kraken2_db from kraken2_db_ch

    output:
    tuple val(sample_id), file('*_kraken2.report')
    tuple val(sample_id), file('*_bracken.tsv')
    tuple val(sample_id), file('*.filtered.bracken') into kraken_biom_ch
    file '*_kraken2.output'

    script:
    """
    tar -xvf $kraken2_db
    kraken2 \
        --db . \
        --report ${sample_id}_kraken2.report \
        --output ${sample_id}_kraken2.output \
        --use-names \
        --threads ${task.cpus} \
        ${reads_file[0]}

    bracken \
        -d . \
        -i ${sample_id}_kraken2.report \
        -o ${sample_id}_bracken.tsv
    
    python /krakentools/filter_bracken.out.py -i ${sample_id}_bracken.tsv -o ${sample_id}.filtered.bracken --exclude 9606
    """
}

process convert_kraken_to_biom {
    publishDir "$params.outdir/kraken2/biom/", mode: 'copy'
    container 'shaunchuah/kraken_biom'

    input:
    file(kraken2_report_files) from kraken_biom_ch.collect()

    output:
    file 'collated_kraken2_direct.biom'

    script:
    """
    kraken-biom ${kraken2_report_files} -o collated_kraken2_direct.biom
    """
}

