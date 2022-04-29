#!/usr/bin/env nextflow

/*
==================
CFDNA_NEXTFLOW PIPELINE STRUCTURE
==================
Illumina reads
    |- Fastp --> kraken2 --> bracken --> krakentools to filter host reads --> kraken_biom
            |- MultiQC

==================
PIPELINE INSTRUCTIONS
==================
This pipeline has been constructured for illumina sequencing reads
Sample input folder structure <path for reads here>/<sample_id>/<contains fastq.gz read1 and read2>

nextflow run pipeline.nf -resume -profile test
nextflow run pipeline.nf -resume -profile production

Azure VM Reference:
Standard_E8a_v4 8cpus 64gb ram
Standard_D16_v3 16cpus 64gb ram

Top tip:
Azure storage does not support transferring of folders into containers in the traditional sense.
You have to put the file directly in the root of the container otherwise azure throws errors.
So I have zipped up the bowtie2 reference index and placed it at az://<container>/<tar.gz bowtie2 index>
Unzip the index once inside the container.
*/

// ==================
// PIPELINE PARAMETERS
// ==================
params.reads = "$baseDir/data/*/*_{R1,R2}_*.fastq.gz"
params.outdir = 'reports'
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
FASTP into MULTIQC
==================
Removes Illumina adapters and performs QC checks.
*/
process fastp {
    publishDir "$params.outdir/quality_control/fastp/json/", mode: 'copy', pattern: '*.fastp.json'
    publishDir "$params.outdir/quality_control/fastp/html/", mode: 'copy', pattern: '*.fastp.html'
    container 'biocontainers/fastp:v0.20.1_cv1'
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), file(reads_file) from fastp_reads

    output:
    tuple val(sample_id), file('*.fq.gz') into kraken2_direct_ch, reads_for_GRCh38, reads_for_mito
    file '*.fastp.json' into multiqc_ch
    file '*.fastp.html'

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

process multiqc {
    publishDir "$params.outdir/quality_control/", mode: 'copy'
    container 'shaunchuah/multiqc'

    input:
    file(fastp_files) from multiqc_ch.collect()

    output:
    path ('multiqc_report.html')

    script:
    """
    multiqc .
    """
}

/*
==================
KRAKEN2 DIRECT
==================
Runs all sequences into kraken2 to allow for seqeunces to classify against human
Kraken2 output then run into bracken to perform abundance estimation (default species level)
Bracken output file is then filtered to remove human reads to allow for alpha diversity analysis
*/
process kraken2_bracken_direct {
    publishDir "$params.outdir/kraken2/report/", mode: 'copy', pattern: '*_kraken2.report'
    publishDir "$params.outdir/kraken2/bracken/", mode: 'copy', pattern: '*_bracken.tsv'
    publishDir "$params.outdir/kraken2/filtered_bracken/", mode: 'copy', pattern: '*.filtered.bracken.tsv'
    container 'shaunchuah/kraken_bracken'
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), file(reads_file) from kraken2_direct_ch
    file kraken2_db from kraken2_db_ch

    output:
    tuple val(sample_id), file('*_kraken2.report')
    tuple val(sample_id), file('*_bracken.tsv')
    tuple val(sample_id), file('*.filtered.bracken.tsv')
    tuple val(sample_id), file('*_bracken_kraken2.report') into kraken_biom_ch
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
        --paired ${reads_file[0]} ${reads_file[1]}

    bracken \
        -d . \
        -i ${sample_id}_kraken2.report \
        -o ${sample_id}_bracken.tsv \
        -w ${sample_id}_bracken_kraken2.report \
        -r 100 \
        -l S \
        -t 5
    
    python /krakentools/filter_bracken.out.py -i ${sample_id}_bracken.tsv -o ${sample_id}.filtered.bracken.tsv --exclude 9606
    """
}

process convert_to_biom {
    publishDir "$params.outdir/kraken2/biom/", mode: 'copy'
    container 'shaunchuah/kraken_biom'

    input:
    file(kraken2_report_files) from kraken_biom_ch.collect()

    output:
    file 'collated_kraken_bracken.biom'

    script:
    """
    kraken-biom ${kraken2_report_files} -o collated_kraken_bracken.biom
    """
}

/*
==================
BOWTIE2 ALIGNMENT AGAINST GRCH38 & MITO
takes fastp output then aligns
deduplication using samblaster
sorted sam file as output
==================
*/

process bowtie2_grch38 {
    publishDir "$params.outdir/grch38/samtools_flagstat/", mode: 'copy', pattern: '*_flagstat.txt'
    container 'shaunchuah/bowtie2_samblaster_samtools'
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), path(reads_file) from reads_for_GRCh38
    path db from bowtie2_db_ch

    output:
    tuple val(sample_id), path('*.bam') into deeptools_ch
    path('*_flagstat.txt')

    script:
    """
    tar -xvf $db
    bowtie2 \
    -t -p ${task.cpus} \
    -x GRCh38_noalt_as/GRCh38_noalt_as \
    -1 ${reads_file[0]} \
    -2 ${reads_file[1]} | \
    samblaster | \
    samtools view -@ ${task.cpus} -b | \
    samtools sort -@ ${task.cpus} > ${sample_id}.bam

    samtools flagstat -@ ${task.cpus} ${sample_id}.bam > ${sample_id}_flagstat.txt
    """
}

process bowtie2_mito {
    publishDir "$params.outdir/mito/samtools_flagstat/", mode: 'copy', pattern: '*_flagstat.txt'
    container 'shaunchuah/bowtie2_samblaster_samtools'
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), path(reads_file) from reads_for_mito
    path db from bowtie2_mitodb_ch

    output:
    tuple val(sample_id), path('*.bam')
    path('*_flagstat.txt')

    script:
    """
    tar -xvf $db
    bowtie2 \
    -t -p ${task.cpus} \
    -x human_mito_db/human_mito_db \
    -1 ${reads_file[0]} \
    -2 ${reads_file[1]} | \
    samblaster | \
    samtools view -@ ${task.cpus} -b | \
    samtools sort -@ ${task.cpus} > ${sample_id}.bam

    samtools flagstat -@ ${task.cpus} ${sample_id}.bam > ${sample_id}_flagstat_mito.txt
    """
}

/*
==================
POST ALIGNMENT ANALYSIS STEPS - GRCh38
==================
*/

/*
==================
POST ALIGNMENT ANALYSIS STEPS - Mito
==================
*/
