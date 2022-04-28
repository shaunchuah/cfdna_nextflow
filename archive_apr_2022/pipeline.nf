#!/usr/bin/env nextflow

/*
==================
PIPELINE STRUCTURE
==================
Illumina reads
    |- Fastqc --> MultiQC
    |- Kraken2 --> Kraken_biom
    |- Bowtie2 against GRCh38 no-alt https://benlangmead.github.io/aws-indexes/bowtie
        |- samtools_chr_counts
        |- samtools flagstat
        |- get_unmapped_reads --> Kraken2_unmapped_only --> Kraken_biom_unmapped
    |- Bowtie2 against different mitochondrial DB https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1?report=fasta --> samtools_flagstat_mito

==================
PIPELINE TODO
==================
1. Add deduplication step with GATK Picard

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

params.bowtie2_mito_index = "$baseDir/reference_db/bowtie2/bt2_index.tar.gz"
bowtie2_mitodb_ch = Channel.value(file("${params.bowtie2_mito_index}"))

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
// bowtie2_db_ch = Channel.value(file("${params.bowtie2_reference_index}"))
kraken2_db_ch = Channel.value(file("${params.kraken2_db}"))
reads = Channel.fromFilePairs(params.reads)

reads.into { fastqc_reads; reads_for_direct_kraken2; reads_for_alignment; reads_for_mito_db; reads_for_fastp}

/*
==================
FASTQC
==================
*/
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

/*
==================
RAW READ FILES AGAINST KRAKEN2
==================
*/
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

/*
==================
BOWTIE2 ALIGNMENT AGAINST GRCH38
into chromosome counts
then split into unmapped reads before running Kraken2 against the unmapped reads 
==================
*/
process bowtie2 {
    container 'biocontainers/bowtie2:v2.4.1_cv1'
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), file(reads_file) from reads_for_alignment
    file db from bowtie2_db_ch

    output:
    tuple val(sample_id), file('*.sam') into stats_ch, chr_counts_ch, aligned_ch

    script:
    """
    tar -xvf $db
    bowtie2 -t -p ${task.cpus} -x GRCh38_noalt_as/GRCh38_noalt_as -1 ${reads_file[0]} -2 ${reads_file[1]} -S ${sample_id}.sam 
    """
}

process samtools_chr_counts {
    publishDir "$params.outdir/samtools_chr_counts/", mode: 'copy', pattern: "*_chr_counts.txt"
    publishDir "$params.outdir/samtools_idxstats/", mode: 'copy', pattern: "*_idxstats.tsv"
    publishDir "$params.outdir/mitochondrial_bam_files/", mode: 'copy', pattern: "*_chrM.bam"
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus "$params.cpus".toInteger()
    tag "$sample_id"

    input:
    tuple val(sample_id), file(sam_file) from chr_counts_ch

    output:
    path "${sample_id}_chr_counts.txt"
    path "${sample_id}_idxstats.tsv"
    path "${sample_id}_chrM.bam"

    script:
    """
    samtools view -@ ${task.cpus} -S -b $sam_file | samtools sort - > temp.bam
    samtools index -@ ${task.cpus} temp.bam

    echo "${sample_id} - Counts by Chromosome" > ${sample_id}_chr_counts.txt

    echo "total" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam >> ${sample_id}_chr_counts.txt
    echo "chrM" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chrM >> ${sample_id}_chr_counts.txt
    echo "chr1" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr1 >> ${sample_id}_chr_counts.txt
    echo "chr2" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr2 >> ${sample_id}_chr_counts.txt
    echo "chr3" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr3 >> ${sample_id}_chr_counts.txt
    echo "chr4" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr4 >> ${sample_id}_chr_counts.txt
    echo "chr5" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr5 >> ${sample_id}_chr_counts.txt
    echo "chr6" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr6 >> ${sample_id}_chr_counts.txt
    echo "chr7" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr7 >> ${sample_id}_chr_counts.txt
    echo "chr8" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr8 >> ${sample_id}_chr_counts.txt
    echo "chr9" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr9 >> ${sample_id}_chr_counts.txt
    echo "chr10" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr10 >> ${sample_id}_chr_counts.txt
    echo "chr11" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr11 >> ${sample_id}_chr_counts.txt
    echo "chr12" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr12 >> ${sample_id}_chr_counts.txt
    echo "chr13" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr13 >> ${sample_id}_chr_counts.txt
    echo "chr14" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr14 >> ${sample_id}_chr_counts.txt
    echo "chr15" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr15 >> ${sample_id}_chr_counts.txt
    echo "chr16" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr16 >> ${sample_id}_chr_counts.txt
    echo "chr17" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr17 >> ${sample_id}_chr_counts.txt
    echo "chr18" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr18 >> ${sample_id}_chr_counts.txt
    echo "chr19" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr19 >> ${sample_id}_chr_counts.txt
    echo "chr20" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr20 >> ${sample_id}_chr_counts.txt
    echo "chr21" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr21 >> ${sample_id}_chr_counts.txt
    echo "chr22" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chr22 >> ${sample_id}_chr_counts.txt
    echo "chrX" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chrX >> ${sample_id}_chr_counts.txt
    echo "chrY" >> ${sample_id}_chr_counts.txt
    samtools view -@ ${task.cpus} -c temp.bam chrY >> ${sample_id}_chr_counts.txt
    
    samtools idxstats -@ ${task.cpus} temp.bam > ${sample_id}_idxstats.tsv

    samtools view -@ ${task.cpus} -b temp.bam chrM > ${sample_id}_chrM.bam
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

*/

process get_unmapped_reads {
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus "$params.cpus".toInteger()
    tag "$sample_id"

    input:
    tuple val(sample_id), file(sam_file) from aligned_ch

    output:
    tuple val(sample_id), path("${sample_id}_unmapped_reads.fastq") into kraken2_unmapped_ch

    script:
    """
    samtools view -@ ${task.cpus} -bf 4 $sam_file | samtools fastq -@ ${task.cpus} - > ${sample_id}_unmapped_reads.fastq
    """
}

process kraken2_unmapped {
    publishDir "$params.outdir/kraken2_unmapped_only/report/", mode: 'copy', pattern: "*_kraken2report.txt"
    publishDir "$params.outdir/kraken2_unmapped_only/output/", pattern: "*_kraken2output.txt"
    container 'staphb/kraken2:2.1.2-no-db'
    cpus "$params.cpus".toInteger()
    tag "$sample_id - kraken2_direct"

    input:
    tuple val(sample_id), file(reads_file) from kraken2_unmapped_ch
    file kraken2_db from kraken2_db_ch

    output:
    file '*_unmapped_kraken2report.txt' into kraken_biom_unmapped_ch
    file '*_unmapped_kraken2output.txt'

    script:
    """
    tar -xvf $kraken2_db
    kraken2 \
        --db . \
        --report ${sample_id}_unmapped_kraken2report.txt \
        --output ${sample_id}_unmapped_kraken2output.txt \
        --use-names \
        --threads ${task.cpus} \
        ${reads_file[0]}
    """
}

// For use in PhyloSeq
process kraken_biom_unmapped {
    publishDir "$params.outdir/kraken2_unmapped_only/biom/", mode: 'copy'
    container 'shaunchuah/kraken_biom'

    input:
    file(kraken2_report_files) from kraken_biom_unmapped_ch.collect()

    output:
    file "collated_kraken2_unmapped_only.biom"

    script:
    """
    kraken-biom ${kraken2_report_files} -o collated_kraken2_unmapped_only.biom
    """
}

/*
==================
BOWTIE2 ALIGNMENT AGAINST MITO GENOME DIRECTLY
to check that using GRCH38 is accurate enough
into flagstat mito
==================
*/
process bowtie2_mitodb {
    container 'biocontainers/bowtie2:v2.4.1_cv1'
    cpus "$params.cpus".toInteger()

    input:
    tuple val(sample_id), file(reads_file) from reads_for_mito_db
    file db from bowtie2_mitodb_ch

    output:
    tuple val(sample_id), file('*.sam') into mito_stats_ch

    script:
    """
    tar -xvf $db
    bowtie2 -t -p ${task.cpus} -x human_mito_db/human_mito_db -1 ${reads_file[0]} -2 ${reads_file[1]} -S ${sample_id}.sam 
    """
}

process samtools_flagstat_mito {
    publishDir "$params.outdir/samtools_flagstat_mito/"
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    cpus "$params.cpus".toInteger()
    tag "$sample_id"

    input:
    tuple val(sample_id), file(sam_file) from mito_stats_ch

    output:
    path "${sample_id}_flagstat_mito.txt"

    script:
    """
    samtools flagstat -@ ${task.cpus} $sam_file > ${sample_id}_flagstat_mito.txt
    """
}
