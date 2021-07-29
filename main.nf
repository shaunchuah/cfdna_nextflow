reads = Channel.fromPath( 's3://scgenomics.reads/test/*.fastq.gz' )
reads.into { fastqc_reads; reads_for_alignment }

params.outdir = 'reports'
params.minimap_reference_db = "$baseDir/databases/hg38.mmi"
params.scratchdir = 'intermediate'
params.bowtie2_human_db = "$baseDir/databases/bowtie2_index/"
params.kraken2_db = "$baseDir/databases/k2_standard_8gb_20210517/"

process fastqc_run {
    publishDir "$params.outdir/fastqc/", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'
    tag "$sample"

    input:
    file sample from fastqc_reads

    output:
    file '*' into multiqc_ch

    script:
    """
    fastqc $sample -o .
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

/*
process bowtie2 {
    container 'biocontainers/bowtie2:v2.4.1_cv1'

    input:
    file sample from reads_for_alignment
    path bowtie2_index from params.bowtie2_human_db

    output:
    file "${sample}_bowtie2.sam"

    script:
    """
    bowtie2 -x $bowtie2_index -U $sample -S ${sample}_bowtie2.sam
    """
}
*/

process minimap {
    container 'staphb/minimap2'
    memory '16 GB'

    input:
    file sample from reads_for_alignment
    path minimap_reference_db from params.minimap_reference_db

    output:
    file "${sample}_minimap2.sam" into align_ch

    script:
    """
    gzip -d $sample --stdout > temp.fastq
    minimap2 -ax sr $minimap_reference_db temp.fastq > ${sample}_minimap2.sam
    """
}

process get_unmapped_reads {
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    publishDir "$params.scratchdir/unmapped/"

    input:
    file sam_file from align_ch

    output:
    file "${sam_file}_unmapped.fastq" into kraken2_ch

    script:
    """
    samtools view -bf 4 ${sam_file} > aln.bam
    samtools fastq aln.bam > "${sam_file}_unmapped.fastq"
    """
}

process kraken2 {
    container 'staphb/kraken2'
    publishDir "$params.outdir/kraken2/", mode: 'copy'
    memory '16 GB'

    input: 
    file unmapped_reads from kraken2_ch
    path kraken2_db from params.kraken2_db

    output:
    file "${unmapped_reads}_k2report.txt"
    file "${unmapped_reads}_k2.out"

    script:
    """
    kraken2 --db $kraken2_db --report ${unmapped_reads}_k2report.txt --output ${unmapped_reads}_k2.out --threads $task.cpus --use-names $unmapped_reads
    """

}
