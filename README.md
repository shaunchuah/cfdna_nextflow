# Nextflow pipeline for cfDNA analysis

Test illumina data included for pipeline development in /data folder credits to https://github.com/hartwigmedical/testdata  

## Installation

You need Docker and Nextflow for this to work well.

### 1. Install Nextflow

Full documentation here: https://www.nextflow.io/docs/latest/getstarted.html#installation

```
curl https://get.nextflow.io | bash
```

### 2. Install Docker

Please google for instructions.

## Pipeline Settings

Quick settings overview (full settings in nextflow.config file):

```
profiles {

    standard {
        process.executor = 'local'
    }


    gls {
        params.reads = 'gs://scgenomics/data/*/*_{R1,R2}_*.fastq.gz'
        params.outdir = 'gs://scgenomics/reports/'
        process.executor = 'google-lifesciences'
        workDir = 'gs://scgenomics/work'
        google.location = 'europe-west2'
        google.region  = 'europe-west2'
        google.project = 'nextflow-321415'
        google.lifeSciences.preemptible = true
        process.errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
        process.maxRetries = 3
    }

    az {
        params.reads = 'az://scgenomics/test_data/*/*_{R1,R2}_*.fastq.gz'
        params.outdir = 'az://scgenomics/test_reports'
        process.executor = 'azurebatch'
        workDir = 'az://scgenomics/work/'
        azure {
            storage {
                accountName = "<account name here>"
                accountKey = "<account key here>"
            }
            batch {
                location = 'uksouth'
                accountName = '<account name here>'
                accountKey = '<account key here>'
                autoPoolMode = true
                deletePoolsOnCompletion = true
            }
        }    

    }
```

## Main Architectural Concept

Develop and test pipeline with test data locally then orchestrate your cloud execution from your local machine. You will need to let your local machine run until all execution in the cloud is finished.

All bioinformatic tools in docker containers. Nextflow handles the start and stopping of containers automatically. You may need to run Docker pull manually on your local machine but all cloud container orchestration will be handled automatically by GCP or Azure.

1. Develop your pipeline using test data and execute locally
2. Upload your real fastq files to a cloud storage bucket
3. Edit nextflow.config to point to your cloud storage bucket for input and output locations
4. Run the pipeline with the correct flag (-gls or -az for Google or Azure)

## Sample Process

Conceptually simple. Name your process. Specify the container, the inputs, the outputs and the actual script to run the tool. You can additionally specify cpus and memory as well as a file output directory if you need access to the files.

```
process fastqc_run {
    publishDir "$params.outdir/fastqc/$sample_id/", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'
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
```

## Local Execution for Testing/Development
```
nextflow pipeline.nf -resume
```

## Cloud Execution Notes

Difficulty of set up (easiest to most difficult): Google Cloud > Azure > AWS.

### Why GCP > Azure > AWS?

Google Cloud simplifies the IAM management by allowing you to download a json key and setting up the service role only requires adding 4 policies.

Azure is slightly more complicated in terms of bucket/container set up. Note azure requires your credentials to be in nextflow.config. Beware! Do not commit the credentials to your git repo.

AWS is the most difficult because you need to create a custom Amazon Machine Image for the EC2 instances that can execute the awscli. This proved to be impossible as I ran into errors and troubleshooting was pretty much a nightmare and I did not succeed in the end.

### Google Cloud Execution


Basic tutorial here: https://cloud.google.com/life-sciences/docs/tutorials/nextflow

* To edit cloud configuration open nextflow.config and edit the desired settings under the profile section
* You will need to create a json key and download it

#### 1. Set the environment variables
```
export NXF_MODE=google
export GOOGLE_APPLICATION_CREDENTIALS=<path to key_filename.json>
```

Tip: Consider adding the above environment variables to the .bashrc file

#### 2. Running the pipeline
```
nextflow pipeline.nf -resume -profile gls
```

### Azure Cloud Execution

Slightly easier nextflow setup but setting up the azure services is slightly less user-friendly as Azure Portal is quite overwhelming compared to Google Cloud Console.

#### 1. Set azure credentials

Open up nextflow.config and fill in the storage and batch account names and keys.

#### 2. Running the pipeline
```
nextflow pipeline.nf -resume -profile az
```