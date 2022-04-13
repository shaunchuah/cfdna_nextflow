# Nextflow pipeline for cfDNA analysis

Welcome to this Nextflow pipeline for cfDNA analysis. After much trial and error, I have had the greatest overall success getting this pipeline up and running with Azure Batch.

For context, I develop this pipeline locally on my desktop machine (Core i5 6600K - 4 cores + 16GB RAM) running ubuntu 20.04 + Docker.
Then I have the same test data in Azure storage and I check that the pipeline orchestrates well in the cloud.
Finally I upload the real data and switch the directories to the real input and output directories and let it run.

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
        azure {
            storage {
                accountName = azure_config["storageAccountName"]
                accountKey = azure_config["storageAccountKey"]
            }
        }
    }
    az {
        params.reads = 'az://scgenomics/test_data/*/*_{R1,R2}_*.fastq.gz'
        params.outdir = 'az://scgenomics/test_reports'
        params.bowtie2_reference_index = "az://scgenomics/bt2_index.tar.gz"
        process.executor = 'azurebatch'
        workDir = 'az://scgenomics/work/'
        azure {
            storage {
                accountName = azure_config["storageAccountName"]
                accountKey = azure_config["storageAccountKey"]
            }
            batch {
                location = 'uksouth'
                accountName = azure_config["batchAccountName"]
                accountKey = azure_config["batchAccountKey"]
                autoPoolMode = true
                pools {
                    auto {
                        vmType = 'Standard_D4_v3'
                        vmCount = 4
                        maxVmCount = 12
                        autoScale = true
                    }
                }
            }
        }
    }
}
```



## Main Architectural Concept

Develop and test pipeline with test data locally then orchestrate your cloud execution from your local machine. You will need to let your local machine run until all execution in the cloud is finished. 

To take this to the next level, you may consider deploying a VM in the cloud as the head server for true portability but that is beyond the scope of this pipeline tutorial.

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

** Addendum 05/08/2021:

- Turns out it is easier to get a higher vCPU quota on Azure than GCP (minutes vs days)
- Note for reference databases in azure, you need to place **files at the root of the container**. Azure doesn't work well with subdirectories and a folder of files. I zipped them up in a .tar.gz and unzipped the reference folders in the process. Took me 2 days to figure out.



### Azure Cloud Execution

Slightly easier nextflow setup but setting up the azure services is slightly less user-friendly as Azure Portal is quite overwhelming compared to Google Cloud Console.

#### 1. Set azure settings

Open up nextflow.config and edit the path variables.

#### Pipeline Credentials

Azure credentials stored and read from credentials.json. Copy sample_credentials.json to a file named 'credentials.json' in the root of this repo and fill in your Azure login details. This helps split sensitive variables out of the git repository/nextflow config file.

#### 2. Running the pipeline
```
nextflow pipeline.nf -resume -profile az
```

## Pipeline Updates

13 April 2022: Moved old version into `archive_2021`. New pipeline streamlined. Input reads are split into QC and kraken2 steps.
