# funQCD - Quality Control and Contamination Detection workflow.

funQCD is a quality control pipeline for short read fungal data, with a focus on Candida auris sequencing data. It is a modified version of QCD (https://github.com/Snitkin-Lab-Umich/QCD), and is intended to be used on a HPC cluster. The workflow is split into three sections (assembly, prediction, and annotation). Basic QC is checked after the assembly and prediction steps, with failed samples removed from the pipeline.

### Summary

This pipeline performs the following steps:

* [Fastqc](https://github.com/s-andrews/FastQC) is used to generate HTML reports to asses quality of sequencing reads before and after trimming reads.
* [QUAST](https://quast.sourceforge.net/) is used for quality assesment of the assembly.
* [Trimmomatic](https://github.com/usadellab/Trimmomatic) is used to trim and filter low-quality bases and adapter sequences from raw FASTQ reads.
* [fastq-scan](https://github.com/rpetit3/fastq-scan) is used to estimate genome coverage of FASTQ files.
* [SPAdes](https://github.com/ablab/spades) is used to assemble trimmed reads into contigs.
* [Kraken2](https://github.com/DerrickWood/kraken2) is used to provide detailed reports on the taxonomic composition of the trimmed raw reads.
* [Funannotate](https://github.com/nextgenusfs/funannotate) is used for structural and functional annotation of the assembly from SPAdes.
* [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) is used to soft-mask the genome.
* [GeneMark](https://exon.gatech.edu/) is used for to for structural annotation of the genome, as part of Funannotate.
* [InterProScan](https://github.com/ebi-pf-team/interproscan) is used for functional annotation of the genome, with its results used by Funannotate.
* [EggNOG-Mapper](https://github.com/eggnogdb/eggnog-mapper) is used for functional annotation of the genome, with its results used by Funannotate.
* [AuriClass](https://github.com/RIVM-bioinformatics/auriclass) is used to determine the clade of Candida auris.
* [BUSCO](https://busco.ezlab.org/) is used to evaluate the completeness of the assembly.
* [Multiqc](https://github.com/MultiQC/MultiQC) aggregates multiple QC results across all samples to produce a single HTML report.


Steps from QCD not currently used:

* The assembled contigs from [SPAdes](https://github.com/ablab/spades) are passed through [Prokka](https://github.com/tseemann/prokka) for annotation, [QUAST](https://quast.sourceforge.net/) for assembly statistics, [MLST](https://github.com/tseemann/mlst) for determining sequence type based on sequences of housekeeping genes, [AMRFinderPlus](https://github.com/ncbi/amr) to identify antimicrobial resistance genes, [skani](https://github.com/bluenote-1577/skani) to identify closest reference genome, and [BUSCO](https://busco.ezlab.org/) for assembly completeness statistics.


## Installation 


> Clone the github directory. 

```

git clone https://github.com/Snitkin-Lab-Umich/funQCD
cd funQCD

```

> Load Bioinformatics, snakemake, and singularity modules.

```

module load Bioinformatics snakemake singularity

```


## Setup

This pipeline requires several files that specify where your raw reads are located, where external databases can be found, how much memory is available, and several other parameters. Once funQCD is installed, you'll need to update the following files as well. Usually, only a few edits are needed. If you are just testing the pipeline, these files already contain example paths you can use.

### 1) Run setup.py

The pipeline requires a list of sample names to process, in the form of a CSV with `sample_id` as the only column. You can make this yourself, or run the setup script to do so automatically. Just provide the script with the path to the directory containing your Illumina reads:

```

python setup.py [path_to_reads]

```
The output file will be `config/samples.csv`. This will only work if all of your read names are unique, and end with `_R1.fastq.gz` or `_R1.fastq.gz`. For example, `TO315_R1.fastq.gz` would work, but `TO315_R1_001.fastq.gz` would not. It is recommended to check samples.csv after running the script to make sure it worked as expected.

This script will also attempt to copy a GeneMark license (.gm_key) from your home directory into the funQCD directory. This file is required for the pipeline to work, and can be downloaded [here](https://genemark.bme.gatech.edu/license_download.cgi). Place the key either in your home directory (~/) or the funQCD working directory.

### 2) config/config.yaml

This file contains information about where the sequencing reads and databases can be found. Replace the `short_reads` section with the same file path you used for the setup script in step 1.

(Optional) You can also use this file to change the databases and references that funQCD uses. Detailed instructions for each section are provided in the file itself, with a short description of the most important ones here.

* short_reads: This must contain the full path to a directory containing your raw sequencing reads, as .fastq.gz files with standardized names.
* samples: This is the path to a .csv file containing the names of the raw reads you want to process (see next section).
* prefix: This is the name of the directory containing all outputs, which will appear under `results/` after you run funQCD. This should be a unique name for each run of funQCD.
* adapter_file: This is the path to the adapter file used for Trimmomatic.
* skani_db: This is the path to the database used while running skani.
* funqcd_lib: This is the path to a directory containing several databases required for InterProScan, EggNOG-mapper, and funannotate, as well as RNA-seq data for funannotate. This is also where GeneMark should be installed.


### 3) profile/config.yaml

This file contains the options for running snakemake. Replace the `slurm_account` section with your own account. 

(Optional) If you changed the databases and references used above, you'll need to change the directories bound via singularity in this file as well. You can also set limits on the resources funQCD will attempt to use here. For example, adding the line `resources: mem_mb=40000` will limit all jobs to ~40G of memory.


## Running funQCD

> First, perform a dry run of funQCD's assembly step by running the following command. This will show the steps and commands that funQCD will execute in the real run, without actually executing anything. Remove the --quiet flag for a more detailed view.

```

snakemake -s workflow/funQCD_ASSEMBLY.smk -p --configfile config/config.yaml --profile ./profile/ -n --quiet

```

> The snakemake options present in profile/config.yaml should be visible in the detailed dry run (such as memory and runtime for each rule). By default, --slurm is enabled in these options, and snakemake will submit jobs to the cluster using the account in your profile. If everything looks correct, start the run using a job script with minimal CPUs, moderate memory, and a long runtime. An example job script is provided in `run_fqcd.job`.

```
#!/bin/bash

#SBATCH --mail-user=[your_email]
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --partition=standard
#SBATCH --account=[your_account]
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=20:00:00

module load Bioinformatics snakemake singularity
snakemake -s workflow/funQCD_ASSEMBLY.smk -p --configfile config/config.yaml --profile ./profile/
snakemake -s workflow/funQCD_PREDICTION.smk -p --configfile config/config.yaml --profile ./profile/
snakemake -s workflow/funQCD_ANNOTATION.smk -p --configfile config/config.yaml --profile ./profile/

```

> It's recommended to check the output of each intermediate step (assembly and prediction) before running the next one. funQCD is currently being update to automate this process. 
