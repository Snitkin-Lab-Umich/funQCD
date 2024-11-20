# funQCD - Quality Control and Contamination Detection workflow.

funQCD is a quality control pipeline for short read fungal data, with a focus on Candida auris sequencing data. It is a modified version of QCD, a workflow for microbial Illumina sequencing quality control and contamination detection implemented in Snakemake (https://github.com/Snitkin-Lab-Umich/QCD).

### Summary

This pipeline performs the following steps:

* [Fastqc](https://github.com/s-andrews/FastQC) is used to generate HTML reports to asses quality of sequencing reads before and after trimming reads.
* [QUAST](https://quast.sourceforge.net/) is used for quality assesment of the assembly.
* [Trimmomatic](https://github.com/usadellab/Trimmomatic) is used to trim and filter low-quality bases and adapter sequences from raw FASTQ reads.
* [fastq-scan](https://github.com/rpetit3/fastq-scan) is used to estimate genome coverage of FASTQ files.
* [SPAdes](https://github.com/ablab/spades) is used to assemble trimmed reads into contigs.
* [Kraken2](https://github.com/DerrickWood/kraken2) is used to provide detailed reports on the taxonomic composition of the trimmed raw reads.
* [Funannotate](https://github.com/nextgenusfs/funannotate) is used for structural and functional annotation of the assembly from SPAdes.
* [InterProScan](https://github.com/ebi-pf-team/interproscan) is used for functional annotation of the genome, with its results used by Funannotate.
* [EggNOG-Mapper](https://github.com/eggnogdb/eggnog-mapper) is used for functional annotation of the genome, with its results used by Funannotate.
* [AuriClass](https://github.com/RIVM-bioinformatics/auriclass) is used to determine the clade of Candida auris.
* [BUSCO](https://busco.ezlab.org/) is used to evaluate the completeness of the assembly.
* [Multiqc](https://github.com/MultiQC/MultiQC) aggregates multiple QC results across all samples to produce a single HTML report.


Steps from QCD not currently used:

* The assembled contigs from [SPAdes](https://github.com/ablab/spades) are passed through [Prokka](https://github.com/tseemann/prokka) for annotation, [QUAST](https://quast.sourceforge.net/) for assembly statistics, [MLST](https://github.com/tseemann/mlst) for determining sequence type based on sequences of housekeeping genes, [AMRFinderPlus](https://github.com/ncbi/amr) to identify antimicrobial resistance genes, [skani](https://github.com/bluenote-1577/skani) to identify closest reference genome, and [BUSCO](https://busco.ezlab.org/) for assembly completeness statistics.

The workflow generates all the output in the output prefix folder set in the config file (instructions on setup found [below](#config)). Each workflow step gets its own individual folder as shown:

```
results
└─—— run_name
   ├── auriclass
   ├── busco
   ├── downsample
   ├── multiqc
   ├── quality_aftertrim
   ├── quality_raw
   ├── quast
   ├── raw_coverage
   ├── spades
   ├── trimmomatic
   └── funannotate
		├──interproscan
		└──eggnog

```


## Installation 


> If you are using Great Lakes HPC, ensure you are cloning the repository in your scratch directory. Change `your_uniqname` to your uniqname. 

```

cd /scratch/esnitkin_root/esnitkin1/your_uniqname/

```

> Clone the github directory onto your system. 

```

git clone https://github.com/Snitkin-Lab-Umich/funQCD

```

> Ensure you have successfully cloned QCD. Type `ls` and you should see the newly created directory **_funQCD_**. Move to the newly created directory.

```

cd funQCD

```

> Load Bioinformatics, snakemake and singularity modules from Great Lakes modules.

```

module load Bioinformatics snakemake singularity

```
<!--
```

module load snakemake singularity

```
-->

This workflow makes use of singularity containers available through [State Public Health Bioinformatics group](https://github.com/StaPH-B/docker-builds). If you are working on Great Lakes (umich's HPC cluster), you can load snakemake and singularity modules as shown above. However, if you are running it on your local or other computing platform, ensure you have snakemake and singularity installed.


## Setting up configuration files

This pipeline requires several files that specify where your raw reads are located, where external databases can be found, how much memory is available, and several other parameters. Once funQCD is installed, you'll need to update the following files as well. If you are just testing the pipeline, these files already contain example paths you can use.

### 1) config/config.yaml
This file contains information about where the sequencing reads and databases can be found. Detailed instructions for each section are provided in the file itself, with a short description of the most important ones here.

* short_reads: This is the most important section to edit. It must contain the full path to a directory containing your raw sequencing reads, as .fastq.gz files with standardized names (see below).
* samples: This is the path to a .csv file containing the names of the raw reads you want to process (see next section).
* prefix: This is the name of the directory containing all outputs, which will appear under `results/` after you run funQCD. This should be a unique name for each run of funQCD.
* adapter_file: This is the path to the adapter file used for Trimmomatic.
* skani_db: This is the path to the database used while running skani.
* funqcd_lib: This is the path to a directory containing several databases required for InterProScan, EggNOG-mapper, and funannotate, as well as RNA-seq data for funannotate.

Currently, funQCD expects your sequencing reads to be named in a specific format:
[SampleName]_R[number].fastq.gz

If your reads are named differently, then the script will fail to recognize them. For example, `TO315_R1.fastq.gz` would work, but `TO315_R1_001.fastq.gz` would not. funQCD is currently being updated to fix this issue. 

### 2) config/samples.csv
This file contains the names of your raw read files. It should be a comma-seperated file consisting of two columns: `sample_id` and `illumina_r1`. These name in these columns should exactly match the names used in the directory specified by `short_reads`.

* `sample_id` is the prefix that should be extracted from your FASTQ reads. For example, if you have files called `Rush_KPC_110_R1.fastq.gz` and `Rush_KPC_110_R2.fastq.gz` in your raw FASTQ directory, your sample_id would be `Rush_KPC_110`.

* `illumina_r1` is the full name of the corresponding raw FASTQ file containing the forward reads. Using the same example as above, your sample_id would be `Rush_KPC_110_R1.fastq.gz`. **_Only include the forward reads here._**

You can create sample.csv file using the following for loop. Replace *path_to_your_raw_reads* below with the actual path to your raw sequencing reads.

```

echo "sample_id,illumina_r1" > config/sample.tsv

for read1 in path_to_your_raw_reads/*_R1.fastq.gz; do
    sample_id=$(basename $read1 | sed 's/_R1.fastq.gz//g')
    read1_basename=$(basename $read1)
    echo $sample_id,$read1_basename
done >> config/sample.tsv

```

### 3) profile/config.yaml
This file contains the options for running snakemake. The most important section to change is `slurm_account`, which needs to match an account that can submit jobs to the Great Lakes HPC cluster. If you changed the location of the funqcd_lib above, you'll need to change the directories bound via singularity in this file as well. You can also set limits on the resources funQCD will attempt to use here. For example, adding the line `resources: mem_mb=40000` will limit all jobs to ~40G of memory. 

## Running funQCD

> First, perform a dry run of funQCD by running the following command. This will show the steps and commands that funQCD will execute in the real run, without actually executing anything.

```

snakemake -s workflow/fQCD.smk -p --configfile config/config.yaml --profile ./profile/ -n

```

> The snakemake options present in profile/config.yaml should be visible in the dry run (such as memory and runtime for each rule). By default, --slurm is enabled in these options, and snakemake will submit jobs to the cluster using the account in your profile. If everything looks correct, start the run using a job script with minimal CPUs, minimal memory, and a long runtime.  

```
#!/bin/bash

#SBATCH --mail-user=[your_email]
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --partition=standard
#SBATCH --account=[your_account]
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500m
#SBATCH --time=20:00:00

module load Bioinformatics snakemake singularity
snakemake -s workflow/fQCD.smk -p --configfile config/config.yaml --profile ./profile/

```

> If you want to run the snakemake pipeline without cluster execution mode, change the `slurm: True` line to `slurm: False` in your profile and run the same command.

```
module load Bioinformatics snakemake singularity
snakemake -s workflow/fQCD.smk -p --configfile config/config.yaml --profile ./profile/
 
```

## Dependencies

### Near Essential
* [Snakemake>=7.32.4](https://snakemake.readthedocs.io/en/stable/#)
* [Conda](https://docs.conda.io/en/latest/)

<!--All the necessary software stack required for the workflow will be installed using conda package manager.-->

### Tool stack used in workflow

* [fastq-scan](https://github.com/rpetit3/fastq-scan)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [SPades](https://github.com/ablab/spades)
* [AMRFinderPlus](https://github.com/ncbi/amr)
* [bioawk](https://github.com/lh3/bioawk)
* [Prokka](https://github.com/tseemann/prokka)
* [mlst](https://github.com/tseemann/mlst)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
