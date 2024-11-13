# Author: Ali Pirani and Dhatri Badri
configfile: "config/config.yaml"

#include: "fQCD_report.smk"

import pandas as pd
import os

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

PREFIX = config["prefix"]

SHORTREADS = list(samples_df['sample_id'])

if not os.path.exists("results/"):
    os.system("mkdir %s" % "results/")

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)

# This is only needed for a local installation of the funannotate databases
# if not os.path.exists("lib/"):
#     os.system("mkdir %s" % "lib/")
#
# This directory can be empty, but it cannot be missing (otherwise the singularity bind will fail)
# For IPS, this needs to end up in a specific directory as well ('/opt/interproscan/data/')
# if not os.path.exists("lib/interproscan_data/"):
#     os.system("mkdir %s" % "lib/interproscan_data/")

# if not os.path.exists("lib/interproscan_data/data/"):
#     os.system("mkdir %s" % "lib/interproscan_data/data/")

# A similar directory is needed for BUSCO
# This needs to end up in '/opt/databases/saccharomycetales' for funannotate to see it
# if not os.path.exists("lib/busco/"):
#     os.system("mkdir %s" % "lib/busco/")

# if not os.path.exists("lib/busco/saccharomycetales/"):
#     os.system("mkdir %s" % "lib/busco/saccharomycetales/")


def downsample_reads(file, file2, out1, out2, genome_size):
    file = file.pop()
    file2 = file2.pop()
    out1 = out1.pop()
    out2 = out2.pop()

    # Extract basic fastq reads stats with seqtk

    gsize = genome_size.pop()

    print("Using Genome Size: %s to calculate coverage" % gsize)
    
    seqtk_check = "/nfs/esnitkin/bin_group/seqtk/seqtk fqchk -q3 %s > %s_fastqchk.txt" % (file, file)

    print(seqtk_check)

    try:
        os.system(seqtk_check)
    except sp.CalledProcessError:
        print('Error running seqtk for extracting fastq statistics.')
        sys.exit(1)

    with open("%s_fastqchk.txt" % file, 'rU') as file_open:
        for line in file_open:
            if line.startswith('min_len'):
                line_split = line.split(';')
                min_len = line_split[0].split(': ')[1]
                max_len = line_split[1].split(': ')[1]
                avg_len = line_split[2].split(': ')[1]
            if line.startswith('ALL'):
                line_split = line.split('\t')
                total_bases = int(line_split[1]) * 2
    file_open.close()

    print('Average Read Length: %s' % avg_len)

    print('Total number of bases in fastq: %s' % total_bases)

    # Calculate original depth and check if it needs to be downsampled to a default coverage.
    ori_coverage_depth = int(total_bases / gsize)

    print('Original Covarage Depth: %s x' % ori_coverage_depth)

    # proc = sp.Popen(["nproc"], stdout=sp.PIPE, shell=True)
    # (nproc, err) = proc.communicate()
    # nproc = nproc.strip()

    if ori_coverage_depth >= 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        # r1_sub = "/tmp/%s" % os.path.basename(file)
        r1_sub = out1

        # Downsample using seqtk
        try:
            #print("Generating seqtk Downsampling command")
            print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file, factor, r1_sub))

            seqtk_downsample = "/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (
                file, factor, r1_sub)
            os.system(seqtk_downsample)
            #call(seqtk_downsample, logger)
        except sp.CalledProcessError:
            print('Error running seqtk for downsampling raw fastq reads.')
            sys.exit(1)

        if file2:
            # r2_sub = "/tmp/%s" % os.path.basename(file2)
            r2_sub = out2
            try:
                print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file2, factor, r2_sub))
                os.system("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file2, factor, r2_sub))
                #call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (file2, factor, nproc, os.path.basename(file2)), logger)
            except sp.CalledProcessError:
                print('Error running seqtk for downsampling raw fastq reads.')
                sys.exit(1)
        else:
            r2_sub = "None"

    elif ori_coverage_depth < 100:
        r1_sub = file
        r2_sub = file2
        os.system("cp %s %s" % (file, out1))
        os.system("cp %s %s" % (file2, out2))

rule all:
    input:
        coverage = expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", sample=SAMPLE, prefix=PREFIX),
        # fastqc_raw = expand("results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_001_fastqc.html", sample=SAMPLE, prefix=PREFIX),
        fastqc_raw = expand("results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html", sample=SAMPLE, prefix=PREFIX),
        trim = expand("results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE, prefix=PREFIX),
        fastqc_aftertrim = expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html", sample=SAMPLE, prefix=PREFIX),
        downsample_read = expand("results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE, prefix=PREFIX),
        spades_assembly = expand("results/{prefix}/spades/{sample}/contigs.fasta", sample=SAMPLE, prefix=PREFIX),
        spades_l1000_assembly = expand("results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta", sample=SAMPLE, prefix=PREFIX),
        ##prokka_gff = expand("results/{prefix}/prokka/{sample}/{sample}.gff", sample=SAMPLE, prefix=PREFIX),
        quast_report = expand("results/{prefix}/quast/{sample}/report.tsv", sample=SAMPLE, prefix=PREFIX),
        auriclass_report = expand("results/{prefix}/auriclass/{sample}/{sample}_report.tsv", sample=SAMPLE, prefix=PREFIX),
        ##mlst_report = expand("results/{prefix}/mlst/{sample}/report.tsv", sample=SAMPLE, prefix=PREFIX),
        #skani_ref_genome_results = expand("results/{prefix}/skani/{sample}/{sample}_skani_output.txt", sample=SAMPLE, prefix=PREFIX),
        #coverage_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Final_Coverage.txt", prefix=PREFIX),
        #skani_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Skani_report_final.csv", prefix=PREFIX),
        multiqc_report = expand("results/{prefix}/multiqc/{prefix}_QC_report.html", prefix=PREFIX, sample = SAMPLE),
        ##mlst_final_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_MLST_results.csv", prefix=PREFIX),
        #QC_summary = expand("results/{prefix}/{prefix}_Report/data/{prefix}_QC_summary.csv", prefix=PREFIX),
        #QC_plot = expand("results/{prefix}/{prefix}_Report/fig/{prefix}_Coverage_distribution.png", prefix=PREFIX)
        #funannotate = expand("results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),
        #funannotate_setup = "lib/busco/saccharomycetales/dataset.cfg",
        funannotate_mask = expand("results/{prefix}/funannotate/{sample}/{sample}_masked.fa", sample=SAMPLE, prefix=PREFIX),
        funannotate_train = expand("results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam", sample=SAMPLE, prefix=PREFIX),
        funannotate_predict = expand("results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),
        funannotate_update = expand("results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),        
        interproscan_data_dl = config["funqcd_lib"] + "interproscan_data/data/antifam/",
        #interproscan_data_init = config["funqcd_lib"] + "interproscan_data/test_log.txt",
        interproscan = expand("results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml", sample=SAMPLE, prefix=PREFIX),
        eggnog_data_dl = config["funqcd_lib"] + "eggnog_data/eggnog.db",
        eggnog = expand("results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations", sample=SAMPLE, prefix=PREFIX),
        funannotate_annotate = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),
        busco_final = expand("results/{prefix}/busco/busco_output/batch_summary.txt", prefix = PREFIX)

rule coverage:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
        #r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1_001.fastq.gz")),
        #r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2_001.fastq.gz")),
    output:
        coverage = f"results/{{prefix}}/raw_coverage/{{sample}}/{{sample}}_coverage.json",
    params:
        size=config["genome_size"]
    #conda:
    #    "envs/fastq-scan.yaml"
    singularity:
        "docker://staphb/fastq-scan:1.0.1"
    shell:
        "zcat {input.r1} {input.r2} | fastq-scan -g {params.size} > {output.coverage}"

rule quality_raw:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
        #r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1_001.fastq.gz")),
        #r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2_001.fastq.gz")),
    output:
        raw_fastqc_report_fwd = f"results/{{prefix}}/quality_raw/{{sample}}/{{sample}}_Forward/{{sample}}_R1_fastqc.html",
        raw_fastqc_report_rev = f"results/{{prefix}}/quality_raw/{{sample}}/{{sample}}_Reverse/{{sample}}_R2_fastqc.html",
    log:
        "logs/{prefix}/quality_raw/{sample}/{sample}.log"
    params:
        outdir="results/{prefix}/quality_raw/{sample}/{sample}"
    #conda:
    #    "envs/fastqc.yaml"
    singularity:
        "docker://staphb/fastqc:0.12.1"
    #envmodules:
    #    "Bioinformatics",
    #    "fastqc"
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        """

rule trimmomatic_pe:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz"))
        #r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1_001.fastq.gz")),
        #r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2_001.fastq.gz")),  
    output:
        r1 = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R1_trim_paired.fastq.gz",
        r2 = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R2_trim_paired.fastq.gz", 
        # reads where trimming entirely removed the mate
        r1_unpaired = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R2_trim_unpaired.fastq.gz",
    params:
        adapter_filepath=config["adapter_file"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],
        #threads = config["ncores"],
    log:
        "logs/{prefix}/trimmomatic/{sample}/{sample}.log"
    threads: 8
    # threads: workflow.cores
    # This has been changed to specifically use 8 cores (the previous default number of cores used)
    # Before, this always allocated the maximum number of cores
    #conda:
    #    "envs/trimmomatic.yaml"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    #envmodules:
    #    "Bioinformatics",
    #    "trimmomatic"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {threads} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}"

rule quality_aftertrim:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        aftertrim_fastqc_report_fwd = f"results/{{prefix}}/quality_aftertrim/{{sample}}/{{sample}}_Forward/{{sample}}_R1_trim_paired_fastqc.html",
        aftertrim_fastqc_report_rev = f"results/{{prefix}}/quality_aftertrim/{{sample}}/{{sample}}_Reverse/{{sample}}_R2_trim_paired_fastqc.html",
    log:
        "logs/{prefix}/{sample}/quality_aftertrim/{sample}.log"
    params:
        outdir="results/{prefix}/quality_aftertrim/{sample}/{sample}"
    #conda:
    #    "envs/fastqc.yaml"
    singularity:
        "docker://staphb/fastqc:0.12.1"
    #envmodules:
    #    "Bioinformatics",
    #    "fastqc"
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        """

rule downsample:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        outr1 = f"results/{{prefix}}/downsample/{{sample}}/{{sample}}_R1_trim_paired.fastq.gz",
        outr2 = f"results/{{prefix}}/downsample/{{sample}}/{{sample}}_R2_trim_paired.fastq.gz",
    params:
        gsize = config["genome_size"],
    run:
        downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize})

# rule kraken:
#     input:
#         r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
#         r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
#     output:
#         kraken_out = f"results/{{prefix}}/kraken/{{sample}}/{{sample}}_kraken_out",
#         kraken_report = f"results/{{prefix}}/kraken/{{sample}}/{{sample}}_kraken_report.tsv",
#     params:
#         db = config["kraken_db"],
#         threads = 12
#         # threads = config["threads"]
#     #conda:
#     #    "envs/kraken.yaml"
#     singularity:
#         "docker://staphb/kraken2:2.1.3"
#     shell:
#         "kraken2 --db {params.db} --threads {params.threads} --paired --gzip-compressed --quick --output {output.kraken_out} --report {output.kraken_report} {input.r1} {input.r2}"

rule assembly:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        spades_assembly = f"results/{{prefix}}/spades/{{sample}}/contigs.fasta",
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        #db = config["kraken_db"],
    #conda:
    #    "envs/spades.yaml"
    singularity:
        "docker://staphb/spades:4.0.0"
    #envmodules:
    #    "Bioinformatics",
    #    "spades/4.0.0"
    shell:
        "spades.py --isolate --pe1-1 {input.r1} --pe1-2 {input.r2} -o {params.out_dir}"

#rule bioawk:
#    input:
#        spades_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/contigs.fasta"),
#    output:
#        spades_l1000_assembly = f"results/{{prefix}}/{{sample}}/spades/{{sample}}_contigs_l1000.fasta",
#    params:
#        out_dir = "results/{prefix}/{sample}/spades/",
#        bioawk_params = config["bioawk"],
#        prefix = "{sample}",
#    conda:
#        "envs/bioawk.yaml"
#    shell:
#        """
#        {params.bioawk_params} {input.spades_assembly} > {output.spades_l1000_assembly} && grep '>' {output.spades_l1000_assembly} > {params.out_dir}/spades_assembly_header_info.txt && sed -i 's/>NODE_/>{params.prefix}_/g' {output.spades_l1000_assembly} && sed -i 's/_length_.*_cov_.*//g' {output.spades_l1000_assembly}
#        """

rule bioawk:
    input:
        spades_assembly = lambda wildcards: f"results/{wildcards.prefix}/spades/{wildcards.sample}/contigs.fasta"
    output:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta"
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        bioawk_params = config["bioawk"],
        prefix = "{sample}"
    #conda:
    #    "envs/bioawk.yaml"
    singularity:
        "docker://lbmc/bioawk:1.0"
    shell:
        """
        ./bioawk.sh {input.spades_assembly} {output.spades_l1000_assembly} {params.out_dir} {params.prefix}
        """

# rule funannotate:
#     input:
#         spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/{wildcards.sample}_contigs_l1000.fasta"),
#     output:
    #     prokka_gff = f"results/{{prefix}}/prokka/{{sample}}/{{sample}}.gff",
    # params: 
    #     prokka_params = config["prokka"],
    #     outdir = "results/{prefix}/prokka/{sample}",
    #     prefix = "{sample}",
    #conda:
    #    "envs/prokka.yaml"
#     singularity:
#         "docker://nextgenusfs/funannotate"
#     #envmodules:
#     #    "Bioinformatics",
#     #    "funannotate"
#     shell:
#         "prokka -outdir {params.outdir} --strain {params.prefix} --prefix {params.prefix} {params.prokka_params} {input.spades_l1000_assembly}"
#     #cd /nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2023_09_21_Funannotate_References
# #funannotate predict -i /nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2023_09_21_Funannotate_References/B11205_mask.fna -o /nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2023_09_21_Funannotate_References/B11205_funannotate -s "Candida auris" --augustus_species "candida_albicans" --AUGUSTUS_CONFIG_PATH=/nfs/turbo/umms-esnitkin/conda/funannotate/config/ --cpus 8 

rule quast:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        quast_report = f"results/{{prefix}}/quast/{{sample}}/report.tsv",
    params: 
        outdir = "results/{prefix}/quast/{sample}/",
        prefix = "{sample}",
    #conda:
    #    "envs/quast.yaml"
    singularity:
        "docker://staphb/quast:5.0.2"
    #envmodules:
    #    "Bioinformatics",
    #    "quast"
    shell:
        """
        ./quast.sh {input.spades_l1000_assembly} {params.outdir} 
        """
        #"quast.py {input.spades_l1000_assembly} -o {params.outdir} --contig-thresholds 0,1000,5000,10000,25000,50000"

rule auriclass:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        auriclass_report = f"results/{{prefix}}/auriclass/{{sample}}/{{sample}}_report.tsv",
    params: 
        #outdir = "results/{prefix}/mlst/{sample}/",
        sample = "{sample}",
    singularity:
        "docker://quay.io/biocontainers/auriclass:0.5.4--pyhdfd78af_0"
    shell:
        "auriclass --name {params.sample} -o {output.auriclass_report} {input.spades_l1000_assembly}"

#rule amrfinder:
#    input:
#        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/{wildcards.sample}_contigs_l1000.fasta"),
#    output:
#        amrfinder = f"results/{{prefix}}/{{sample}}/amrfinder/{{sample}}_amrfinder.tsv",
#    params: 
#        outdir = "results/{prefix}/{sample}/amrfinder",
#        prefix = "{sample}",
#        organism = config['amrfinder_organism']
    #conda:
    #    "envs/amrfinder.yaml"
#    singularity:
#        "docker://staphb/ncbi-amrfinderplus:latest"
#    shell:
#        "amrfinder --plus --output {output.amrfinder} --debug --log {params.outdir}/{params.prefix}.log --nucleotide_output {params.outdir}/{params.prefix}_reported_nucl.fna -n {input.spades_l1000_assembly} -O {params.organism}"

rule busco:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        busco_out = f"results/{{prefix}}/{{sample}}/busco/busco.txt",
    params: 
        outdir = "results/{prefix}/busco/{sample}/",
        prefix = "{sample}",
        # threads = config["ncores"],
    #conda:
    #    "envs/busco.yaml"
    singularity:
        #"docker://staphb/busco:5.7.1-prok-bacteria_odb10_2024-01-08"
        "docker://ezlabgva/busco:v5.7.1_cv1"
    threads: 8
    # threads: workflow.cores
    #envmodules:
    #    "Bioinformatics",
    #    "busco"
    shell:
        "busco -f -i {input.spades_l1000_assembly} -m genome -l bacteria_odb10 -o {params.outdir}"

rule skani:
    input:
        spades_contigs_file = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/contigs.fasta")
    output:
        skani_output = f"results/{{prefix}}/skani/{{sample}}/{{sample}}_skani_output.txt"
    params:
        skani_ani_db = config["skani_db"],
        #threads = 4
    threads: 4
    #conda:
    #    "envs/skani.yaml"
    singularity:
        "docker://staphb/skani:0.2.1"
    shell:
        "skani search {input.spades_contigs_file} -d {params.skani_ani_db} -o {output.skani_output} -t {threads}"
        
rule multiqc:
    input:
        quast_report = expand("results/{prefix}/quast/{sample}/report.tsv", sample = SAMPLE, prefix = PREFIX),
        #prokka_gff = "results/{prefix}/{sample}/prokka/{sample}.gff",
        #spades_assembly = "results/{prefix}/spades/{sample}/contigs.fasta",
        #kraken_report = "results/{prefix}/{sample}/kraken/{sample}_kraken_report.tsv",
        #coverage = "results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json",
        aftertrim_fastqc_report_fwd = expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html", sample = SAMPLE, prefix = PREFIX),
        aftertrim_fastqc_report_rev = expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Reverse/{sample}_R2_trim_paired_fastqc.html", sample = SAMPLE, prefix = PREFIX),
        raw_fastqc_report_fwd = expand("results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html", sample = SAMPLE, prefix = PREFIX),
        raw_fastqc_report_rev = expand("results/{prefix}/quality_raw/{sample}/{sample}_Reverse/{sample}_R2_fastqc.html", sample = SAMPLE, prefix = PREFIX)
    output:
        multiqc_report = "results/{prefix}/multiqc/{prefix}_QC_report.html",
    params:
        resultsoutdir = "results/{prefix}",
        outdir = "results/{prefix}/multiqc",
        prefix = "{prefix}",
        quast_dir = "results/{prefix}/quast/",
        raw_fastqc_dir = "results/{prefix}/quality_raw/",
        aftertrim_fastqc_dir = "results/{prefix}/quality_aftertrim/",
    singularity:
        "docker://multiqc/multiqc:v1.25.1"
    shell:
        """
        multiqc -f --outdir {params.outdir} -n {params.prefix}_QC_report -i {params.prefix}_QC_report \
        {params.quast_dir} {params.raw_fastqc_dir} {params.aftertrim_fastqc_dir}
        """

# This adds another BUSCO database for funannotate to use
# This is bound via singularity to add it to the existing databases
# rule funannotate_setup:
#     output:
#         busco_db = "lib/busco/saccharomycetales/dataset.cfg"
#     params:
#         out_dir = "lib/busco/",
#     singularity:
#         "docker://nextgenusfs/funannotate:v1.8.17"
#     shell:
#         """
#         funannotate setup --busco_db saccharomycetales --install busco --database {params.out_dir} 
#         """

# This needs to be split into individual steps for clean/sort/mask, predict, and annotate
# I've split out the predict and annotate steps
# The new order is setup, mask, train, predict, update, interproscan, eggnong, annotate
rule funannotate_mask:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        masked_assembly = "results/{prefix}/funannotate/{sample}/{sample}_masked.fa"
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}"
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate clean -i {input.spades_l1000_assembly} -o {params.out_dir}{params.sample}_cleaned.fa
        funannotate sort -i {params.out_dir}{params.sample}_cleaned.fa -o {params.out_dir}{params.sample}_sorted.fa --minlen 0
        funannotate mask -i {params.out_dir}{params.sample}_sorted.fa -o {params.out_dir}{params.sample}_masked.fa
        """
# removed: funannotate annotate -i {output.funannotate_predict_out} -o {params.out_dir} --cpus {threads}
# removed: funannotate predict -i {params.out_dir}{params.sample}_masked.fa -o {params.out_dir} --species '{params.sample}' --augustus_species candida_albicans --cpus {threads}

# in progress
# this requires the RNA-seq data from Teresa, and assumes that the files are in a specific format
# specifically, that read 1 contains '_R1_' in the file name, that all files are in fastq.gz format and in RF order for stranded RNA-seq, and that they can be ordered alphabetically
rule funannotate_train:
    input:
        masked_assembly = "results/{prefix}/funannotate/{sample}/{sample}_masked.fa"
    output:
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
        funannotate_training_pasa_gff = "results/{prefix}/funannotate/{sample}/training/funannotate_train.pasa.gff3",
        funannotate_training_stringtie = "results/{prefix}/funannotate/{sample}/training/funannotate_train.stringtie.gtf",
        funannotate_training_transc_align = "results/{prefix}/funannotate/{sample}/training/funannotate_train.transcripts.gff3",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}",
        rna_data_r1 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R1_' in x and 'fastq.gz' in x])),
        rna_data_r2 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R2_' in x and 'fastq.gz' in x])),
        mem_g = str(int(config["mem_mb"]/1000)) + "G",
        test = vars(workflow.resources)
    # threads: workflow.cores
    threads: 8
    resources:
        mem_mb = config["mem_mb"]
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate train --input {input.masked_assembly} --out {params.out_dir} \
        --left {params.rna_data_r1} --right {params.rna_data_r2} --stranded RF \
        --jaccard_clip --species "Candida auris" --isolate {params.sample} --cpus {threads} --memory {params.mem_g}
        """

# in progress
# This should automatically detect the four training files generated previously, even without explicit input
# All steps should run with 'pasa' or 'rna-bam' under Training-Methods. Nothing should run with 'busco'.
rule funannotate_predict:
    input:
        masked_assembly = "results/{prefix}/funannotate/{sample}/{sample}_masked.fa",
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
        busco_db = config["funqcd_lib"] + "busco/saccharomycetales/dataset.cfg"
    output:
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa",
        funannotate_predict_out = directory("results/{prefix}/funannotate/{sample}/predict_results/")
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}"
    threads: 8
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        # removed --species "Candida auris" --isolate {params.sample}
        # consider changing output names and re-adding
        """
        funannotate predict --input {input.masked_assembly} --out {params.out_dir} \
        --species {params.sample} \
        --busco_seed_species candida_albicans --busco_db saccharomycetales --cpus {threads}
        """

# in progress
rule funannotate_update:
    input:
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa"
    output:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
        funannotate_update_out = directory("results/{prefix}/funannotate/{sample}/update_results/")
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}"
    threads: 8
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate update --input {params.out_dir} --cpus {threads}
        """


# This downloads the databases needed for InterProScan and moves them to the directory bound via singularity
# in progress for local download (replace path with config["funqcd_lib"])
rule interproscan_data_dl:
    output:
        #interproscan_data_dl_dir = directory('interproscan-5.71-102.0/data/')
        interproscan_data = directory(config["funqcd_lib"] + "interproscan_data/data/antifam/")
    shell:
        """
        curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/alt/interproscan-data-5.71-102.0.tar.gz
        curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/alt/interproscan-data-5.71-102.0.tar.gz.md5
        if md5sum -c --quiet interproscan-data-5.71-102.0.tar.gz.md5; then
        tar -pxzf interproscan-data-5.71-102.0.tar.gz
        mv interproscan-5.71-102.0/data/* lib/interproscan_data/data/
        rm interproscan-data-5.71-102.0.tar.gz
        rm interproscan-data-5.71-102.0.tar.gz.md5
        rm -r interproscan-5.71-102.0
        else
        echo 'Error downloading InterProScan data'
        fi
        """

# This initializes and tests the InterProScan data, using a script called ips_setup.sh
# rule interproscan_data_init:
#     input: 
#         interproscan_data = config["funqcd_lib"] + "interproscan_data/data/antifam/"
#     output:
#         test_log = config["funqcd_lib"] + "interproscan_data/test_log.txt"
#     singularity:
#         "docker://interpro/interproscan:5.71-102.0"
#     shell:
#         """
#         bash ips_setup.sh > lib/interproscan_data/test_log.txt
#         """

rule interproscan:
    input:
        interproscan_data = config["funqcd_lib"] + "interproscan_data/data/antifam/",
        test_log = config["funqcd_lib"] + "interproscan_data/test_log.txt",
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa"
    output:
        interproscan_out = "results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml"
    singularity:
        "docker://interpro/interproscan:5.71-102.0"
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/interproscan/",
        sample = "{sample}",
        # cpus = config["ncores"]
    threads: 8
    shell:
        """
        bash /opt/interproscan/interproscan.sh --input {input.funannotate_predict_proteins} --output-dir {params.out_dir} \
        --disable-precalc --cpu {threads}
        """

# in progress
# this downloads the databases needed for eggnog to run
rule eggnog_data_dl:
    output:
        eggnog_data = config["funqcd_lib"] + "eggnog_data/eggnog.db"
    singularity:
        "docker://nanozoo/eggnog-mapper:2.1.9--4f2b6c0"
    shell:
        """
        download_eggnog_data.py -y --data_dir lib/eggnog_data/
        """

# in progress
rule eggnog:
    input:
        eggnog_data = config["funqcd_lib"] + "eggnog_data/eggnog.db",
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa"
    output:
        eggnog_out = "results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations"
    singularity:
        "docker://nanozoo/eggnog-mapper:2.1.9--4f2b6c0"
    params:
        eggnog_data_dir = config["funqcd_lib"] + "eggnog_data/",
        out_dir = "results/{prefix}/funannotate/{sample}/eggnog/",
        sample = "{sample}"
    threads: 8
    shell:
        """
        emapper.py -i {input.funannotate_predict_proteins} --itype proteins --data_dir {params.eggnog_data_dir} -m diamond \
        --output {params.sample} --output_dir {params.out_dir} --cpu {threads}
        """

# in progress
rule funannotate_annotate:
    input:
        funannotate_predict_out = "results/{prefix}/funannotate/{sample}/update_results/",
        interproscan_out = "results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml",
        eggnog_out = "results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations",
        busco_db = config["funqcd_lib"] + "busco/lineages/saccharomycetes_odb10/dataset.cfg"
    output:
        funannotate_annotate_proteins = "results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa"
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}",
        # cpus = config["ncores"]
    threads: 8
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate annotate -i {input.funannotate_predict_out} -o {params.out_dir} --cpus {threads} \
        --iprscan {input.interproscan_out} --eggnog {input.eggnog_out} --busco_db saccharomycetes_odb10
        """

rule busco_final:
    input:
        funannotate_annotate_proteins = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa", prefix = PREFIX, sample = SAMPLE)
    output:
        busco_out = "results/{prefix}/busco/busco_output/batch_summary.txt"
    params:
        prefix = "{prefix}",
        busco_db = config["funqcd_lib"] + "busco/"
    threads: 8
    singularity:
        "docker://ezlabgva/busco:v5.7.0_cv1"
    shell:
        """
        mkdir -p results/{params.prefix}/busco/input/
        cp results/{params.prefix}/funannotate/*/annotate_results/*.proteins.fa results/{params.prefix}/busco/input/
        busco -f --in results/{params.prefix}/busco/input/ --mode protein --lineage_dataset saccharomycetes_odb10 --out_path results/{params.prefix}/busco/ -c {threads} --out busco_output --offline --download_path {params.busco_db}
        """
