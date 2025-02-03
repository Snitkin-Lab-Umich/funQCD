# Author: Ali Pirani, Dhatri Badri, and Joseph Hale

configfile: "config/config.yaml"

#include: "fQCD_report.smk"

import pandas as pd
import os

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

PREFIX = config["prefix"]

SHORTREADS = list(samples_df['sample_id'])

# make all directories
# if not os.path.exists(f'results/{PREFIX}/{PREFIX}_Report'):
#     os.system(f'mkdir -p results/{PREFIX}/{PREFIX}_Report')
# for d in ['scripts','data','multiqc','fig']:
#     if not os.path.exists(f'results/{PREFIX}/{PREFIX}_Report/{d}'):
#         os.system(f'mkdir results/{PREFIX}/{PREFIX}_Report/{d}')

if not os.path.exists("results/"):
    os.system("mkdir %s" % "results/")

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)










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
        trimmomatic_pe = expand("results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE, prefix=PREFIX),
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
        funannotate_sort = expand("results/{prefix}/funannotate/{sample}/{sample}_sorted.fa", sample=SAMPLE, prefix=PREFIX),
        repeatmasker = expand("results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa", sample = SAMPLE, prefix = PREFIX),
        funannotate_train = expand("results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam", sample=SAMPLE, prefix=PREFIX),
        funannotate_predict = expand("results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),
        funannotate_update = expand("results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),        
        #interproscan_data_dl = config["funqcd_lib"] + "interproscan_data/data/antifam/",
        interproscan = expand("results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml", sample=SAMPLE, prefix=PREFIX),
        #eggnog_data_dl = config["funqcd_lib"] + "eggnog_data/eggnog.db",
        eggnog = expand("results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations", sample=SAMPLE, prefix=PREFIX),
        funannotate_annotate = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),
        busco_final = expand("results/{prefix}/busco/busco_output_prot/batch_summary.txt", prefix = PREFIX),
        qc_report_final = expand("results/{prefix}/multiqc/{prefix}_final_qc_summary.tsv", prefix = PREFIX),

rule coverage:
    input:
        r1 = config["short_reads"] + "/" + "{sample}_R1.fastq.gz",
        r2 = config["short_reads"] + "/" + "{sample}_R2.fastq.gz",
    output:
        coverage = "results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json",
    params:
        size=config["genome_size"]
    resources:
        mem_mb=2000,
        runtime=20
    singularity:
        "docker://staphb/fastq-scan:1.0.1"
    shell:
        "zcat {input.r1} {input.r2} | fastq-scan -g {params.size} > {output.coverage}"

rule quality_raw:
    input:
        r1 = config["short_reads"] + "/" + "{sample}_R1.fastq.gz",
        r2 = config["short_reads"] + "/" + "{sample}_R2.fastq.gz",
    output:
        raw_fastqc_report_fwd = "results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html",
        raw_fastqc_report_rev = "results/{prefix}/quality_raw/{sample}/{sample}_Reverse/{sample}_R2_fastqc.html",
    log:
        "logs/{prefix}/quality_raw/{sample}/{sample}.log"
    params:
        outdir="results/{prefix}/quality_raw/{sample}/{sample}"
    resources:
        mem_mb = 1000,
        runtime = 10
    singularity:
        "docker://staphb/fastqc:0.12.1"
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        """

rule trimmomatic_pe:
    input:
        r1 = config["short_reads"] + "/" + "{sample}_R1.fastq.gz",
        r2 = config["short_reads"] + "/" + "{sample}_R2.fastq.gz",
    output:
        r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_paired.fastq.gz", 
        # reads where trimming entirely removed the mate
        r1_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_unpaired.fastq.gz",
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
        lead_trail_qual=config["lead_trail_qual"],
    log:
        "logs/{prefix}/trimmomatic/{sample}/{sample}.log"
    threads: 8
    singularity:
        "docker://staphb/trimmomatic:0.39"
    resources:
        mem_mb = 5000,
        runtime = 30
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {threads} \
        ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}
        """

rule quality_aftertrim:
    input:
        r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_paired.fastq.gz",
    output:
        aftertrim_fastqc_report_fwd = "results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html",
        aftertrim_fastqc_report_rev = "results/{prefix}/quality_aftertrim/{sample}/{sample}_Reverse/{sample}_R2_trim_paired_fastqc.html",
    log:
        "logs/{prefix}/{sample}/quality_aftertrim/{sample}.log"
    params:
        outdir="results/{prefix}/quality_aftertrim/{sample}/{sample}"
    singularity:
        "docker://staphb/fastqc:0.12.1"
    resources:
        mem_mb = 3000,
        runtime = 30
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        """

rule downsample:
    input:
        r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_paired.fastq.gz",
    output:
        outr1 = "results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz",
        outr2 = "results/{prefix}/downsample/{sample}/{sample}_R2_trim_paired.fastq.gz",
    params:
        gsize = config["genome_size"],
    resources:
        mem_mb = 3000,
        runtime = 30
    run:
        downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize})


rule assembly:
    input:
        r1 = "results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/downsample/{sample}/{sample}_R2_trim_paired.fastq.gz",
    output:
        spades_assembly = "results/{prefix}/spades/{sample}/contigs.fasta",
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        mem_g = 15
    singularity:
        "docker://staphb/spades:4.0.0"
    threads: 8
    resources:
        mem_mb = 15000,
        runtime=30
    shell:
        "spades.py --isolate --pe1-1 {input.r1} --pe1-2 {input.r2} -o {params.out_dir} --threads {threads} --memory {params.mem_g}"


rule bioawk:
    input:
        spades_assembly = "results/{prefix}/spades/{sample}/contigs.fasta",
    output:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta"
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        bioawk_params = config["bioawk"],
        prefix = "{sample}"
    resources:
        mem_mb = 5000,
        runtime = 30
    singularity:
        "docker://lbmc/bioawk:1.0"
    shell:
        """
        ./bioawk.sh {input.spades_assembly} {output.spades_l1000_assembly} {params.out_dir} {params.prefix}
        """


rule quast:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        quast_report = "results/{prefix}/quast/{sample}/report.tsv",
    params: 
        outdir = "results/{prefix}/quast/{sample}/",
        prefix = "{sample}",
    resources:
        mem_mb = 5000,
        runtime = 30
    singularity:
        "docker://staphb/quast:5.0.2"
    shell:
        """
        ./quast.sh {input.spades_l1000_assembly} {params.outdir} 
        """

rule auriclass:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        auriclass_report = "results/{prefix}/auriclass/{sample}/{sample}_report.tsv",
    params: 
        sample = "{sample}",
    resources:
        mem_mb = 5000,
        runtime = 30
    singularity:
        "docker://quay.io/biocontainers/auriclass:0.5.4--pyhdfd78af_0"
    shell:
        "auriclass --name {params.sample} -o {output.auriclass_report} {input.spades_l1000_assembly}"


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

# The order is setup, mask, train, predict, update, interproscan, eggnong, annotate
rule funannotate_sort:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        sorted_assembly = "results/{prefix}/funannotate/{sample}/{sample}_sorted.fa"
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}"
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    resources:
        mem_mb = 5000,
        runtime=60
    shell:
        """
        funannotate clean -i {input.spades_l1000_assembly} -o {params.out_dir}{params.sample}_cleaned.fa
        funannotate sort -i {params.out_dir}{params.sample}_cleaned.fa -o {params.out_dir}{params.sample}_sorted.fa --minlen 0
        """

# This will generate a directory named RM_* with a large number of temporary files. 
# There does not appear to be a way to change this output from RepeatMasker itself, so these files are deleted in a later rule
# (Deleting RM_* directories in this rule breaks RepeatMasker runs in progress)
rule repeatmasker:
    input:
        sorted_assembly = "results/{prefix}/funannotate/{sample}/{sample}_sorted.fa"
    output:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa"
    params:
        out_dir = "results/{prefix}/repeatmasker/{sample}/",
        repeat_lib = config["funqcd_lib"] + "repeat_libraries/fungi_b8441/b8441_fungi_repeatlib.fa",
        prefix = "{prefix}",
        sample = "{sample}",
    threads: 8
    resources:
        mem_mb = 5000,
        runtime = 120,
    singularity:
        "docker://dfam/tetools:1.89.2"
    shell:
        """
        RepeatMasker -xsmall -dir {params.out_dir} -lib {params.repeat_lib} \
        results/{params.prefix}/funannotate/{params.sample}/{params.sample}_sorted.fa -pa {threads}
        mv {params.out_dir}/{params.sample}_sorted.fa.masked {params.out_dir}/{params.sample}_masked.fa
        """

# this requires the RNA-seq data from Teresa, and assumes that the files are in a specific format
# specifically, that read 1 contains 'R1' in the file name, that all files are in fastq.gz format and in RF order for stranded RNA-seq, and that they can be ordered alphabetically
# This is also currently hard-coded to use 30G of memory
rule funannotate_train:
    input:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa"
    output:
        funannotate_training_dir = directory("results/{prefix}/funannotate/{sample}/training/"),
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
        funannotate_training_pasa_gff = "results/{prefix}/funannotate/{sample}/training/funannotate_train.pasa.gff3",
        funannotate_training_stringtie = "results/{prefix}/funannotate/{sample}/training/funannotate_train.stringtie.gtf",
        funannotate_training_transc_align = "results/{prefix}/funannotate/{sample}/training/funannotate_train.transcripts.gff3",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}",
        rna_data_r1 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R1_' in x and 'fastq.gz' in x])),
        rna_data_r2 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R2_' in x and 'fastq.gz' in x])),
        mem_g = "30G",
        test = vars(workflow.resources)
    threads: 8
    resources:
        mem_mb = 32000,
        runtime = 780
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate train --input {input.masked_assembly} --out {params.out_dir} \
        --left {params.rna_data_r1} --right {params.rna_data_r2} --stranded RF \
        --jaccard_clip --species "Candida auris" --isolate {params.sample} --cpus {threads} --memory {params.mem_g}
        """

# This should automatically detect the four training files generated previously, even without explicit input
# All steps should run with 'pasa' or 'rna-bam' under Training-Methods. Nothing should run with 'busco'.
rule funannotate_predict:
    input:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa",
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
        busco_db = config["funqcd_lib"] + "busco/lineages/saccharomycetes_odb10/dataset.cfg"
    output:
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa",
        funannotate_predict_out = directory("results/{prefix}/funannotate/{sample}/predict_results/")
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}",
        genemark_path = config["funqcd_lib"] + "genemark/gmes_linux_64_4/",
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 500
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        # removed --species "Candida auris" --isolate {params.sample}
        # consider changing output names and re-adding
        """
        funannotate predict --input {input.masked_assembly} --out {params.out_dir} \
        --species {params.sample} --force \
        --busco_seed_species candida_albicans --busco_db saccharomycetes_odb10 --cpus {threads} \
        --GENEMARK_PATH {params.genemark_path}
        """

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
    resources:
        mem_mb = 15000,
        runtime = 600
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate update --input {params.out_dir} --cpus {threads}
        """


# This downloads the databases needed for InterProScan and moves them to the directory bound via singularity
# rule interproscan_data_dl:
#     output:
#         #interproscan_data_dl_dir = directory('interproscan-5.71-102.0/data/')
#         interproscan_data = directory(config["funqcd_lib"] + "interproscan_data/data/antifam/")
#     shell:
#         """
#         curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/alt/interproscan-data-5.71-102.0.tar.gz
#         curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/alt/interproscan-data-5.71-102.0.tar.gz.md5
#         if md5sum -c --quiet interproscan-data-5.71-102.0.tar.gz.md5; then
#         tar -pxzf interproscan-data-5.71-102.0.tar.gz
#         mv interproscan-5.71-102.0/data/* lib/interproscan_data/data/
#         rm interproscan-data-5.71-102.0.tar.gz
#         rm interproscan-data-5.71-102.0.tar.gz.md5
#         rm -r interproscan-5.71-102.0
#         else
#         echo 'Error downloading InterProScan data'
#         fi
#         """

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
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 360
    shell:
        """
        bash /opt/interproscan/interproscan.sh --input {input.funannotate_predict_proteins} --output-dir {params.out_dir} \
        --disable-precalc --cpu {threads}
        """

# this downloads the databases needed for eggnog to run
# rule eggnog_data_dl:
#     output:
#         eggnog_data = config["funqcd_lib"] + "eggnog_data/eggnog.db"
#     singularity:
#         "docker://nanozoo/eggnog-mapper:2.1.9--4f2b6c0"
#     shell:
#         """
#         download_eggnog_data.py -y --data_dir lib/eggnog_data/
#         """

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
    resources:
        mem_mb = 10000,
        runtime = 300
    shell:
        """
        emapper.py -i {input.funannotate_predict_proteins} --itype proteins --data_dir {params.eggnog_data_dir} -m diamond \
        --output {params.sample} --output_dir {params.out_dir} --cpu {threads} --override
        """

rule funannotate_annotate:
    input:
        funannotate_predict_out = "results/{prefix}/funannotate/{sample}/update_results/",
        interproscan_out = "results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml",
        eggnog_out = "results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations",
        busco_db = config["funqcd_lib"] + "busco/lineages/saccharomycetes_odb10/dataset.cfg"
    output:
        funannotate_annotate_proteins = "results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa",
        funannotate_annotate_assembly = "results/{prefix}/funannotate/{sample}/annotate_results/{sample}.scaffolds.fa",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}",
    threads: 8
    resources:
        mem_mb = 3000,
        runtime = 80
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate annotate -i {input.funannotate_predict_out} -o {params.out_dir} --cpus {threads} \
        --iprscan {input.interproscan_out} --eggnog {input.eggnog_out} --busco_db saccharomycetes_odb10
        """

# The line 'rm -rf RM_*' removes the directories that RepeatMasker generates in the working directory
rule busco_final:
    input:
        funannotate_annotate_proteins = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa", prefix = PREFIX, sample = SAMPLE),
        funannotate_annotate_nucleotides = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.scaffolds.fa", prefix = PREFIX, sample = SAMPLE),       
    output:
        busco_out_p = "results/{prefix}/busco/busco_output_prot/batch_summary.txt",
        busco_out_n = "results/{prefix}/busco/busco_output_nucl/batch_summary.txt",
    params:
        prefix = "{prefix}",
        busco_db = config["funqcd_lib"] + "busco/"
    threads: 8
    resources:
        mem_mb = 10000,
        runtime = 240
    singularity:
        "docker://ezlabgva/busco:v5.7.0_cv1"
    shell:
        """
        mkdir -p results/{params.prefix}/busco/input/prot/
        mkdir -p results/{params.prefix}/busco/input/nucl/
        cp results/{params.prefix}/funannotate/*/annotate_results/*.proteins.fa results/{params.prefix}/busco/input/prot
        cp results/{params.prefix}/funannotate/*/annotate_results/*.scaffolds.fa results/{params.prefix}/busco/input/nucl
        busco -f --in results/{params.prefix}/busco/input/prot --mode protein --lineage_dataset saccharomycetes_odb10 --out_path results/{params.prefix}/busco/ -c {threads} --out busco_output_prot --offline --download_path {params.busco_db}
        busco -f --in results/{params.prefix}/busco/input/nucl --mode genome --lineage_dataset saccharomycetes_odb10 --out_path results/{params.prefix}/busco/ -c {threads} --out busco_output_nucl --offline --download_path {params.busco_db}
        rm -rf RM_*
        """

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
        raw_fastqc_report_rev = expand("results/{prefix}/quality_raw/{sample}/{sample}_Reverse/{sample}_R2_fastqc.html", sample = SAMPLE, prefix = PREFIX),
        busco_out_p = "results/{prefix}/busco/busco_output_prot/batch_summary.txt",
        busco_out_n = "results/{prefix}/busco/busco_output_nucl/batch_summary.txt",
        #busco_out = "results/{prefix}/busco/busco_output/batch_summary.txt",
    output:
        multiqc_report = "results/{prefix}/multiqc/{prefix}_QC_report.html",
    params:
        #resultsoutdir = "results/{prefix}",
        outdir = "results/{prefix}/multiqc",
        prefix = "{prefix}",
        quast_dir = "results/{prefix}/quast/",
        raw_fastqc_dir = "results/{prefix}/quality_raw/",
        aftertrim_fastqc_dir = "results/{prefix}/quality_aftertrim/",
        busco_dir_p = "results/{prefix}/busco/busco_output_prot/",
        busco_dir_n = "results/{prefix}/busco/busco_output_nucl/",
    resources:
        mem_mb = 2000,
        runtime = 100
    singularity:
        "docker://multiqc/multiqc:v1.25.1"
    shell:
        """
        multiqc -f --outdir {params.outdir} -n {params.prefix}_QC_report -i {params.prefix}_QC_report \
        {params.quast_dir} {params.raw_fastqc_dir} {params.aftertrim_fastqc_dir} {params.busco_dir_p} {params.busco_dir_n}
        """

rule qc_report_final:
    input:
        multiqc_report = "results/{prefix}/multiqc/{prefix}_QC_report.html",
        auriclass_report = expand("results/{prefix}/auriclass/{sample}/{sample}_report.tsv", sample = SAMPLE, prefix = PREFIX), 
        raw_coverage_report = expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", sample = SAMPLE, prefix = PREFIX),
    output:
        summary_output = "results/{prefix}/multiqc/{prefix}_final_qc_summary.tsv"
    params:
        multiqc_dir = "results/{prefix}/multiqc/{prefix}_QC_report_data/",
        auriclass_dir = "results/{prefix}/auriclass/", 
        raw_coverage_dir = "results/{prefix}/raw_coverage/",
    resources:
        mem_mb = 2000,
        runtime = 60,
    script:
        "QC_summary.py"
    
    

        
