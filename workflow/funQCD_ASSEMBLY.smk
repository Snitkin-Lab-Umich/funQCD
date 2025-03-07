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

if not os.path.exists("results/"):
    os.system("mkdir %s" % "results/")

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)


# define function for the downsample rule

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






# main pipeline

rule all:
    input:
        check_assembly = expand("results/{prefix}/quast/{sample}/assembly_check.tsv",sample=SAMPLE,prefix=PREFIX),

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



# define a function for adjusting runtime based on attempt number
def get_assembly_runtime(wildcards,attempt):
    return(30 + (60 * (attempt - 1)))

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
        #runtime = 30,
        runtime = get_assembly_runtime,
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
    resources:
        mem_mb = 5000,
        runtime = 30
    singularity:
        "docker://quay.io/biocontainers/auriclass:0.5.4--pyhdfd78af_0"
    shell:
        "auriclass --name {wildcards.sample} -o {output.auriclass_report} {input.spades_l1000_assembly}"


def check_assembly_fun(
    quast_report,intermediate_qc,example_qc,sample_name,output_file,new_sample_file,
    max_contigs,min_n50,min_assembly_length,max_assembly_length):
    # determine the format of the QC file from the example
    with open(example_qc,'r') as fh_example:
        header_line = fh_example.readline()
        fail_line = fh_example.readline()
    # if the intermediate qc file is empty, make it with the header
    if not os.path.isfile(intermediate_qc):
        with open(intermediate_qc,'w') as fh_inter:
            _ = fh_inter.write(header_line)
    # do the same with the new sample file
    if not os.path.isfile(new_sample_file):
        with open(new_sample_file,'w') as fh_new:
            _ = fh_new.write('sample_id\n')
    # read the quast report for basic assembly stats
    with open(quast_report,'r') as fh_quast:
        qc_data = {}
        for line in fh_quast:
            d = line.strip().split('\t')
            if d[0] in ['# contigs','Total length','N50']:
                qc_data[d[0]] = int(d[1])
    # based on this data, write to the intermediate qc file and the output file (which is essentially just for snakemake to see completion)
    # this does assume that the file names are at fixed positions (0 and 5)
    with open(intermediate_qc,'a') as fh_inter, open(output_file,'w') as fh_out,open(new_sample_file, 'a') as fh_new:
        if any([
            qc_data['# contigs'] > max_contigs, qc_data['N50'] < min_n50, 
            qc_data['Total length'] < min_assembly_length, qc_data['Total length'] > max_assembly_length]):
            _ = fh_out.write('ASSEMBLY QC FAILED')
            fail_line_list = fail_line.strip().split('\t')
            fail_line_list[0] = sample_name
            fail_line_list[5] = sample_name + '_R1_fastq.gz'
            _ = fh_inter.write('\t'.join(fail_line_list) + '\n')
        else:
            _ = fh_out.write('ASSEMBLY QC PASSED')
            _ = fh_new.write(sample_name + '\n')

rule check_assembly:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
        quast_report = "results/{prefix}/quast/{sample}/report.tsv",
        auriclass_report = "results/{prefix}/auriclass/{sample}/{sample}_report.tsv",
        aftertrim_fastqc_report_fwd = "results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html",
        raw_fastqc_report_fwd = "results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html",
        coverage = "results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json",
    output:
        assembly_check = "results/{prefix}/quast/{sample}/assembly_check.tsv",
    resources:
        mem_mb = 5000,
        runtime = 20,
    params:
        example_qc = config["example_qc_summary_file"],
        #new_sample_file = "config/assembly_pass_samples.csv",
        new_sample_file = "results/{prefix}/assembly_pass_samples.csv",
        intermediate_qc = "results/{prefix}/quast/failed_assembly_qc_summary.tsv",
        max_contigs = config["max_contigs"],
        min_n50 = config["min_n50"],
        min_alen = config["min_assembly_length"],
        max_alen = config["max_assembly_length"],
    run:
        check_assembly_fun(
            input.quast_report,params.intermediate_qc,params.example_qc,wildcards.sample,output.assembly_check,
            params.new_sample_file,params.max_contigs,params.min_n50,params.min_alen,params.max_alen)

