# Author: Ali Pirani, Dhatri Badri, and Joseph Hale

configfile: "config/config.yaml"

#include: "fQCD_report.smk"

import pandas as pd
import os

PREFIX = config["prefix"]

samples_df = pd.read_csv(f'results/{PREFIX}/predict_pass_samples.csv')
# this is the new sample file created by the previous step in the pipeline (funQCD_PREDICTION)
SAMPLE = list(samples_df['sample_id'])


rule all:
    input:
        qc_report_final = expand("results/{prefix}/multiqc/{prefix}_final_qc_summary.tsv",prefix=PREFIX),

rule interproscan:
    input:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
       #funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa",
    output:
        interproscan_out = "results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml"
    singularity:
        "docker://interpro/interproscan:5.71-102.0"
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/interproscan/",
        #interproscan_data = config["funqcd_lib"] + "interproscan_data/data/antifam/",
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 360
    shell:
        """
        bash /opt/interproscan/interproscan.sh --input {input.funannotate_update_proteins} --output-dir {params.out_dir} \
        --disable-precalc --cpu {threads}
        """

rule eggnog:
    input:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
        #funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa"
    output:
        eggnog_out = "results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations"
    singularity:
        "docker://nanozoo/eggnog-mapper:2.1.9--4f2b6c0"
    params:
        eggnog_data_dir = config["funqcd_lib"] + "eggnog_data/",
        out_dir = "results/{prefix}/funannotate/{sample}/eggnog/",
    threads: 8
    resources:
        mem_mb = 10000,
        runtime = 300
    shell:
        """
        emapper.py -i {input.funannotate_update_proteins} --itype proteins --data_dir {params.eggnog_data_dir} -m diamond \
        --output {wildcards.sample} --output_dir {params.out_dir} --cpu {threads} --override
        """

rule funannotate_annotate:
    input:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
        #funannotate_predict_out = "results/{prefix}/funannotate/{sample}/update_results/",
        interproscan_out = "results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml",
        eggnog_out = "results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations",
    output:
        funannotate_annotate_proteins = "results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa",
        funannotate_annotate_assembly = "results/{prefix}/funannotate/{sample}/annotate_results/{sample}.scaffolds.fa",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        funannotate_update_dir = "results/{prefix}/funannotate/{sample}/update_results/",
    threads: 8
    resources:
        mem_mb = 3000,
        runtime = 80
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate annotate -i {params.funannotate_update_dir} -o {params.out_dir} --cpus {threads} \
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
        busco_db = config["funqcd_lib"] + "busco/",
    threads: 8
    resources:
        mem_mb = 20000,
        runtime = 2800,
    singularity:
        "docker://ezlabgva/busco:v5.7.0_cv1"
    shell:
        """
        mkdir -p results/{wildcards.prefix}/busco/input/prot/
        mkdir -p results/{wildcards.prefix}/busco/input/nucl/
        cp results/{wildcards.prefix}/funannotate/*/annotate_results/*.proteins.fa results/{wildcards.prefix}/busco/input/prot
        cp results/{wildcards.prefix}/funannotate/*/annotate_results/*.scaffolds.fa results/{wildcards.prefix}/busco/input/nucl
        busco -f --in results/{wildcards.prefix}/busco/input/prot --mode protein --lineage_dataset saccharomycetes_odb10 --out_path results/{wildcards.prefix}/busco/ -c {threads} --out busco_output_prot --offline --download_path {params.busco_db}
        busco -f --in results/{wildcards.prefix}/busco/input/nucl --mode genome --lineage_dataset saccharomycetes_odb10 --out_path results/{wildcards.prefix}/busco/ -c {threads} --out busco_output_nucl --offline --download_path {params.busco_db}
        rm -rf RM_*
        """

rule multiqc:
    input:
        quast_report = expand("results/{prefix}/quast/{sample}/report.tsv", sample = SAMPLE, prefix = PREFIX),
        aftertrim_fastqc_report_fwd = expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html", sample = SAMPLE, prefix = PREFIX),
        aftertrim_fastqc_report_rev = expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Reverse/{sample}_R2_trim_paired_fastqc.html", sample = SAMPLE, prefix = PREFIX),
        raw_fastqc_report_fwd = expand("results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html", sample = SAMPLE, prefix = PREFIX),
        raw_fastqc_report_rev = expand("results/{prefix}/quality_raw/{sample}/{sample}_Reverse/{sample}_R2_fastqc.html", sample = SAMPLE, prefix = PREFIX),
        busco_out_p = "results/{prefix}/busco/busco_output_prot/batch_summary.txt",
        busco_out_n = "results/{prefix}/busco/busco_output_nucl/batch_summary.txt",
    output:
        multiqc_report = "results/{prefix}/multiqc/{prefix}_QC_report.html",
    params:
        outdir = "results/{prefix}/multiqc",
        quast_dir = "results/{prefix}/quast/",
        raw_fastqc_dir = "results/{prefix}/quality_raw/",
        aftertrim_fastqc_dir = "results/{prefix}/quality_aftertrim/",
        busco_dir_p = "results/{prefix}/busco/busco_output_prot/",
        busco_dir_n = "results/{prefix}/busco/busco_output_nucl/",
    resources:
        mem_mb = 5000,
        runtime = 600
    threads: 4
    singularity:
        "docker://multiqc/multiqc:v1.25.1"
    shell:
        """
        multiqc -f --outdir {params.outdir} -n {wildcards.prefix}_QC_report -i {wildcards.prefix}_QC_report \
        {params.quast_dir} {params.raw_fastqc_dir} {params.aftertrim_fastqc_dir} {params.busco_dir_p} {params.busco_dir_n}
        """

rule qc_report_final:
    input:
        multiqc_report = "results/{prefix}/multiqc/{prefix}_QC_report.html",
        auriclass_report = expand("results/{prefix}/auriclass/{sample}/{sample}_report.tsv", sample = SAMPLE, prefix = PREFIX), 
        raw_coverage_report = expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", sample = SAMPLE, prefix = PREFIX),
    output:
        summary_output = "results/{prefix}/multiqc/{prefix}_final_qc_summary.tsv",
    params:
        multiqc_dir = "results/{prefix}/multiqc/{prefix}_QC_report_data/",
        auriclass_dir = "results/{prefix}/auriclass/", 
        raw_coverage_dir = "results/{prefix}/raw_coverage/",
        intermediate_qc_assembly = "results/{prefix}/quast/failed_assembly_qc_summary.tsv",
        intermediate_qc_prediction = "results/{prefix}/funannotate/failed_prediction_qc_summary.tsv",
    resources:
        mem_mb = 2000,
        runtime = 60,
    script:
        "QC_summary.py"

