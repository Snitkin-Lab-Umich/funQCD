# Author: Ali Pirani, Dhatri Badri, and Joseph Hale

configfile: "config/config.yaml"


import pandas as pd
import os

PREFIX = config["prefix"]

#samples_df = pd.read_csv(f'results/{PREFIX}/assembly_pass_samples.csv')
# this is the new sample file created by the previous step in the pipeline (funQCD_ASSEMBLY)
samples_df = pd.read_csv(f'results/{PREFIX}/assembly_pass_samples.csv')
SAMPLE = list(samples_df['sample_id'])


rule all:
    input:
        check_update = expand("results/{prefix}/funannotate/{sample}/update_results/annotation_check.txt",prefix=PREFIX,sample=SAMPLE),

# The order is sort, mask, train, predict, update, interproscan, eggnong, annotate
rule funannotate_sort:
    input:
        #assembly_check = "results/{prefix}/quast/{sample}/assembly_check.tsv",
        # consider adding the quast and auriclass reports here as inputs when running downloaded assemblies through the pipeline
        # fastqc and coverage steps cannot be run on the assemblies themselves, since it requires the reads
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        sorted_assembly = "results/{prefix}/funannotate/{sample}/{sample}_sorted.fa",
    params:
        #spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
        out_dir = "results/{prefix}/funannotate/{sample}/",
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    resources:
        mem_mb = 5000,
        runtime=60,
    shell:
        """
        funannotate clean -i {input.spades_l1000_assembly} -o {params.out_dir}{wildcards.sample}_cleaned.fa
        funannotate sort -i {params.out_dir}{wildcards.sample}_cleaned.fa -o {params.out_dir}{wildcards.sample}_sorted.fa --minlen 0
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
    threads: 8
    resources:
        mem_mb = 10000,
        runtime = 960,
    singularity:
        "docker://dfam/tetools:1.89.2"
    shell:
        """
        RepeatMasker -xsmall -dir {params.out_dir} -lib {params.repeat_lib} \
        results/{wildcards.prefix}/funannotate/{wildcards.sample}/{wildcards.sample}_sorted.fa -pa {threads}
        mv {params.out_dir}/{wildcards.sample}_sorted.fa.masked {params.out_dir}/{wildcards.sample}_masked.fa
        """


# This time, I want to attempt to use RNA-seq data that has already been concatenated, trimmed, and normalized

# this requires the RNA-seq data from Teresa, and assumes that the files are in a specific format
# specifically, that read 1 contains 'R1' in the file name, that all files are in fastq.gz format and in RF order for stranded RNA-seq, and that they can be ordered alphabetically
# if this rule fails, empty training files will be created
rule funannotate_train:
    input:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa"
    output:
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
        funannotate_training_pasa_gff = "results/{prefix}/funannotate/{sample}/training/funannotate_train.pasa.gff3",
        funannotate_training_stringtie = "results/{prefix}/funannotate/{sample}/training/funannotate_train.stringtie.gtf",
        funannotate_training_transc_align = "results/{prefix}/funannotate/{sample}/training/funannotate_train.transcripts.gff3",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        #rna_data_r1 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R1_' in x and 'fastq.gz' in x])),
        #rna_data_r2 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R2_' in x and 'fastq.gz' in x])),
        #rna_data_r1 = sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R1_' in x and 'fastq.gz' in x])[0],
        #rna_data_r2 = sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R2_' in x and 'fastq.gz' in x])[0],
        rna_data_r1 = config["funqcd_lib"] + 'training_data/rna_r1_cat_reads_trimmed.fastq.gz',
        rna_data_r2 = config["funqcd_lib"] + 'training_data/rna_r2_cat_reads_trimmed.fastq.gz',
        rna_data_r1_norm = config["funqcd_lib"] + 'training_data/normalized/rna_r1_cat_reads_trimmed.fastq.gz.normalized_K25_maxC50_minC5_maxCV10000.fq',
        rna_data_r2_norm = config["funqcd_lib"] + 'training_data/normalized/rna_r2_cat_reads_trimmed.fastq.gz.normalized_K25_maxC50_minC5_maxCV10000.fq',
        mem_g = "30G",
        funannotate_training_dir = "results/{prefix}/funannotate/{sample}/training/",
    threads: 8
    resources:
        mem_mb = 32000,
        runtime = 480,
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        set +e
        timeout 7h funannotate train --input {input.masked_assembly} --out {params.out_dir} \
        --left {params.rna_data_r1} --left_norm {params.rna_data_r1_norm} \
        --right {params.rna_data_r2} --right_norm {params.rna_data_r2_norm} --stranded RF \
        --jaccard_clip --species "Candida auris" --isolate {wildcards.sample} --cpus {threads} --memory {params.mem_g}
        exitcode=$?
        if [ $exitcode != 0 ];
        then
            rm -r {params.funannotate_training_dir}
            mkdir -p {params.funannotate_training_dir}
            echo -n > {output.funannotate_training_rna_bam}
            echo -n > {output.funannotate_training_pasa_gff}
            echo -n > {output.funannotate_training_stringtie}
            echo -n > {output.funannotate_training_transc_align}
            exit 0
        fi
        """

# This should automatically detect the four training files generated previously, even without explicit input
# All steps should run with 'pasa' or 'rna-bam' under Training-Methods. Nothing should run with 'busco'.
rule funannotate_predict:
    input:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa",
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
        #busco_db = config["funqcd_lib"] + "busco/lineages/saccharomycetes_odb10/dataset.cfg"
    output:
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        predict_out_dir = "results/{prefix}/funannotate/{sample}/predict_results/",
        #predict_tmp_dir = "results/{prefix}/funannotate/{sample}/predict_misc/{sample}_predict_tmp/",
        predict_tmp_dir = "/tmp/{sample}_predict_tmp/",
        genemark_path = config["funqcd_lib"] + "genemark/gmes_linux_64_4/",
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 360,
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    priority: 1
    shell:
        # this needs a special check to skip prediction if the training files are empty
        # (by default, prediction can still run even if training failed)
        """
        set +e
        trainfile=$(wc -l < {input.funannotate_training_rna_bam})
        if [ $trainfile != 0 ];
        then
            timeout 5h funannotate predict --input {input.masked_assembly} --out {params.out_dir} \
            --species {wildcards.sample} --force \
            --busco_seed_species candida_albicans --busco_db saccharomycetes_odb10 --cpus {threads} \
            --GENEMARK_PATH {params.genemark_path} --tmpdir {params.predict_tmp_dir}
        fi
        exitcode=$?
        if [ $exitcode != 0 ] || [ $trainfile == 0 ];
        then
            rm -r {params.predict_out_dir}
            mkdir -p {params.predict_out_dir}
            echo -n > {output.funannotate_predict_proteins}
            exit 0
        fi
        """

rule funannotate_update:
    input:
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa",
    output:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 360,
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    priority: 2
    shell:
        """
        set +e
        predictfile=$(wc -l < {input.funannotate_predict_proteins})
        if [ $predictfile != 0 ];
        then
            timeout 5h funannotate update --input {params.out_dir} --cpus {threads}
        fi
        exitcode=$?
        if [ $exitcode != 0 ] || [ $predictfile == 0 ];
        then
            echo -n > {output.funannotate_update_proteins}
        fi
        """


# this step removes the entire contents of training/, predict_misc/, and update_misc/
# these files are needed for funannotate_update, but are too large for permanent storage
rule more_cleanup:
    input:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
    output:
        cleanup_check = "results/{prefix}/funannotate/{sample}/training/{sample}_cleanup_complete.txt",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
    threads: 1
    resources:
        mem_mb = 5000,
        runtime = 120,
    priority: 5
    shell:
        """
        rm -r -f {params.out_dir}training/*
        rm -r -f {params.out_dir}predict_misc/*
        rm -r -f {params.out_dir}update_misc/*
        touch {output.cleanup_check}
        """


# this step removes several large files generated by funannotate_train
# these files are needed for funannotate_update, but are too large for permanent storage
# funannotate_update can still be run with the remaining files present in the normalize/ directory
# rule cleanup:
#     input:
#         funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
#     output:
#         cleanup_check = "results/{prefix}/funannotate/{sample}/training/{sample}_cleanup_complete.txt",
#     params:
#         out_dir = "results/{prefix}/funannotate/{sample}/",
#     threads: 1
#     resources:
#         mem_mb = 5000,
#         runtime = 60,
#     priority: 5
#     shell:
#         """
#         rm -f {params.out_dir}training/left.fq.gz
#         rm -f {params.out_dir}training/right.fq.gz
#         rm -r -f {params.out_dir}training/trimmomatic/
#         rm -r -f {params.out_dir}training/trinity_gg/
#         rm -f {params.out_dir}training/normalize/*CV10000.fq
#         touch {output.cleanup_check}
#         """


def check_update_fun(predict_file,output_file,example_qc,intermediate_qc,new_sample_file,sample_name):
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
    with open(intermediate_qc,'a') as fh_inter, open(output_file,'w') as fh_out,open(new_sample_file, 'a') as fh_new:
        if os.path.getsize(predict_file) == 0:
            _ = fh_out.write('FUNANNOTATE PREDICTION FAILED')
            fail_line_list = fail_line.strip().split('\t')
            fail_line_list[0] = sample_name
            fail_line_list[5] = sample_name + '_R1.fastq.gz'
            _ = fh_inter.write('\t'.join(fail_line_list) + '\n')
        else:
            _ = fh_out.write('FUNANNOTATE PREDICTION PASSED')
            _ = fh_new.write(sample_name + '\n')

rule check_update:
    input:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
        cleanup_check = "results/{prefix}/funannotate/{sample}/training/{sample}_cleanup_complete.txt",
    output:
        outp = "results/{prefix}/funannotate/{sample}/update_results/annotation_check.txt",
    params:
        example_qc = config["example_qc_summary_file"],
        #new_sample_file = "config/predict_pass_samples.csv",
        new_sample_file = "results/{prefix}/predict_pass_samples.csv",
        intermediate_qc = "results/{prefix}/funannotate/failed_prediction_qc_summary.tsv",
    resources:
        mem_mb = 2000,
        runtime = 20,
    run:
        check_update_fun(
            input.funannotate_update_proteins,output.outp,params.example_qc,
            params.intermediate_qc,params.new_sample_file,wildcards.sample)
