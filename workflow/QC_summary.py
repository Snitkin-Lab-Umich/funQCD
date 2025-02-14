
import pandas as pd
import os
import json
from functools import reduce
import numpy as np

def make_summary_report(input_path, output_path, report, type = 'tsv'):
    input_path,output_path = [x+'/' if x[-1] != '/' else x for x in [input_path,output_path]]
    # this is the path to the folder with each individual auriclass report
    input_dict = {}
    for input_dir in [x for x in os.listdir(input_path) if os.path.isdir(input_path + x)]:
        for input_file in os.listdir(input_path + input_dir):
            if '.' + type in input_file:
                if type == 'tsv':
                    fdata = pd.read_csv(input_path + input_dir + '/' + input_file, sep = '\t', header = 0)
                    if fdata['Sample'][0] in input_dict:
                        print('Error: overlapping data in input directories')
                    input_dict[fdata['Sample'][0]] = [fdata['QC_decision'][0], fdata['Clade'][0]]
                if type == 'json':
                    with open(input_path + input_dir + '/' + input_file, 'r') as fh:
                        fdata = json.load(fh)
                    sample = input_file.split('_coverage.json')[0]
                    if sample in input_dict:
                        print('Error: overlapping data in input directories')
                    fdata = [fdata['qc_stats']['read_total'], fdata['qc_stats']['total_bp'], fdata['qc_stats']['read_mean'], fdata['qc_stats']['coverage']]
                    input_dict[sample] = [str(x) for x in fdata]
    with open(output_path + report + '_summary.tsv','w') as fhout:
        if report == 'auriclass':
            _ = fhout.write('Sample\tauriclass_QC_decision\tauriclass_clade\n')
        if report == 'raw_coverage':
            _ = fhout.write('Sample\ttotal_reads\ttotal_bp\tmean_read_length\tcoverage\n')
        for s in input_dict:
            _ = fhout.write('\t'.join([s] + input_dict[s]) + '\n')


def final_qc_summary(multiqc_path,auriclass_path,raw_coverage_path):
    # QC reports to combine:
    # 'final coverage' table (at "results/{PREFIX}/{PREFIX}_Report/data/{PREFIX}_Final_Coverage.txt")
    # this file is made via the coverage_report rule, which pulls from raw_coverage outputs
    # this should now be handled by make_summary report
    # fastqc before trimmomatic
    # fastqc after trimmomatic
    # quast
    # contig distribution (this is covered by quast)
    # busco (proteins and nucleotides)
    # auriclass (this should be handled by make_summary_report)

    if multiqc_path[-1] != '/':
        multiqc_path+='/'
    fastqc_path = multiqc_path + 'multiqc_fastqc.txt'
    quast_path = multiqc_path + 'multiqc_quast.txt'
    busco_path = multiqc_path + 'multiqc_busco.txt'

    # load raw_coverage report
    data_raw_coverage = pd.read_csv(raw_coverage_path, sep = '\t', header = 0)

    # load fastqc report (multiqc_fastqc.txt)
    data_fastqc = pd.read_csv(fastqc_path, sep = '\t', header = 0)
    # fastqc for raw data
    data_raw_fastqc = data_fastqc.copy()[data_fastqc['Filename'].str.contains('_R1.fastq.gz')]
    # data_raw_fastqc['Sample'].replace('_R1','',regex=True, inplace=True)
    data_raw_fastqc['Sample'] = data_raw_fastqc['Sample'].replace('_R1','',regex=True)
    # fastqc after trimmomatic
    data_trim_fastqc = data_fastqc.copy()[data_fastqc['Sample'].str.contains('_R1_trim_paired')]
    data_trim_fastqc = data_trim_fastqc.drop(['Filename','File type','Encoding'],axis=1)
    # this removes redundant columns
    # data_trim_fastqc['Sample'].replace('_R1_trim_paired','',regex=True, inplace=True)
    data_trim_fastqc['Sample'] = data_trim_fastqc['Sample'].replace('_R1_trim_paired','',regex=True)
    data_trim_fastqc = data_trim_fastqc.add_prefix('after_trim_')
    data_trim_fastqc.rename(columns = {'after_trim_Sample':'Sample'}, inplace = True)

    # load quast report
    data_quast = pd.read_csv(quast_path, sep = '\t', header = 0)
    # data_quast['Sample'].replace('_contigs_l1000','',regex=True, inplace=True)
    data_quast['Sample'] = data_quast['Sample'].replace('_contigs_l1000','',regex=True)
    data_quast = data_quast[['Sample','N50','Total length','# contigs (>= 0 bp)']]

    # load busco report
    data_busco = pd.read_csv(busco_path, sep = '\t', header = 0)
    data_busco = data_busco[~(data_busco['Sample'] == 'short')]
    # removes the 'short' row
    data_busco['Sample'] = [x[-1] for x in data_busco['Sample'].str.split('odb10.')]
    # sample names should now be either {sample} or {sample}.proteins
    data_busco['Percent Complete BUSCOs'] = round((data_busco['complete'] / data_busco['total']) * 100,2)
    data_busco = data_busco[['Sample','Percent Complete BUSCOs']]
    # calculates and rounds the percentage of complete buscos, removes other columns
    data_busco_prot = data_busco.copy()[data_busco['Sample'].str.contains('.proteins')]
    data_busco_nucl = data_busco.copy()[~(data_busco['Sample'].str.contains('.proteins'))]
    # splits busco into prot and nucl results
    # data_busco_prot['Sample'].replace('.proteins','',regex=True, inplace=True)
    data_busco_prot['Sample'] = data_busco_prot['Sample'].replace('.proteins','',regex=True)
    data_busco_prot.rename(columns = {'Percent Complete BUSCOs':'Percent Complete Protein BUSCOs'}, inplace = True)
    data_busco_nucl.rename(columns = {'Percent Complete BUSCOs':'Percent Complete Nucleotide BUSCOs'}, inplace = True)
    # renames columns to distinguish nucl from prot
    data_busco_merge = pd.merge(left=data_busco_nucl,right=data_busco_prot,how='inner',on='Sample')
    # combine tables

    # load auriclass report
    data_auriclass = pd.read_csv(auriclass_path, sep = '\t', header = 0)

    # merge everything on the sample column
    data_list = [data_raw_coverage, data_raw_fastqc, data_trim_fastqc, data_quast, data_busco_merge, data_auriclass]
    data_merge = reduce(lambda left,right: pd.merge(left,right,on='Sample',how='inner'), data_list)
    
    # choose which tests to use for fastqc 
    fastqc_tests = [
        'after_trim_basic_statistics','after_trim_per_base_sequence_quality','after_trim_per_tile_sequence_quality',
        'after_trim_per_sequence_quality_scores','after_trim_per_base_sequence_content','after_trim_per_sequence_gc_content',
        'after_trim_per_base_n_content','after_trim_sequence_length_distribution','after_trim_sequence_duplication_levels',
        'after_trim_overrepresented_sequences','after_trim_adapter_content'
        ]
    data_merge['fastqc_tests_passed'] = [sum(data_merge.loc[x,fastqc_tests] == 'pass') for x in range(len(data_merge.index))]
    
    return(data_merge)


def qc_evaluate(input_df, n50_min, contig_number_max, contig_number_min, assembly_length_max, assembly_length_min, average_coverage_min, fastqc_tests_passed_min, busco_n_score_min):
    qc_check = [
        input_df['N50'] < n50_min,
        input_df['# contigs (>= 0 bp)'] > contig_number_max,
        input_df['# contigs (>= 0 bp)'] < contig_number_min,
        input_df['Total length'] > assembly_length_max,
        input_df['Total length'] < assembly_length_min,
        input_df['coverage'] < average_coverage_min,
        input_df['fastqc_tests_passed'] < fastqc_tests_passed_min,
        input_df['Percent Complete Nucleotide BUSCOs'] < busco_n_score_min,
    ]
    
    status = ['FAIL'] * len(qc_check)
    
    input_df['QC_EVALUATION'] = np.select(qc_check, status, default='PASS')
    return(input_df)


def main(multiqc_path, auriclass_path, raw_coverage_path, output_path, n50_min, contig_number_max, contig_number_min, assembly_length_max, assembly_length_min, average_coverage_min, fastqc_tests_passed_min, busco_n_score_min):
    auriclass_path,raw_coverage_path = [x+'/' if x[-1] != '/' else x for x in [auriclass_path,raw_coverage_path]]
    # generate auriclass summary
    make_summary_report(input_path = auriclass_path, output_path = auriclass_path, report = 'auriclass', type = 'tsv')
    # generate raw coverage summary
    make_summary_report(input_path = raw_coverage_path, output_path = raw_coverage_path, report = 'raw_coverage', type = 'json')
    mult = multiqc_path
    auri = auriclass_path + 'auriclass_summary.tsv'
    rcov = raw_coverage_path + 'raw_coverage_summary.tsv'
    final = final_qc_summary(multiqc_path = mult, auriclass_path = auri, raw_coverage_path = rcov)
    final_eval = qc_evaluate(
        final,n50_min=n50_min,contig_number_max=contig_number_max,contig_number_min=contig_number_min,assembly_length_max=assembly_length_max,assembly_length_min=assembly_length_min,
        average_coverage_min=average_coverage_min,fastqc_tests_passed_min=fastqc_tests_passed_min,busco_n_score_min=busco_n_score_min
        )
    final_eval.to_csv(output_path, sep='\t',index=False)
    

if __name__ == "__main__":
    mult = snakemake.params["multiqc_dir"]
    auri = snakemake.params["auriclass_dir"]
    rcov = snakemake.params["raw_coverage_dir"]
    outp = snakemake.output[0]
    n50_min = snakemake.config["min_n50"]
    contig_number_min = snakemake.config["min_contigs"]
    contig_number_max = snakemake.config["max_contigs"]
    assembly_length_min = snakemake.config["min_assembly_length"]
    assembly_length_max = snakemake.config["max_assembly_length"]
    average_coverage_min = snakemake.config["min_avg_coverage"]
    fastqc_tests_passed_min = snakemake.config["min_fastqc_tests_passed"]
    busco_n_score_min = snakemake.config["min_busco_nucl_score"]
    main(
        multiqc_path=mult, auriclass_path=auri, raw_coverage_path=rcov, output_path=outp, n50_min=n50_min, contig_number_max=contig_number_max, contig_number_min=contig_number_min, 
        assembly_length_max=assembly_length_max, assembly_length_min=assembly_length_min, average_coverage_min=average_coverage_min, 
        fastqc_tests_passed_min=fastqc_tests_passed_min, busco_n_score_min=busco_n_score_min
        )

