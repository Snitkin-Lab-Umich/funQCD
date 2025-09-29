import pandas as pd
import sys
import os
import subprocess
import argparse
from collections import Counter

def find_fail_point(sample_name,funqcd_dir):
    # search through all of the funQCD outputs to find the earliest step that failed
    # failed steps should have placeholder output files with zero lines
    # return the name of the step that failed
    # step 0: check if the sample actually exists
    sample_dir = os.path.join(funqcd_dir,'spades',sample_name)
    if not os.path.exists(sample_dir):
        return('sample_not_found')
    # step 1: assembly check
    assembly_check_file = os.path.join(funqcd_dir,'quast',sample_name,'assembly_check.tsv')
    # this file will say "ASSEMBLY QC FAILED" if it did not pass this step
    with open(assembly_check_file,'r') as fh:
        lines = fh.readlines()
    if lines != ['ASSEMBLY QC PASSED']:
        return('failed_assembly_qc')
    # step 2: funannotate train
    # if this step failed, then all four of the output files will have zero lines
    # however, the contents of this directory are unfortunately deleted as part of the cleanup step
    # this means I cannot tell which samples failed this step
    #train_file = os.path.join(funqcd_dir,'funannotate',sample_name,'training/funannotate_train.coordSorted.bam')
    # step 3: funannotate predict
    predict_file = os.path.join(funqcd_dir,'funannotate',sample_name,f'predict_results/{sample_name}.proteins.fa')
    if not os.path.exists(predict_file):
        return('failed_funannotate_predict')
    if count_lines(predict_file) == 0:
        return('failed_funannotate_predict')
    # step 4: funannotate update
    update_file = os.path.join(funqcd_dir,'funannotate',sample_name,f'update_results/{sample_name}.proteins.gff3')
    if not os.path.exists(update_file):
        return('failed_funannotate_update')
    if count_lines(update_file) == 0:
        return('failed_funannotate_update')
    # step 5: funannotate annotate
    # if isolates reached this step, then they will only fail due to the qc results
    return('no_failure_detected')


def count_lines(fname):
    # count the number of lines in a file
    with open(fname,'r') as fh:
        lines = fh.readlines()
    return(len(lines))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide the path to a funQCD output directory.''',
        default=None
        )
    parser.add_argument(
        '--qc-file','-qc',type=str,
        help='''Provide the path to a qc summary file that contains the names of failed samples.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a path for output table.''',
        default=None
        )    
    args = parser.parse_args()
    qc = pd.read_csv(args.qc_file,sep='\t')
    failed_samples = qc[qc['QC_EVALUATION'] == 'FAIL']['Sample'].tolist()
    print(f'Found {len(failed_samples)} failed samples in the QC file.')
    # get the status for each failed sample
    failed_sample_status = [find_fail_point(x,args.input) for x in failed_samples]
    print(f'Identified failure points for {len(failed_sample_status)} samples.')
    print(Counter(failed_sample_status))
    # convert to pandas dataframe and write to file
    df = pd.DataFrame({'Sample':failed_samples,'Failure_Point':failed_sample_status})
    # remove any rows where the sample was not found
    df = df[df['Failure_Point'] != 'sample_not_found']
    df.to_csv(args.output,sep='\t',index=False)


if __name__ == '__main__':
    main()


