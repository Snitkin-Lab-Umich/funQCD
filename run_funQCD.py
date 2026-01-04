
import sys
import os
import subprocess
import argparse

# attempt to run the snakemake commands for funQCD in the proper order, based on the provided arguments
# note: remember to add the remaining steps from the setup script to this process, to start directly from raw reads

def run_snakemake_commands(start_step, skip_train):
    command_format = ['snakemake','-s','placeholder','-p','--configfile','config/config.yaml','--profile','profile/']
    workflow_list = ['ASSEMBLY','PREDICTION','ANNOTATION']
    if skip_train:
        workflow_list[1] = 'PREDICTION_SKIP_TRAIN'
    if start_step == 'prediction':
        workflow_list = workflow_list[1:]
    elif start_step == 'annotation':
        workflow_list = workflow_list[2:]
    for workflow in workflow_list:
        workflow_file = f'workflow/funQCD_{workflow}.smk'
        command = command_format.copy()
        command[2] = workflow_file
        print(f"Running snakemake workflow: {workflow_file}")
        result = subprocess.run(command)
        if result.returncode != 0:
            print(f"Error: snakemake workflow {workflow_file} failed with return code {result.returncode}")
            sys.exit(result.returncode)
    print('Finished running all funQCD commands')

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--samples','-s',type=str,
        help='''Provide a path to the samples table for this batch. Note that any existing config/samples.csv file will be overwritten.''',
        default=None
        )
    parser.add_argument(
        '--batch','-b',type=str,
        help='''Provide a name for the run. If resuming a previous run, this must match the previous batch name.''',
        default=None
        )
    parser.add_argument(
        '--reads','-r',type=str,
        help='''Provide a path to the directory containing the raw reads for this run. The file names must match those in the samples table.
        The files also must end in _R1.fastq.gz or _R2.fastq.gz.''',
        default=None
        )
    parser.add_argument(
        '--start',type=str, choices=['assembly','prediction','annotation'],
        help='''Provide the step you want to start the pipeline at. In order, the steps are assembly, prediction, and annotation.
        Downstream steps will be performed as well. Use prediction or annotation to resume a previous run.''',
        default='assembly'
        )
    parser.add_argument(
        '--skip_train','-st',action='store_true',
        help='''Indicate if you want to skip the RNA-seq training step during gene prediction. This will greatly speed up the pipeline,
        but will also significantly reduce the quality of gene predictions.''',
        default=False
        )
    parser.add_argument(
        '--skip_config','-sc',action='store_true',
        help='''Indicate if you want to keep the config file as-is. No changes will be made to config/samples.csv. The samples.csv file
        will also not be modified.''',
        default=False
        )
    args = parser.parse_args()
    if not args.skip_config:
        # overwrite config/samples.csv with args.samples
        if not os.path.isfile(args.samples):
            print(f"Error: could not locate samples table at {args.samples}")
            quit(1)
        if not os.path.isdir(args.reads):
            print(f"Error: could not locate reads directory at {args.reads}")
            quit(1)
        subprocess.run(['cp',args.samples,'config/samples.csv'])
        # replace short_reads and prefix fields in config.yaml
        with open('config/config.yaml','r') as f:
            config_lines = f.readlines()
        with open('config/config.yaml','w') as f:
            for line in config_lines:
                if line.startswith('short_reads:'):
                    line = f'short_reads: {args.reads}\n'
                    if args.reads.endswith('/'):
                        line = f'short_reads: {args.reads[:-1]}\n'
                if line.startswith('prefix:'):
                    line = f'prefix: {args.batch}\n'
                f.write(line)
    # determine which snakemake command to run based on args.start
    run_snakemake_commands(args.start, args.skip_train)
    

if __name__ == "__main__":
    main()


