
import sys
import os
import subprocess
import argparse

def copy_gm_key():
    if not os.path.exists('./.gm_key'):
        homedirkey = os.path.expanduser('~') + '/.gm_key'
        if not os.path.exists(homedirkey):
            print('GeneMark license not found in home directory. Please download .gm_key from https://genemark.bme.gatech.edu/license_download.cgi')
        else:
            subprocess.call(['cp',homedirkey,'./.gm_key'])
            print('Copied GeneMark license to working directory')
    else:
        print('GeneMark license found in funQCD directory')

def make_samples_csv(path):
    flist = os.listdir(path)
    sample_id_set = set()
    if os.path.exists('config/samples.csv'):
        print('Overwriting config/samples.csv with new version')
    for f in flist:
        if '_R1.fastq.gz' in f or '_R2.fastq.gz' in f:
            sample_id = '_R'.join(f.split('_R')[:-1])
            # this should always return the full text to the left of '_R1.fastq.gz' or '_R2.fastq.gz', even if it contains '_R' somewhere other than the end
            sample_id_set.add(sample_id)
    with open('config/samples.csv','w') as fhout:
        _ = fhout.write('sample_id\n')
        for sid in sorted(list(sample_id_set)):
            _ = fhout.write(sid + '\n')

def add_path_to_config(path,prefix,config_file = 'config/config.yaml'):
    lines = []
    with open(config_file,'r') as fh:
        for line in fh:
            if line.startswith('short_reads:'):
                line = f'short_reads: {path}\n'
            if line.startswith('prefix:'):
                line = f'prefix: {prefix}\n'
            lines.append(line)
    with open(config_file,'w') as fh_out:
        for line in lines:
            _ = fh_out.write(line)

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--path','-p',type=str,
        help='''Provide a path to the raw data to process with the pipeline. This path will be added to config/config.yaml, and 
        a config/samples.csv file will be created from the directory contents. Files to process must end in _R1.fastq.gz or _R2.fastq.gz''',
        default=None,required=True
        )
    parser.add_argument(
        '--prefix','-f',type=str,
        help='''Provide a name for this batch. This will be added to config/config.yaml.''',
        default=None,required=True
        )
    args = parser.parse_args()
    if not os.path.isdir(args.path):
        print(f'Could not locate directory at {args.path}')
        quit(1)
    # ensure the provided path does not end with '/'
    if args.path.endswith('/'):
        args.path = args.path[:-1]
    copy_gm_key()
    make_samples_csv(args.path)
    add_path_to_config(args.path,args.prefix)

if __name__ == "__main__":
    main()


