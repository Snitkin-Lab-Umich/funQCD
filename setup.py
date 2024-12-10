
import sys
import os
import subprocess

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
        if '.fastq.gz' in f:
            sample_id = '_'.join(f.split('_')[:-1])
            # this should always return the full text to the left of '_R1.fastq.gz' or '_R2.fastq.gz', even if it contains '_' characters
            sample_id_set.add(sample_id)
    with open('config/samples.csv','w') as fhout:
        _ = fhout.write('sample_id\n')
        for sid in sorted(list(sample_id_set)):
            _ = fhout.write(sid + '\n')

if __name__ == "__main__":
    copy_gm_key()
    make_samples_csv(sys.argv[1])


