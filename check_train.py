
import sys
import os
import subprocess
import argparse

def make_new_csv(input_path,output_path,funannotate_dir):
    with open(input_path,'r') as fh_in, open(output_path,'w') as fh_out:
        _ = next(fh_in)
        fh_out.write('sample_id\n')
        for line in fh_in:
            sample = line.strip()
            checkpaths = []
            checkpaths.append(f'{funannotate_dir}{sample}/training/left.fq.gz')
            checkpaths.append(f'{funannotate_dir}{sample}/training/right.fq.gz')
            checkpaths.append(f'{funannotate_dir}{sample}/training/trinity_gg')
            checkpaths.append(f'{funannotate_dir}{sample}/training/trimmomatic')
            if all([os.path.exists(x) for x in checkpaths]):
                _ = fh_out.write(line)

if __name__ == "__main__":
    inp = 'results/2025-05-06_ncbiPathogen_t30_v1/assembly_pass_samples.csv'
    outp = 'results/2025-05-06_ncbiPathogen_t30_v1/assembly_pass_samples_PART1.csv'
    funannotate_dir = 'results/2025-05-06_ncbiPathogen_t30_v1/funannotate/'
    make_new_csv(inp,outp,funannotate_dir)


