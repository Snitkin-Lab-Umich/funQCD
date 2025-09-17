
import os
import argparse

def remake_predict_pass_file(dir_path):
    funannotate_dir = dir_path + 'funannotate/'
    passlist = []
    if not os.path.isdir(funannotate_dir):
        print(f'Could not locate funannotate directory at {funannotate_dir}')
        quit(1)
    for sample_name in os.listdir(funannotate_dir):
        fpath = funannotate_dir + sample_name
        if os.path.isdir(fpath):
            checkpath = fpath + '/update_results/annotation_check.txt'
            if os.path.isfile(checkpath):
                with open(checkpath,'r') as fh:
                    d = fh.readlines()
                    if d == ['FUNANNOTATE PREDICTION PASSED']:
                        passlist.append(sample_name)        
    with open(dir_path + 'predict_pass_samples.csv','w') as fhout:
        _ = fhout.write('sample_id\n')
        for s in passlist:
            _ = fhout.write(s + '\n')


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--dir','-d',type=str,
        help='''Provide a path to the results directory for this batch. The predict_pass_samples.csv file will be overwritten.''',
        default=None,required=True
        )
    args = parser.parse_args()
    if not os.path.isdir(args.dir):
        print(f'Could not locate directory at {args.dir}')
        quit(1)
    # ensure the provided path does not end with '/'
    if not args.dir.endswith('/'):
        args.dir = args.dir + '/'
    remake_predict_pass_file(args.dir)

if __name__ == "__main__":
    main()


