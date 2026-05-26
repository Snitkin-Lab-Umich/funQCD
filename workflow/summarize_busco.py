
#import argparse
import pandas as pd
import os

def summarize_busco(busco_output_dir, output_file):
    # take a directory containing BUSCO output directories
    # return a single tab-separated file summarizing all BUSCO outputs
    busco_dict_main = {
        'Input_file': [], 'Dataset': [], 'Complete': [], 'Single': [], 'Duplicated': [], 'Fragmented': [], 'Missing': [], 'n_markers': [],
        'Internal stop codon percent': [], 'Scaffold N50': [], 'Contigs N50': [], 'Percent gaps': [], 'Number of scaffolds': []
        }
    for input_file_name in os.listdir(busco_output_dir):
        busco_path = os.path.join(busco_output_dir, input_file_name)
        busco_dict_single = {
            'Input_file': input_file_name, 'Dataset': None, 'Complete': None, 'Single': None, 'Duplicated': None, 'Fragmented': None, 'Missing': None, 'n_markers': None,
            'Internal stop codon percent': None, 'Scaffold N50': None, 'Contigs N50': None, 'Percent gaps': None, 'Number of scaffolds': None
        }
        if not os.path.isdir(busco_path):
            continue
        # find the .txt summary file in this directory
        summary_file = None
        for file_name in os.listdir(busco_path):
            if file_name.startswith('short_summary') and file_name.endswith('.txt') and input_file_name in file_name:
                summary_file = os.path.join(busco_path, file_name)
                break
        if summary_file is None:
            print(f'No summary file found for {input_file_name}, this file has been skipped')
            continue
        with open(summary_file, 'r') as fh:
            for line in fh:
                if 'The lineage dataset is:' in line:
                    dataset = line.strip().split(': ')[1]
                    dataset = dataset.split(' ')[0]
                    busco_dict_single['Dataset'] = dataset
                elif 'C:' in line and '[S:' in line:
                    line2 = line.strip()
                    if 'E:' not in line2:
                        line2 = line2 + ',E:0.0%'
                    busco_dict_single['Complete'] = line2.split('%[S:')[0].split('C:')[1].strip()
                    busco_dict_single['Single'] = line2.split('%,D:')[0].split('%[S:')[1].strip()
                    busco_dict_single['Duplicated'] = line2.split('%],F:')[0].split('%,D:')[1].strip()
                    busco_dict_single['Fragmented'] = line2.split('%,M:')[0].split('%],F:')[1].strip()
                    busco_dict_single['Missing'] = line2.split('%,n:')[0].split('%,M:')[1].strip()
                    busco_dict_single['n_markers'] = line2.split(',E:')[0].split('%,n:')[1].strip()
                    busco_dict_single['Internal stop codon percent'] = line2.split(',E:')[1].split('%')[0].strip()
                elif 'Scaffold N50' in line:
                    busco_dict_single['Scaffold N50'] = line.strip().split('KB')[0].strip()
                elif 'Contigs N50' in line:
                    busco_dict_single['Contigs N50'] = line.strip().split('KB')[0].strip()
                elif 'Percent gaps' in line:
                    busco_dict_single['Percent gaps'] = line.strip().split('Percent gaps')[0].strip()
                elif 'Number of scaffolds' in line:
                    busco_dict_single['Number of scaffolds'] = line.strip().split('Number of scaffolds')[0].strip()
        for key in busco_dict_single.keys():
            busco_dict_main[key].append(busco_dict_single[key])
    busco_df = pd.DataFrame(busco_dict_main)
    busco_df.to_csv(output_file, sep='\t', index=False)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--nucl_input','-ni',type=str,
        help='''Provide the path to the nucleotide BUSCO outputs to summarize.''',
        required=True
        )
    parser.add_argument(
        '--prot_input','-pi',type=str,
        help='''Provide the path to the protein BUSCO outputs to summarize.''',
        required=True
        )
    parser.add_argument(
        '--nucl_output','-no',type=str,
        help='''Provide a path to a directory for the nucleotide batch summary file.''',
        required=True
        )
    parser.add_argument(
        '--prot_output','-po',type=str,
        help='''Provide a path to a directory for the protein batch summary file.''',
        required=True
        )
    args = parser.parse_args()
    for p in [args.nucl_input, args.prot_input]:
        if not os.path.isdir(p):
            print(f'Could not locate directory at {p}')
            quit(1)
    for p in [args.nucl_output, args.prot_output]:
        if not os.path.isdir(p):
            os.makedirs(p)
    summarize_busco(args.nucl_input, os.path.join(args.nucl_output, 'busco_nucl_batch_summary.txt'))
    summarize_busco(args.prot_input, os.path.join(args.prot_output, 'busco_prot_batch_summary.txt'))

if __name__ == '__main__':
    #main()
    nucl_input = snakemake.params['busco_output_n_dir']
    prot_input = snakemake.params['busco_output_p_dir']
    nucl_output = snakemake.params['busco_summary_n_dir']
    prot_output = snakemake.params['busco_summary_p_dir']
    summarize_busco(nucl_input, os.path.join(nucl_output, 'busco_nucl_batch_summary.txt'))
    summarize_busco(prot_input, os.path.join(prot_output, 'busco_prot_batch_summary.txt'))

