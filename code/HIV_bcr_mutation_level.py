#!/usr/bin/env python
import pandas as pd
import numpy as np
from scipy.stats import skew
import argparse
import sys

def percentie99(x):
    return np.percentile(x, 99)

def main():
    parser = argparse.ArgumentParser(description='Get per subject V segment mutation levels',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('data', metavar='read_data', help='file with read data and V segment mutation levels')
    
    args = parser.parse_args()
    column_types = {'subject': str,'sample': str, 'isotype': str, 'lineage': str, 'v_j_in_frame': str, 'has_stop_codon': str, 'cdr3_sequence':str, 'v_segment':str, 'j_segment':str , 'mutation_level': float}
    read_data = pd.read_csv(args.data, sep = ',', usecols=column_types.keys(), dtype=column_types)

    read_data = read_data[read_data['v_j_in_frame'].notna()]
    read_data = read_data[read_data['has_stop_codon'].notna()]
    read_data = read_data[read_data['isotype'].notna()]

    status = {'True': True, 'False': False}
    read_data['has_stop_codon'] = read_data['has_stop_codon'].map(status)
    read_data['v_j_in_frame'] = read_data['v_j_in_frame'].map(status)

    read_data = read_data[(read_data['v_j_in_frame'] == True) & (read_data['has_stop_codon'] == False)]

    del read_data['v_j_in_frame']
    del read_data['has_stop_codon']

    clone_data   = read_data.groupby(['subject', 'sample', 'isotype', 'lineage'], dropna=False).\
                aggregate({'mutation_level': [np.size, np.mean, np.median ]})
                
    subject_data = clone_data.groupby(['subject', 'sample', 'isotype'], dropna=False).\
                aggregate({('mutation_level', 'size'):   [np.size, np.mean, np.median, np.sum],
                           ('mutation_level', 'mean'):   [np.size, np.mean, np.median, skew],
                           ('mutation_level', 'median'): [np.size, np.mean, np.median, percentie99, skew]})
                
                
    subject_data.columns = ['.'.join(c) for c in subject_data.columns] 
    subject_data.to_csv(sys.stdout, float_format='%.5f')
if __name__ == '__main__':
    sys.exit(main())
