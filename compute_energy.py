#! /usr/bin/env python

# This script computes energies for all sequences or a given subset of a fasta file using a given position count matrix (PCM)
# The PCM is converted into a PWM using the overall nucleotide frequencies of all sequences in the fasta file as background frequencies.
# The computed energies are pd.Series objects containing window-energies for windows of the size of the PCM length and include energies for binding on either strand (incl. reverse complement).

import numpy as np, pandas as pd, pybedtools as pb, multiprocessing as mp, sys, os.path, feather
sys.path.insert(0, '/project/wig/tobias/reg_evo/jupyter/')
from functions_PWM import *
from functions import *
from Bio import SeqIO

def compute_energy(seq):
    e = sum([pwm.loc[tup] for tup in zip(range(len(seq)),list(seq))])
    e_rc = sum([pwm_RC.loc[tup] for tup in zip(range(len(seq)),list(seq))])
    return(np.log(np.exp(e)+np.exp(e_rc)))

def main():
    if len(sys.argv) < 4:
        print 'Usage: compute_energy.py <fasta_file> <PCM_file> <factor_name> <record_ids (optional, comma-separated)>'
        sys.exit(0)
        
    _, fasta_file, pcm_file, factor_name = sys.argv[:4]
    if len(sys.argv) == 5:
        record_ids = str.split(sys.argv[4], ',')
    output_file = factor_name + '.feather'

    print 'compute background frequencies'
    bg_counts = np.zeros(4)
    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            bg_counts += np.array([record.seq.count('A'), record.seq.count('C'), record.seq.count('G'), record.seq.count('T')])
    bg_freqs = bg_counts / bg_counts.sum()
    bg_freqs

    print 'compute PWM'
    global pwm, pwm_RC
    pwm = get_PWM(pcm_file, bg_freqs=bg_freqs, pseudo_count=1)
    pwm_RC = pwm.loc[::-1,::-1] # reverse complement
    pwm_RC.rename(index=dict(zip(pwm_RC.index,pwm.index)),
                  columns = dict(zip(pwm_RC.columns,pwm.columns)),
                  inplace=True)
    pwm['N'], pwm_RC['N'] = 0, 0

    # runtime ~30-45s per GRB --> ~7h for ~600 GRBs
    print 'compute energies'
    l = pwm.shape[0]
    with open('stdout_%s.txt' %factor_name, 'w') as o:
        o.write('record.id\n')
    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            try:
                if not record.id in record_ids:
                    continue
            except NameError:
                pass
            with open('stdout_%s.txt' %factor_name, 'a') as o:
                o.write(record.id+'\n')
                output_file_i = '%s_%s.feather' %(factor_name, record.id)
            if not os.path.exists(output_file_i):
                s = record.seq
                p = mp.Pool(80)
                e = p.map(compute_energy, [s[i:i+l] for i in np.arange(0,len(s)-l+1,1)])
                p.close()
                p.join()
                df = pd.DataFrame(e, columns=[record.id])
                feather.write_dataframe(df, output_file_i)
                
    # write results to file using `feather` (compressed data frame readable in python and R)
    # only if energies were calculated for all sequences in fasta file, otherwise leave single feather files
    if 'record_ids' not in locals():
        print 'write results to %s' %output_file
        files = [f for f in os.listdir('.') if f.startswith(factor_name) & f.endswith('feather') & (f != output_file)]
        df_list = [feather.read_dataframe(f) for f in files]
        df = pd.concat(df_list, axis=1)
        df = df[sorted(df.columns, key=int)]
        feather.write_dataframe(df, output_file) # write combined file
        _ = [os.remove(f) for f in files] # delete single files
    print 'Done'
    
if __name__ == '__main__':
    main()
