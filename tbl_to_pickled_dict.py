#! /usr/bin/env python

# This script is used before mass-projecting genomic coordinates, e.g. using project_grb.sh
# It takes a species list from a line-break delimited file (i.e. 1 column), reads all corresponding pairwise alignment .tbl files from the alignment data folder and saves a serialized and pickled dict to the same directory.

import numpy as np, pandas as pd, sys, pyarrow as pa, pickle

def main():
    if not len(sys.argv) == 3:
        print('Usage: ./tbl_to_pickled_dict.py species_file outfile')
        sys.exit(0)
    _, species_file, outfile = sys.argv
    
    # read pairwise aln files for a defined list of species
    # read pwaln for all available pairs is not advisable because the resulting pickled pyarrow file gets way bigger because of the huge teleost pairs and takes 5 times as long to read
    pwaln_dir = '/project/wig/tobias/reg_evo/data/alignment/'
    with open (species_file, 'r') as f:
        species = np.array([x.strip() for x in f.readlines()])
    cols = ['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']
    pwaln = {ref: {qry: pd.read_csv(pwaln_dir + 'tbl/%s.%s.tbl' %(ref,qry), sep='\t', header=None, names=cols) for qry in species if not qry==ref} for ref in species}

    # serialize every data frame in the dict using pyarrow, then pickle the dict to file
    pwaln_pyarrow = {x : {y : pa.serialize(v).to_buffer() for y,v in pwaln[x].items()} for x in pwaln.keys()}
    pickle_out = open(outfile, 'wb')
    pickle.dump(pwaln_pyarrow, pickle_out)
    pickle_out.close()

    return

if __name__=='__main__':
    main()
