#! /usr/bin/env python

### This script takes the states from the enhancer and promoter foreground modules, keeps the emisP, but shuffles / unitizes the transPs
### such that the foreground is only different to the background in terms of transitions, i.e. emphasizing on the molecular structure.

import numpy as np, pandas as pd, sys

def main():
    if not len(sys.argv) == 4:
        print('Usage: python remove_states_from_ehmm_model.py <enhancer_model.txt> <promoter_model.txt> <output_model.txt>')
        sys.exit(0)

    _, enhancer_model, promoter_model, output_model = sys.argv

    # MAYBE USE EHMM PACKAGE IN R FOR THIS SCRIPT?
    
    with open(input_model, 'r') as f:
        lines = f.readlines()
    lines = np.array([l.strip() for l in lines])
    header = lines[:8]
    nstates = int(lines[1])
    states = np.array(header[3].split(','))
    emisP = lines[9:(9+nstates)]
    transP = pd.DataFrame([x.split('\t') for x in lines[(9+nstates+1):(9+2*nstates+1)]], dtype=np.float) # convert to pandas DataFrame for manipulation
    transP = transP.divide(transP.sum(axis=1), axis=0) # normalize back to row-sums of 1
    transP = np.array(transP.apply(lambda x : '\t'.join([str(y) for y in x]), axis=1)) # convert back into string in array format
    initP = [np.float(x) for x in lines[-nstates:]]
    
    outlines = np.concatenate((header, np.array(['emisP']), emisP, np.array(['transP']), transP, np.array(['initP']), initP))
    with open(output_model, 'w') as f:
        for l in outlines:
            f.write(l + '\n')

    return
            
if __name__ == '__main__':
    main()
