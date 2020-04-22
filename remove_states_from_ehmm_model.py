#! /usr/bin/env python

import numpy as np, pandas as pd, sys

def main():
    if not len(sys.argv) == 4:
        print('Usage: python remove_states_from_ehmm_model.py <input_model.txt> <states,to,remove> <output_model.txt>')
        sys.exit(0)

    _, input_model, rmstates, output_model = sys.argv

    rmstates = np.array(rmstates.split(','))
    with open(input_model, 'r') as f:
        lines = f.readlines()
    lines = np.array([l.strip() for l in lines])
    header = lines[:8]
    nstates = int(lines[1])
    header[1] = str(nstates - len(rmstates)) # number of states
    states = np.array(header[3].split(','))
    idx_rmstates = np.intersect1d(states, rmstates, return_indices=True)[1]
    idx_keepstates = np.delete(np.arange(nstates), idx_rmstates)
    keepstates = np.delete(states, idx_rmstates)
    header[3] = ','.join(np.delete(np.array(header[3].split(',')), idx_rmstates))
    header[5] = ','.join(np.delete(np.array(header[5].split(',')), idx_rmstates))
    emisP = lines[9:(9+nstates)]
    emisP = np.delete(emisP, idx_rmstates) # remove states from emisP
    transP = pd.DataFrame([x.split('\t') for x in lines[(9+nstates+1):(9+2*nstates+1)]], dtype=np.float) # convert to pandas DataFrame for manipulation
    transP = transP.iloc[idx_keepstates, idx_keepstates] # remove states
    transP = transP.divide(transP.sum(axis=1), axis=0) # normalize back to row-sums of 1
    transP = np.array(transP.apply(lambda x : '\t'.join([str(y) for y in x]), axis=1)) # convert back into string in array format
    initP = [np.float(x) for x in lines[-nstates:]]
    initP = np.delete(initP, idx_rmstates) # remove states
    initP /= initP.sum() # normalize back to sum of 1
    initP = [str(x) for x in initP]

    outlines = np.concatenate((header, np.array(['emisP']), emisP, np.array(['transP']), transP, np.array(['initP']), initP))
    with open(output_model, 'w') as f:
        for l in outlines:
            f.write(l + '\n')

    return
            
if __name__ == '__main__':
    main()
