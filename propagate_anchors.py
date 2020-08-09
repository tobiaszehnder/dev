#! /usr/bin/env python

import numpy as np, pandas as pd, pyarrow as pa
from functions_project_genomic_coordinates import *

# this script takes a region in the format chrN:xxxxxx and returns a table with all used anchors that lead to a minimized anchor span in the target species.
# run this script in parallel for multiple regions stored in a bed-file using propagate_anchors_wrapper.sh

def main():
     if not len(sys.argv) == 7: 
          print('Usage: python propagate_anchors.py reference_species target_species coordinate<chrN:xxxxxx> id path_pwaln_pkl output_directory')
          sys.exit(0)
     _, reference, target, coord, coord_id, path_pwaln_pkl, outdir = sys.argv
 
     # read pwaln and genome size files
     with open(path_pwaln_pkl, 'rb') as pkl: # you can also pass the ce_pyarrow.pkl path to path_pwaln_pkl, so this works for both pairwise alignments and C(N)Es.
          pwaln_pyarrow = pickle.load(pkl)
     pwaln = {x : {y : pa.deserialize(v) for y,v in pwaln_pyarrow[x].items()} for x in pwaln_pyarrow.keys()} # unserialize msgpack
      
     # propagate anchors
     ob, path = propagate_anchors(reference, target, coord, coord_id, pwaln, silent=True)
     if not path.shape == (0,0): # only write to file if both direct up- and downstream anchors are found from refernce to target. if not, an empty path will have been returned. end script without writing to file but print a message.
          path.to_csv('%s/%s.aprop' %(outdir,coord_id), sep='\t', header=True)
     else:
          print('%s %s: Not both up- and downstream anchors were found between the reference and target species. The job was aborted.' %(coord_id,coord))
     
     return

if __name__ == '__main__':
     main()
