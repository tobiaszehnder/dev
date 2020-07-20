#! /usr/bin/env python

import numpy as np, pandas as pd, pyarrow as pa, multiprocessing, sys, pickle
from joblib import Parallel, delayed
from functions_project_genomic_coordinates_pwaln import *

def project(coord, ref, qry, species, pwaln, genome_size):
     direct_projection = project_genomic_location(ref, qry, coord, 1.0, pwaln, genome_size)
     shortest_path_to_qry, shortest_path, orange = get_shortest_path(ref, qry, coord, species, pwaln, genome_size, verbose=False)
     columns = pd.MultiIndex.from_tuples([('coords','ref'),('coords','direct'),('coords','dijkstra'),('score','direct'),('score','dijkstra'),
                                          ('ref_anchors','direct'),('ref_anchors','dijkstra'),('qry_anchors','direct'),('qry_anchors','dijkstra'),
                                          ('ref_anchorspan','direct'),('ref_anchorspan','dijkstra'),('qry_anchorspan','direct'),('qry_anchorspan','dijkstra'),
                                          ('bridging_species','dijkstra')]) 
     # Note: shortest_path_to_qry['ref_anchors'][1] gives the second species in the path with the anchors from the ref species before
     values = np.array([coord, direct_projection[1], shortest_path_to_qry.loc[qry,'coords'],
                        direct_projection[0], shortest_path_to_qry.loc[qry,'score'],
                        direct_projection[2], shortest_path_to_qry['ref_anchors'][1], # ref anchors of first species in path
                        direct_projection[3], shortest_path_to_qry['qry_anchors'][-1], # qry anchors of last species in path
                        abs(np.diff([int(x.split(':')[1]) for x in direct_projection[2]])[0]), abs(np.diff([int(x.split(':')[1]) for x in shortest_path_to_qry['ref_anchors'][1]])[0]), # ref anchors width
                        abs(np.diff([int(x.split(':')[1]) for x in direct_projection[3]])[0]), abs(np.diff([int(x.split(':')[1]) for x in shortest_path_to_qry['qry_anchors'][-1]])[0]), # qry anchors width
                        ','.join(shortest_path_to_qry.index.values[1:-1])])
     df = pd.DataFrame(values, index=columns).T
     return df

def main(ref, qry, coord, coord_id):
     # define paths
     data_dir = '/project/wig/tobias/reg_evo/data'
     pwaln_dir = data_dir + '/alignment/'
     assembly_dir = data_dir + '/assembly/'
     
     # read pwaln and genome size files
     with open(pwaln_dir + 'pwaln_pyarrow.pkl', 'rb') as pkl:
          pwaln_pyarrow = pickle.load(pkl)
     pwaln = {x : {y : pa.deserialize(v) for y,v in pwaln_pyarrow[x].items()} for x in pwaln_pyarrow.keys()} # unserialize msgpack
     species = np.array(list(pwaln.keys())) # The list of species is determined based on the keys of the supplied pairwise aln (pwaln) dict
     genome_size = {s : read_genome_size(assembly_dir + s + '.sizes') for s in species}
     
     # project coordinates
     df = project(coord, ref, qry, species, pwaln, genome_size)
     df.index = [coord_id]
     df.to_csv('tmp/%s.proj.tmp' %coord_id, sep='\t', header=True)
     
     return

if __name__ == '__main__':
     main()
