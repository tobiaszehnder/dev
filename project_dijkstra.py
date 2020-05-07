#! /usr/bin/env python

import numpy as np, pandas as pd, pyarrow as pa, multiprocessing, sys, pickle
from joblib import Parallel, delayed
if 'functions_project_genomic_coordinates' in sys.modules:
     del sys.modules['functions_project_genomic_coordinates'] # deletes previously imported module so that potential changes will be loaded
from functions_project_genomic_coordinates import *

if not len(sys.argv) == 5:
     print('Usage: python project_dijkstra.py reference_species query_species coord<chrN:xxxxxx> id')
     sys.exit(0)

_, ref, qry, coord, coord_id = sys.argv

def project(coord, ref, qry, species, ce, genome_size):
     direct_projection = project_genomic_location(ref, qry, coord, 1.0, ce, genome_size)
     shortest_path_to_qry, shortest_path, orange = get_shortest_path(ref, qry, coord, species, ce, genome_size, verbose=False)
     columns = pd.MultiIndex.from_tuples([('coords','ref'),('coords','direct'),('coords','dijkstra'),('score','direct'),('score','dijkstra'),
                                                    ('ref_anchors','direct'),('ref_anchors','dijkstra'),('bridging_species','dijkstra')]) 
     # Note: shortest_path_to_qry['qry_anchors'][1] gives the second species in the path and the anchors from the ref species before
     values = np.array([coord, direct_projection[1], shortest_path_to_qry.loc[qry,'coords'],
                        direct_projection[0], shortest_path_to_qry.loc[qry,'score'],
                        direct_projection[2], shortest_path_to_qry['ref_anchors'][1],
                        shortest_path_to_qry.index.values[1:-1]])
     df = pd.DataFrame(values, index=columns).T
     return df

def main():
     # define paths and variables
     species = species = np.array(['hg38','mm10','galGal5','xenTro9','danRer10','lepOcu1','calMil1','oryLat2','gasAcu1','tetNig2','fr3','latCha1','cteIde1','cypCar1','carAur01'])
     data_dir = '/project/wig/tobias/reg_evo/data'
     cne_dir = data_dir + '/CNEs/CNEr/'
     assembly_dir = data_dir + '/assembly/'
     ce_files = {s1 : {s2 : get_ce_path(cne_dir, s1, s2) for s2 in species[species != s1]} for s1 in species}
     genome_size = {s : read_genome_size(assembly_dir + s + '.sizes') for s in species}

     # read CE files (~30 sec)
     pkl = open('/project/wig/tobias/reg_evo/data/CNEs/CNEr/ce_pyarrow.pkl', 'rb')
     ce_pyarrow = pickle.load(pkl)
     ce = {x : {y : pa.deserialize(v) for y,v in ce_pyarrow[x].items()} for x in ce_pyarrow.keys()} # unserialize msgpack
     
     df = project(coord, ref, qry, species, ce, genome_size)
     df.index = [coord_id]
     df.to_csv('tmp/%s.proj.tmp' %coord_id, sep='\t', header=True)

     return

if __name__ == '__main__':
     main()
