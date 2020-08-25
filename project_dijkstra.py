#! /usr/bin/env python

import numpy as np, pandas as pd, pyarrow as pa, multiprocessing, sys, pickle
from joblib import Parallel, delayed
from functions_project_genomic_coordinates import *

def project(coord, ref, qry, species, pwaln, genome_size, scaling_factor):
     direct_projection = project_genomic_location(ref, qry, coord, 1.0, pwaln, genome_size, scaling_factor)
     # print(direct_projection)
     shortest_path_to_qry, shortest_path, orange = get_shortest_path(ref, qry, coord, species, pwaln, genome_size, scaling_factor, verbose=False)
     # print(shortest_path_to_qry)
     
     # # get anchors for (dijkstra-) projected location. they can be pointing to any species as long as they are:
     # # a) collinear (get_anchors() checks that) and
     # # b) closest to the projected location among the closest anchors from all species
     # proj_chrom, proj_loc = shortest_path_to_qry.loc[qry,'coords'].split(':')
     # qry_anchors_allspecies = pd.concat([get_anchors(pwaln[qry][sp], proj_chrom, int(proj_loc)).assign(sp=sp) for sp in pwaln[qry]], axis=0).rename_axis('direction').reset_index()
     # idx_closest_upstream = qry_anchors_allspecies.loc[qry_anchors_allspecies.direction=='upstream','ref_coord'].idxmax()
     # id
     # x_closest_downstream = qry_anchors_allspecies.loc[qry_anchors_allspecies.direction=='downstream','ref_coord'].idxmin()
     # qry_anchors_narrow = qry_anchors_allspecies.loc[[idx_closest_upstream, idx_closest_downstream],:].set_index('direction')

     columns = pd.MultiIndex.from_tuples([('coords','ref'),('coords','direct'),('coords','dijkstra'),('score','direct'),('score','dijkstra'),
                                          ('ref_anchors','direct'),('ref_anchors','dijkstra'),('qry_anchors','direct'),('qry_anchors','dijkstra'),
                                          ('ref_dist_closest_anchor','direct'),('ref_dist_closest_anchor','dijkstra'),('qry_dist_closest_anchor','direct'),('qry_dist_closest_anchor','dijkstra'),
                                          ('bridging_species','dijkstra')])
     #                                      ('narrow_anchors','coords'),('narrow_anchors','span'),('narrow_anchors','species')])
     # Note: shortest_path_to_qry['ref_anchors'][1] gives the second species in the path with the anchors from the ref species before
     ref_coord_int = int(coord.split(':')[1])
     qry_coord_int_dijkstra = int(shortest_path_to_qry.loc[qry,'coords'].split(':')[1])
     try:
          qry_coord_int_direct = int(direct_projection[1].split(':')[1])
     except IndexError:
          qry_coord_int_direct = np.inf
     values = np.array([coord, direct_projection[1], shortest_path_to_qry.loc[qry,'coords'],
                        direct_projection[0], shortest_path_to_qry.loc[qry,'score'],
                        direct_projection[2], shortest_path_to_qry['ref_anchors'][1], # ref anchors of first species in path
                        direct_projection[3], shortest_path_to_qry['qry_anchors'][-1], # qry anchors of last species in path
                        min([np.inf]+[abs(int(x.split(':')[1]) - ref_coord_int) for x in direct_projection[2]]), min([np.inf]+[abs(int(x.split(':')[1]) - ref_coord_int) for x in shortest_path_to_qry['ref_anchors'][1]]), # dist to closest ref anchor
                        min([np.inf]+[abs(int(x.split(':')[1]) - qry_coord_int_direct) for x in direct_projection[3]]), min([np.inf]+[abs(int(x.split(':')[1]) - qry_coord_int_dijkstra) for x in shortest_path_to_qry['qry_anchors'][-1]]), # dist to closest qry anchor
                        # abs(np.diff([int(x.split(':')[1]) for x in direct_projection[2]])[0]), abs(np.diff([int(x.split(':')[1]) for x in shortest_path_to_qry['ref_anchors'][1]])[0]), # ref anchors width
                        # abs(np.diff([int(x.split(':')[1]) for x in direct_projection[3]])[0]), abs(np.diff([int(x.split(':')[1]) for x in shortest_path_to_qry['qry_anchors'][-1]])[0]), # qry anchors width
                        ','.join(shortest_path_to_qry.index.values[1:-1])])
     #                    tuple(qry_anchors_narrow[['ref_chrom','ref_coord']].astype(str).agg(':'.join, axis=1)), np.diff(qry_anchors_narrow.ref_coord.values)[0], tuple(qry_anchors_narrow.sp)])
     df = pd.DataFrame(values, index=columns).T
     return df

def main():
     if not len(sys.argv) == 7:
          print('Usage: python project_dijkstra.py reference_species query_species coord<chrN:xxxxxx> id half_life_distance path_pwaln_pkl')
          sys.exit(0)
     _, ref, qry, coord, coord_id, half_life_distance, path_pwaln_pkl = sys.argv
     
     # define paths
     data_dir = '/project/wig/tobias/reg_evo/data'
     assembly_dir = data_dir + '/assembly/'
     
     # read pwaln and genome size files
     with open(path_pwaln_pkl, 'rb') as pkl: # you can also pass the ce_pyarrow.pkl path to path_pwaln_pkl, so this works for both pairwise alignments and C(N)Es.
          pwaln_pyarrow = pickle.load(pkl)
     pwaln = {x : {y : pa.deserialize(v) for y,v in pwaln_pyarrow[x].items()} for x in pwaln_pyarrow.keys()} # unserialize msgpack
     species = np.array(list(pwaln.keys())) # The list of species is determined based on the keys of the supplied pairwise aln (pwaln) dict
     genome_size = {s : read_genome_size(assembly_dir + s + '.sizes') for s in species}

     # determine scaling factor based on desired distance_half_life (at which distance to an anchor in the reference species is the score supposed to be 0.5)
     scaling_factor = get_scaling_factor(genome_size[ref], int(half_life_distance))
     
     # project coordinates
     df = project(coord, ref, qry, species, pwaln, genome_size, scaling_factor)
     df.index = [coord_id]
     df.to_csv('tmp/%s.proj.tmp' %coord_id, sep='\t', header=True)
     
     return

if __name__ == '__main__':
     main()
