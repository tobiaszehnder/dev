#! /usr/bin/env python

'''
This script computes sequence similarities based on the overlap of k-mers (normalized with the expected number of overlapping k-mers)

-----------------

There are two functions to be chosen from:

1: pairwise
This function computes pairwise sequence similarities between pairs of sequences from a given bed-file.

2: find_homolog
This function takes a single sequence in species X, projects it genomic location to species Y based on interpolated CNE coordinates
and returns sequence similarities in a 500 Kbp window around this projected center in a bigwig file.

-----------------

Count k-mers of a reference region and its projected target (i.e. 500 Kbp region centered on the projected target coordinate)
Combine rev / comp / rev-comp in kmer-counts

Introduce a dict with IDs such that d['AAAT'] = d['TAAA'] = d['TTTA'] = d['ATTT']
Elongate every k-mer with '#' (and add this character to the alphabet) to prevent k-mers with different k to have the same ID

ID_kmer# = sum(a(s_i) len(a)^i)
with alphabet a={A:0,C:1,G:2,T:3,#:4}
such that e.g. IDTATC#=3x50+0x51+3x52+1x53+4x54=3+0+75+125+2500=2703
For every possible k-mer, assign the minimum of the IDs of the k-mer, its reverse, complement and reverse complement to all of those k-mers.
'''

import numpy as np, pandas as pd, pybedtools as pb, scipy, sys, os.path, itertools, subprocess, json, datetime
from Bio import SeqIO
from multiprocessing.pool import ThreadPool
sys.path.insert(1, '/project/wig/tobias/reg_evo/jupyter/')
if 'functions' in sys.modules:
  del sys.modules['functions'] # deletes previously imported module so that potential changes will be loaded
from functions import *

def get_BSgenome_name(assembly):
  bsgenome_dict = {'hg19': 'BSgenome.Hsapiens.UCSC.hg19', 'mm10' : 'BSgenome.Mmusculus.UCSC.mm10'}
  if not assembly in bsgenome_dict.keys():
    print 'Error: specified assembly "%s" is not yet supported in this script. Feel free to add.' %assembly
    sys.exit(1)
  return bsgenome_dict[assembly]

def approximate_expected_overlap(kmer_ids, ks):
  os_approx = {}
  ns = [float(len(set([v for k,v in kmer_ids.items() if len(k) == K]))) for K in ks]
  for i,K in enumerate(ks):
    n = ns[i]
    k = 500.
    expected_m = n - n*(1-1/n)**k
    expected_o = expected_m**2 / n
    os_approx[K] = expected_o
  return os_approx

def pairwise(regions_bed, assembly, ks):
  print 'Fetch sequences'
  seqs = {}
  subprocess.check_call(['Rscript', 'get_sequence.R', regions_bed, get_BSgenome_name(assembly)], shell=False)
  with open(regions_bed.replace('bed', 'fasta'), 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
      seqs[record.id] = str(record.seq)

  # create dict with IDs and initiate id_counter dict (~20 sec)
  print 'Assign k-mer IDs'
  possible_kmers = [''.join(p) for k in ks for p in itertools.product(list('ACGT'), repeat=k)]
  kmer_ids = {kmer : get_kmer_id(kmer) for kmer in possible_kmers}

  # count kmers
  print 'Count k-mers'
  kmer_counts_dict = {seq_id : {k : count_kmers(seq, kmer_ids, ks=[k]) for k in ks} for seq_id, seq in seqs.items()}

  # approximate E[o] = the expected number of unique k-mer ids overlapping between two random sequence windows of size 500 bp (~25 sec)
  print 'Approximate expected number of shared k-mers'
  os_approx = approximate_expected_overlap(kmer_ids, ks)
  
  # compute similarity scores for every K (~6 min for 1 Mbp)
  print 'Compute pairwise similarities'
  m_dict = {k : np.full([len(seqs), len(seqs)], np.nan) for k in ks}
  seq_ids = sorted(seqs.keys())
  for k in ks:
    for i,id_i in enumerate(seq_ids):
      for j in range(i,len(seqs)):
        id_j = seq_ids[j]
        m_dict[k][i,j] = compute_similarity(kmer_counts_dict[id_i][k], kmer_counts_dict[id_j][k], os_approx[k])

  # write table to file
  for k in ks:
    outfile = regions_bed.replace('.bed', '_%s-mer_pairwise_similarities.table' %k).replace('bed/','output/')
    pd.DataFrame(m_dict[k]).to_csv(outfile, sep='\t', index=False, header=False)
  return

def find_homolog(reference_bed, reference_assembly, target_assembly, cnefile, ks, array_size=5e5, window_size=500):
  bsgenome_dict = {'hg19': 'BSgenome.Hsapiens.UCSC.hg19', 'mm10' : 'BSgenome.Mmusculus.UCSC.mm10'}
  if not ((reference_assembly in bsgenome_dict.keys()) & (target_assembly in bsgenome_dict)):
    print 'Error: specified assemblies are not yet supported in this script. Feel free to add.'
    sys.exit(1)
  
  # read reference bed and fetch sequence from genome fasta by running the get_sequence.R (writes it to bedfile.fasta)
  print 'Fetch reference sequence'
  subprocess.check_call(['Rscript', 'get_sequence.R', reference_bed, bsgenome_dict[reference_assembly]], shell=False)
  with open(reference_bed.replace('bed', 'fasta'), 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
      ref_id = record.id
      ref_seq = str(record.seq)
      break # ignore any subsequent sequences for now

  # read CNE file and project center of reference sequence to the target genome
  # extend location of target genome to 1 MB and write to target.bed. Use it with `get_sequence.R` to fetch the sequence.
  print 'Project to target genome'
  cnes = pd.read_csv(cnefile, sep='\t', header=None).sort_values([0,1,2])
  ref_region = pd.read_csv(reference_bed, sep='\t', header=None, usecols=range(4), index_col=3, names=['chrom','start','end','name']).iloc[0]
  ref_center = (ref_region['start']+ref_region['end'])/2
  cne_upstream = cnes.loc[(cnes[0] == ref_region['chrom']) & ((cnes[1]+cnes[2])/2 < (ref_center)),].tail(1)
  cne_downstream = cnes.loc[(cnes[0] == ref_region['chrom']) & ((cnes[1]+cnes[2])/2 > ref_center),].head(1)
  cne_neighbours_center = [float((cne_upstream[1] + cne_upstream[2]) / 2), float((cne_downstream[1] + cne_downstream[2]) / 2)]
  ref_relative_position = (ref_center - cne_neighbours_center[0]) / (cne_neighbours_center[1] - cne_neighbours_center[0])
  cnes_target_center = sorted([float((cne_upstream[4] + cne_upstream[5])) / 2, float((cne_downstream[4] + cne_downstream[5]) / 2)]) # use sort as we don't know target orientation
  ref_target_position = int(cnes_target_center[0] + (cnes_target_center[1] - cnes_target_center[0]) * ref_relative_position)
  target_bed = reference_bed.replace('.bed', '_target.bed')
  pd.DataFrame(np.array([[cne_upstream.values[0,3], int(ref_target_position-array_size/2), int(ref_target_position+array_size/2), ref_region.name]])).to_csv(target_bed, sep='\t', header=None, index=False)

  # write fasta and fetch sequence using get_sequence.R
  print 'Fetch target sequence'
  subprocess.check_call(['Rscript', 'get_sequence.R', target_bed, bsgenome_dict[target_assembly]], shell=False)
  with open(target_bed.replace('bed', 'fasta'), 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
      target_seq = str(record.seq)
  
  # create dict with IDs and initiate id_counter dict (~20 sec)
  print 'Assign k-mer IDs'
  possible_kmers = [''.join(p) for k in ks for p in itertools.product(list('ACGT'), repeat=k)]
  kmer_ids = {kmer : get_kmer_id(kmer) for kmer in possible_kmers}

  # compute kmers for reference
  print 'Count reference k-mers'
  kmer_counts_ref_dict = {k : count_kmers(ref_seq, kmer_ids, ks=[k]) for k in ks}

  # compute kmers for target (~8 min for 1 Mbp)
  print 'Count target k-mers'
  kmer_counts_target_dict = {k : count_kmers_wrapper(target_seq, kmer_ids, ks=[k], window_size=window_size) for k in ks}

  # read target bed and create data frame with 1 bp resolution
  target_regions = pd.read_csv(target_bed, sep='\t', header=None).loc[0]
  start = np.arange(target_regions[1], target_regions[1] + len(kmer_counts_target_dict.values()[0])) # the bigwig's last entry starts at 500 bp before the target region's end coordinate.
  end = start + 1
  chrs = np.repeat(target_regions[0], len(start))   

  # approximate E[o] = the expected number of unique k-mer ids overlapping between two random sequence windows of size 500 bp (~25 sec)
  print 'Approximate expected number of shared k-mers between two random sequences'
  os_approx = approximate_expected_overlap(kmer_ids, ks)
  
  # compute similarity scores for every K (~6 min for 1 Mbp)
  print 'Compute similarities'
  L = len(target_seq)
  N = L - window_size + 1
  sims = {k : [compute_similarity(kmer_counts_ref_dict[k], kmer_counts_target_dict[k][j], os_approx[k]) for j in range(N)] for k in ks}

  # compute max / mean similarity ratio and write to file
  ratios = [np.max(sims[k]) / np.mean(sims[k]) for k in ks]
  pd.Series(ratios).to_csv(target_bed.replace('.bed', '_ratios_max_mean_similarity').replace('bed/','output/'), index=False)
  
  # write bigwig with similarities (~7 sec per file)
  print 'Write bigwig'
  for k in ks:
    bigwig_file = target_bed.replace('.bed', '_%s-mer_similarity_score.bw' %k).replace('bed/','output/')
    write_bigwig(chrs, start, end, sims[k], bigwig_file, 'mm10')
  print '%s: Done' %ref_region.name
  
  return

def main():
  usage_msg = '''Usage:
python compute_kmer_similarity.py find_homolog reference.bed reference_assembly target_assembly CNEs_reference_target.bed Ks_comma_separated
python compute_kmer_similarity.py pairwise regions.bed assembly Ks_comma_separated

pairwise:
This function computes pairwise sequence similarities between pairs of sequences from a given bed-file.

find_homolog:
This function takes a single sequence in species X, projects it genomic location to species Y based on interpolated CNE coordinates
and returns sequence similarities in a 1 Mbp window around this projected center in a bigwig file.'''

  window_size = 500
  array_size = 5e5
  try:
    if sys.argv[1] == 'pairwise':
      _, fnc, regions_bed, assembly, ks = sys.argv
      ks = map(int, ks.split(','))
      pairwise(regions_bed, assembly, ks)
    elif sys.argv[1] == 'find_homolog':
      _, fnc, reference_bed, reference_assembly, target_assembly, cnefile, ks = sys.argv
      ks = map(int, ks.split(','))
      find_homolog(reference_bed, reference_assembly, target_assembly, cnefile, ks, array_size, window_size)
    else:
      print usage_msg
      sys.exit(1)
  except (IndexError, ValueError):
    print usage_msg
    sys.exit(1)

  print 'Done'
  
if __name__ == '__main__':
  main()
