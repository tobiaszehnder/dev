#! /usr/bin/env python

import numpy as np, pandas as pd, pybedtools as pb, scipy, sys, os.path, itertools, subprocess, json, datetime
from Bio import SeqIO
from multiprocessing.pool import ThreadPool
sys.path.insert(1, '/project/wig/tobias/reg_evo/jupyter/')
if 'functions' in sys.modules:
  del sys.modules['functions'] # deletes previously imported module so that potential changes will be loaded
from functions import *

'''
Count k-mers of a reference region and its projected target (i.e. 1 MB region centered on the projected target coordinate)
Combine rev / comp / rev-comp in kmer-counts

Introduce a dict with IDs such that d['AAAT'] = d['TAAA'] = d['TTTA'] = d['ATTT']
Elongate every k-mer with '#' (and add this character to the alphabet) to prevent k-mers with different k to have the same ID

ID_kmer# = sum(a(s_i) len(a)^i)
with alphabet a={A:0,C:1,G:2,T:3,#:4}
such that e.g. IDTATC#=3x50+0x51+3x52+1x53+4x54=3+0+75+125+2500=2703
For every possible k-mer, assign the minimum of the IDs of the k-mer, its reverse, complement and reverse complement to all of those k-mers.
'''

def count_kmers(seq, kmer_ids, ks=[4,6,8,10]):
  id_counter = {}
  kmer_counts = {key : val for d in [collections.Counter(seq[start:start+K] for start in xrange(len(seq)-K+1)) for K in ks] for key, val in d.items()}
  kmer_counts_items = np.delete(kmer_counts.items(), np.where(['N' in x for x in kmer_counts.keys()]), axis=0) # skips k-mers with 'N'
  for k,v in kmer_counts_items:
    try:
      id_counter[kmer_ids[k]] += int(v)
    except KeyError:
      id_counter[kmer_ids[k]] = int(v)
  return id_counter

def compute_similarity(kmer_counts_ref, kmer_counts_target, expected_overlap_approx):
  # this function computes a similarity score that is normalized with the expected number of overlapping k-mer IDs between two sequences.
  # note that this function only works for one K at a time.
  common_ids = set(kmer_counts_ref).intersection(set(kmer_counts_target))
  similarity_score = len(common_ids) / expected_overlap_approx
  return similarity_score

def write_bigwig(chrs, start, end, scores, outfile, assembly):
  genome_size = read_genome_size(assembly)
  bw = pyBigWig.open(outfile, 'w')
  bw.addHeader(genome_size)
  bw.addEntries(chrs, start, end, scores)
  bw.close()

def main():
  if not len(sys.argv) == 5:
    print 'Usage: python count_kmers.py reference.bed reference_assembly target_assembly CNEs_reference_target.bed'
    sys.exit(1)
  _, reference_bed, reference_assembly, target_assembly, cnefile = sys.argv

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
  pd.DataFrame(np.array([[cne_upstream.values[0,3], int(ref_target_position-5e5), int(ref_target_position+5e5), ref_region.name]])).to_csv(target_bed, sep='\t', header=None, index=False)

  # write fasta and fetch sequence using get_sequence.R
  print 'Fetch target sequence'
  subprocess.check_call(['Rscript', 'get_sequence.R', target_bed, bsgenome_dict[target_assembly]], shell=False)
  with open(target_bed.replace('bed', 'fasta'), 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
      target_seq = str(record.seq)
  
  # create dict with IDs and initiate id_counter dict
  # runtime ~ 20 sec
  print 'Assign k-mer IDs'
  ks = [4,6,8,10]
  possible_kmers = [''.join(p) for k in ks for p in itertools.product(list('ACGT'), repeat=k)]
  kmer_ids = {kmer : get_kmer_id(kmer) for kmer in possible_kmers}

  # compute kmers for reference
  print 'Count reference k-mers'
  global kmer_counts_ref
  kmer_counts_ref_dict = {k : count_kmers(ref_seq, kmer_ids, ks=[k]) for k in ks}

  # compute kmers for target
  print 'Count target k-mers'
  L = len(target_seq)
  n = 1000
  N = L - n + 1
  t1 = datetime.datetime.now()
  kmer_counts_target_dict = {}
  for k in ks:
    print k
    t1 = datetime.datetime.now()
    kmer_counts_target_dict[k] = [count_kmers(target_seq[i:i+n], kmer_ids, ks=[k]) for i in range(N)]
    print datetime.datetime.now() - t1
  # kmer_counts_target_dict = {k : [count_kmers(target_seq[i:i+n], kmer_ids, ks=[k]) for i in range(N)] for k in ks}
  with open('time_for_counting_target_kmers', 'w') as f:
    f.write(str(datetime.datetime.now() - t1))
# read target bed and create data frame with 1 bp resolution                                                                                                            
  target_regions = pd.read_csv(target_bed, sep='\t', header=None).loc[0]                                                                                                  
  start = np.arange(target_regions[1], target_regions[1] + len(kmer_counts_target_dict[4])) # the bigwig's last entry starts at 1000 bp before the target region's end co\
ordinate.                                                                                                                                                                 
  end = start + 1                                                                                                                                                         
  chrs = np.repeat(target_regions[0], len(start))   
  # read target bed and create data frame with 1 bp resolution
  target_regions = pd.read_csv(target_bed, sep='\t', header=None).loc[0]
  start = np.arange(target_regions[1], target_regions[1] + len(kmer_counts_target_dict[4])) # the bigwig's last entry starts at 1000 bp before the target region's end coordinate.
  end = start + 1
  chrs = np.repeat(target_regions[0], len(start))

  # approximate E[o] = the expected number of unique k-mer ids overlapping between two random sequence windows of size 1000 bp
  print 'Approximate expected number of shared k-mers between two random sequences'
  os_approx = {}
  ns = [float(len(set([v for k,v in kmer_ids.items() if len(k) == K]))) for K in ks]
  t1 = datetime.datetime.now()
  for i,K in enumerate(ks):
    n = ns[i]
    k = 1000.
    expected_m = n - n*(1-1/n)**k
    expected_o = expected_m**2 / n
    os_approx[K] = expected_o
  with open('time_for_counting_E_o', 'w') as f:
    f.write(str(datetime.datetime.now() - t1))
  # compute similarity scores for every K
  print 'Compute binarized cosine similarities'
  sims = {k : [compute_similarity(kmer_counts_ref_list_dict[k], kmer_counts_target_dict[k][j], os_approx[k]) for j in range(N)] for k in ks}
    
  # compute binarized cosine similarities
  # sim = [compare_to_enhancer(kmer_counts_target[i]) for i in range(N)]
  
  # write bigwig with similarities
  print 'Write bigwig'
  for k in ks:
    bigwig_file = target_bed.replace('.bed', '_%s-mer_similarity_score.bw' %k).replace('bed/','output/')
    write_bigwig(chrs, start, end, sims[k], bigwig_file, 'mm10')
  print '%s: Done' %ref_region.name
  
if __name__ == '__main__':
  main()
