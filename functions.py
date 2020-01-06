# This script contains functions that can be loaded into other scripts / notebooks using the following line
# from functions import *

import numpy as np, pandas as pd, pickle, pybedtools as pb, os, pyBigWig, collections

# function to set scores at given indices (mask_coords) to zero. this is used for the RepeatMasker regions that overlap the target region.
def mask_filter(scores, mask_coords):
  scores[mask_coords] = np.nan
  return scores.tolist()

def get_kmer_id(x):
  # compute a unique ID for a given kmer as the sum of the products of the letter at position i (A=0, C=1, ...) times alphabet size to the power of the position number
  # (e.g. 5^0 at the first position). include a character `x` to prevent k-mers with different k to have the same ID.
  alph = dict(zip(list('ACGT#'), range(5)))
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
  rev = x[::-1]
  comp = ''.join([complement[x_i] for x_i in x])
  rev_comp = comp[::-1]
  kmrs = [y + '#' for y in [x, rev, comp, rev_comp]]
  l = len(alph)
  return min([sum([alph[kmr[i]] * l**i for i in range(len(kmr))]) for kmr in kmrs])

def count_kmers(seq, kmer_ids, ks=[4,6,8,10]):
  id_counter = {}
  kmer_counts = {key : val for d in [collections.Counter(seq[start:start+K] for start in range(len(seq)-K+1)) for K in ks] for key, val in d.items()}
  kmer_counts_items = np.delete(np.array(list(kmer_counts.items())), np.where(['N' in x for x in kmer_counts.keys()]), axis=0) # skips k-mers with 'N'
  for k,v in kmer_counts_items:
    try:
      id_counter[kmer_ids[k]] += int(v)
    except KeyError:
      id_counter[kmer_ids[k]] = int(v)
  return id_counter

def count_kmers_wrapper(seq, kmer_ids, ks=[4,6,8,10], window_size=500):
  # this function takes a sequence and computes kmer counts for moving windows of 1000 bp
  L = len(seq)
  N = L - window_size + 1
  kmer_counts = [{} for i in range(N)]
  kmer_counts[0] = count_kmers(str(seq[:window_size]), kmer_ids, ks)
  for i in range(1,N):
    kmer_counts[i] = kmer_counts[i-1].copy()
    kmers_to_remove = collections.Counter(seq[(i-1):(i-1+K)] for K in ks)
    kmers_to_remove = {k:v for k,v in kmers_to_remove.items() if not 'N' in k}
    kmers_to_add = collections.Counter(seq[:(i+window_size)][-K:] for K in ks)
    kmers_to_add = {k:v for k,v in kmers_to_add.items() if not 'N' in k}
    for k,v in kmers_to_remove.items():
      if kmer_counts[i][kmer_ids[k]] == v:
        del kmer_counts[i][kmer_ids[k]] # remove the id from the dict in case the count was set to zero in order to keep the dicts sparse
      else:
        kmer_counts[i][kmer_ids[k]] -= v
    for k,v in kmers_to_add.items():
      try:
        kmer_counts[i][kmer_ids[k]] += v
      except KeyError:
        kmer_counts[i][kmer_ids[k]] = v
  return kmer_counts

def compute_similarity(kmer_counts_ref, kmer_counts_target, method, expected_overlap_approx=None):
  # This function computes a similarity score between two sequences. Note that this function only works for one K and method at a time.
  sim = {}
  ref = kmer_counts_ref.copy()
  target = kmer_counts_target.copy()
  if len(target) == 0: # if the whole target window is 'NNNNNN...', then there will be no k-mer counts and the L2 norm will get zero and a division by zero will occur. return a similarity of zero.
    return 0
  common_ids = set(ref).intersection(set(target))
  union_keys = {key : val for d in [ref, target] for key, val in d.items()}.keys()
  if method == 'overlap':
    return len(common_ids) / expected_overlap_approx
  if 'tanimoto' in method:
    ref.update({key : 0 for key in list(set(union_keys) - set(ref.keys()))})
    target.update({key : 0 for key in list(set(union_keys) - set(target.keys()))})
    if method == 'tanimoto':
      ref_vector = [ref[key] for key in union_keys]
      target_vector = [target[key] for key in union_keys]
      return np.dot(ref_vector, target_vector) / float(np.dot(ref_vector, ref_vector) + np.dot(target_vector, target_vector) - np.dot(ref_vector, target_vector))
    elif method == 'tanimoto_binarized':
      ref_binarized = np.where(np.array(ref_vector) > 0, 1, 0)
      target_binarized = np.where(np.array(target_vector) > 0, 1, 0)
      return np.dot(ref_binarized, target_binarized) / float(np.dot(ref_binarized, ref_binarized) + np.dot(target_binarized, target_binarized) - np.dot(ref_binarized, target_binarized))
  if method in ['dot', 'dot_binarized', 'cosine', 'cosine_binarized']:
    ref_common = np.array([ref[cid] for cid in common_ids], dtype=float)
    target_common = np.array([target[cid] for cid in common_ids], dtype=float)
  if method == 'dot':
    return np.dot(ref_common, target_common) / float(len(union_keys))
  if 'cosine' in method:
    ref_L2_norm = np.linalg.norm(list(ref.values())) # important: norm is not from sparse but from full values, not just common
    target_L2_norm = np.linalg.norm(list(target.values()))
  if method == 'cosine':
    return np.dot(ref_common, target_common) / (ref_L2_norm * target_L2_norm)
  if 'binarized' in method:
    ref_common_binarized = np.where(ref_common > 0, 1, 0)
    target_common_binarized = np.where(target_common > 0, 1, 0)
  if method == 'dot_binarized':
    return np.dot(ref_common_binarized, target_common_binarized) / float(len(union_keys))
  if method == 'cosine_binarized':
    return np.dot(ref_common_binarized, target_common_binarized) / (ref_L2_norm * target_L2_norm)
  return

# function to write bigwig files with the pyBigWig module
def write_bigwig(chrs, start, end, scores, outfile, assembly):
  genome_size = read_genome_size(assembly)
  bw = pyBigWig.open(outfile, 'w')
  bw.addHeader(genome_size)
  bw.addEntries(chrs, start, end, scores)
  bw.close()
  return

# function to read genome size of given assembly
def read_genome_size(assembly):
  command = os.popen('find /project/wig-data/*/genome/ -name mm10.genome 2>&1 | grep -v "Permission denied"')
  filename = command.read().strip('\n')
  _ = command.close()
  if filename == '':
    print('Assembly not found. Add file to /project/wig-data/<organism>/genome/<assembly>.genome')
    sys.exit(0)
  with open(filename, 'r') as f:
    genome_size = [tuple(str.split(x, '\t')) for x in f.read().splitlines()]
    genome_size = [(x[0], int(x[1])) for x in genome_size] # convert size variable to integer
  return genome_size

# function to compute cosine similarity of two vectors
def cosine_similarity(a,b):
  
  a_sparse
  return np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b)) # np.linalg.norm is the L2 norm of a vector (magnitude)

# this function converts the dict with aggregated energies into a long dataframe (melted)
def energies_to_long_dataframe(energies_agg, spcs, tf):
  df = pd.concat([pd.DataFrame(energies_agg[spcs][k][tf], columns=[k]) for k in energies_agg[spcs].keys()], axis=1).melt(var_name='label', value_name='energy')
  df.insert(0, 'TF', tf)
  df.insert(0,'species',spcs)
  return df

# fuction to put energy values from dict to wide dataframe
def energies_to_wide_dataframe(energies_agg, spcs, tf):
  df = pd.concat([pd.DataFrame(energies_agg[spcs][tf][k], columns=[k]) for k in energies_agg[spcs][tf].keys()], axis=1).melt(var_name='label', value_name='energy')
  df.insert(0, 'TF', tf)
  df.insert(0,'species',spcs)
  return df

# this function takes precalculated energies for the grb and allocates them to a set of enhancers (E = log(sum_i(exp(E_i))) with i being the enhancers' windows
# or, if no enhancer object is passed, for the whole GRB (background)
def aggregate_energies(energies, grb, enhancers=None, window_size=500):
  grb = grb.to_dataframe()
  l = int((grb.end - grb.start) - energies.shape[0] + 1) # motif length
  if enhancers == None:
    idx_start = np.arange(0, len(energies), window_size) # compute energies for whole grb (background)
  else:
    idx_start = np.array(enhancers.to_dataframe().start) - np.array(grb.start) # compute energies for enhancers only
  E =  [np.log(np.exp(energies.loc[x : (x + window_size-l),'89']).sum()) for x in idx_start]
  return E

# Return a table with one region obtained by combining the chromosome-name and start coordinate of the first entry and the end coordinate of the last entry of a given table.
# (I apply this function in get_Astates() to the df grouped by enhancer-ID.)
def collapse(df):
  return pd.DataFrame({'chr':[df.iloc[0].chr], 'start':[df.iloc[0].start], 'end':[df.iloc[-1].end]}, columns=['chr', 'start', 'end'])

# function to resize a pb.BedTool interval relative to the regions' centers.
def resize_interval(x, size):
  center = (x.start + x.end) / 2
  x.start = max(center - size/2., 0)
  x.end = center + size/2.
  return x

# function to split a pb.BedTool interval into regions of a given size.
# first, resize the whole interval to a multiple of the given size, then split.
def split_interval(x, size):
  n = max(int(np.round(float(x.length) / size)), 1)
  l = n * size
  x_resized = resize_interval(x, l)
  chrom = np.repeat(x.chrom, n)
  start = np.array([x_resized.start + size*i for i in range(n)])
  end = start + 500
  x_split_list = [pb.BedTool(' '.join(x), from_string=True) for x in zip(chrom, [str(x) for x in start], [str(x) for x in end])]
  if len(x_split_list) == 1:
    return x_split_list[0]
  x_split = x_split_list[0].cat(*x_split_list[1:], postmerge=False)
  return x_split

# Load a given eHMM segmentation file, extract the A-states and reduce neighboring regions. Center and resize to 200 bp to get a standardized measure of conservation.
def get_Astates(path, size=None):
  # load segmentation
  seg = pd.read_table(path, header=None, usecols=range(4), names=['chr', 'start', 'end', 'state'])
  # take only A-states
  acc = seg.loc[seg['state'].str.startswith('E_A')]
  # assign common IDs to neighbors
  # compare the difference of the start coordinate and the end coordinate of the previous entry. if it is zero, the regions are neighbours and get the same ID that can later be used to group the regions.
  diff = np.array(acc.start[1:]) - np.array(acc.end[:-1])
  id = np.array([0])
  for i in np.arange(0,len(diff)):
    if diff[i] == 0:
      j = id[-1]
    else:
      j = id[-1] + 1
    id = np.append(id, j)
  acc.index = id
  # group regions by ID and resize to 200 bp
  acc_reduced = acc.groupby(by=acc.index).apply(collapse).reset_index(drop=True)
  center = np.array((acc_reduced['start'] + acc_reduced['end']) / 2, dtype=int)
  acc_resized = acc_reduced
  if size:
    acc_resized['start'] = max(center - size/2, 0)
    acc_resized['end'] = center + size/2
  return acc_resized

# write a 'pickled' object to file
def pickle_obj(obj, name):
  with open(name, 'wb') as f:
    pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
  return

# load a 'pickled' object from file
def unpickle_obj(name):
  with open(name, 'rb') as f:
    return pickle.load(f)
