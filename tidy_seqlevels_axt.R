#! /usr/bin/env Rscript

### This script takes all *.net.axt* files in the directory and checks the seqlevels of the ranges so that they start with 'chr'.
### It also corrects lepOcu1-specific things (e.g. AHA12345.1 to chrUn_AHA12345)

# args <- commandArgs(trailingOnly=T)
# if (length(args) != 1) stop('Usage: ./tidy_seqlevels_axt.R alignment.net.axt.gz')

library(CNEr); library(parallel); suppressMessages(library(rtracklayer))

tidy_seqlevels <- function(s) {
  gsub_pairs_start <- list(c('AHA', 'chrUn_AHA'), c('JH', 'chrUn_JH'))
  for (pair in gsub_pairs_start) {
    s[startsWith(s,pair[1])] <- gsub(pair[1], pair[2], s[startsWith(s,pair[1])])
  }
  gsub_pairs <- list(c('MT', 'M'), c('\\.1', ''))
  for (pair in gsub_pairs) {
    s <- gsub(pair[1], pair[2], s)
  }
  # prepend 'chr' for all seqlevels except the calMil1 ones ('KI...' and 'AAV...')
  s[!startsWith(s, 'chr') & !startsWith(s, 'KI') & !startsWith(s, 'AAV')] <- paste0('chr', s[!startsWith(s, 'chr')])
  return(s)
}

files <- list.files(path='.', pattern='*.net.axt*')
axts <- mclapply(files, function(filename) {
  # read axt file and tidy seqlevels
  axt <- readAxt(filename)
  l <- lapply(list(axt@first, axt@second), function(x) tidy_seqlevels(seqlevels(x)))

  # compare if tidy seqlevels are all in assembly files
  spcs <- strsplit(filename, '\\.')[[1]][1:2]
  sl <- sapply(spcs, function(sp) seqlevels(TwoBitFile(sprintf('/project/wig/tobias/reg_evo/data/assembly/%s.2bit', sp))))
  diff <- sapply(seq_along(l), function(i) setdiff(l[[i]], sl[[spcs[i]]]))
  for (i in seq_along(diff)) {
    if (length(diff[[i]])>0) stop(sprintf('The following seqlevels of %s are not in %s.2bit: %s', spcs[i], spcs[i], diff[i]))
  }

  # only rewrite files if there was any change in seqlevels
  change <- F
  if (any(seqlevels(axt@first) != l[[1]])) {
    seqlevels(axt@first) <- l[[1]]
    change <- T
  }
  if (any(seqlevels(axt@second) != l[[2]])) {
    seqlevels(axt@second) <- l[[2]]
    change <- T
  }
  if (change == T) {
	filename_tmp <- gsub('.gz', '.newseqlevels', filename)
	writeAxt(axt, filename_tmp)
	system(sprintf('gzip %s', filename_tmp))
	system(sprintf('mv %s.gz %s', filename_tmp, filename))
    print(sprintf('%s: Seqlevels tidied up', filename))
  } else {
    print(sprintf('%s: Seqlevels already tidy. Nothing done.', filename))
  }
  return(axt)
}, mc.cores=min(length(files), 40))
