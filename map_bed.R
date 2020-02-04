#! /usr/bin/env Rscript

### This script maps all the regions in a given query bed from a query to a reference genome using GRB and CNE information
### I implemented a 2-step algorithm including a bridging species with the goal to increase projection accuracy by increasing the number of CNEs as anchor points for the local regression
### Local regression models are fitted separately for every GRB that overlaps with a query region.
### Hence, query-reference GRBs as well as CNEs from the query and reference to the bridging species are required.

args <- commandArgs(trailingOnly=T)
if (length(args) != 7) stop('Usage: ./map_bed.R query_regions_to_project.bed> <grbs_ref_query.bed> <grbs_query_ref.bed> <cnes_qry_ref> <cnes_ref_bridge.bed> <cnes_qry_bridge.bed> <out.bed>')

cat('Load packages\n')
suppressMessages(library(rtracklayer)); suppressMessages(library(ggplot2)); suppressMessages(library(CNEr)); suppressMessages(library(magrittr))
source('functions.R')

bed_qry_file <- args[1]
grbs_ref_qry_file <- args[2]
grbs_qry_ref_file <- args[3]
cnes_qry_ref_file <- args[4]
cnes_ref_bridge_file <- args[5]
cnes_qry_bridge_file <- args[6]
bed_outfile <- args[7]

bed_outfile_direct <- gsub('.bed', '_direct.bed', bed_outfile)
bed_outfile_bridge <- gsub('.bed', '_bridge.bed', bed_outfile)

cat('Load data\n')
bed_qry <- import.bed(bed_qry_file)
grbs <- GRangePairs(import.bed(grbs_qry_ref_file), import.bed(grbs_ref_qry_file)) %>% .[overlapsAny(.@first,bed_qry)]
grbs@first$id <- grbs@second$id <- seq_along(grbs)
# allocate a GRB ID. This only works for regions that can be allocated to one unique GRB
ov <- findOverlaps(bed_qry, grbs@first)
bed_qry$grb_id <- NA
bed_qry$grb_id[queryHits(ov)] <- subjectHits(ov)
get_cne_pairs <- function(filename) import.bed(filename) %>% GRangePairs(., GRanges(data.frame(.)$name))
cnes <- list(ref=list(bridge=get_cne_pairs(cnes_ref_bridge_file)),
             qry=list(ref=get_cne_pairs(cnes_qry_ref_file),
			 bridge=get_cne_pairs(cnes_qry_bridge_file)))

cat('Fit direct model and project regions\n')
# fit models for every GRB
bed_ref_qry_direct <- do.call('c', mclapply(seq_along(grbs), function(i) {
  cnes_i <- get_GRangePairs(grbs@first[i], cnes$qry$ref@first)
  # remove CNE outliers, i.e. CNEs that are unusually far from the median position of CNEs in GRB in only one of both species
  cnes_i <- remove_CNE_outliers(cnes_i, grbs@first[i])
  model_i <- fit_loess(cnes_i)
  gr <- bed_qry[bed_qry$grb_id == i]
  s <- sapply(start(gr), function(x) predict(model_i, data.frame(x=x), se=F))
  e <- sapply(end(gr), function(x) predict(model_i, data.frame(x=x), se=F))
  # coordinates outside the CNE range won't be extrapolated and return NA. remove them from all ranges and print message
  idx <- which(!is.na(s) & !is.na(e))
  if (length(idx) < length(gr)) {
    print(sprintf("GRB %s: The following regions lie outside the GRB's CNE range and were not mapped:", i))
    print(gr[setdiff(seq_along(gr), idx)])
  }
  if (length(idx) == 0) {
    gr_ref_qry_i <- GRangePairs()
  } else {
    s <- s[idx]
    e <- e[idx]
    gr <- gr[idx]
    coords <- order_coords(s,e)
    gr_out <- GRanges(seqnames=unique(seqnames(cnes_i@second)), IRanges(coords['start',], coords['end',]))
    gr_ref_qry_i <- GRangePairs(gr, gr_out)
  }
  return(gr_ref_qry_i)
}, mc.cores=length(grbs)))

cat('Fit bridging species model and project regions\n')
# fit models for every GRB
bed_ref_qry_bridge <- do.call('c', mclapply(seq_along(grbs), function(i) {
  cnes_i <- list(xb=get_GRangePairs(grbs@first[i], cnes$qry$bridge@first),
                 by=invert_GRangePairs(get_GRangePairs(grbs@second[i], cnes$ref$bridge@first)))
  # remove CNE outliers, i.e. CNEs that are unusually far from the median position of CNEs in GRB in only one of both species
  cnes_i$xb <- remove_CNE_outliers(cnes_i$xb, grbs@first[i])
  model_i <- lapply(list(xb=cnes_i$xb,by=cnes_i$by), function(grp) fit_loess(grp))
  gr <- bed_qry[bed_qry$grb_id == i]
  s1 <- sapply(start(gr), function(x) predict(model_i$xb, data.frame(x=x), se=F))
  s2 <- sapply(s1, function(y) predict(model_i$by, data.frame(x=y), se=F))
  e1 <- sapply(end(gr), function(x) predict(model_i$xb, data.frame(x=x), se=F))
  e2 <- sapply(e1, function(y) predict(model_i$by, data.frame(x=y), se=F))
  # coordinates outside the CNE range won't be extrapolated and return NA. remove them from all ranges and print message
  idx <- which(!is.na(s1) & !is.na(s2) & !is.na(e1) & !is.na(e2))
  if (length(idx) < length(gr)) {
    print(sprintf("GRB %s: The following regions lie outside the GRB's CNE range and were not mapped:", i))
    print(gr[setdiff(seq_along(gr), idx)])
  }
  if (length(idx) == 0) {
    gr_ref_qry_i <- GRangePairs()
  } else {
    s1 <- s1[idx]
    s2 <- s2[idx]
    e1 <- e1[idx]
    e2 <- e2[idx]
    gr <- gr[idx]
    coords1 <- order_coords(s1,e1)
    coords2 <- order_coords(s2,e2)
    gr_out <- GRanges(seqnames=unique(seqnames(cnes_i$by@second)), IRanges(coords2['start',],coords2['end',]))
    gr_ref_qry_i <- GRangePairs(gr, gr_out)
  }
  return(gr_ref_qry_i)
}, mc.cores=length(grbs)))

cat('Write bed file\n')
export(bed_ref_qry_direct@second, bed_outfile_direct)
export(bed_ref_qry_bridge@second, bed_outfile_bridge)

cat('Done\n')