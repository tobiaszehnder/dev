#! /usr/bin/env Rscript-4

### This script masks the regions overlapping a repeat from a given input bed-file (i.e. removes them)

# read arguments
args <- commandArgs(trailingOnly=T)
if (length(args) != 2) stop('Usage: Rscript mask_repeats.R input.bed repeats.bed')

# load packages
suppressMessages(library(rtracklayer))

# read bed-files
gr <- import.bed(args[1])
rpt <- import.bed(args[2])

# remove regions that overlap repeats (full removal if only part of it overlaps a repeat)
gr_masked <- suppressWarnings(sort(gr[!gr %over% rpt,]))

# write to file
export.bed(gr_masked, gsub('.bed', '.repeatMasked.bed', args[1]))
