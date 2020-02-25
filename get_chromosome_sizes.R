#! /usr/bin/env Rscript

# This script takes a 2bit, computes the chromosome sizes and writes it into the same directory (x.2bit -> x.sizes)
# This is faster than UCSC's `faSize`

args <- commandArgs(trailingOnly=T)
if (length(args) != 1) stop('Usage: ./get_chromosome_sizes.R <assembly.2bit>')
filename <- args[1]

suppressMessages(library(CNEr)); suppressMessages(library(rtracklayer))
info <- seqinfo(TwoBitFile(filename))
write.table(data.frame(seqnames(info), seqlengths(info)), gsub('.2bit', '.sizes', filename), sep='\t', quote=F, row.names=F, col.names=F)
