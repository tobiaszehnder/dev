#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if (length(args) != 1) stop('Usage: Rscript axtToBed.R <axtFile>')

suppressMessages(library(rtracklayer)); suppressMessages(library(CNEr)); suppressMessages(library(ehmm))

# read genome size file for target species
# reference <- strsplit(basename(args[1]), '\\.')[[1]][1]
target <- strsplit(basename(args[1]), '\\.')[[1]][2]
sizes_target <- ehmm:::readGenomeSize(sprintf('/project/wig/tobias/reg_evo/data/assembly/%s.sizes', target), dropNonStandardLevels=F)

# read axt file
axt <- readAxt(args[1])

# CNEr shows coordinates relative to the chromosome END if region is on '-' strand. correct that by subtracting the values from the respective chromosome sizes
# the result for those regions is that the start coordinate is higher than the end coordinate
# CNEr shows coordinates relative to the chromosome END if region is on '-' strand. correct that by subtracting the values from the respective chromosome sizes
# the result for those regions is that the start coordinate is higher than the end coordinate
chrsize <- sizes_target[as.vector(seqnames(axt@second[strand(axt@second)=='-']))]
gr1 <- axt@first
gr2 <- axt@second
seqnames_gr2_negativeStrand <- seqnames(axt@second[strand(axt@second)=='-'])
start_gr2_negativeStrand <- as.numeric(chrsize - end(axt@second[strand(axt@second)=='-']))
end_gr2_negativeStrand <- as.numeric(chrsize - start(axt@second[strand(axt@second)=='-']))
gr2_negativeStrand <- GRanges(seqnames_gr2_negativeStrand, IRanges(start=start_gr2_negativeStrand, end=end_gr2_negativeStrand), strand='-')
gr2_negativeStrand$seqs <- gr2[which(strand(gr2)=='-')]$seqs
gr2[which(strand(gr2)=='-')] <- gr2_negativeStrand
gr1$name <- paste0(seqnames(gr2), ':', start(gr2)-1, '-', end(gr2))
gr2$name <- paste0(seqnames(gr1), ':', start(gr1)-1, '-', end(gr2))
outfile1 <- gsub('.net.axt', '.bed', gsub('.net.axt.gz', '.bed', args[1]))
outfile2 <- gsub('.bed', '.reversed.bed', outfile1)
export.bed(gr1, outfile1)
export.bed(gr2, outfile2)
