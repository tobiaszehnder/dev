#! /usr/bin/env Rscript

# This script takes a bed file and returns a fasta with the sequence from the BSgenome package.

args <- commandArgs(trailingOnly=T)
if (length(args) != 2) {
   stop('Usage: Rscript get_sequence.R bed_file BSgenome_package')
}

bedfile <- args[1]
package <- args[2]

suppressMessages(library(rtracklayer)); suppressMessages(library(package, character.only=T))

gr <- import.bed(bedfile)
if ((length(names(mcols(gr))) == 0) || (names(mcols(gr)) != c('name'))) {
   stop('Error: bed_file must contain 4 columns with the 4th containing the name of the region(s)')
}
if(all(is.na(gr$name))) gr$name <- seq_len(length(gr))

# this horrible snippet is an hommage to R
# because BSgenome::getSeq requires either the loaded BSgenome.xxx.UCSC.xx object or the short version (e.g. Hsapiens), it is not possible to use the passed package character object.
if (grepl('Hsapiens', package)) {
   seq <- getSeq(Hsapiens, gr)
} else if (grepl('Mmusculus', package)) {
   seq <- getSeq(Mmusculus, gr)
} else stop('Error: script only includes Hsapiens and Mmusculus. Feel free to add other assemblies if necessary.')

fastafile <- gsub('bed', 'fasta', bedfile)
lines <- sapply(seq_len(length(gr)), function(i) c(paste0('>',gr$name[i]), as.character(seq[i])))
writeLines(lines, fastafile)
