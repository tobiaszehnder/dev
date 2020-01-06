#! /usr/bin/env Rscript

# This script takes a bed file and returns a fasta with the sequence from the BSgenome package.
# In addition, you can specify a filter bedfile that is used to mask the overlapping base pairs in the sequence with N's.
# For example, if the bedfile is chr1:101-110 and the filter bedfile is chr1:103-106, the resulting sequence may look like this: 'ACNNNNGACT'

args <- commandArgs(trailingOnly=T)
if (!(length(args) %in% c(2,3))) {
   stop('Usage: Rscript get_sequence.R bed_file BSgenome_package filter_bed_file')
}

bedfile <- args[1]
package <- args[2]

suppressMessages(library(rtracklayer)); suppressMessages(library(package, character.only=T))

gr <- import.bed(bedfile)
if ((length(names(mcols(gr))) == 0) || (names(mcols(gr)) != c('name'))) {
   stop('Error: bed_file must contain 4 columns with the 4th containing the name of the region(s)')
}
if(all(is.na(gr$name))) gr$name <- seq_len(length(gr))

seq <- getSeq(get(package), gr)

### I moved the following masking step to after the similarity score calculation.
### If I'd do it here (replace nucleotides of repeats with 'N'), sequences with mostly N and only a few k-mer counts can obtain high similarity scores because their L2 norm gets really small
# # load filter regions
# if (length(args) == 3) {
#    filter <- import.bed(args[3], which=gr)
# } else {
#    filter <- GRanges()
# }
#
# # mask filter regions
# coords <- do.call('c', sapply(seq_along(filter), function(i) {
#   start <- max(1, start(filter)[i] - start(gr))
#   end <- min(end(filter)[i] - start(gr), width(gr))
#   return(seq(start, end, 1))
# }))
# seq_masked <- as.character(seq)
# seq_masked_split <- strsplit(seq_masked, '')[[1]]
# seq_masked_split[coords] <- 'N'
# seq_masked <- do.call('paste0', as.list(seq_masked_split))

fastafile <- gsub('bed', 'fasta', bedfile)
lines <- sapply(seq_len(length(gr)), function(i) c(paste0('>',gr$name[i]), as.character(seq[i])))
writeLines(lines, fastafile)
