#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if (length(args) != 4) stop('Usage: ./call_conserved_elements.R <species1> <species2> <type[CE/CNE]> <processing[merged/final]')

suppressMessages(library(rtracklayer)); suppressMessages(library(CNEr)); suppressMessages(library(Organism.dplyr))
data_dir <- '/project/wig/tobias/reg_evo/data'

# ------------------------------
# args
# ------------------------------

spcs1 <- args[1]
spcs2 <- args[2]
type <- tolower(args[3])
if (!(type %in% c('ce', 'cne'))) stop('type must be either CE or CNE')
processing <- args[4]
if (!(processing %in% c('merged', 'final'))) stop('processing must be either merged or final')
window <- 50L
identity <- 35L
cne_dir <- '/project/wig/tobias/reg_evo/data/CNEs/CNEr'

# ------------------------------
# functions
# ------------------------------

get_assemblyFn <- function(spcs) file.path(data_dir, sprintf('assembly/%s.2bit', spcs))
get_axtFn <- function(spcs1, spcs2) file.path(data_dir, sprintf('alignment/axtNet/%s.%s.net.axt.gz', spcs1, spcs2))

# function to load repeats from '/project/wig/tobias/reg_evo/data/repeats'
get_repeats <- function(spcs, repeat_dir='/project/wig/tobias/reg_evo/data/repeats') {
  if (!dir.exists(repeat_dir)) stop('repeat_dir does not exist')
  gr <- GRanges()
  repeats_file <- file.path(repeat_dir, sprintf('repeats_%s.bed', spcs))
  if (file.exists(repeats_file)) {
    gr <- import.bed(repeats_file)
  }
  gr <- tryCatch(keepStandardChromosomes(gr, pruning.mode='coarse'),
                 error=function(cond) return(gr))
  return(gr)
}

# function to retrieve exons from TxDb if available
get_exons <- function(spcs) {
  gr <- tryCatch({
    txdb <- supportedOrganisms()$TxDb[which(grepl(spcs, supportedOrganisms()$TxDb))]
	if (length(txdb) > 1) txdb <- txdb[grepl('knownGene', txdb)] # because mm10 has two, .ensGene and .knownGene
    library(txdb, character.only=T)
    exons(eval(as.name(txdb)))
    }, error=function(cond) return(GRanges()))
  gr <- tryCatch(keepStandardChromosomes(gr, pruning.mode='coarse'),
                 error=function(cond) {
                   message(cond)
                   return(gr)
                   })
  return(gr)
}

write_conserved_elements <- function(x, processing=c('CNEFinal', 'CNEMerged'), elementType=c('cne', 'ce'), cne_dir='/project/wig/tobias/reg_evo/data/CNEs/CNEr') {
  # function to write CEs / CNEs.
  x <- sort(get(processing)(x)) # equivalent to CNEFinal(x) or CNEMerged(x) given the passed argument processing
  processing <- tolower(gsub('CNE', '', processing))
  gr1 <- x@first
  gr1$name <- paste0(seqnames(x@second), ':', start(x@second), '-', end(x@second))
  gr2 <- x@second
  gr2$name <- paste0(seqnames(x@first), ':', start(x@first), '-', end(x@first))
  strand(gr1) <- strand(gr2)
  export.bed(gr1, file.path(cne_dir, sprintf('%s_%s_%s_%s_%s_%s.bed', elementType, processing, spcs1, spcs2, identity, window)))
  # data1 <- data.frame(x)[,c(1,2,3,6,7,8,10,11,12)]
  # write.table(data1, file.path(cne_dir, sprintf('%s_%s_%s_%s_%s_%s', elementType, processing, spcs1, spcs2, identity, window)), sep='\t', quote=F, row.names=F, col.names=F)
  export.bed(gr2, file.path(cne_dir, sprintf('%s_%s_%s_%s_%s_%s.bed', elementType, processing, spcs2, spcs1, identity, window)))
  # data2 <- data.frame(x)[,c(6,7,8,1,2,3,10,11,12)]
  # write.table(data2, file.path(cne_dir, sprintf('%s_%s_%s_%s_%s_%s', elementType, processing, spcs2, spcs1, identity, window)), sep='\t', quote=F, row.names=F, col.names=F)
  return(invisible())
}

# ------------------------------
# ------------------------------

# create CNE objects
cne_obj <- CNE(assembly1Fn=get_assemblyFn(spcs1),
               assembly2Fn=get_assemblyFn(spcs2),
               window=window,
               identity=identity,
               axt12Fn=get_axtFn(spcs1,spcs2),
               axt21Fn=get_axtFn(spcs2,spcs1))

# load filter regions
cat('Loading filter regions\n')
r <- sapply(c(spcs1, spcs2), function(x) get_repeats(x))
if (type == 'ce') {
  fltr <- r
} else if (type == 'cne') {
  e <- sapply(c(spcs1, spcs2), function(x) get_exons(x))
  fltr <- sapply(c(spcs1, spcs2), function(i) c(r[[i]], e[[i]]))
}

# if no exons are found for both species, CEs = CNEs --> create symbolic link to previously calculated CEs
if (type == 'cne' & all(fltr[[spcs1]] == r[[spcs1]]) & all(fltr[[spcs2]] == r[[spcs2]])) {
   outfile12 <- file.path(cne_dir, sprintf('cne_%s_%s_%s_%s.bed', spcs1, spcs2, identity, window))
   outfile21 <- file.path(cne_dir, sprintf('cne_%s_%s_%s_%s.bed', spcs2, spcs1, identity, window))
   system(sprintf('ln -s %s %s', gsub('cne', 'ce', outfile12), outfile12))
   system(sprintf('ln -s %s %s', gsub('cne', 'ce', outfile21), outfile21))
} else {
  # call and merge CEs/CNEs and remove unannotated repeats (i.e. CNEs that map to > 8 locations)
  cat(sprintf('scanning %ss\n', toupper(type)))
  y <- ceScan(x=cne_obj, tFilter=fltr[[1]], qFilter=fltr[[2]], window=window, identity=identity)[[sprintf('%s_%s', identity, window)]]

  cat(sprintf('merging %ss\n', toupper(type)))
  y <- cneMerge(y) # ~ 10 sec
  if (processing == 'merged') {
    cat(sprintf('writing merged %ss to file\n', toupper(type)))
  	write_conserved_elements(y, 'CNEMerged', type, cne_dir)
  }
  
  if (processing == 'final') {
    cat('mapping back to genome, removing repeats\n')
  	y <- blatCNE(y) # 1.5 hours
  	cat(sprintf('writing final %ss to file\n', toupper(type)))
  	write_conserved_elements(y, 'CNEFinal', type, cne_dir)
  }
}

cat('Done\n')
# ------------------------------
