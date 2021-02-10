#! /usr/bin/make -f

### This script is for downloading, decompressing and indexing genome fasta files

### Usage: ./download_genome_fasta.sh outdir=/project/MDL_ChIPseq/data/genome/fasta build=mm10 link=ftp.ncbi.nlm.nih.gov/*.fna.gz

### IMPORTANT: omit leading 'ftp://' or 'http://' in link. make has problems with colons...

TARGETS = $(outdir)/$(build).fa.fai

# ------------------------------------------------------------------------------

all: $(TARGETS)

# ------------------------------------------------------------------------------

.DELETE_ON_ERROR:
.SECONDARY:

# ------------------------------------------------------------------------------

$(outdir)/$(build).fa.gz:
	curl -o $@ $(link)

%.fa: %.fa.gz
	gzip -d $<

%.fa.fai: %.fa
	samtools faidx $<

# ------------------------------------------------------------------------------
