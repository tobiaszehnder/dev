#! /usr/bin/make -f

### This script is the analogue to process_fastq_pairedend.sh specific for ATAC-seq data.
### First, ATAC data needs to be treated as single-end during bigwig computation with a shift of zero because we are interested in the fragment cut sites, not the fragment centers.
### Second, ATAC sam files need to be filtered for non-mitochondrial DNA because chrM is highly accessible and skews the distribution.

### This script only works for one biological replicate as bwa mem takes the two fastqs for the paired ends and puts them into a single .full.sam file.
### Use it for one biological replicate at a time and manually merge them later.

###
# Usage: ./process_fastq_atac.sh root=/project/wig/tobias/reg_evo/data assembly=mm10 sample=ATAC_FL_E10.5_R1 srr='SRR111111_A SRR111111_B'
###

# ------------------------------------------------------------------------------

GENOME_FASTA=$(root)/fasta/$(assembly).fa
GENOME=$(root)/assembly/$(assembly).sizes

DATA_DIR=.

FASTQDIR=$(DATA_DIR)/fastq
BAMDIR=$(DATA_DIR)/bam
BWDIR=$(DATA_DIR)/bigwig
SAMPLE=$(sample)

TARGETS = $(BAMDIR)/$(SAMPLE).sort.rmdup.bam.csi # $(BWDIR)/$(SAMPLE).cpm.bw

# ------------------------------------------------------------------------------

all: $(TARGETS)

# track files
# ------------------------------------------------------------------------------

# %$(BWDIR)/$(SAMPLE).cpm.bw: $(BAMDIR)/$(SAMPLE).sort.rmdup.bam
# 	bamToBigWig -v --paired-as-single-end --normalize-track=cpm --binning-method "mean overlap" --bin-size 2 $^ $@

# keep intermediate files and delete files whenever an error occurs
# ------------------------------------------------------------------------------

.DELETE_ON_ERROR:
.SECONDARY:

# process sra archives
# ------------------------------------------------------------------------------

# extract archives
%.fastq:
	fastq-dump -O $(dir $@) $(notdir $(basename $@))

# compress fastq
%.fastq.gz: %.fastq
	gzip -c $^

# compute alignment
$(BAMDIR)/$(SAMPLE).full.sam: $(addprefix $(FASTQDIR)/, $(addsuffix .fastq.gz, $(srr)))
	bwa mem -t 64 $(GENOME_FASTA) $^ > $@

# filter reads (minimum quality of 30) and convert to bam
%.full.bam: %.full.sam
	samtools view -hSq 30 $< | grep -v 'XA:Z:' | grep -v 'chrM' | samtools view -Sb - > $@
	samtools sort -n $@ > $@.tmp
	mv $@.tmp $@
	samtools fixmate $@ $@.tmp
	mv $@.tmp $@

# sort bam reads
%.sort.bam: %.full.bam
	samtools sort $< > $@

# remove PCR duplicates
%.sort.rmdup.bam: %.sort.bam
	samtools rmdup -s $< $@

# create index
%.sort.rmdup.bam.csi: %.sort.rmdup.bam
	samtools index -c $< # -c flag creates a CSI index. BAI has a chromosome size limit of 512 Mbp.

# ------------------------------------------------------------------------------

clean:
	$(RM) $(DATA_DIR)/*.bam.csi
	$(RM) $(DATA_DIR)/*.bam
	$(RM) $(DATA_DIR)/*.sam
    $(RM) $(DATA_DIR)/*.sai

distclean: clean
	$(RM) $(DATA_DIR)/*.fastq  
