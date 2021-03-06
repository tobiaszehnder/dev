#! /usr/bin/make -f

### THIS SCRIPT PROCESSES PAIRED-END READS (use bwa aln -0 for single-end only) ###

### ONLY FOR ONE BIOLOGICAL REPLICATE ###

### This script only works for one biological replicate as bwa mem takes the two fastqs for the paired ends (R1 and R2) and puts them into a single .full.sam file.
### Use it for one biological replicate at a time and manually merge them later.

###
# Usage: ./process_fastq_pairedend.sh root=/project/wig/tobias/reg_evo/data assembly=mm10 sample=H3K27ac_FL_E10.5_R1 srr='SRR111111_R1 SRR111111_R2'
###

GENOME_FASTA=$(root)/fasta/$(assembly).fa
GENOME=$(root)/assembly/$(assembly).sizes

DATA_DIR=.

FASTQDIR=$(DATA_DIR)/fastq
BAMDIR=$(DATA_DIR)/bam
BWDIR=$(DATA_DIR)/bigwig
SAMPLE=$(sample)

TARGETS = $(BAMDIR)/$(SAMPLE).sort.rmdup.bam.csi

# ------------------------------------------------------------------------------

all: $(TARGETS)


# # track files
# # ------------------------------------------------------------------------------

%$(BWDIR)/$(SAMPLE).cpm.bw: $(BAMDIR)/$(SAMPLE).sort.rmdup.bam
	bamToBigWig -v --normalize-track=cpm --binning-method "mean overlap" --bin-size 2 $^ $@ 

# $(BWDIR)/$(SAMPLE).cpm.200bpBinsMeanOverlap.bw: $(BAMDIR)/$(SAMPLE).merged.bam
# 	bamToBigWig -v --normalize-track=cpm --binning-method "mean overlap" --bin-size 200 $^ $@

# $(BAMDIR)/$(SAMPLE).merged.bam: $(addprefix $(BAMDIR)/, $(addsuffix .sort.rmdup.bam, $(srr)))
# 	samtools merge $@ $^

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

# compute alignment from both ends (R1 and R2)
$(BAMDIR)/$(SAMPLE).full.sam: $(addprefix $(FASTQDIR)/, $(addsuffix .fastq.gz, $(srr)))
	bwa mem -t 64 $(GENOME_FASTA) $^ > $@

# filter reads (minimum quality of 30) and convert to bam
%.full.bam: %.full.sam
	samtools view -Sb -q 30 $< > $@

# sort bam reads
%.sort.bam: %.full.bam
	samtools sort $< > $@

# remove PCR duplicates
%.sort.rmdup.bam: %.sort.bam
	samtools rmdup $< $@

# create index
%.bam.csi: %.bam
	samtools index -c $< # -c flag creates a CSI index. BAI has a chromosome size limit of 512 Mbp.

# ------------------------------------------------------------------------------

clean:
	$(RM) $(DATA_DIR)/*.bam.csi
	$(RM) $(DATA_DIR)/*.bam
	$(RM) $(DATA_DIR)/*.sam
	$(RM) $(DATA_DIR)/*.sai

distclean: clean
	$(RM) $(DATA_DIR)/*.fastq
