#! /usr/bin/make -f

### THIS SCRIPT PROCESSES SINGLE-END READS (use bwa mem for paired-end) ###

###
# Usage: ./process_fastq_singleend.sh sample=H3KXyZ_stageXY srr='SRR111111 SRR222222 ...'
###

ROOT=/project/wig/tobias/reg_evo/data
GENOME_FASTA=$(ROOT)/fasta/$(assembly).fa
GENOME=$(ROOT)/assembly/$(assembly).sizes

DATA_DIR=.

BAMDIR=$(DATA_DIR)/bam
BWDIR=$(DATA_DIR)/bigwig
SAMPLE=$(sample)

TARGETS = $(BWDIR)/$(SAMPLE).cpm.200bpBinsMeanOverlap.bw $(BAMDIR)/$(SAMPLE).merged.bam.bai

# ------------------------------------------------------------------------------

all: $(TARGETS)


# track files
# ------------------------------------------------------------------------------
$(BWDIR)/$(SAMPLE).cpm.200bpBinsMeanOverlap.bw: $(BAMDIR)/$(SAMPLE).merged.bam
	bamToBigWig -v --normalize-track=cpm --binning-method "mean overlap" --bin-size 200 $^ $@

$(BAMDIR)/$(SAMPLE).merged.bam: $(addprefix $(BAMDIR)/, $(addsuffix .sort.rmdup.bam, $(srr)))
	samtools merge $@ $^

# keep intermediate files and delete files whenever an error occurs
# ------------------------------------------------------------------------------

.DELETE_ON_ERROR:
.SECONDARY:

# process sra archives
# ------------------------------------------------------------------------------

# extract archives
%.fastq:
	fastq-dump -O $(dir $@) $(notdir $(basename $@))

# compute alignment
%.sai: %.fastq
	# bwa aln -t 64 -k 2 -n 3 -0 $(GENOME_FASTA) $^ > $@
	bwa aln -t 64 -k 2 -n 3 $(patsubst %.fa,%,$(GENOME_FASTA)) $^ > $@

# convert to sam
%.full.sam: %.sai %.fastq
	# bwa samse $(GENOME_FASTA) $^ > $@
	bwa samse $(patsubst %.fa,%,$(GENOME_FASTA)) $^ > $@

# filter reads (minimum quality of 30) and convert to bam
%.full.bam: %.full.sam
	samtools view -Sb $< > $@ # samtools view -Sb -q 30 $< > $@

# sort bam reads
%.sort.bam: %.full.bam
	samtools sort $< > $@

# remove PCR duplicates
%.sort.rmdup.bam: %.sort.bam
	samtools rmdup -s $< $@

# create index
%.bam.bai: %.bam
	samtools index $<

# ------------------------------------------------------------------------------

clean:
	$(RM) $(DATA_DIR)/*.bam.bai
	$(RM) $(DATA_DIR)/*.bam
	$(RM) $(DATA_DIR)/*.sam
	$(RM) $(DATA_DIR)/*.sai

distclean: clean
	$(RM) $(DATA_DIR)/*.fastq
