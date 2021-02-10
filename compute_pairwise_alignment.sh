#! /usr/bin/make -f

# usage: ./compute_pairwise_alignment.sh data_dir=/path/to/data S1=species1 S2=species2
# $data_dir must be the folder where all assembly, alignment and fasta folders are stored.

# DATA_DIR = $(realpath $(data_dir))
DATA_DIR = $(data_dir)
ASSEMBLY_DIR = $(DATA_DIR)/assembly
LASTDB_DIR = $(DATA_DIR)/alignment/lastdb
FASTA_DIR = $(DATA_DIR)/fasta
ALIGNMENT_DIR = $(DATA_DIR)/alignment
OUTDIR = $(ALIGNMENT_DIR)/axtNet

TARGETS = $(OUTDIR)/$(S1).$(S2).net.axt.gz

## -----------------------------------------------------------------------------

all: $(TARGETS)

## generate rules
## -----------------------------------------------------------------------------

$(ASSEMBLY_DIR)%.2bit: # $(FASTA_DIR)/%.fa
	get_assembly.sh $(ASSEMBLY_DIR) $(notdir $(patsubst %.2bit,%,$@))

# # ifeq ("$(wildcard $@)", "")
# # 	@echo "Error: $@ does not exist"
# # 	false
# # endif

%.sizes: %.2bit
	get_chromosome_sizes.R $<

# $(FASTA_DIR)/%.fa: $(ASSEMBLY_DIR)/%.2bit
# 	twoBitToFa $< $@

$(FASTA_DIR)/%.fa:
	echo "Download $(notdir $(patsubst %.fa,%,$@)) genome fasta file and save it to $@"

$(LASTDB_DIR)/%.prj: $(FASTA_DIR)/%.fa
	lastdb -c $(patsubst %.prj,%,$@) $<

$(ALIGNMENT_DIR)/$(S1).$(S2).maf: $(LASTDB_DIR)/$(S1).prj $(FASTA_DIR)/$(S2).fa
	lastal $(patsubst %.prj,%,$(word 1,$^)) $(word 2,$^) > $@

%.psl: %.maf
	maf-convert psl $< > $@

%.axtChain: %.psl $(ASSEMBLY_DIR)/$(S1).2bit $(ASSEMBLY_DIR)/$(S2).2bit
	axtChain -psl -linearGap=loose $< $(ASSEMBLY_DIR)/$(S1).2bit $(ASSEMBLY_DIR)/$(S2).2bit $@

%.all.chain: %.axtChain
	chainMergeSort $< > $@

%.all.pre.chain: %.all.chain $(ASSEMBLY_DIR)/$(S1).sizes $(ASSEMBLY_DIR)/$(S2).sizes
	chainPreNet $< $(ASSEMBLY_DIR)/$(S1).sizes $(ASSEMBLY_DIR)/$(S2).sizes $@

%.noClass.net: %.all.pre.chain $(ASSEMBLY_DIR)/$(S1).sizes $(ASSEMBLY_DIR)/$(S2).sizes
	chainNet $< $(ASSEMBLY_DIR)/$(S1).sizes $(ASSEMBLY_DIR)/$(S2).sizes stdout /dev/null | netSyntenic stdin $@

%.net.axt: %.noClass.net %.all.pre.chain $(ASSEMBLY_DIR)/$(S1).2bit $(ASSEMBLY_DIR)/$(S2).2bit
	netToAxt $(word 1,$^) $(word 2,$^) $(ASSEMBLY_DIR)/$(S1).2bit $(ASSEMBLY_DIR)/$(S2).2bit $@

$(OUTDIR)/%.net.axt.gz: $(ALIGNMENT_DIR)/%.net.axt
	gzip -c $< > $@

# .PRECIOUS: %.net.axt %.noClass.net %.all.pre.chain %.all.chain %.axtChain %.maf $(LASTDB_DIR)/%.prj %.sizes %.fa %.psl
