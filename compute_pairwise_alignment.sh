#! /usr/bin/make -f

DATA_DIR = /project/wig/tobias/reg_evo/data
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

%.2bit:
ifeq ("$(wildcard $@)", "")
	@echo "Error: $@ does not exist"
	false
endif

%.sizes: %.2bit
	get_chromosome_sizes.R $<

$(FASTA_DIR)/%.fa: $(ASSEMBLY_DIR)/%.2bit
	twoBitToFa $< $@

$(LASTDB_DIR)/%.prj: $(FASTA_DIR)/%.fa
	lastdb -c $(patsubst %.prj,%,$@) $<

$(ALIGNMENT_DIR)/$(S1).$(S2).maf: $(LASTDB_DIR)/$(S1).prj $(FASTA_DIR)/$(S2).fa
	lastal $(patsubst %.prj,%,$(word 1,$^)) $(word 2,$^) > $@

%.psl: %.maf
	maf-convert psl $< > $@

%.axtChain: %.psl
	axtChain -psl -linearGap=loose $< $(ASSEMBLY_DIR)/$(S1).2bit $(ASSEMBLY_DIR)/$(S2).2bit $@

%.all.chain: %.axtChain
	chainMergeSort $< > $@

%.all.pre.chain: %.all.chain
	chainPreNet $< $(ASSEMBLY_DIR)/$(S1).sizes $(ASSEMBLY_DIR)/$(S2).sizes $@

%.noClass.net: %.all.pre.chain
	chainNet $< $(ASSEMBLY_DIR)/$(S1).sizes $(ASSEMBLY_DIR)/$(S2).sizes stdout /dev/null | netSyntenic stdin $@

%.net.axt: %.noClass.net %.all.pre.chain
	netToAxt $^ $(ASSEMBLY_DIR)/$(S1).2bit $(ASSEMBLY_DIR)/$(S2).2bit $@

$(OUTDIR)/%.net.axt.gz: $(ALIGNMENT_DIR)/%.net.axt
	gzip -c $< > $@

.SECONDARY: %.net.axt %.noClass.net %.all.pre.chain %.all.chain %.axtChain %.psl $(S1).$(S2).maf $(LASTDB_DIR)/%.prj %.sizes %.fa
