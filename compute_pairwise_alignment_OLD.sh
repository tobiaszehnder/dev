#! /usr/bin/make -f

DATA_DIR = /project/wig/tobias/reg_evo/data
ASSEMBLY_DIR = $(DATA_DIR)/assembly
FASTA_DIR = $(DATA_DIR)/fasta
OUTDIR = $(DATA_DIR)/alignment/axtNet

SPECIES = hg38 mm10 #galGal5 #xenTro9 danRer10 lepOcu1 calMil1 oryLat2 gasAcu1 tetNig2 fr3 latCha1 cteIde cypCar carAur
TARGETS = $(foreach S1,$(SPECIES),$(foreach S2,$(filter-out $(S1),$(SPECIES)), $(OUTDIR)/$(S1).$(S2).net.axt.gz))

## -----------------------------------------------------------------------------

all: $(TARGETS)

## generate rules
## -----------------------------------------------------------------------------

%.sizes: %.2bit
	get_chromosome_sizes.R $<
	touch $@

%.fa : $(ASSEMBLY_DIR)/%.2bit
	twoBitToFa $< > $@
	touch $@

lastdb/%.prj: $(FASTA_DIR)/%.fa
	lastdb -c $(patsubst %.prj,%,$@) $<
	touch $@

define rule =
$(1).$(2).maf: lastdb/$(1).prj $(FASTA_DIR)/$(2).fa
	lastal $$(patsubst %.prj,%,$$(word 1,$$^)) $$(word 2,$$^) > $$@
	touch $$@

$(1).$(2).all.pre.chain: $(1).$(2).all.chain
	chainPreNet $$< $(ASSEMBLY_DIR)/$(1).sizes $(ASSEMBLY_DIR)/$(2).sizes $$@
	touch $$@

$(1).$(2).noClass.net: $(1).$(2).all.pre.chain
	chainNet $$< $(ASSEMBLY_DIR)/$(1).sizes $(ASSEMBLY_DIR)/$(2).sizes stdout /dev/null | netSyntenic stdin $$@
	touch $$@

$(OUTDIR)/$(1).$(2).net.axt: $(1).$(2).noClass.net $(1).$(2).all.pre.chain
	netToAxt $$^ $(ASSEMBLY_DIR)/$(1).2bit $(ASSEMBLY_DIR)/$(2).2bit $$@
	touch $$@

$(OUTDIR)/$(1).$(2).net.axt.gz: $(OUTDIR)/$(1).$(2).net.axt
	gzip $$<
	touch $$@
endef

.INTERMEDIATE: $(foreach S1,$(SPECIES),$(foreach S2,$(filter-out $(S1),$(SPECIES)), $(OUTDIR)/$(S1).$(S2).net.axt))

$(foreach S1,$(SPECIES),$(foreach S2,$(filter-out $(S1),$(SPECIES)),$(eval $(call rule,$(S1),$(S2)))))

%.psl: %.maf
	maf-convert psl $< > $@
	touch $@

%.axtChain: %.psl
	axtChain -psl -linearGap=loose $< $(1).2bit $(2).2bit $@
	touch $@

%.all.chain: %.axtChain
	chainMergeSort $< > $@
	touch $@

# %.net.axt.gz: %.net.axt
# 	gzip $<
