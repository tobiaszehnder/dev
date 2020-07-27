#! /usr/bin/make -f

### USAGE: ./call_CNEs.sh nthreads=<int>
ifndef nthreads
$(error nthreads is not set)
endif
ifndef cneType
$(error cneType is not set. It must be either "merged" or "final")
endif

DATA_DIR = /project/wig/tobias/reg_evo/data
OUTDIR = $(DATA_DIR)/CNEs/CNEr
SEM_ID = "call_CEs_$(shell hostname)"

SPECIES = hg38 mm10 galGal5 xenTro9 danRer10 lepOcu1 calMil1 oryLat2 gasAcu1 tetNig2 fr3 latCha1 cteIde1 cypCar1 carAur01
# SPECIES = hg38 mm10 galGal5 monDom5 ornAna3 xenTro9 danRer10 lepOcu1 calMil1
# SPECIES = danRer10 mm10 calMil1
# SPECIES = hg38 mm10 galGal5 monDom5 ornAna3 xenTro9 danRer10 lepOcu1 calMil1 oryLat2 gasAcu1 tetNig2 fr3 latCha1 cteIde1 cypCar1 carAur01
TARGETS = $(foreach S1,$(SPECIES),$(foreach S2,$(filter-out $(S1),$(SPECIES)), $(OUTDIR)/ce_$(cneType)_$(S1)_$(S2)_35_50.bed))

## -----------------------------------------------------------------------------

all: $(TARGETS)

## generate rules
## -----------------------------------------------------------------------------

define rule =
$(OUTDIR)/ce_$(cneType)_$(S1)_$(S2)_35_50.bed:
	sem --id $(SEM_ID) -j$(nthreads) call_conserved_elements.R $(1) $(2) CE $(cneType)
endef

$(foreach S1,$(SPECIES),$(foreach S2,$(filter-out $(S1),$(SPECIES)),$(eval $(call rule,$(S1),$(S2)))))
