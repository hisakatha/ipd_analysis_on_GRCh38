SHELL := bash
.DELETE_ON_ERROR:
.SECONDEXPANSION:
all: $$(targets)
.PHONY: all

DSLOP := 20
FA := /glusterfs/hisakatha/ensembl_GRCh38.p13/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
FAI := $(FA:=.fai)
# https://github.com/hisakatha/collect_ipd_per_kmer
COLLECT_IPD_PER_KMER := ../../../svr_ipd/collect_ipd_per_kmer/collect_ipd

fst_p6 := load_and_save.shi2016_P6C4.R.done
$(fst_p6):
	Rscript load_and_save.shi2016_P6C4.R && touch $@

targets += $(fst_p6)

deep_kinetics_p6 := deep_kinetics_region.shi2016_P6C4.bed
DEEP_KINETICS := $(deep_kinetics_p6)

$(deep_kinetics_p6): $(fst_p6)
	Rscript convert_kinetics_deep_region_into_bed.shi2016_P6C4.R

DEEP_KINETICS_MERGED_SORTED := $(DEEP_KINETICS:=.merged.sorted)
$(DEEP_KINETICS_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS_SLOP := $(DEEP_KINETICS_MERGED_SORTED:=.slop$(DSLOP).bed)
$(DEEP_KINETICS_SLOP): %.slop$(DSLOP).bed: %
	bedtools slop -i $< -g $(FAI) -b $(DSLOP) > $@ && touch $@
DEEP_KINETICS_SLOP_FA := $(DEEP_KINETICS_SLOP:=.fa)
$(DEEP_KINETICS_SLOP_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) > $@ && touch $@
# For counting the size of the regions
DEEP_KINETICS_SLOP_MERGED := $(DEEP_KINETICS_SLOP:=.merged)
$(DEEP_KINETICS_SLOP_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 4,5,6 -o distinct,mean,distinct > $@ && touch $@
targets += $(DEEP_KINETICS_SLOP_FA) $(DEEP_KINETICS_SLOP_MERGED)

m6A_gff_p6 := ipd_summary.m6A_cov25or12.shi2016_P6C4.gff
$(m6A_gff_p6): ../shi2016_P6C4/ipd_summary.m6A_cov25or12.gff
	ln -s $< $@
m6A_context_p6 := ipd_summary.m6A_cov25or12.shi2016_P6C4.gff.fa
$(m6A_context_p6): ../shi2016_P6C4/ipd_summary.m6A_cov25or12.gff.fa
	ln -s $< $@
targets += $(m6A_gff_p6) $(m6A_context_p6)

MOD := $(m6A_gff_p6)
MSLOP := 20

MOD_SLOP := $(MOD:.gff=.slop$(MSLOP).gff)
$(MOD_SLOP): %.slop$(MSLOP).gff: %.gff
	bedtools slop -header -i $< -g $(FAI) -b $(MSLOP) > $@
MOD_FULL := $(MOD_SLOP:.gff=.fullLength.gff)
$(MOD_FULL): %.fullLength.gff: %.gff
	cat $< | grep -v '^#' | awk -v EXPECTED=$$(($(MSLOP) * 2 + 1)) '{if($$5 - $$4 + 1 == EXPECTED){print}}' > $@
MOD_FULL_FA := $(MOD_FULL:=.fa)
$(MOD_FULL_FA): %.fullLength.gff.fa: %.fullLength.gff
	cat $< | bedtools getfasta -bed - -fi $(FA) -s > $@
MOD_FULL_MERGED := $(MOD_FULL:=.merged)
$(MOD_FULL_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 3,6,7 -o distinct,mean,distinct > $@

targets += $(MOD_FULL_FA) $(MOD_FULL_MERGED)

monomer_ipd_stats.shi2016_P6C4.csv: ../shi2016_P6C4/output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5
	$(COLLECT_IPD_PER_KMER) -k 1 -l 0 -t 25 -c ACGTN -o $@ $^

targets += monomer_ipd_stats.shi2016_P6C4.csv

#monomer_ipd_stats2.shi2016_P6C4.csv: ../shi2016_P6C4/output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5
#	$(COLLECT_IPD_PER_KMER) -k 1 -l 10 -t 25 -c ACGTN -o $@ $^
#targets += monomer_ipd_stats2.shi2016_P6C4.csv
