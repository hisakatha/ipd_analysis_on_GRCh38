SHELL := bash

PATTERN := $(shell cat PATTERN)
ifndef PATTERN
$(error PATTERN is undefined)
endif

### Variables you have to set
REF := /glusterfs/hisakatha/ensembl_GRCh38.p13/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
FAI := $(REF:=.fai)
###

occ_prefix := motif_occ

motif_ipds := motif_ipd.shi2016_P6C4.csv
motif_wga_ipdratios := motif_ipdratio.shi2016_P6C4.csv
motif_kinetics_csv := $(motif_ipds)
motif_kinetics_csv_celegans := $(motif_kinetics_csv:.csv=.c_elegans.csv)
motif_kinetics_csv_celegans += $(motif_wga_ipdratios:.csv=.c_elegans.csv)
motif_kinetics_csv_ecoli := $(motif_kinetics_csv_celegans:.c_elegans.csv=.e_coli.csv)
motif_modelPrediction := motif_modelPrediction.shi2016_P6C4.csv
motif_modelPrediction_celegans := $(motif_modelPrediction:.csv=.c_elegans.csv)
motif_modelPrediction_ecoli := $(motif_modelPrediction:.csv=.e_coli.csv)
motif_summary_input := $(motif_kinetics_csv_celegans) $(motif_kinetics_csv_ecoli) $(motif_modelPrediction_celegans) $(motif_modelPrediction_ecoli)
motif_summary_target := $(motif_summary_input:.csv=.summary.csv)
motif_kinetics_pdf := plot_motif_kinetics.pdf
motif_kinetics_stderr_pdf := plot_motif_kinetics.stderr.pdf
motif_kinetics_stderr_widey_pdf := plot_motif_kinetics.stderr.wide_y.pdf
cor_pdf := plot_kinetics_correlation.pdf
cor_pdf2 := plot_ipdratio_correlation.pdf
cor_pdf_celegans := plot_kinetics_correlation.c_elegans.pdf
#cor_pdf2_celegans := plot_ipdratio_correlation.c_elegans.pdf
cor_pdf_celegans_ecoli := plot_kinetics_correlation.celegans_ecoli.pdf
cor_pdf_celegans_ecoli_subset1 := plot_kinetics_correlation.celegans_ecoli_subset1.pdf
cor_pdf_celegans_ecoli_subset1_positive_strand := plot_kinetics_correlation.celegans_ecoli_subset1_positive_strand.pdf

.SECONDEXPANSION:
#target := $(motif_kinetics_pdf) $(motif_kinetics_stderr_pdf) $(motif_kinetics_stderr_widey_pdf) $(cor_pdf) $(motif_wga_ipdratios) $(cor_pdf2)
#target += $(motif_kinetics_csv_celegans) $(cor_pdf_celegans)
#target += $(motif_modelPrediction) $(motif_modelPrediction_celegans)
#target += $(motif_kinetics_csv_ecoli) $(motif_modelPrediction_ecoli) $(cor_pdf_celegans_ecoli) $(cor_pdf_celegans_ecoli_subset1) $(cor_pdf_celegans_ecoli_subset1_positive_strand)
#target += $(motif_summary_target)
all: $$(target)

.DELETE_ON_ERROR:

# NOTE: It is important to use non-merged files after slopping/extension (*.slop20.bed.fa),
# rather than merged files (*.slop20.merged.bed.fa) for locating motifs near feature loci.
# I want to exclude occurrences of motifs at the edges of extended regions:
# --xxxxxXxxxxx------------
# ------------xxxxxXxxxxx--
# ---YYY-----YYY-----------
# X: feature locus
# x: extended region
# YYY: target motif
# The first motif should be captured, but the second should be ignored.

# TODO?: replace deep_region_cov25, which is not high mapq, into deep_kinetics_region
# pattern occurrence for background
#occ_background := $(occ_prefix).shi2016_P6C4.bed $(occ_prefix).GRCh38.p13.bed
occ_background := $(occ_prefix).GRCh38.p13.bed
$(occ_prefix).GRCh38.p13.bed: $(REF)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
# Root directory for sample analyses
DIR := ../../..
FILE := mapped.alignmentset.merged.high_mapq.bam.cov25.slop20.bed.fa
$(occ_prefix).shi2016_P6C4.bed: $(DIR)/shi2016_P6C4/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@

deep_kinetics_beds := $(wildcard ../../deep_kinetics_region.*.slop20.bed.fa)
occ_deep_kinetics := $(deep_kinetics_beds:../../%=$(occ_prefix).%.bed)
$(occ_deep_kinetics): $(occ_prefix).%.bed: ../../%
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@

# pattern occurrence for features (such as DNA modification and high IPD)
occ_feature_fa := $(wildcard ../../*.slop20.fullLength.gff.fa)
occ_feature := $(occ_feature_fa:../../%=$(occ_prefix).%.bed)
$(occ_feature): $(occ_prefix).%.bed: ../../%
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@

occ := $(occ_background) $(occ_feature) $(occ_deep_kinetics)
occ := $(occ_background) $(occ_feature) $(occ_deep_kinetics)
merged_occ := $(occ:=.merged_occ)
$(merged_occ): %.merged_occ: %
	cat $< | ../merge_occurrence.sh > $@ && touch $@

target += $(merged_occ)

# These lines require a lot of RAM
$(motif_ipds): get_motif_ipd.R.log
	touch $@
get_motif_ipd.R.log: $(occ_prefix).deep_kinetics_region.shi2016_P6C4.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipd.v2.R $(PATTERN) &> $@ && touch $@
$(motif_wga_ipdratios): get_motif_wga_ipdratio.R.log
	touch $@
get_motif_wga_ipdratio.R.log: $(occ_prefix).deep_kinetics_region.shi2016_P6C4.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_wga_ipdratio.R $(PATTERN) &> $@ && touch $@
$(motif_ipdratios): get_motif_ipdratio.R.log
	touch $@

$(motif_kinetics_csv_celegans): %.c_elegans.csv: %.csv
	sample=$$(echo $< | sed -E "s/motif[_0-9a-zA-Z]+\.([_0-9a-zA-Z]+)\.csv/\1/"); \
	Rscript ../extract_celegans_data.R $< $(occ_prefix).deep_kinetics_region.$${sample}.bed*.slop20.bed.fa.bed.merged_occ $@

$(motif_kinetics_csv_ecoli): %.e_coli.csv: %.csv
	sample=$$(echo $< | sed -E "s/motif[_0-9a-zA-Z]+\.([_0-9a-zA-Z]+)\.csv/\1/"); \
	Rscript ../extract_ecoli_data.R $< $(occ_prefix).deep_kinetics_region.$${sample}.bed*.slop20.bed.fa.bed.merged_occ $@

$(motif_summary_target): %.summary.csv: %.csv
	Rscript ../summarize_motif_kinetics.R $< $@

$(motif_kinetics_pdf): $(motif_kinetics_csv)
	Rscript ../plot_motif_kinetics.R $(PATTERN) && touch $@

$(motif_kinetics_stderr_pdf): $(motif_kinetics_csv)
	Rscript ../plot_motif_kinetics.stderr.R $(PATTERN) && touch $@

$(motif_kinetics_stderr_widey_pdf): $(motif_kinetics_csv)
	Rscript ../plot_motif_kinetics.stderr.wide_y.R $(PATTERN) && touch $@

$(cor_pdf): $(motif_kinetics_csv)
	Rscript ../plot_kinetics_correlation.R

$(cor_pdf_celegans): $(motif_kinetics_csv_celegans)
	Rscript ../plot_kinetics_correlation.c_elegans.R

$(cor_pdf2): $(motif_wga_ipdratios) $(motif_hawaiian_ipdratios)
	Rscript ../plot_ipdratio_correlation.R

$(cor_pdf_celegans_ecoli): $(motif_kinetics_csv_celegans) $(motif_kinetics_csv_ecoli)
	Rscript ../plot_kinetics_correlation.celegans_ecoli.R

$(cor_pdf_celegans_ecoli_subset1): $(motif_kinetics_csv_celegans) $(motif_kinetics_csv_ecoli)
	Rscript ../plot_kinetics_correlation.celegans_ecoli_subset1.R

$(cor_pdf_celegans_ecoli_subset1_positive_strand): $(motif_kinetics_csv_celegans) $(motif_kinetics_csv_ecoli)
	Rscript ../plot_kinetics_correlation.celegans_ecoli_subset1_positive_strand.R

$(motif_modelPrediction): motif_modelPrediction.%.csv: motif_ipd.%.csv motif_ipdratio.%.csv
	Rscript ../get_motif_modelPrediction_from_ipdRatio.R $< $(word 2,$^) $@ && touch $@
$(motif_modelPrediction_celegans): motif_modelPrediction.%.c_elegans.csv: motif_ipd.%.c_elegans.csv motif_ipdratio.%.c_elegans.csv
	Rscript ../get_motif_modelPrediction_from_ipdRatio.R $< $(word 2,$^) $@ && touch $@
$(motif_modelPrediction_ecoli): motif_modelPrediction.%.e_coli.csv: motif_ipd.%.e_coli.csv motif_ipdratio.%.e_coli.csv
	Rscript ../get_motif_modelPrediction_from_ipdRatio.R $< $(word 2,$^) $@ && touch $@

count_occ_simple := count_occ_simple_per_motif.csv
$(count_occ_simple): $(merged_occ)
	../count_occ_simple_per_motif.sh > $@

target += $(count_occ_simple)
