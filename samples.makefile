SHELL := bash
.DELETE_ON_ERROR:
.SECONDEXPANSION:
all: $$(targets)
.PHONY: all

ipd_summary_gff := ipd_summary.gff
ipd_summary_gff_src1 := output/tasks/pbcoretools.tasks.gather_gff-1/file.gff
ipd_summary_gff_src2 := output/tasks/kinetics_tools.tasks.ipd_summary-0/basemods.gff
$(ipd_summary_gff):
	if [[ -e $(ipd_summary_gff_src1) ]]; then ln -s $(ipd_summary_gff_src1) $(ipd_summary_gff); \
		elif [[ -e $(ipd_summary_gff_src2) ]]; then ln -s $(ipd_summary_gff_src2) $(ipd_summary_gff); \
		else exit 1; fi

ipd_summary_gff_m6A_high_cov := $(ipd_summary_gff:.gff=.m6A_cov25or12.gff)
$(ipd_summary_gff_m6A_high_cov): $(ipd_summary_gff)
	cat $< | ../extract_m6A_cov25or12.sh > $@
ipd_summary_gff_m6A_high_cov_fa := $(ipd_summary_gff_m6A_high_cov:=.fa)
$(ipd_summary_gff_m6A_high_cov_fa): $(ipd_summary_gff_m6A_high_cov)
	cat $< | grep -v "^#" | awk '{context = gensub(/.*context=([a-zA-Z]+).*/, "\\1", 1, $$9); print ">"$$1":"$$4"-"$$5"\n"context}' > $@
targets += $(ipd_summary_gff_m6A_high_cov_fa)
