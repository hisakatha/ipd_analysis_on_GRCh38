#!/usr/bin/env bash
# Variant of count_occ_per_motif.sh: binary region (feature region vs non-feature region) and Fisher's exact test
# TODO?: Fisher's exact test assumes fixed merginal counts, which is not the case. We need to find another test?
# This script tests the difference between the two sets in the following contingency table:
# (the number of occurrence of a motif in a feature region (high IPD or low IPD)) vs. (the number of occurrence of a motif in deep_kinetics_region), and
# (the number of bases in a feature region) vs. (the number of bases in deep_kinetics_region)
#
# deep_kinetics_region: region with a valid IPD count more than a threshold (default: count >= 25)

MOTIF=$(cat PATTERN)
# sample id (e.g. ab, cd)
ID=$1
major_pattern="^[1-9XY]\|^MT"
mod_test_num=$2
output_name_all="all_chr"
output_name_major="major_chr"

bg=$(cat motif_occ.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | wc -l)
mod=$(cat motif_occ.ipd_summary.m6A_cov25or12.$ID.slop20.fullLength.gff.fa.bed.merged_occ | wc -l)
nonmod=$((bg - mod))
# all: all the kinds of bases (motif N)
all_bg=$(cat ../../deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.merged | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_mod=$(cat ../../ipd_summary.m6A_cov25or12.$ID.slop20.fullLength.gff.merged | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_nonmod=$((all_bg - all_mod))
bg_majorchr=$(cat motif_occ.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | (grep -e "$major_pattern" || [[ $? == 1 ]]) | wc -l)
mod_majorchr=$(cat motif_occ.ipd_summary.m6A_cov25or12.$ID.slop20.fullLength.gff.fa.bed.merged_occ | (grep -e "$major_pattern" || [[ $? == 1 ]]) | wc -l)
nonmod_majorchr=$((bg_majorchr - mod_majorchr))
all_bg_majorchr=$(cat ../../deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.merged | (grep -e "$major_pattern" || [[ $? == 1 ]]) | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_mod_majorchr=$(cat ../../ipd_summary.m6A_cov25or12.$ID.slop20.fullLength.gff.merged | (grep -e "$major_pattern" || [[ $? == 1 ]]) | awk 'BEGIN{sum = 0}{sum += $3 - $2}END{print sum}')
all_nonmod_majorchr=$((all_bg_majorchr - all_mod_majorchr))
echo "# all_bg_allchr = $all_bg, all_mod_allchr = $all_mod"
echo "# all_bg_majorchr = $all_bg_majorchr, all_mod_majorchr = $all_mod_majorchr"

mod_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $mod / $bg))")
all_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $all_mod / $all_bg))")
odds_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $mod / $nonmod / $all_mod * $all_nonmod))")
mod_p=$(python3 -c "import scipy.stats as ss; oddsratio,pvalue = ss.fisher_exact([[$all_mod, $all_nonmod], [$mod, $nonmod]]); print(pvalue)")
mod_adjusted_p=$(Rscript -e "cat(sprintf('%.3g\n', $mod_p * $mod_test_num))")
echo "sample_id,chr,motif,feature_type,occ_feature,occ_nonfeature,occ_ratio,base_feature,base_nonfeature,base_ratio,odds_ratio,fisher.pvalue,adjusted.pvalue"
echo "$ID,$output_name_all,$MOTIF,mod_6mA,$mod,$nonmod,$mod_ratio,$all_mod,$all_nonmod,$all_ratio,$odds_ratio,$mod_p,$mod_adjusted_p"

mod_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $mod_majorchr / $bg_majorchr))")
all_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $all_mod_majorchr / $all_bg_majorchr))")
odds_ratio=$(Rscript -e "cat(sprintf('%.3g\n', $mod_majorchr / $nonmod_majorchr / $all_mod_majorchr * $all_nonmod_majorchr))")
mod_p=$(python3 -c "import scipy.stats as ss; oddsratio,pvalue = ss.fisher_exact([[$all_mod_majorchr, $all_nonmod_majorchr], [$mod_majorchr, $nonmod_majorchr]]); print(pvalue)")
mod_adjusted_p=$(Rscript -e "cat(sprintf('%.3g\n', $mod_p * $mod_test_num))")
echo "$ID,$output_name_major,$MOTIF,mod_6mA,$mod_majorchr,$nonmod_majorchr,$mod_ratio,$all_mod_majorchr,$all_nonmod_majorchr,$all_ratio,$odds_ratio,$mod_p,$mod_adjusted_p"
