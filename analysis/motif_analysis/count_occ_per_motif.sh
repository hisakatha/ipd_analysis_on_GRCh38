#!/usr/bin/env bash
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
all_bg=$(cat ../../deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.merged | awk '{sum += $3 - $2}END{print sum}')
all_mod=$(cat ../../ipd_summary.m6A_cov25or12.$ID.slop20.fullLength.gff.merged | awk '{sum += $3 - $2}END{print sum}')
bg_majorchr=$(cat motif_occ.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | (grep -e "$major_pattern" || [[ $? == 1 ]]) | wc -l)
mod_majorchr=$(cat motif_occ.ipd_summary.m6A_cov25or12.$ID.slop20.fullLength.gff.fa.bed.merged_occ | (grep -e "$major_pattern" || [[ $? == 1 ]]) | wc -l)
all_bg_majorchr=$(cat ../../deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.merged | (grep -e "$major_pattern" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
all_mod_majorchr=$(cat ../../ipd_summary.m6A_cov25or12.$ID.slop20.fullLength.gff.merged | (grep -e "$major_pattern" || [[ $? == 1 ]]) | awk '{sum += $3 - $2}END{print sum}')
echo "# all_bg_allchr = $all_bg, all_mod_allchr = $all_mod"
echo "# all_bg_majorchr = $all_bg_majorchr, all_mod_majorchr = $all_mod_majorchr"

mod_ratio=$(awk -v feature=$mod -v bg=$bg 'BEGIN{printf "%.3g\n", feature / bg}')
mod_p=$(Rscript -e "cat(sprintf('%.3g\n', chisq.test(matrix(c($all_mod, $all_bg, $mod, $bg), 2))\$p.value))")
mod_adjusted_p=$(awk -v pvalue=$mod_p -v test_num=$mod_test_num 'BEGIN{printf "%.3g\n", pvalue * test_num}')
echo "sample_id,chr,motif,feature_type,occ_feature,occ_background,ratio,chisq.pvalue,adjusted.pvalue"
echo "$ID,$output_name_all,$MOTIF,mod_6mA,$mod,$bg,$mod_ratio,$mod_p,$mod_adjusted_p"

mod_ratio=$(awk -v feature=$mod_majorchr -v bg=$bg_majorchr 'BEGIN{printf "%.3g\n", feature / bg}')
mod_p=$(Rscript -e "cat(sprintf('%.3g\n', chisq.test(matrix(c($all_mod_majorchr, $all_bg_majorchr, $mod_majorchr, $bg_majorchr), 2))\$p.value))")
mod_adjusted_p=$(awk -v pvalue=$mod_p -v test_num=$mod_test_num 'BEGIN{printf "%.3g\n", pvalue * test_num}')
#echo "sample_id,chr,motif,feature_type,occ_feature,occ_background,ratio,chisq.pvalue,adjusted.pvalue"
echo "$ID,$output_name_major,$MOTIF,mod_6mA,$mod_majorchr,$bg_majorchr,$mod_ratio,$mod_p,$mod_adjusted_p"
