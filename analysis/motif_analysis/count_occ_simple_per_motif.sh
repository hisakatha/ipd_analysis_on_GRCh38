#!/usr/bin/env bash
MOTIF=$(cat PATTERN)
major_pattern="^[1-9XY]\|^MT"

output_name_all="all_chr"
output_name_major="major_chr"

# header
echo "motif,sample_id,feature,chromosome_set,num_occ"

ref_id="GRCh38.p13"
ref_feature_name="whole_genome"
ref_all=$(cat motif_occ.GRCh38.p13.bed.merged_occ | wc -l)
ref_major=$(cat motif_occ.GRCh38.p13.bed.merged_occ | (grep -e "$major_pattern" || [[ $? == 1 ]]) | wc -l)
echo "$MOTIF,$ref_id,$ref_feature_name,$output_name_all,$ref_all"
echo "$MOTIF,$ref_id,$ref_feature_name,$output_name_major,$ref_major"

sample_ids="shi2016_P6C4"
for ID in $sample_ids
do
    deep_feature_name="deep_kinetics_region"
    deep_all=$(cat motif_occ.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | wc -l)
    deep_major=$(cat motif_occ.deep_kinetics_region.$ID.bed.merged.sorted.slop20.bed.fa.bed.merged_occ | (grep -e "$major_pattern" || [[ $? == 1 ]]) | wc -l)
    echo "$MOTIF,$ID,$deep_feature_name,$output_name_all,$deep_all"
    echo "$MOTIF,$ID,$deep_feature_name,$output_name_major,$deep_major"

    mod_feature_name="m6A_cov25or12"
    mod_all=$(cat motif_occ.ipd_summary.m6A_cov25or12.$ID.slop20.fullLength.gff.fa.bed.merged_occ | wc -l)
    mod_major=$(cat motif_occ.ipd_summary.m6A_cov25or12.$ID.slop20.fullLength.gff.fa.bed.merged_occ | (grep -e "$major_pattern" || [[ $? == 1 ]]) | wc -l)
    echo "$MOTIF,$ID,$mod_feature_name,$output_name_all,$mod_all"
    echo "$MOTIF,$ID,$mod_feature_name,$output_name_major,$mod_major"
done
