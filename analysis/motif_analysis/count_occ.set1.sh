#!/usr/bin/env bash

motifs_high_ipd="ACGCRTG ATCAGCTG GGN_4 RGTA AAAABT AGAGAGTA AGCTATAT AGGCAGGC ATGGGAYA ATTGTTAC CAAACTAC CAGYTG CCAATCAG CRACGAS DCGAGACC DGCTTC GAAGGATC GATATRGY GCGACCTA GCGCGCGC GGHGGY GTAGATCA GTATCGTA TGACGTCA TGGTGSA YGGAR"
#n_motifs_high_ipd=26

motifs_low_ipd="ACATMTGG CTGDAR TACTGTAG AATMAATA AGACGCAG CGGYTTGA CTKCAA GCGCGTCA GTGTGTGY TACCCCKA TACCTTGA TGACGTCA"
#n_motifs_low_ipd=12

#motifs_select1="GAGG AGAA"

csv_prefix="count_occ.set1"
pre_csv="$csv_prefix.pre_csv"

#motifs=$(echo $motifs_high_ipd $motifs_low_ipd $motifs_select1 | xargs -n1 | sort | uniq)
motifs=$(echo $motifs_high_ipd $motifs_low_ipd | xargs -n1 | sort | uniq)
n_motifs=$(echo $motifs | wc -w)
for motif in $motifs
do
    cd $motif || exit 1
    for sample in shi2016_P6C4
    do
        ../count_occ_per_motif.sh $sample $n_motifs || exit 1
    done
    cd ../
done > $pre_csv

(cat $pre_csv | grep sample_id | head -n1; cat $pre_csv | grep -v "#" | grep -v sample_id) > $csv_prefix.csv
