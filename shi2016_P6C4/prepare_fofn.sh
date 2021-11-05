#!/usr/bin/env bash
cat /glusterfs/hisakatha/shi2016_human_HX1_data/kit_version.csv | grep P6 | cut -f1 -d, | sed -E -e 's@(.*)@/glusterfs/hisakatha/shi2016_human_HX1_data/\1/@' | xargs -n1 -I{} find {} -name '*.subreads.bam'
