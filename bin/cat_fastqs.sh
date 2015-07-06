#!/bin/bash

file=$1;
project_dir=$2;
output_dir=$3;
mkdir -p $output_dir
seq_run=$4;
fbase=${file##*Sample_}_${seq_run}
echo "Zipping $fbase"
echo "Fastqs: $project_dir"
echo "Datapath: $output_dir"
for R in 1 2; do
      zcat ${project_dir}/$file/*L[0-9][0-9][0-9]_R${R}_*.fastq.gz | gzip -n >${output_dir}/"${fbase}".$R.fastq.gz;
done

