#!/bin/bash

library=$1
species=$2
genome=$3
path2scripts=$(dirname $0)
echo "path to scripts = ${path2scripts}"
bam="${library}_${species}_${genome}.bam"
nodup="${library}_${species}_${genome}_nodup.bam"

samtools view -F0x104 -uh@8 "$bam" | samtools sort -@16 -m5G -t BC -O sam -T /tmp | \
  python $path2scripts/dedup_cDNA.py --num_shift=5 --num_mismatches_UMI=1 | samtools view -hb@8 > "$nodup"
