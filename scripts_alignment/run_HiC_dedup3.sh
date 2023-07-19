#!/bin/bash

library=$1
species=$2
genome=$3
path2scripts=$(dirname $0)
echo "path to scripts = ${path2scripts}"
chrom_size_file="/work/magroup/tianming/Researches/sc-hic/data/genome/hg38_mm10/hg38_mm10.chrom.sizes"
echo "chrom size file = $chrom_size_file"
pairs="${library}_${species}_${genome}.pairs.gz"
pairs_nodup="${library}_${species}_${genome}_nodup.pairs.gz"
pairs_nodup_wo1k="${library}_${species}_${genome}_nodup_wo1k.pairs.gz"

pigz -cdp8 "$pairs" | grep -P 'UU|UR|RU' | cut -f2-7,11- | scripts/dedup_Hi-C_cpp 500 500 9216 47 | \
  grep -vP '\tDD\t' | ( tee >(pigz -6p8 > "${pairs_nodup}") ) | \
  awk '$2!=$4 || $5-$3>=1000' | pigz -6p8 > "${pairs_nodup_wo1k}"
