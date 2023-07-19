#!/bin/bash

library=$1
species=$2
genome=$3
walk_policy=$4
bam=$5
if [[ $bam == "bam" ]]; then
  bam="../alignment/${library}_${species}_${genome}.bam"
fi
echo "bam file = $bam"
path2scripts=$(dirname $0)
echo "path to scripts = ${path2scripts}"
# ----  for pairtools  ----
assembly="hg38_mm10"
echo "assembly = $assembly"
chrom_size_file="/work/magroup/tianming/Researches/sc-hic/data/genome/hg38_mm10/hg38_mm10.chrom.sizes"
echo "chrom size file = $chrom_size_file"
pairs="../alignment_${walk_policy}/${library}_${species}_${genome}.pairs.gz"
# ----  pairtools  ----
# run pairtools
# simplify its output
# flip and mark reads
function simplify {
  awk -v OFS='\t' '{
    NF = 11
    p = match($9,"BC:Z:")+5; $11 = substr($9,p,match(substr($9,p),"\x19")-1)
    $9="."; $10="."
    print
  }'
}
function flip {
  awk -v OFS='\t' -F'[\t]' '
    BEGIN { chra2i["!"] = 0 } ARGIND==1 { chra2i[$1] = FNR ; next }
    { if ( chra2i[$2] > chra2i[$4] || ( $2 == $4 && $3 > $5 ) ) {
      print $1 OFS $4 OFS $5 OFS $2 OFS $3 OFS $7 OFS $6 OFS substr($8,2,1) substr($8,1,1) OFS $10 OFS $9 OFS $11 OFS "F"
    } else { print $0 OFS "O" } }
  ' $chrom_size_file -
}
function keep_walk {
  awk -v OFS='\t' -F'[\t]' '
    {
      if (qname == $1) { cnt += 1; print last_line }
      else { if ( cnt > 1 ) { print last_line } cnt = 1 }
      qname = $1; last_line = $0
    }
    END { if ( cnt > 1 ) print last_line }
  '
}
if [[ $walk_policy == "wowalk" ]]; then
  function f { pairtools parse --walks-policy mask \
    --no-flip --min-mapq=10 --assembly=$assembly --chroms-path=$chrom_size_file | \
    grep -v '^#' | simplify | flip | grep '' --line-buffered ; }
elif [[ $walk_policy == "linear" ]]; then
  function f { pairtools parse --walks-policy all \
    --no-flip --min-mapq=10 --assembly=$assembly --chroms-path=$chrom_size_file | \
    grep -v '^#' | simplify | flip | grep '' --line-buffered ; }
elif [[ $walk_policy == "complete" ]]; then
  function f { pairtools parse --walks-policy all \
    --no-flip --min-mapq=10 --assembly=$assembly --chroms-path=$chrom_size_file | \
    grep -v '^#' | simplify | python $path2scripts/rescue_walk.py | flip | grep '' --line-buffered ; }
elif [[ $walk_policy == "walkonly" ]]; then
  function f { pairtools parse --walks-policy all \
    --no-flip --min-mapq=10 --assembly=$assembly --chroms-path=$chrom_size_file | \
    grep -v '^#' | simplify | keep_walk | python $path2scripts/rescue_walk.py | flip | grep '' --line-buffered ; }
else echo "Not Implemented"; fi
if [[ $bam = "/dev/stdin" ]]; then
  f | gzip -6 > "$pairs"
else
  samtools view -h@8 "$bam" | f | gzip -6 > "$pairs"
fi
