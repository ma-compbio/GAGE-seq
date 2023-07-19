#!/bin/bash

library=$1
species=$2
genome=$3
path2scripts=$(dirname $0)
echo "path to scripts = ${path2scripts}"
# ----  for demultiplexing  ----
# barcodes
adaptor_file_bc1="adaptor_Hi-C_bc1"
adaptor_file_bc2="adaptor_Hi-C_bc2"
adaptor_file_l1="adaptor_Hi-C_l1"
adaptor_file_l2="adaptor_Hi-C_l2"
# simplify read names
function trimSeqName { awk -v OFS='\t' 'NR==1 { idx=1; for(c=0;c<4;idx++) if(substr($0,idx,1)==":") c++; } NR%4==1 { $0 = "@"substr($0,idx) } 1' ; }
barcode2keep=""
echo "barcode2keep = $barcode2keep"
# ----  for bwa alignment  ----
#genome_index="/work/magroup/tianming/Researches/sc-hic/data/genome/$genome/bwa/$genome"
genome_index="/work/magroup/tianming/Researches/sc-hic/data/genome/$genome/bwa_all/$genome"
echo "genome index = $genome_index"
bam="${library}_${species}_${genome}.bam"
echo "bam file = $bam"

# ----  demultiplexing  ----
python scripts/demultiplex.py \
  --barcode=$adaptor_file_bc2 --pos1="[]" --pos2="slice(0,8)"   --mode=1 \
  --barcode=$adaptor_file_l2  --pos1="[]" --pos2="slice(8,23)"  --mode=5 \
  --barcode=$adaptor_file_bc1 --pos1="[]" --pos2="slice(23,31)" --mode=1 \
  --infile1  <(gzip -cd "../raw/${library}_R1.fastq.gz" | trimSeqName) \
  --infile2  <(gzip -cd "../raw/${library}_R2.fastq.gz" | trimSeqName) \
  --outfile1 >(grep -a . > "${library}_demulti_R1.fastq") \
  --outfile2 >(grep -a . > "${library}_demulti_R2.fastq") \
  --barcode2keep="$barcode2keep"
# ----  alignment  ----
bwa mem -SP5M -t36 "$genome_index" \
  <(awk 'NR%2==0{print substr($0,13)} NR%2==1{print}' "${library}_demulti_R1.fastq") \
  <(awk 'NR%2==0{print substr($0,36)} NR%2==1{print}' "${library}_demulti_R2.fastq") | \
awk -v OFS='\t' '$0~/^@/ { print; next } { split($1,bc,"_"); $1 = bc[1]; print $0 "\tBC:Z:" bc[2] "," bc[4] }' \
| samtools view -hb@16 > "$bam"
# ----  pairs  ----
samtools view -h@16 "$bam" | $path2scripts/run_HiC_pair.sh $library $species $genome "complete" /dev/stdin
