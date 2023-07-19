#!/bin/bash

library=$1
species=$2
genome=$3
bam="${library}_${species}_${genome}.bam"
path2scripts=$(dirname $0)
echo "path to scripts = ${path2scripts}"
# ----  demultiplexing  ----
# barcodes
adaptor_file_bc1="adaptor_cDNA_bc1"
adaptor_file_bc2="adaptor_cDNA_bc2"
adaptor_file_l0="adaptor_cDNA_l0"
adaptor_file_l1="adaptor_cDNA_l1"
adaptor_file_l2="adaptor_cDNA_l2"
# simplify read names
function trimSeqName { awk -v OFS='\t' 'NR==1 { idx=1; for(c=0;c<4;idx++) if(substr($0,idx,1)==":") c++; } NR%4==1 { $0 = "@"substr($0,idx) } 1' ; }
barcode2keep=""
echo "barcode2keep = $barcode2keep"
# ----  for STAR alignment  ----
genomeDir="/work/magroup/tianming/Researches/sc-hic/data/genome/${genome}/STAR_anno"
echo $genomeDir
genomeLoad="LoadAndKeep"
#  genomeLoad="NoSharedMemory"

python $path2scripts/demultiplex.py \
  --barcode=$adaptor_file_bc2 --pos1="slice(10, 18)" --pos2="[]" --mode=1 \
  --barcode=$adaptor_file_l2  --pos1="slice(18, 48)" --pos2="[]" --mode=6 \
  --barcode=$adaptor_file_bc1 --pos1="slice(48, 56)" --pos2="[]" --mode=1 \
  --barcode=$adaptor_file_l1  --pos1="slice(56, 71)" --pos2="[]" --mode=6 \
  --infile1  <(gzip -cd "../raw/${library}_R1.fastq.gz" | trimSeqName) \
  --infile2  <(gzip -cd "../raw/${library}_R2.fastq.gz" | trimSeqName) \
  --outfile1 >(grep -a . > "${library}_demulti_R1.fastq") \
  --outfile2 >(grep -a . > "${library}_demulti_R2.fastq") \
  --barcode2keep="$barcode2keep" \
  2>"${file_delog[i]}"
STAR --runThreadN=16 --genomeLoad $genomeLoad --genomeDir="$genomeDir" --outSAMtype SAM --outStd SAM \
  --outSAMunmapped Within --outFileNamePrefix="${library}_${species}_${genome}_" --readFilesIn=<(
    paste -d'\0' "${library}_demulti_R2.fastq" <(awk 'NR%4==2{printf "_" substr($0,1,10) "\n\n\n\n"}' "${library}_demulti_R1.fastq")
) | awk -v OFS='\t' '{split($1,bcumi,"_"); $1=bcumi[1]; print $0 "\tBC:Z:" bcumi[2] "," bcumi[4] "\tRX:Z:" bcumi[6]}' | \
  samtools view -hb@8 > "$bam"
