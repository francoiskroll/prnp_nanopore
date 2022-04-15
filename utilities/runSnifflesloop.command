#!/bin/bash

# run Sniffles on a folder of BAM files

mkdir SVs_sniffles/vcf

for i in promoterbamfilt/*.bam # loop thru files in directory that finish by .bam
do

  BAM="$i"

  PDG="$(echo "$BAM" | cut -d'_' -f 1 | cut -d'/' -f 2)" # everything before '_', then everything after '/'
  # e.g. promoterbamfilt/55052_promoterSnif.vcf >> promoterbamfilt/55052, then 55052

  VCF="$(echo $"SVs_sniffles/vcf/"$PDG$"_promoterSnif.vcf")"


  echo "---- [ SNIFFING PROMOTER "$" $BAM ] ----"

  sniffles --mapped_reads "$BAM" \
      --vcf "$VCF" \
      --genotype \
      --min_support 10 \
      --min_length 10 \
      --max_distance 50

done

######

for i in bamfilt/*.bam # loop thru files in directory that finish by .bam
do

  BAM="$i"

  PDG="$(echo "$BAM" | cut -d'_' -f 1 | cut -d'/' -f 2)" # everything before '_', then everything after '/'
  # e.g. bamfilt/55052_promoterSnif.vcf >> promoterbamfilt/55052, then 55052

  VCF="$(echo $"SVs_sniffles/vcf/"$PDG$"_bodySnif.vcf")"

  echo "---- [ SNIFFING GENEBODY "$" $BAM ] ----"

  sniffles --mapped_reads "$BAM" \
      --vcf "$VCF" \
      --genotype \
      --min_support 300 \
      --min_length 20 \
      --max_distance 50

done
