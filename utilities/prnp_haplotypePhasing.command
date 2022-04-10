#!/bin/bash

# script below is meant to be run being in ~/Dropbox/prnp_nanopore/

# for SNVs phasing/haplotagging, there are 2 possibilities

  # 1- there is only 1 SNV available for that sample, so VCF written by prepareVCFforWhatshap.R is 'as if' it went through whatshap phase already
    # the reason being: whatshap phase does not work with a single SNV as the goal is to phase variants in relation to each other
      # however, whatshap haplotag works fine and that what we need to separate the reads from each haplotype; but the variant need to 'look' phased by whatshap phase

  # 2- there are more than 1 SNV available for that sample, go through whatshap phase, then whatshap haplotag (normal process)


####
# Possibility 1-
  # if only 1 SNV; prepareVCFforWhatshap.R wrote VCF as pdg_whap.vcf (p for phased)
for i in haplotypephasing/vcf_phased/*_whap.vcf # loop through those VCF files (that end with whap.vcf) in the folder
do
  # zip the VCF, index it, run directly whatshap haplotag
  vcf="$i"

  echo
  echo
  echo
  echo "---- [ using VCF"$" $vcf ] ----"
  echo

  gz="$(echo "$vcf"$".gz")" # append .gz to the VCF filename

  bgzip -c $vcf > $gz # zip it (! need to bgzip and not gzip for tabix below)
  tabix -p vcf $gz # index the zipped vcf

  # get the PDG number
  pdg="$(echo "$vcf" | cut -d'_' -f 2 | cut -d'/' -f 2)"
  # vcf above will be eg. 2021haplotypephasing/vcf_phased/55048_whap.vcf
  # above will cut to 55048

  mkdir haplotypephasing/bam_haplotagged
  hptmp="$(echo $"haplotypephasing/bam_haplotagged/hptmp.bam")" # will be haplotagged bam, before samtools sort
  bamp="$(echo $"haplotypephasing/bam_haplotagged/""$pdg"$"_hp.bam")" # final output: bam where reads are haplotagged (added tag to each read to tell which haplotype it is from)
  # hp for haplotype-phased

  bamf="$(echo $"bamfilt/""$pdg"$"_f.bam")" # filtered bam (before phasing, input to whatshap haplotag)
  # should be eg. bamfilt/1906_f.bam

  # run directly whatshap haplotag
  whatshap haplotag \
  --ignore-read-groups \
  -o $hptmp \
  --reference /Users/francoiskroll/Dropbox/rotations/Rotation3/hg38/hg38.fa \
  $gz $bamf

  # sort and index the phased bam
  samtools sort $hptmp > $bamp
  samtools index $bamp

  # delete tmp file
  rm $hptmp

done

###
# Possibility 2-
  # if more than 1 SNV; prepareVCFforWhatshap.R wrote VCF as pdg_wha.vcf
    # need to first run whatshap phase to phase variants in relation to each other
    # then whatshap haplotag to tag the haplotype on each read
for i in haplotypephasing/vcf_unphased/*_wha.vcf # loop through those VCF files (that end with wha.vcf) in the folder
do
  # zip the VCF and index it (same as Possibility 1-)
  vcf="$i"

  echo
  echo
  echo "---- [ using VCF"$" $vcf ] ----"
  echo
  echo

  gz="$(echo "$vcf"$".gz")" # append .gz to the VCF filename

  bgzip -c $vcf > $gz # zip it (! need to bgzip and not gzip for tabix below)
  tabix -p vcf $gz # index the zipped vcf

  # get the PDG number
  pdg="$(echo "$vcf" | cut -d'_' -f 2 | cut -d'/' -f 2)"
  # vcf above will be eg. 2021haplotypephasing/vcf_unphased/1906_wha.vcf
  # above will cut to 1906

  vcfp="$(echo $"haplotypephasing/vcf_phased/""$pdg"$"_whap.vcf")" # vcf phased (output of whatshap phase)

  bamf="$(echo $"bamfilt/""$pdg"$"_f.bam")" # filtered bam (before phasing, input to whatshap haplotag)
  # should be eg. bamfilt/1906_f.bam

  mkdir haplotypephasing/bam_haplotagged # will place some temporary files (that we will delete), then the final haplotagged file
  hptmp="$(echo $"haplotypephasing/bam_haplotagged/hptmp.bam")" # will be haplotagged bam, before samtools sort
  bamp="$(echo $"haplotypephasing/bam_haplotagged/""$pdg"$"_hp.bam")" # final output: bam where reads are haplotagged (added tag to each read to tell which haplotype it is from)
  # hp for haplotype-phased

  # run whatshap phase, using the filtered bam
  whatshap phase \
    --ignore-read-groups \
    -o $vcfp \
    --reference=/Users/francoiskroll/Dropbox/rotations/Rotation3/hg38/hg38.fa \
    $gz \
    $bamf

  # now use phased vcf to tag the reads with the haplotype from which they are from (whatshap haplotag)

  # zip & index the vcf created above
  gzp="$(echo "$vcfp"$".gz")" # append .gz to the phased VCF filename

  bgzip -c $vcfp > $gzp # zip it (! need to bgzip and not gzip for tabix below)
  tabix -p vcf $gzp # index the zipped vcf

  # use it to haplotag reads
  whatshap haplotag \
  --ignore-read-groups \
  -o $hptmp \
  --reference /Users/francoiskroll/Dropbox/rotations/Rotation3/hg38/hg38.fa \
  $gzp $bamf

  # sort and index the phased bam
  samtools sort $hptmp > $bamp
  samtools index $bamp

  # delete tmp file
  rm $hptmp

done
