#!/bin/bash

# Filtering the promoter BAM files
  # Maximum length
    # reads should be chr20:4685060--4688047 = 2988 bp
    # allowing some extra length for adapter/barcode >> say maximum 3500 bp
    # Why care? Longer reads should not be possible and so are suspicious
      # see eg. https://community.nanoporetech.com/posts/a-quick-guide-to-tiling-am, below Barcode mismatches

  # Minimum length
    # this can be slightly more flexible
    # say minimum 1000 bp

  # Soft-clipping
    # from looking at some reads in IGV, soft-clipping is typically can up to ~ 200 bp on each side and read still look good
    # i.e. 400 bp / 2988 bp = 13%
    # longer is potentially suspicious
    # say maximum 15%

  # Keep only primary alignments
    # should guarantee that read IDs are unique

mkdir promoterbamfilt

####
for i in promoterBams/*.bam # loop through the bam files
do

  bam="$i"

  echo
  echo
  echo
  echo "---- [ filtering BAM"$" $bam ] ----"
  echo

  # get the PDG number
  pdg="$(echo "$bam" | cut -d'_' -f 1 | cut -d'/' -f 2)"
  # BAM can be eg. 2021bam/1906_sorted_mdfix.bam
  # above will cut to 2021bam/1906, then 1906

  tmp="$(echo $"promoterbamfilt/tmp.bam")"
  bamfil="$(echo $"promoterbamfilt/""$pdg"$"_f.bam")"

  # filter the bam file
  filterBam.command -i $bam -f 1000 -c 3500 -s 0.15 -p yes -o $tmp

  # sort and index the filtered bam
  samtools sort $tmp > $bamfil
  samtools index $bamfil

  rm $tmp
  rm promoterbamfilt/tmp.bam.bai # not sure where it comes from!

done
