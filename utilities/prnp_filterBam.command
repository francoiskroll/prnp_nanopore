#!/bin/bash

# Filtering the gene-body BAM files before OPR analysis
  # Maximum length
    # reads should be chr20:4687732--4701756 = 14024 bp
    # allowing some extra length for long OPR insertions + adapter/barcode >> say maximum 15000 bp
    # Why care? Longer reads should not be possible and so are suspicious
      # see eg. https://community.nanoporetech.com/posts/a-quick-guide-to-tiling-am, below Barcode mismatches

  # Minimum length
    # this can be slightly more flexible
    # but keeping only reads that span most of region should help haplotype phasing (as most reads should cover all available SNVs)
    # say minimum 10000 bp

  # Soft-clipping
    # from looking at some reads in IGV, soft-clipping is typically ~ 105 left + 50 right / 14,024 = 1.1% of the read
    # longer is potentially suspicious
    # say maximum 5%

  # Keep only primary alignments
    # should guarantee that read IDs are unique

mkdir bamfilt

####
for i in genebodyBams/*.bam # loop through the bam files
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

  tmp="$(echo $"bamfilt/tmp.bam")"
  bamfil="$(echo $"bamfilt/""$pdg"$"_f.bam")"

  # filter the bam file
  filterBam.command -i $bam -f 10000 -c 15000 -s 0.05 -p yes -o $tmp

  # sort and index the filtered bam
  samtools sort $tmp > $bamfil
  samtools index $bamfil

  rm $tmp
  rm bamfilt/tmp.bam.bai # not sure where it comes from!

done
