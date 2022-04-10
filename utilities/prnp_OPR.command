# see somaticOPR_README.md for more complete explanations
# you probably want to be in ~/.../haplotypephasing/bam_haplotagged/

#!/bin/bash

for i in *.bam # loop thru files in directory that finish by .bam
do

  BAM="$i"


####
# 1- trim BAM to keep only OPR
  # genomic windows we do NOT want are in OPRpos.bed

  trim="trim.bam" # temporary file, will delete later
  pdg="$(echo "$BAM" | cut -d'_' -f 1)" # everything before first _
  # names of BAM files eg. 1906_sorted_mdfix.bam, so will give PDG (1906)

  mkdir clipLogs
  log="$(echo $"clipLogs/""$pdg"$"_clipLog.txt")"

  echo
  echo
  echo
  echo "---- [ Sample"$" $BAM ] ----"
  echo
  echo "  Trimming to keep only OPR..."

  samtools ampliconclip --hard-clip --both-ends --clipped \
  -b ~/Dropbox/prnp_nanopore/needleSam/OPRpos.bed \
  $BAM \
  -o $trim \
  -f $log


####
# 2- filter trim.bam
  # minimum read length
    # = 21 bp
    # equivalent to 4 OPRD (only R1 left, 27 bp) - 6 bp for low noise
  # maximum read length
    # = 702 bp
    # equivalent to reference (123 bp) + 24 OPRI (576 bp) + 3 bp for high noise
  # maximum proportion of read soft-clipped = 0.1
  # keep only primary alignments = yes

  echo "  Filtering the trimmed BAM..."

  mkdir OPRtrim
  out="$(echo $"OPRtrim/""$pdg"$"_opr.bam")"

  filterBam.command -i $trim -f 21 -c 702 -s 0.1 -p yes -o $out
  # Note; this will partly overlap with the filtering performed in prnp_haplotypePhasing.command

  # we do not need to keep trim.bam
  rm $trim


####
# 3- count read lengths
  # from https://www.biostars.org/p/65216/

  mkdir readlengths
  readlghts="$(echo $"readlengths/""$pdg"$"_lengths.txt")"


  # readlghts="$(echo "$pdg"$"_lengths.txt")"

  samtools view -F 4 "$out"| # -F 4 = only mapped reads
    cut -f 10 |
    perl -ne 'chomp;print length($_) . "\n"' |
    sort -n |
    uniq -c > $readlghts


####
# 4- convert to SAM so can work on it in R
  mkdir OPRtrim_sam
  tmp1="$(echo $"OPRtrim_sam/tmp1.sam")"
  tmp2="$(echo $"OPRtrim_sam/tmp2.sam")"
  sam="$(echo $"OPRtrim_sam/""$pdg"$"_opr.sam")"

  samtools view -h -o $tmp1 $out # convert to SAM (first a temporary file)

  # remove the SA tag from the temporary SAM
  # (it tags chimeric read)
    # it is a pain in R later because it is not present for all reads, so it shifts the columns for certain reads
  sed 's/\tSA\:Z\:[^\t]*//' $tmp1 > $tmp2

  # drop the QUAL column (in practice, switch it to *)
    # I am not using it and it is full of special characters, including @ which is the comment character when importing in R
      # makes it very challenging to import in R correctly if not removing it
  cat $tmp2 | awk -v OFS='\t' '$11="*"' > $sam
    # above is not a particularly elegant solution because it also affects the header lines, but I am not importing them in R...

  # delete the temporary SAM
  rm $tmp1
  rm $tmp2

####

done
