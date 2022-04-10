#!/bin/bash

# filter reads in a BAM file

# currently handles
  # minimum Phred score
  # minimum reference span (previously read length)
    # Note; reference span is basically read length as you see it in IGV, i.e. how much of the reference it covers
    # it includes deletions/insertions and soft-clipped
    # e.g. read has 150 bp aligned + 10 bp deletion + 20 bp soft-clipped, reference span is 150 + 10 + 20 = 180 bp
    # however, if we were looking at its length (as in number of characters), it would be 150 + 20 = 170 bp, so doing so has the risk to bias against deletions
  # maximum read length (not currently in this version, see in v1 if need to add it back)
  # maximum proportion soft-clipped
  # only primary alignments (or all)

# explanation
  # filterBAM.command -i XXX.bam -e int -f int -s float -p (yes or no) -o XXX.bam
    # -i = input = bam file to process
    # -e = PhrEd score = minimum Phred score
    # -f = floor = minimum read span
    # -s = maximum proportion Soft-clipped
    # -o = output = name of bam file to output

# example
  #  filterBam.command -i F08.bam -e 40 -f 149 -s 0.2 -p no -o clean.bam

# v2: relies more heavily on R script readsToThrow.R to name the reads to remove
# easier as can edit there if want to add filters

######

# read the flags/arguments
while getopts i:e:f:s:p:o: flag
do
    case "${flag}" in
        i) bam=${OPTARG};;
        e) min_phred=${OPTARG};;
        f) min_readspan=${OPTARG};;
        s) max_softprop=${OPTARG};;
        p) primary=${OPTARG};;
        o) out=${OPTARG};;
    esac
done

echo
echo "Input BAM = $bam";
echo "Minimum Phred score = $min_phred";
echo "Minimum reference span = $min_readspan";
echo "Maximum proportion soft-clipped = $max_softprop";
echo "Keep only primary alignments = $primary";
echo "Output BAM = $out";
echo

######

nreads=$(samtools view -c $bam)
echo "Number of alignments: $nreads"

# will write a temporary SAM to remove secondary alignments and/or filter with readsToThrow.R script
tmp1=$"tmp1.sam"

samtools view -h -o $tmp1 $bam # convert to sam for below

######

# remove secondary alignments if needed
  # write temporary SAM = without secondary alignments
  # (or just a copy of tmp2.sam if -primary no)
tmp2=$"tmp2.sam"

if [ $primary = "yes" ]
then

  # counting number of alignments before
  before=$(samtools view -c $tmp1)

  # removing secondary alignments
  samtools view -h -F 256 -F 4 -F 2048 $tmp1 > $tmp2
    # See https://broadinstitute.github.io/picard/explain-flags.html
    # Note I think - F 4 (unmapped) -F 256 (not primary) are either not used or were removed at some point before
      # but safer to keep
    # -F 2048 = supplementary alignment, seems to include -2064 as well
    # but if I explicitly add -F 2064, it removes all Reverse (-F 16) for some reason

  # counting number of alignments after
  after=$(samtools view -c $tmp2)

  # how many did we remove?
  ndel=$( expr $before - $after )
  echo "Number of secondary alignments removed: $ndel"

else # if No (or anything else)
  # just store copy of tmp1 as tmp2 so we can continue below
  cp $tmp1 $tmp2

fi

# in both cases, also copy a BAM copy for below
tmp2b=$"tmp2.bam"
samtools view -bS $tmp2 -o $tmp2b


######

# Note; how to get full path of a file, example $sam
# echo "$(cd "$(dirname "$sam")"; pwd -P)/$(basename "$sam")"

# filter by soft-clipping

sampath="$( echo "$(cd "$(dirname "$tmp3")"; pwd -P)/$(basename "$tmp2")" )" # get full path of SAM file we wrote above

# removereads.txt = list of read names to be removed as too much soft-clipping
touch removereads.txt # write an empty file
# get its full path
removereadspath="$( echo "$(cd "$(dirname removereads.txt)"; pwd -P)/$(basename removereads.txt)" )" # get full path of SAM file we wrote above

Rscript ~/packages/myscripts/readsToThrow.R \
  $sampath \
  $min_phred \
  $min_readspan \
  $max_softprop \
  $removereadspath

# tell user number of reads which will remove as above soft-clipping threshold
ndel=$(wc -l < removereads.txt)
echo "Number of alignments removed by Phred/read span length/soft-clipped filters:$ndel"

# now exclude these reads from tmp2b.bam (the reads listed in removereads.txt)
tmp3=$"tmp3.bam"

# previous solution; but takes ages: samtools view -h $tmp3b | grep -vf removereads.txt | samtools view -bS -o $tmp4
# v2
echo # FilterSamReads still talks even though quiet, so leave some space
echo

# ! skip that step if 0 reads to remove
if [ $ndel -ne 0 ] # if not 0 reads to remove
then

  java -jar ~/packages/picard/picard.jar FilterSamReads \
    I=$tmp2b \
    O=$tmp3 \
    READ_LIST_FILE=removereads.txt FILTER=excludeReadList \
    use_jdk_deflater=true \
    use_jdk_inflater=true \
    quiet=TRUE

else
  # just store copy of tmp3b as tmp4 so we can continue below
  cp $tmp2b $tmp3

fi


# sort tmp2.bam and write final bam file
# index the final bam file
samtools sort $tmp3 > $out
samtools index $out

# tell final number of alignments
fin=$(samtools view -c $out)
echo
echo
echo "Final number of alignments: $fin"
echo
echo
echo

######

# clean after ourselves
rm $tmp1
rm $tmp2
rm $tmp2b
rm $tmp3
rm removereads.txt
