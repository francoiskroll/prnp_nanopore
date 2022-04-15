# 100720
# returns number of reads in a fastq file

# usage: numberReadsFastq.command XXX.fastq

#!/bin/bash

FASTQGZ="$1"

zcat < $FASTQGZ | echo $((`wc -l`/4))
