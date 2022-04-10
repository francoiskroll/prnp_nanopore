# take a SAM file as input, outputs the names of the reads which are soft-clipped


# packages & small functions ----------------------------------------------

library(stringr)

substrEnding <- function(x, n){ # x = string (or vector of strings); n = last n characters
  substr(x, nchar(x)-n+1, nchar(x))
}




# function readNametoThrow() ----------------------------------------------

# given one read (row) from a SAM file, checks whether we need to throw it or not based on filters
# if no need to throw the read; does not do anything
# if need to throw the read; outputs its name

# Note; probably useful in the future: how to split a CIGAR in its component:
# str_extract_all('29S45M633H59S', '[0-9]+[A-Z]')

readNametoThrow <- function(read, min_phred, min_readspan, max_softprop) { # read = one row of a correctly-imported SAM file
  
  ###
  
  # 1- check Phred score
  
  # get the read's Phred score
  phred <- as.numeric(read[5]) # ! assumes it is in column 5
  
  if (phred < min_phred) { # if it is below minimum Phred score, name that read
    cat('\t \t \t \t >>> Read', read[1], 'excluded because its Phred score is', phred, '\n')
    return(as.character(read[1]))
  } # Note; if Phred is below threshold, read name is output and function stops here due to return() call
  # is below only read if Phred is okay;
  
  # for debugging
  # cat('\n \n Phred is ok \n')
  
  ###
  
  # 2- check the read's reference span
  
  # get the read's CIGAR
  cigar <- as.character(read[6]) # ! assumes it is column 6
  
  tmp <- str_extract_all(cigar, '[0-9]+[A-Z]')[[1]] # splits the CIGAR into its componenets; e.g. 59S129M32S >> 59S 129M 32S
  span <- sum(as.numeric(unlist(strsplit(tmp, '[A-Z]')))) # sum the number parts, e.g. 59 + 129 + 32 = 220; this gives the reference span
  
  if (span < min_readspan) { # if read's span is below minimum length, name that read
    cat('\t \t \t \t >>> Read', read[1], 'excluded because its reference span is', span, 'bp \n')
    return(as.character(read[1]))
  }
  
  # for debugging
  # cat('Span is ok \n')
  
  
  ###
  
  # 3- check the read's soft-clipping
  # uses CIGAR from above
  
  if(length(grep(pattern='S', cigar)) != 0) { # if read is soft-clipped; i.e. there is S in CIGAR
    tmp <- str_extract_all(cigar, '[0-9]+S')[[1]] # extracts the S components; eg. 59S129M32S >> 59S 32S
    softbp <- sum(as.numeric(unlist(strsplit(tmp, '[A-Z]'))))
    # split at any letter (here will be only S), so example above: 59 32
    # convert to numeric & sum
    
    # divide by read length to get proportion soft-clipped
    softp <- as.numeric(softbp/nchar(read[10]))
    
    if (softp > max_softprop) { # if proportion soft-clipped is above unwanted threshold; get the read name
      cat('\t \t \t \t >>> Read', read[1], 'excluded because', round(softp * 100), '% of it is soft-clipped \n')
      return(as.character(read[1]))
    }
    
    # to get until here where nothing happens, read has to pass all filters above
    
  }
  
  # for debugging
  # cat('Soft-clipping is ok \n')
  
}


# main function -----------------------------------------------------------


namesOfReadsToThrow <- function(sampath, min_phred, min_readspan, max_softprop, outpath) {
  
  # check the path looks ok
  if(substrEnding(sampath, 4) != '.sam') stop('\t \t \t \t >>> That does not look like the path to a SAM file')
  
  # check the max_softprop argument looks ok
  if(!(max_softprop > 0 & max_softprop < 1)) stop('\t \t \t \t >>> softmax must be 0--1')
  
  # if all good -- import the SAM file
  numcols <- max(count.fields(sampath, sep='\t'), na.rm=TRUE)
  sam <- read.table(file=sampath, fill=TRUE, comment.char='@', sep='\t',
                    col.names=sprintf('col%i', 1:numcols))
  
  
  # apply getSoftReadName() function to all the reads
  # and write txt output
  excludereads <- as.data.frame(unlist(apply(sam, 1,
                                             readNametoThrow, min_phred, min_readspan, max_softprop)))
  
  # only keep unique names
  excludereads <- unique(excludereads)
  
  # say how many we are removing
  cat('\n \n \t \t \t \t >>> Excluding total', nrow(excludereads), 'reads \n \n \n')
  
  write.table(excludereads, file=outpath,
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  
}


# run that function -------------------------------------------------------

# first read the Terminal arguments
args <- commandArgs(trailingOnly=TRUE)

cat('\n \n readsToThrow.R got arguments: \n')

sampath <- args[1] # read first argument from Terminal = path to SAM file
cat('\t Path to SAM file = ', sampath, '\n')

min_phred <- args[2] # read second argument from Terminal = minimum Phred score
cat('\t Minimum Phred score = ', min_phred, '\n')

min_readspan <- args[3] # read third argument from Terminal = minimum reference span
cat('\t Minimum reference span = ', min_readspan, 'bp \n')

max_softprop <- args[4] # read fourth argument from Terminal = maximum proportion soft-clipped
cat('\t Maximum proportion soft-clipped ', max_softprop, '\n')

outpath <- args[5] # read fifth argument from Terminal = name for output (output being list of reads to remove)
cat('\t Output =', outpath, '\n')


# run the main function = namesOfSoftReads()
namesOfReadsToThrow(sampath=sampath,
                    min_phred=min_phred,
                    min_readspan=min_readspan,
                    max_softprop=max_softprop,
                    outpath=outpath)