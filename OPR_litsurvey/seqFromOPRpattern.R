# generate sequence from OPR pattern

library(openxlsx)
library(here)
library(dplyr)
library(tidyr)


# import OPR patterns -----------------------------------------------------

oprs <- read.xlsx(here('OPR_litsurvey', 'OPRLitCatalog.xlsx'))

# any duplicates in pattern?
sum(duplicated(oprs$pattern)) # No, OK

# any duplicates in allele?
sum(duplicated(oprs$allele)) # No, OK


# import reference for each OPR -------------------------------------------

oref <- read.xlsx(here('OPR_litsurvey', 'eachOPR.xlsx'))



# function to generate full sequence from a pattern -----------------------

pattern2Seq <- function(pattern) {
  
  oseqs <- sapply(strsplit(pattern, '/')[[1]],
                  function(r){
                    as.character(subset(oref, Rname==r, sequence))
                  })
  
  fullseq <- paste(oseqs, collapse='')
  
  return(fullseq)
  
}


# apply the function ------------------------------------------------------

oprs$fullsequence <- sapply(oprs$pattern,
                            function(pa){
                              pattern2Seq(pa)
                              })


# add OPRD1_Lee2016 sequence manually -------------------------------------
# allele OPRD1_Lee2016 is a deletion that overlaps two OPRs, add sequence manually

oprs[which(oprs$allele=='OPRD1_Lee2016'), 'fullsequence'] <-
  'CCTCAGGGCGGTGGTGGCTGGGGGCAGCCTCATGGTGGTGGCTGGGGGCAGCCTCATGGTGGTGGCTGGGGGCAGCCCCATGGTGGTGGCTGGGGTCAA'

# check that all the lengths make sense
which((nchar(oprs$fullsequence) == oprs$OPR_length)==FALSE) # all correct


# write the file ----------------------------------------------------------

# will need to copy the $pattern column from original file to keep the colouring
write.xlsx(oprs, here('OPR_litsurvey', 'OPRLitCatalogSeq.xlsx'))
