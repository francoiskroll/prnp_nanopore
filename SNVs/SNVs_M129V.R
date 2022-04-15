# Reviewer's question (paraphrased):
# M129V individuals tend to have a higher number of non-coding SNVs (Table 1).
# Are these SNVs occuring on the same haplotype?

# packages ----------------------------------------------------------------

library(openxlsx)
library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)


# which samples are M129V -------------------------------------------------
# look in samples sheet of AdditionalFile1
meta <- read.xlsx('AdditionalFile1.xlsx', sheet='samples')

mvs <- meta[which(meta$codon129=='M129V'), 'sample'] # mvs for M129V samples


# import phased VCFs for these samples ------------------------------------

# temporary list to store the VCFs
vcfs <- vector(mode='list', length=length(mvs))

# for each M129V sample
# build the path to its phased VCF file
# import it (i.e. add it to the list)
for (sp in 1:length(mvs)) {
  # read.table() below understands it has to skip the header lines starting with #
  # import that one vcf
  vc <- read.table(here('haplotypephasing', 'vcf_phased', paste0(mvs[sp], '_whap.vcf')), sep='\t')
  
  # add the column names back
  colnames(vc) <- c('chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'phasing')
  
  # add the sample ID as first column
  vc <- vc %>%
    add_column(sampleid=mvs[sp], .before='chr')
  
  # add that vcf to the list
  vcfs[[sp]] <- vc
}

# now collapse the list in one dataframe
# below will overwrite the list, we do not need it anymore
vcfs <- as.data.frame(rbindlist(vcfs))


# check that only one phase set in each VCF -------------------------------

# from https://whatshap.readthedocs.io/en/latest/guide.html
# under Phasing represented by PS (“phase set”) tag
# >> it can be that not all SNVs in the VCF are phased together
# below, we will assume that 1 or 0 refers to the same haplotype throughout each VCF, so we need to check this assumption

# create new column phase set
vcfs$phaseset <- unlist(lapply(strsplit(vcfs$phasing, split=':'), function(s) {s[2]}))

# simplify phasing column by removing phasingset info
vcfs$phasing <- unlist(lapply(strsplit(vcfs$phasing, split=':'), function(s) {s[1]}))

# check all sample only have 1 phaseset
# below, counts number of distinct phasesets for each sample
vcfs %>%
  group_by(sampleid) %>%
  summarise(n_distinct(phaseset))
# all 1, so OK to assume all SNVs for each sample are phased together


# get phasing info so M129V is always on the same haplotype ---------------
# within each sample, 1 or 0 is arbitrary
# below will be simpler is M129V is always on the same haplotype

# is M129V always on the same haplotype?
# M129V is position chr20:4699605
vcfs[which(vcfs$pos=='4699605'), 'phasing']
# no, 4/9 are 0|1 and 5/9 are 1|0
# this is essentially what you would expect

# within any sample where M129V is 0|1, flip all the phasing info

# small function to flip phasing info
flipPhase <- function(phase) {
  
  # sapply below so it can be ran on a vector directly
  phaseflip <- sapply(phase, function(pha) {
    
    if (pha!='0|1' & pha!='1|0') {
      stop('\t \t \t \t >>> Error: phasing info is not 1|0 or 0|1 \n')
    }
    
    else if (pha=='0|1') {
      return('1|0')
    }
    
    else if (pha=='1|0') {
      return('0|1')
    }
    
  })
  
  # simplify it
  phaseflip <- as.vector(phaseflip)
  
  # return it
  return(phaseflip)
  
}
###


# which are the samples we need to flip?
# i.e. which are the rows where pos is 4699605 and phasing is 0|1
flipsps <- vcfs[which(vcfs$pos=='4699605' & vcfs$phasing=='0|1') , 'sampleid'] # flip these samples

# now to flip these samples
# take all the SNVs belonging to these samples and flip them
vcfs[which(vcfs$sampleid %in% flipsps), 'phasing'] <-  flipPhase( vcfs[which(vcfs$sampleid %in% flipsps), 'phasing'] )
# Note: do not need to do sample by sample
# for each sample, we need to flip all the phasing info once
# so same to flip all the SNVs at once

###

# now check that M129V always 1|0
vcfs[which(vcfs$pos=='4699605'), 'phasing']
# yes -- now all 1|0


# for each SNVs, how many times is it phased with M129V? ------------------

# go through each unique SNV (by position),
# then for each sample which has this SNV,
# is it on the same haplotype as M129V?
# we set M129V to always be 1|0, so just need to count how many times it is 1|0

# first, for each SNV, count how many samples have it
wom <- vcfs %>%
  group_by(pos) %>%
  tally(name='n_samples')
# which will preallocate wom for with/without M129V

# count by position by phase
n_pha <- vcfs %>%
  group_by(pos, phasing) %>%
  tally()

# only the SNVs which were phased as 1|0:
n10 <- subset(n_pha, phasing=='1|0')
n10$phasing <- NULL # delete phasing info, we know it is 1|0
# merge these samples count to wom
wom <- left_join(x=wom, y=n10, by='pos')

# now repeat for SNVs which were phased as 0|1
n01 <- subset(n_pha, phasing=='0|1')
n01$phasing <- NULL # delete phasing info, we know it is 0|1
# merge these samples count to wom
wom <- left_join(x=wom, y=n01, by='pos')

# now rename the column
colnames(wom) <- c('pos', 'n_samples', 'n_withM129V', 'n_woM129V')

# everything that is still NA is actually 0 sample
wom$n_withM129V[is.na(wom$n_withM129V)] <- 0
wom$n_woM129V[is.na(wom$n_woM129V)] <- 0

# few checks
# for each SNV, sum of number of samples which have it as 1|0 (i.e. with M129V) + number of samples which do not have it as 0|1 (i.e. without M129V)
# should be equal total number of samples with that SNV
unique(wom$n_withM129V + wom$n_woM129V == wom$n_samples)
# OK
# looking at the table: SNVs are always on one haplotype or the other

# ! should delete M129V from the table otherwise we will be counting it
wom <- wom[-which(wom$pos==4699605),]

# in summary
# out of 25 SNVs (not counting M129V), 
nrow(subset(wom, n_withM129V>0 & n_woM129V==0))
# 22 were always found with M129V

nrow(subset(wom, n_woM129V>0 & n_withM129V==0))
# only 3 SNVs were found not with M129V

#####

# in samples without M129V ------------------------------------------------

# you would expect to *never* find the 23 SNVs from above
# and only find the 3 SNVs (plus potentially others)

####

# I think looking at phasing of the SNVs in the no-M129V samples is beyond the question
# will essentially consider them as all haplotypes not carrying M129V

# as we will not use the phasing information, we can simply look at the AdditionalFile1
# easier than deal with the VCFs

# which samples do not carry M129V? ---------------------------------------
# look in samples sheet of AdditionalFile1
meta <- read.xlsx('AdditionalFile1.xlsx', sheet='samples')

ovs <- meta[which(meta$codon129!='M129V'), 'sample'] # ovs for other (not M129V) samples

# import the SNVs
snvs <- read.xlsx('AdditionalFile1.xlsx', sheet='SNVs_filtered')

# keep only the SNVs of the samples without M129V
snvo <- subset(snvs, SAMPLE %in% ovs)


# do any of the samples carry the 'with M129V-SNVs' -----------------------

# in the M129V samples, these SNVs are always found on the same haplotype
snv_with <- as.vector(unlist(wom[which(wom$n_withM129V > 0), 'pos']))
# so you would expect to never find them in them in the no M129V samples

intersect(snvo$POS, snv_with)
# 13 of 22 are found
# i.e. these 14 SNVs are NOT 'M129V specific' **
# the other 9 seem M129V specific
# but difficult to be certain! It may be that we see them co-segregate in a larger cohort

# the 9 that may be M129V specific are:
setdiff(snv_with, snvo$POS)

subset(wom, pos %in% setdiff(snv_with, snvo$POS))

# ** for example, chr20:4690579
# in M129V samples, it was only found on the M129V haplotype (9 of 9)
wom[which(wom$pos=='4690579'),]

# however, it is also found in samples without M129V:
snvs[which(snvs$POS=='4690579' & snvs$SAMPLE %in% ovs),] # found in samples 46345, 52331, 55826 which is not M129V

# !!! 52331 was not tested for M129V by Sanger, so written as not M129V
# but detected as M129V by Nanopore, so very likely M129V


# conversely...
# in the M129V samples, these SNVs are always found on the OTHER haplotype than M129V
snv_without <- as.vector(unlist(wom[which(wom$n_woM129V > 0), 'pos']))
# so you would expect to find them in the no M129V samples

intersect(snvo$POS, snv_without)
# 2 out of 3 are found
# not found: 4698921 but it was only present in one sample w/ M129V so relatively rare it seems
snvs[which(snvs$POS=='4698921'),] # of the whole cohort there is only sample with that SNV
# so roughly as expected: SNVs found on the other haplotype as M129V are also found in samples without M129V