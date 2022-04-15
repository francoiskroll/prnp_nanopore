# Note, below is all taken directly from needleSam_exploration.R
# here, applied to the control PCR on sample #56635
# will do as little modifications as possible, so will look convoluted as original script is applied on all the samples

# after needleSam.R
  # needleSam.R creates the dataset allOPRreads.csv
  # needleSam_exploration.R: exploring the dataset for possible OPR somatic mutations

exportOrNo <- FALSE # Do you want to export the plots (useful to put FALSE when just want to load everything)

# packages ----------------------------------------------------------------

library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(ggbeeswarm)
library(ggpubr)



# small functions ---------------------------------------------------------

# function to pool columns from a dataframe into one column
  # ! assumes that each row, there is only one element from the columns to pool which is not NA
    # eg. at row 5, the columns we want to pool are NA, 9, NA >> gives 9
    # but if 2, 9, NA >> Error! (as more than one element)

# usage example; poolColumns(df, c('col2', 'col5'))
poolColumns <- function(dataframe, vectorOfColNames) {
  apply(dataframe, 1, function(row) {
    tmp <- as.vector(unlist(c(row[vectorOfColNames]))) # concatenate the columns we need to pool at this row
    # check there is only one not NA
    if (sum(!is.na(tmp)) != 1) stop ('\t \t \t \t >>> Error: One or more rows has more than one element in the columns to be pooled \n')
    # if ok, remove the NA and we are left with a single element
    return(tmp[!is.na(tmp)]) # return for within the apply
  })
}

###

# counts & frequencies within bins of a continuous variable (each bin = one range of values) per group
# e.g. data is number of mismatches for each read
# and 2 groups: reference / mutated
# for each group; how many reads in bin 0--10 mismatches; how many reads in bin 11-20 etc.
# both in counts and frequencies (% of total within each group)
# Note; bins represent the right boundary, e.g. a case where we have proportions, 0.1 means the bin that contains the data from 0 to 0.1
# same in the dataframe it returns
# example usage/ freqInBins(df = mydata, datcol = 'mismatchpercentage', grpcol = 'treated', bins=seq(0, 1, by=0.1))
# by default grpcol = NA, meaning it does do it by any group, just takes all the data in column datcol
freqInBins <- function(df, datcol, grpcol=NA, bins) {
  
  ###
  # if no grouping
  if (is.na(grpcol)) {
    dat <- as.numeric(df[, datcol]) # data we need to bin
    cos <- hist(dat, breaks=bins, include.lowest=TRUE, plot=FALSE)$counts # counts
    fre <- cos / sum(cos) # frequencies
    gdf <- data.frame(bin=bins[-1], counts=cos, freq=fre) # small dataframe for this group, two columns counts and frequencies
    
    return(gdf)
    
  }
  
  ###
  # if grouping
  
  else {
    
    ngrps <- length(unique(df[,grpcol])) # number of groups
    grpnms <- unique(df[,grpcol]) # group names
    
    cfpg <- vector(mode='list', length=ngrps) # preallocate list counts + frequencies per group
    
    for(i in 1:length(grpnms)) { # loop through groups
      dat <- as.numeric(df[ which(df[,grpcol] == grpnms[i]) , datcol]) # data for this group (the data we need to bin)
      cos <- hist(dat, breaks=bins, include.lowest=TRUE, plot=FALSE)$counts # counts
      fre <- cos / sum(cos) # frequencies
      
      gdf <- data.frame(grp=grpnms[i], bin=bins[-1], counts=cos, freq=fre) # small dataframe for this group, two columns counts and frequencies
      
      # group column should be called the same as grpcol
      colnames(gdf)[1] <- grpcol
      
      # add it to the list
      cfpg[[i]] <- gdf
    }
    
    # stick the list together and return it
    fbing <- as.data.frame(rbindlist(cfpg)) # frequencies in each bin by group
    return(fbing)
  }
  
}


# plot colours ------------------------------------------------------------

# somatic call TRUE/FALSE
somcols <- c('#f1876b', '#B3BDC4')

# strand
  # will use same colours as IGV
  # forward = pink; reverse = blue
strcols <- c('#ebb0af', '#aeafd8') # ! Forward, Reverse

# haplotypes
hapcols <- c('#f1876b', '#4D4D4D', '#B3BDC4') # ! 1, 2, 0

# cohorts
  # inherited = yellow
  # sCJD = khaki
  # control = grey
cohcols <- c('#fcb505', '#78ac63', '#B3BDC4')




# import ------------------------------------------------------------------

# all the reads
ied <- read.csv(here('controlPCR', 'conPCR_allOPRreads.csv'))

# sample metadata
meta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='samples')

# format meta as in needleSam.R
  # - >> NA
meta$SV[which(meta$SV=='-')] <- NA

  # special case: PDG 46345 is 5 OPRI/1 OPRD
    # make a new column SV2 for it (just for the 1 OPRD)
meta <- meta %>%
  mutate(SV2=NA, .after=SV)
meta[which(meta$sample==46345), 'SV2'] <- unlist(strsplit(meta[which(meta$sample==46345), 'SV'], '/'))[2]
meta[which(meta$sample==46345), 'SV'] <- unlist(strsplit(meta[which(meta$sample==46345), 'SV'], '/'))[1]

  # 2 OPRD >> OPRD2, etc.
    # first lapply takes what is after the space = OPRD or OPRI
    # second lapply takes what is before the space = 2, 1, ...
meta$SV <- paste0(lapply(strsplit(meta$SV, ' '), function(i) {i[2]}),
                  lapply(strsplit(meta$SV, ' '), function(i) {i[1]}))
    # the NA above become 'NANA', put them back to NA
meta$SV[which(meta$SV=='NANA')] <- NA

  # same for SV2
meta$SV2 <- paste0(lapply(strsplit(meta$SV2, ' '), function(i) {i[2]}),
                   lapply(strsplit(meta$SV2, ' '), function(i) {i[1]}))
meta$SV2[which(meta$SV2=='NANA')] <- NA

  # sample: 52331 >> pdg52331, etc.
meta$sample <- sprintf('pdg%s', meta$sample)


# haplotype phasing metadata
  # eg. how many SNVs available for haplotagging
hpmeta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='haplotypephasing')
    # sample: 52331 >> pdg52331, etc.
hpmeta$sample <- sprintf('pdg%s', hpmeta$sample)

  # from this:
    # pdgs of samples that we could haplotag
pdgs_tag <- as.character(unlist(subset(hpmeta, haplotagged, sample)))
    # pdgs of samples we could not haplotag
pdgs_notag <- as.character(unlist(subset(hpmeta, !haplotagged, sample)))

# sequencing results summary
seqmeta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='sequencing_summary')
# sample: 52331 >> pdg52331, etc.
seqmeta$sample <- sprintf('pdg%s', seqmeta$sample)

# we will also need all pdgs
pdgs <- meta$sample
# do not take the ones we do not have loaded (only applies if running script on subset of samples)
pdgs <- intersect(unique(ied$sample), pdgs)




# add sample info to ied for simplicity -----------------------------------

ied <- left_join(ied, meta, by='sample')

ied <- hpmeta %>%
  select('sample', 'haplotagged', 'genebody_heterozygousSNVs') %>%
  left_join(ied, ., by='sample')

ied <- seqmeta %>%
  select('sample', 'genebody_coverage', 'regulatoryregion_coverage') %>%
  left_join(ied, ., by='sample')




# Possible somatic mutations: approach 1 ----------------------------------

# there is no read that has both an OPRD and a OPRI that pass mismatch threshold, which makes things easier
subset(ied,
       (opri_pass & oprd_pass))

# interesting reads are:
  # a- there is a mutation (not reference)
    # ! except for PDG 46345 (compound homozygous 1 OPRD / 5 OPRI), for this sample 'reference' may be a somatic mutation
  # b- Sanger SV is NA (eg. sCJD cases) or is different than the mutation in the read (not interesting if matches the Sanger genotype)
  # c- One of OPRI or OPRD is below mismatch threshold

# in practice:
  # 1- all reads that match what we know from Sanger (so uninteresting here) have a true_nummis number, so can take only true_nummis = NA
  # (true_nummis = number of mismatches vs reference, when it matches what is expected)
    # by doing so, no need to worry about PDG 46345 (compound homozygous 1 OPRD / 5 OPRI), true_nummis is present for both OPRD1 or OPRI5 reads
      # (and is NA if read says reference, as it should)
  # 2- ! true_nummis = NA also includes unassigned reads, so need to exclude those
    # (! remember to explicitly include othercategory = NA otherwise != excludes NA)
  # 3- make use of mismatch threshold, opri_pass or oprd_pass

iedSom1 <- subset(ied,
                 is.na(true_nummis) & # 1- 
                   (is.na(othercategory) | othercategory!='unassigned') & # 2-
                   (opri_pass | oprd_pass | ref_pass)) # 3-, one of opri_pass, oprd_pass, ref_pass is TRUE


# Possible somatic mutations: approach 2 ----------------------------------
# Another approach...
  # allows to confirm the above if we reach the same conclusion with two different logics

# tag the 'expected' reads, i.e. those that do not look like a somatic mutation
ied$expected <- NA # preallocate all as NA

  # all reads whose OPRI matches their sample's SV are 'expected'
ied[which( ied$opri == ied$SV ), 'expected'] <- TRUE

  # all reads whose OPRI matches their sample's SV2 are 'expected' (only applies to PDG 46345)
ied[which( ied$opri == ied$SV2 ), 'expected'] <- TRUE

  # all reads whose OPRD matches their sample's SV are 'expected'
ied[which( ied$oprd == ied$SV ), 'expected'] <- TRUE

  # all reads whose OPRD matches their sample's SV2 are 'expected' (only applies to PDG 46345)
ied[which( ied$oprd == ied$SV2 ), 'expected'] <- TRUE

  # all reads that are reference are 'expected' (! except for PDG 46345, read above) 
ied[which( ied$sample != 'pdg46345' & ied$othercategory == 'reference' ), 'expected'] <- TRUE

# >>> at this stage, expected should tag the same rows as trueread
identical(subset(ied, trueread, read), subset(ied, expected, read))

  # lots of reads are unassigned to an OPR category
    # >> will tag them as expected for simplicity; this will be the difference vs trueread
ied[which( ied$othercategory == 'unassigned' ), 'expected'] <- TRUE

# great majority of reads should be 'expected'
cat('\t \t \t \t >>> ', sum(ied$expected, na.rm=TRUE), 'out of', nrow(ied), 'considered expected \n')
cat('\t \t \t \t \t i.e.', round((sum(ied$expected, na.rm=TRUE) / nrow(ied))*100, 1), '% \n')

  # tag the rest as 'not expected' (FALSE), i.e. potentially somatic mutation
ied[which(is.na(ied$expected)), 'expected'] <- FALSE

  # check there is no NA left
sum(is.na(ied$expected))

# possible somatic = 'not expected', and pass mismatch threshold
iedSom <- subset(ied,
                 !expected &
                   (opri_pass | oprd_pass | ref_pass)) # one of opri_pass, oprd_pass, ref_pass is TRUE


# >>> does it give the same reads as first logic?
identical(iedSom$read , iedSom1$read)
# yes, same result
# will only keep iedSom
rm(iedSom1)

# will be useful for later;
  # add column 'isSomatic' (TRUE or FALSE) to ied, i.e. this read is a somatic mutation call
ied$isSomatic <- FALSE # preallocate as FALSE (most of them)
ied[match(iedSom$read, ied$read), 'isSomatic'] <- TRUE # for the reads that are in iedSom, replace by TRUE
# put isSomatic in order I want
ied$isSomatic <- factor(ied$isSomatic, levels=c('TRUE', 'FALSE'))





# Illustrate mismatch threshold -------------------------------------------

# in read length vs number of mismatches;
  # plot all 'potentially somatic mutations' reads and colour the ones that pass the mismatch threshold
    # the reads that have both OPRI1 and OPRD1 make this annoyingly complicated
    # for those, I think will simply plot the minimum mismatch of OPRI/OPRD (read length does not change)
    # (to be honest it does not really matter as all of these reads are above the mismatch threshold anyways)

# all somatic mutations candidates
  # (i.e. iedSom from above + all not expected but which does not pass the mismatch threshold)
iedsomall <- subset(ied,
                    !expected)

iedsomall$conmis <- NA # number of mismatches vs consensus

# for reference reads, fill in conmis with reference_nummis
iedsomall[which(!is.na(iedsomall$reference_nummis)) , 'conmis'] <- iedsomall[which(!is.na(iedsomall$reference_nummis)) , 'reference_nummis']

# for OPRI reads, fill in conmis with opri_nummis
iedsomall[which(!is.na(iedsomall$opri_nummis)) , 'conmis'] <- iedsomall[which(!is.na(iedsomall$opri_nummis)) , 'opri_nummis']

# for OPRD reads, fill in conmis with oprd_nummis
iedsomall[which(!is.na(iedsomall$oprd_nummis)) , 'conmis'] <- iedsomall[which(!is.na(iedsomall$oprd_nummis)) , 'oprd_nummis']

# check no more NA
sum(is.na(iedsomall$conmis))

# ! now, any reads that has both OPRD & OPRI got the number of mismatches vs OPRD consensus (as OPRD was last)
# for these reads, go back and take the minimum number of mismatches between OPRD and OPRI (the best match)
rows2edit <- which( !is.na(iedsomall$opri) & !is.na(iedsomall$oprd) )
for (r in rows2edit) {
  iedsomall[r,'conmis'] <- min(c(iedsomall[r,'opri_nummis'] , iedsomall[r,'oprd_nummis']))
}
# >>> we have number of mismatches vs consensus for all candidate somatic mutations 

# mismatch threshold
  # mean of % mismatch of all true reads (those that match what we know is there) + sd
thrMis <- mean(ied$true_misbp, na.rm=TRUE) + sd(ied$true_misbp, na.rm=TRUE) # as calculated in needleSam.R

# control PCR: better to use same threshold as before for consistency
thrMis <- 0.058

# more intuitively, the threshold excludes the
1 - (nrow(subset(ied, trueread & true_misbp < thrMis)) / nrow(subset(ied, trueread))) # least accurate reads


# on plot, will add markers for reference OPR lengths (eg. reference 123 bp, 1 OPRI = 123 + 24 bp, etc.)
ref <- 9 * 1 * 3 + 8 * 4 * 3 # reference = 1 nonapeptide + 4 octapeptide
# OPR I / D are simply +/- full Rs, so +/- multiples of 8 * 3 = 24 bp
oct <- 8 * 3

oprlengths <- c(ref-4*oct,
                ref-3*oct,
                ref-2*oct,
                ref-1*oct,
                #ref,
                ref+1*oct,
                ref+2*oct,
                ref+3*oct,
                ref+4*oct,
                ref+5*oct,
                ref+6*oct,
                ref+7*oct,
                ref+8*oct,
                ref+9*oct,
                ref+10*oct,
                ref+11*oct,
                ref+12*oct,
                ref+13*oct)

somMisThr <- ggplot(iedsomall, aes(x=read_length, y=conmis, colour=isSomatic)) +
  geom_vline(xintercept=ref, linetype=2, colour='#b0b0b0', size=0.8) +
  geom_vline(xintercept=oprlengths, linetype=2, colour='#ebebeb', size=0.8) +
  geom_segment(x=0, xend=450,
               y=0, yend=450*thrMis,
               size=1.0, linetype=2, colour='#595E60') +
  geom_point(size=0.8) +
  scale_colour_manual(values=somcols) +
  coord_cartesian(xlim=c(30,400), ylim=c(0, 100.1)) + # when cropping to 450 bp (see below), no read above 150 mismatches
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=-2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  xlab('read length') + ylab('mismatches vs template')
somMisThr

if (exportOrNo) {
  ggsave(filename=here('controlPCR', 'somMisThr.pdf'), somMisThr, width=90, height=70, units='mm')
}

# Note; the few reads above 450 bp are crazy repetitions eg. GCGCGCGCGCGC, I think fair to crop the plot at 450 bp
# (and not somatic calls as they have terrible % mismatch)




# What are the somatic calls? ---------------------------------------------
# plot calls vs counts
# calls as in eg. OPRD1, OPRI2, etc.

# somatic OPR calls are separated in 3 columns mostly to allow cases where a read is both OPRD & OPRI

# but in the calls after mismatch threshold...
subset(iedSom, !is.na(opri) & !is.na(oprd)) # no row is both OPRI/OPRD
# and by definition any read that is reference in $othercategory has to be NA in both OPRI & OPRD
# >>> so we can merge all three columns in one column
iedSom$somOPR <- as.character(poolColumns(iedSom, c('othercategory', 'opri', 'oprd')))

# check factor order as I want it
iedSom$somOPR <- factor(iedSom$somOPR, levels=c(sprintf('OPRD%i', 4:1), 'reference', sprintf('OPRI%i', 1:24)))

# now count number of calls per new OPR
sowp <- iedSom %>% # number of somatic calls by new OPR 
  group_by(somOPR) %>%
  tally(name='nSom')

# for plot, add a normal OPR bar
sowp2 <- sowp %>%
  add_row(somOPR='OPR0', nSom=NA, .after=4)

sowp2$somOPR <- factor(sowp2$somOPR, levels=c('OPRD4', 'OPRD3', 'OPRD2', 'OPRD1', 'OPR0', 'OPRI1', 'OPRI2', 'OPRI3', 'OPRI4'))

callsvsnewOPR <- ggplot(sowp2, aes(x=somOPR, y=nSom)) +
  geom_col(fill=NA, colour='#4D4D4D', width=0.8) +
  theme_minimal() +
  theme(
    panel.grid.minor.x=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t=2, r=0, b=0, l=0)),
    axis.text.x=element_text(angle=90, size=7, margin = margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  xlab('') + ylab('somatic mutation calls') +
  scale_x_discrete(labels=c('4 OPRD', '3 OPRD', '2 OPRD', '1 OPRD', '', '1 OPRI', '2 OPRI', '3 OPRI', '4 OPRI'))
callsvsnewOPR

if (exportOrNo) {
  ggsave(filename=here('controlPCR', 'callsvsnewOPR.pdf'), callsvsnewOPR, width=150, height=100, units='mm')
}

# as percentage of all reads

sowp2$somOPRper <- (sowp2$nSom / nrow(ied)) * 100


callsvsnewOPRPer <- ggplot(sowp2, aes(x=somOPR, y=somOPRper)) +
  geom_col(fill='#4D4D4D', colour='#4D4D4D', width=0.8, alpha=0.5) +
  theme_minimal() +
  theme(
    panel.grid.minor.x=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t=2, r=0, b=0, l=0)),
    axis.text.x=element_text(angle=90, size=7, margin = margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  xlab('') + ylab('somatic mutation calls (% of reads)') +
  scale_x_discrete(labels=c('4 OPRD', '3 OPRD', '2 OPRD', '1 OPRD', '', '1 OPRI', '2 OPRI', '3 OPRI', '4 OPRI'))
callsvsnewOPRPer

if (exportOrNo) {
  ggsave(filename=here('controlPCR', 'callsvsnewOPRPer.pdf'), callsvsnewOPRPer, width=120, height=70, units='mm')
}

# I do not think I can do in percentages here
# I do not see percentage of 'what' I would do...



# percentage of somatic calls ---------------------------------------------

100 * (sum(as.logical(ied$isSomatic)) / nrow(ied))

# how many 1 OPRD
nrow(iedSom)

length(which(iedSom$somOPR=='OPRD1'))



# write the somatic calls -------------------------------------------------

# export iedSom
write.xlsx(iedSom, file=here('controlPCR', 'conPCR_somaticcalls.xlsx'))