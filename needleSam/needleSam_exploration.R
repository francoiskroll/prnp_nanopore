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
ied <- read.csv(here('needleSam', 'allOPRreads.csv'))

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
  coord_cartesian(xlim=c(30,377), ylim=c(0, 100.1)) + # when cropping to 450 bp (see below), no read above 150 mismatches
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
  ggsave(filename=here('needleSam', 'somMisThr.pdf'), somMisThr, width=77, height=72, units='mm')
}

# Note; the few reads above 450 bp are crazy repetitions eg. GCGCGCGCGCGC, I think fair to crop the plot at 450 bp
# (and not somatic calls as they have terrible % mismatch)

# Note; now lowered max X to 377, it is excluding 1 read extra





# Do somatic calls have worse % mismatch? ---------------------------------
# or; are somatic reads less accurate than 'true' reads?

# calculate mismatch per bp for all somatic candidates
iedsomall$som_misbp <- iedsomall$conmis / iedsomall$read_length

# mismatch per bp in iedSom currently in three columns (ref / opri / oprd), need to pool them in one
iedSom$som_misbp <- as.numeric(unlist(poolColumns(iedSom, c('ref_misbp', 'opri_misbp', 'oprd_misbp'))))

# check reads that reads which are both in iedsomall and iedSom have the same misbp
round(iedsomall[match(iedSom$read, iedsomall$read) , 'som_misbp'], 2) == round(iedSom$som_misbp, 2)
# (need to round or it says FALSE for numbers which are clearly equal, numbers are just stored with different levels of precision)


# binning for mismatch (will use for histograms below)
misbins <- seq(0, 1, by=0.01) # bins of 1% mismatch

###
###
# first do all somatic reads (i.e. before mismatch threshold) vs trueread

# all somatic reads before mismatch filtering are in iedsomall
# counts / frequencies in bins
som <- freqInBins(df=iedsomall, datcol='som_misbp', grpcol=NA, bins=misbins)
# add column source
som$source <- 'somatic'

# same for all truereads
tru <- freqInBins(df=subset(ied, trueread), datcol='true_misbp', grpcol=NA, bins=misbins)
# add column source
tru$source <- 'trueread'

# stick the two together
somtru <- rbind(som, tru)

# plot
misTruevsSomAll <- ggplot(somtru, aes(x=bin, y=freq, fill=source)) +
  geom_col(position='identity', alpha=0.5) +
  scale_fill_manual(values=somcols) +
  theme_minimal() +
  theme(
    legend.position='none',
  ) +
  coord_cartesian(xlim=c(0, 1)) +
  xlab('mismatch per bp') + ylab('frequency')
misTruevsSomAll

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'misTruevsSomAll.pdf'), misTruevsSomAll, width=120, height=100, units='mm')
}

###
###
# now only selected somatic reads (i.e. after mismatch threshold) vs trueread

# ! need to also crop trueread by the same mismatch threshold if we want to compare frequencies
trure <- subset(ied, trueread & true_misbp < thrMis) # true reads below mismatch threshold
tru <- freqInBins(df=trure, datcol='true_misbp', grpcol=NA, bins=misbins)
# add column source
tru$source <- 'trueread'

# for somatic, use iedSom (i.e. iedsomall + mismatch threshold)
som <- freqInBins(df=iedSom, datcol='som_misbp', grpcol=NA, bins=misbins)
# add column source
som$source <- 'somatic'

# stick the two together
somtru <- rbind(som, tru)

# plot
misTruevsSomThr <- ggplot(somtru, aes(x=bin, y=freq, fill=source)) +
  geom_col(position='identity', alpha=0.5) +
  scale_fill_manual(values=somcols) +
  theme_minimal() +
  theme(
    legend.position='none',
  ) +
  coord_cartesian(xlim=c(0, thrMis+0.3*thrMis)) +
  xlab('mismatch per bp') + ylab('frequency')
misTruevsSomThr

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'misTruevsSomThr.pdf'), misTruevsSomThr, width=120, height=100, units='mm')
}


###
# in this comparison, do somatic reads have significantly worse mismatch per bp?
# compare with a t-test

# gather the reads we plotted above
trure <- select(trure, true_misbp)
colnames(trure) <- 'misbp'
trure$source <- 'trueread'

somre <- select(iedSom, som_misbp)
colnames(somre) <- 'misbp'
somre$source <- 'somatic'

misre <- rbind(trure, somre) # mismatches for reads to compare

wilcox.test(misre$misbp ~ misre$source)
t.test(misre$misbp ~ misre$source)

###
###

# within truereads, are % mismatch different reference vs mutated?
  # if yes -- above is not a fair comparison as somatic reads are almost all mutated
  # (only PDG 46345 can have somatic calls that are reference)

# take all the 'true' reads + their 'othercategory' column
  # within truereads, in 'othercategory' column, should only be reference or NA
  # can assume all NA are OPRI or OPRD in the other columns

trumut <- subset(ied, trueread & is.na(othercategory)) # true reads mutated
truref <- subset(ied, trueread & othercategory=='reference') # true reads reference

trumutb <- freqInBins(df=trumut, datcol='true_misbp', grpcol=NA, bins=misbins)
# add column source
trumutb$source <- 'trueread_mutated'

trurefb <- freqInBins(df=truref, datcol='true_misbp', grpcol=NA, bins=misbins)
# add column source
trurefb$source <- 'trueread_reference'

# stick the two together
somtru <- rbind(trurefb, trumutb)

# plot
misTrueRefvsTrueMut <- ggplot(somtru, aes(x=bin, y=freq, fill=source)) +
  geom_col(position='identity', alpha=0.5) +
  scale_fill_manual(values=somcols) +
  theme_minimal() +
  theme(
    legend.position='none',
  ) +
  coord_cartesian(xlim=c(0, 0.5)) +
  xlab('mismatch per bp') + ylab('frequency')
misTrueRefvsTrueMut

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'misTrueRefvsTrueMut.pdf'), misTrueRefvsTrueMut, width=120, height=100, units='mm')
}


# compare statistically

# gather the reads we plotted
trumut <- select(trumut, true_misbp)
trumut$source <- 'trueread_mutated'

truref <- select(truref, true_misbp)
truref$source <- 'trueread_reference'

trumis <- rbind(trumut, truref)

wilcox.test(trumis$true_misbp ~ trumis$source)
t.test(trumis$true_misbp ~ trumis$source)
# yes - it is different

###
###
# consequence: comparing somatic reads vs all true reads is unfair
  # as most true reads are reference, but most somatic reads are mutated
# to tell if somatic reads really have more mismatches, need to compare with mutated truereads

# take the mutated truereads as above + apply mismatch threshold
trumut <- subset(ied, trueread & is.na(othercategory) & true_misbp < thrMis) # true reads mutated
trumutb <- freqInBins(df=trumut, datcol='true_misbp', grpcol=NA, bins=misbins)
# add column source
trumutb$source <- 'trueread_mutated'

# take the somatic reads as above
somb <- freqInBins(df=iedSom, datcol='som_misbp', grpcol=NA, bins=misbins)
# add column source
somb$source <- 'somatic'

# stick the two together
somtru <- rbind(trumutb, somb)

# plot
misTrueMutvsSomThr <- ggplot(somtru, aes(x=bin, y=freq, fill=source)) +
  geom_col(position='identity', alpha=0.5) +
  scale_fill_manual(values=somcols) +
  theme_minimal() +
  theme(
    legend.position='none',
  ) +
  coord_cartesian(xlim=c(0, thrMis+0.3*thrMis)) +
  xlab('mismatch per bp') + ylab('frequency')
misTrueMutvsSomThr

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'misTrueMutvsSomThr.pdf'), misTrueMutvsSomThr, width=120, height=100, units='mm')
}


# compare statistically

# gather the reads we plotted
trumut <- select(trumut, true_misbp)
colnames(trumut) <- 'misbp'
trumut$source <- 'trueread_mutated'

som <- select(iedSom, som_misbp)
colnames(som) <- 'misbp'
som$source <- 'somatic'

mis <- rbind(trumut, som)

wilcox.test(mis$misbp ~ mis$source)
t.test(mis$misbp ~ mis$source)
# yes - it is different

# it did reduce the difference but excess around 5% mismatch still present
# difference in means now is 0.5%


###
###

# Can it be a coincidence?
# as we have thousands of reference reads but 'only' ~ 100--150 somatic reads

# draw from mutated truereads the same number of somatic reads multiple times and test with Mann-Whitney U Test
ndraws <- 1000
mwu <- c() # Mann-Whitney U Test p-values

for (i in 1:ndraws) {
  
  trudraw <- trumut[sample(nrow(trumut), nrow(iedSom)),]
  mis <- rbind(trudraw, som)
  
  pval <- wilcox.test(mis$misbp ~ mis$source)$p.value
  mwu <- c(mwu, pval)
  
}

length(which(mwu>0.05)) / ndraws
# it is not significant ~ 18% of the time

###
###
# In summary; do somatic calls have worse % mismatch?
# in our dataset, the answer is yes; precisely there seems to be a surplus of reads in the 5% mismatch range
# however, it is possible that we were simply 'unlucky' with the batch of somatic reads we got
# indeed, comparing our batch of somatic reads against subsets of mutated truereads only give a significant difference 82% of the time

# we could be more stringent with the threshold mismatch (lower it to mean for example) and this difference would probably go away
# however, it is important to note 5% mismatch is not surprising for reads that are correct ('trueread'), so would likely exclude some reads that are ok
# and the consequence would be to have low number of somatic calls, which would make all the comparisons between samples below difficult/impossible
# there is also the chance that it represents something real that is not allowed by the OPR consensus sequences
# I will keep the mismatch threshold as is

# Update: I think I found where the excess of ~ 5% mismatch reads come from
# it might be mostly reads where the R1 is affected
# to my knowledge this has never been reported in patients, but I do not know why strand slippage would not affect R1
# (I guess it could only be template strand slippage?)
# from http://www.mad-cow.org/prion_repeat_insertions.html#References; should now allow extra R1 or R4
# this holds in the reads I have inspected closely
# in contrast, I found at least one clean R1 deletion




# Do somatic calls have worse haplotype scores? ---------------------------

# haplotype scores is my rough attempt at putting a score for confidence of assignment to a haplotype
# see needleSam.R for comments

# if somatic calls all have bad haplotype assignment scores, would be worrying

# bin the hap scores for somatic calls
bns <- seq(0, 1, by=0.1)
hscns <- hist(as.numeric(unlist(subset(ied, isSomatic==TRUE, hapscore))), # haplotype score counts
              breaks=bns, include.lowest=TRUE, plot=FALSE)$counts
somhs <- data.frame(source='somatic', binend=bns[-1], counts=hscns, freq=hscns/sum(hscns)) # d for dataframe 

# bin the hap scores for trueread
hscns <- hist(as.numeric(unlist(subset(ied, trueread, hapscore))), # haplotype score counts
              breaks=bns, include.lowest=TRUE, plot=FALSE)$counts
truhs <- data.frame(source='trueread', binend=bns[-1], counts=hscns, freq=hscns/sum(hscns)) # d for dataframe 

# stick them together
hsf <- rbind(somhs, truhs) # haplotype score frequencies


# plot haplotype quality score vs somatic call or not
hapScoreSomatic <- ggplot(hsf, aes(x=binend, y=freq, fill=source)) +
  geom_col(position='identity', alpha=0.4) +
  theme_minimal() +
  theme(
    legend.position='none',
    axis.text.y=element_blank()
  ) +
  scale_fill_manual(values=somcols) +
  xlab('haplotype assignment score') + ylab('frequency')
hapScoreSomatic

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'hapScoreSomatic.pdf'), plot=hapScoreSomatic, width=100, height=100, units='mm')
}


# for test; select the reads we are plotting
iedhs <- ied %>% # ied to compare haplotype scores
  subset(isSomatic==TRUE | trueread)
t.test(iedhs$hapscore ~ iedhs$isSomatic)
# somatic calls do not have worse haplotype assignment scores




# Any evidence of strand bias? --------------------------------------------

# make sure it is in same order in ied & iedSom
ied$strand <- factor(ied$strand, levels=c('+', '-'))
iedSom$strand <- factor(iedSom$strand, levels=c('+', '-'))

# number of calls by strand
cbstr <- iedSom %>%
  group_by(strand) %>%
  tally(name='countStr')

# add total number of Forward/Reverse reads
cbstr <- ied %>%
  group_by(strand) %>%
  tally(name='totalStr') %>%
  left_join(cbstr, ., by='strand')

# ratio of Forward/Reverse of all reads serves as baseline
cbstr <- cbstr %>%
  mutate(ratioTot = totalStr / sum(totalStr))

# add as percentage of Forward/Reverse reads
cbstr$somPer <- (cbstr$countStr / cbstr$totalStr) * 100 # somatic percentage

# ratio of Forward/Reverse reads * total number of somatic calls will give expected counts
cbstr$expCount <- cbstr$ratioTot * sum(cbstr$countStr)

# add a dummy category for x axis below
cbstr$dummy <- factor(1)

# annoyingly, need to reverse the factor levels
cbstr$strand <- factor(cbstr$strand, levels=c('-', '+'))
# ! reverse as well the colours
strcols2 <- rev(strcols)

# plot as horizontal bar with ratios
  # essentially horizontal stacked barplot
strandPro <- ggplot(cbstr, aes(x=dummy, y=countStr, fill=strand)) +
  geom_col(position='stack') +
  coord_flip() +
  scale_fill_manual(values=strcols2) +
  geom_hline(yintercept=cbstr$expCount[2], linetype=2, colour='#ebebeb', size=0.5) +
  theme_minimal() +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    panel.grid=element_blank(),
    legend.position='none'
  )
strandPro

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'strandProportion.pdf'), strandPro, width=60, height=15, units='mm')
}


# plot number of calls by strand

# put back order
cbstr$strand <- factor(cbstr$strand, levels=c('+', '-'))

callsByStrand <- ggplot(cbstr, aes(x=strand, y=countStr, fill=strand)) +
  geom_col() +
  scale_fill_manual(values=strcols) +
  geom_segment(aes(x=0.5, xend=1.5, y=as.numeric(cbstr[1, 'expCount']), yend=as.numeric(cbstr[1, 'expCount'])), # add expected counts markers
               linetype=2, colour='#ebebeb', size=2) +
  geom_segment(aes(x=1.5, xend=2.5, y=as.numeric(cbstr[2, 'expCount']), yend=as.numeric(cbstr[2, 'expCount'])),
               linetype=2, colour='#ebebeb', size=2) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.y=element_blank()
  ) +
  ylab('somatic mutation calls')
callsByStrand

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsByStrand.pdf'), callsByStrand, width=75, height=100, units='mm')
}


# anything unexpected (would be sign of artefact) or within what we expect?
# statistical approach: Chi2 test
# build contingency table
cont <- table(ied$isSomatic, ied$strand)
  # should not do Chi2 test if any expected count below 5 (from http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r)
  # in all the reads, there is some slight surplus of Reverse reads, but we pretty much expect 50/50 Forward/Reverse in the somatic calls
  # we have 139 calls, so around 70/70 Forward/Reverse, so no expected count will be below 5
# now Chi2 test on contingency table
chi <- chisq.test(cont, correct=FALSE)
chi
chi$expected
chi$observed


# another approach I like to visualise this:
  # draw reads at random & see how many Forward/Reverse we get
    # ! better to do the draws properly sample by sample (if any sample has some strand bias)
    # eg. PDG 58398 has 24 calls >> draw at random 24 reads from all of its reads

# how many somatic calls per sample?
cps <- iedSom %>% # calls per sample
  group_by(sample) %>%
  tally(name='nSom')

# function to run the random draws and output the data ready to plot
strandDraws <- function(nsim) {
  
  # preallocate dataframe
  simres <- as.data.frame(matrix(nrow=nsim, ncol=3)) # simulation results;
    # rows = simulations, columns = simulation # / number of forward draws / number of reverse draws
  colnames(simres) <- c('sim', 'plus', 'minus')
  simres$sim <- 1:nrow(simres)
  
  # for each simulation;
  # for each sample, how many reads do we need to draw? (i.e. number of somatic calls for this sample)
  # draw that many reads (! from that sample)
  # pool all the draws together
  # count how many Forward, Reverse we got
  # add the results to the simres
  
  for (i in 1:nsim) {
    
    cat ('\t \t \t \t Simulation', i, 'out of', nsim, '\n')
    
    stradra <- c() # strand drawn (will pool all the draws there)
    
    for (p in 1:length(cps$sample)) { # for each sample (or row of cps)
      spl <- as.character(cps[p, 'sample']) # sample we are talking about
      ndraw <- as.integer(cps[p, 'nSom']) # how many reads do we need to draw?
      drawfrom <- as.character(unlist(subset(ied, sample==spl, strand))) # ! only draw from that sample's reads
      stradra <- c(stradra, sample ( drawfrom , ndraw)) # draw that many reads from that sample's reads & add them to hapdra
    }
    
    # length of stradra now should be = total number of somatic calls
    if (length(stradra) != nrow(iedSom)) stop ('\t \t \t \t >>> Error: something wrong with the draws \n')
    
    # count how many of each strand we got, and add that to our simulation results
    simres[i, 'plus'] <- sum(stradra == '+') # how many forward strand?
    simres[i, 'minus'] <- sum(stradra == '-') # how many reverse strand?
  }
  
  return(simres)
}

nsim <- 100 # ideally 1000, but takes a couple of minutes
stradraws <- strandDraws(nsim=nsim)

# real count of somatic calls from Forward strand
nplus <- nrow(subset(iedSom, strand=='+'))

# real count of somatic calls from Reverse strand
nminus <- nrow(subset(iedSom, strand=='-'))

# plot for each strand
plusDraws <- ggplot(stradraws, aes(plus)) +
  geom_histogram(binwidth=1) +
  geom_vline(xintercept=nplus, linetype=2) +
  theme_minimal() +
  theme(
    axis.text.y=element_blank()
  ) +
  coord_cartesian(xlim=c(25, 100), ylim=c(0, 80)) +
  ylab('frequency') + xlab('somatic calls from forward strand')
plusDraws

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'plusDraws.pdf'), plusDraws, width=100, height=100, units='mm')
}


# plot for each strand
minusDraws <- ggplot(stradraws, aes(minus)) +
  geom_histogram(binwidth=1) +
  geom_vline(xintercept=nplus, linetype=2) +
  theme_minimal() +
  theme(
    axis.text.y=element_blank()
  ) +
  coord_cartesian(xlim=c(25, 100), ylim=c(0, 80)) +
  ylab('frequency') + xlab('somatic calls from reverse strand')
minusDraws
ggsave(filename=here('needleSam', 'minusDraws.pdf'), plusDraws, width=100, height=100, units='mm')




# Somatic mutations calls -- from which cohort? ---------------------------

# make sure it is in this order in ied / iedSom
ied$prion_disease <- factor(ied$prion_disease, levels=c('inherited', 'sCJD', 'control'))
iedSom$prion_disease <- factor(iedSom$prion_disease, levels=c('inherited', 'sCJD', 'control'))

# v1: removing calls from brain samples here in order to keep one sample per individual
# but does not match number in text (129 calls) if do this

# calls by cohort
cbyco <- iedSom %>%
  group_by(prion_disease) %>%
  tally(name='nSom')

# first observation: majority from Inherited cases
nrow(subset(iedSom, prion_disease=='inherited')) / nrow(iedSom) # all samples

cat('\t \t \t \t >>> ', nrow(subset(iedSom, prion_disease=='inherited')), '/', nrow(iedSom), 'calls = ',
    round(nrow(subset(iedSom, prion_disease=='inherited')) / nrow(iedSom) * 100, 2),
    '% from Inherited \n')

# add total number of reads from each cohort
cbyco <- ied %>%
  group_by(prion_disease) %>%
  tally(name='totalCoh') %>%
  left_join(cbyco, ., by='prion_disease')

# somatic calls as percentage of reads from that cohort
cbyco$somPer <- (cbyco$nSom / cbyco$totalCoh) * 100 # somatic percentage

# add ratio each cohort of total number of reads
cbyco <- cbyco %>%
  mutate(proCoh=totalCoh/sum(totalCoh))

# ratio each cohort * total number of somatic calls will give expected counts
cbyco$expCount <- cbyco$proCoh * sum(cbyco$nSom)

# can also add expected percentage
# Note, this is same for all cohorts (as expected counts basically assumes always same percentage)
cbyco$expPer <- (cbyco$expCount / cbyco$totalCoh) * 100

# plot histogram of counts
callsbyCoh <- ggplot(cbyco, aes(x=prion_disease, y=nSom, fill=prion_disease)) +
  geom_col(width=0.7) +
  # add expected counts markers
  # geom_segment(aes(x=0.5, xend=1.5, y=as.numeric(cbyco[1, 'expCount']), yend=as.numeric(cbyco[1, 'expCount'])), linetype=2, colour='#ebebeb', size=2) +
  # geom_segment(aes(x=1.5, xend=2.5, y=as.numeric(cbyco[2, 'expCount']), yend=as.numeric(cbyco[2, 'expCount'])), linetype=2, colour='#ebebeb', size=2) +
  # geom_segment(aes(x=2.5, xend=3.5, y=as.numeric(cbyco[3, 'expCount']), yend=as.numeric(cbyco[3, 'expCount'])), linetype=2, colour='#ebebeb', size=2) +
  scale_fill_manual(values=cohcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    axis.title.x=element_blank(), 
    axis.text.x=element_text(size=9, margin = margin(t=0, r=0, b=0, l=0), angle=90),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  scale_x_discrete(labels=c('inherited', 'sporadic', 'control')) +
  ylab('somatic mutation calls')
callsbyCoh

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsbyCohort.pdf'), callsbyCoh, width=40, height=75, units='mm')
}



# same but percentages
callsbyCohPer <- ggplot(cbyco, aes(x=prion_disease, y=somPer, fill=prion_disease)) +
  geom_col() +
  # add expected percentage
  geom_hline(yintercept=mean(cbyco$expPer), linetype=2, size=2, colour='#ebebeb') +
  scale_fill_manual(values=cohcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    axis.title.x=element_blank(),
  ) +
  scale_x_discrete(labels=c('inherited', 'sporadic', 'control')) +
  ylab('somatic mutation calls (% of reads)')
callsbyCohPer

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsbyCohortPercentage.pdf'), callsbyCohPer, width=100, height=100, units='mm')
}


# Is there really a surplus of Inherited calls?
# Chi2 test
  # build contingency table
# v1: removing brain samples from ied
# v2: keeping all samples
cont <- table(ied$prion_disease, ied$isSomatic)
chi <- chisq.test(cont, correct=FALSE)
chi
# ! check no expected count below 5
chi$expected
chi$observed

# alternative approach to visualise this
  # by drawing reads at random and counting how many we have from each cohort

# function to run the random draws and output the data ready to plot
cohortDraws <- function(nsim) {
  
  # preallocate dataframe
  simres <- as.data.frame(matrix(nrow=nsim, ncol=4)) # simulation results;
  # rows = simulations, columns = simulation # / number of inherited draws / number of sCJD draws / number of control draws
  colnames(simres) <- c('sim', 'inherited', 'sCJD', 'control')
  simres$sim <- 1:nrow(simres)
  
  # for each simulation;
  # for each sample, how many reads do we need to draw? (i.e. number of somatic calls for this sample)
  # draw that many reads (! from that sample)
  # pool all the draws together
  # count how many inherited, sCJD, control
  # add the results to the simres
  
  for (i in 1:nsim) {
    
    cat ('\t \t \t \t Simulation', i, 'out of', nsim, '\n')
    cohdra <- sample(ied$prion_disease, nrow(iedSom)) # cohorts drawn
    
    # count how many of each cohort we got, and add that to our simulation results
    simres[i, 'inherited'] <- sum(cohdra == 'inherited') # how many inherited?
    simres[i, 'sCJD'] <- sum(cohdra == 'sCJD') # how many sCJD?
    simres[i, 'control'] <- sum(cohdra == 'control') # how many sCJD?
  }
  return(simres)
}

nsim <- 1000
cohdraws <- cohortDraws(nsim=nsim)

# real count of somatic calls from inherited
ninh <- nrow(subset(iedSom, prion_disease=='inherited'))

# real count of somatic calls from sCJD
nspo <- nrow(subset(iedSom, prion_disease=='sCJD'))

# real count of somatic calls from control
ncon <- nrow(subset(iedSom, prion_disease=='control'))

# plot for inherited
inhDraws <- ggplot(cohdraws, aes(inherited)) +
  geom_histogram(binwidth=1) +
  geom_vline(xintercept=ninh, linetype=2) +
  theme_minimal() +
  theme(
    axis.text.y=element_blank()
  ) +
  coord_cartesian(xlim=c(0, 150), ylim=c(0, 100)) +
  ylab('frequency') + xlab('somatic calls from inherited samples')
inhDraws

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'inheritedDraws.pdf'), inhDraws, width=100, height=100, units='mm')
}


# plot for sporadic
spoDraws <- ggplot(cohdraws, aes(sCJD)) +
  geom_histogram(binwidth=1) +
  geom_vline(xintercept=nspo, linetype=2) +
  theme_minimal() +
  theme(
    axis.text.y=element_blank()
  ) +
  coord_cartesian(xlim=c(0, 150), ylim=c(0, 100)) +
  ylab('frequency') + xlab('somatic calls from sporadic samples')
spoDraws

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'sporadicDraws.pdf'), spoDraws, width=100, height=100, units='mm')
}


# plot for control
conDraws <- ggplot(cohdraws, aes(control)) +
  geom_histogram(binwidth=1) +
  geom_vline(xintercept=ncon, linetype=2) +
  theme_minimal() +
  theme(
    axis.text.y=element_blank()
  ) +
  coord_cartesian(xlim=c(0, 150), ylim=c(0, 100)) +
  ylab('frequency') + xlab('somatic calls from control samples')
conDraws

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'controlDraws.pdf'), conDraws, width=100, height=100, units='mm')
}





# Somatic mutations calls -- from which samples? --------------------------

# will prepare a dataframe with lots of info for each sample
  # first compute coverage (will use below)
cabs <- ied %>% # number of CAlls By Sample
  group_by(sample) %>%
  tally(name='iedCov') # coverage = simply number of rows per sample in ied

# add all the metadata info
cabs <- left_join(cabs, meta, by='sample')

# check one row per sample
nrow(cabs) == length(unique(ied$sample))

nsomps <- iedSom %>% # number of somatic mutation calls per sample
  group_by(sample) %>%
  tally(name='nSom')

# join this to cabs
cabs <- left_join(x=cabs, y=nsomps, by='sample')

# if sample not in iedSom it means it had 0 calls
  # it should be NA for column nSom in cabs, but check to be safe
nocp <- unique(ied$sample) [ !unique(ied$sample) %in% unique(iedSom$sample) ] # no-call PDGs

as.character(unlist(subset(cabs, is.na(nSom), sample))) %in% nocp
nocp %in% as.character(unlist(subset(cabs, is.na(nSom), sample)))
# ok

# replace cabs column nSom NA in cabs by 0
cabs[which(is.na(cabs$nSom)), 'nSom'] <- 0
# check no more NA in cabs nSom
sum(is.na(cabs$nSom))

# compute % of reads somatic
  # i.e. number of somatic calls / coverage
cabs$somPro <- cabs$nSom / cabs$iedCov
cabs$somPer <- cabs$somPro * 100

# make sure levels of prion_disease are in order of cohcols above
cabs$prion_disease <- factor(cabs$prion_disease, levels=c('inherited', 'sCJD', 'control'))

# order from sample with largest number of somatic calls to lowest
cabs <- cabs[order(-cabs$nSom),]
cabs$sample <- factor(cabs$sample, levels=cabs$sample)

# plot histogram sample vs number of somatic calls

# temporarily replace 0 by a low number (0.1) so it plots the colour of the sample in plot
cabs2 <- cabs
cabs2$nSom <- as.numeric(cabs2$nSom)
cabs2[which(cabs2$nSom==0) , 'nSom'] <- 0.1

cabsple <- ggplot(cabs2, aes(x=sample, y=nSom, fill=prion_disease)) +
  geom_col() +
  scale_fill_manual(values=cohcols, labels=c('inherited', 'sporadic', 'control')) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(angle=90),
    legend.title=element_blank()
  ) +
  ylab('somatic mutation calls')
cabsple

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsbySample.pdf'), cabsple, width=200, height=150, units='mm')
}



# plot with percentages

# order from sample with largest percentage of somatic calls to lowest
cabs <- cabs[order(-cabs$somPer),]
cabs$sample <- factor(cabs$sample, levels=cabs$sample)

# as above, raise a little bit the samples with 0 somatic calls so can see them
cabs2 <- cabs
cabs2$nSom <- as.numeric(cabs2$nSom)
cabs2[which(cabs2$somPer==0) , 'somPer'] <- 0.002

spllabels <- paste0('#', substr(cabs2$sample, 4, 99))

cabsple2 <- ggplot(cabs2, aes(x=sample, y=somPer, fill=prion_disease, alpha=gDNA_tissue)) +
  geom_col(width=0.8) +
  scale_alpha_discrete(range=c(1.0, 0.5)) +
  scale_fill_manual(values=cohcols, labels=c('inherited', 'sporadic', 'control')) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    legend.title=element_blank(),
    legend.position='none',
    axis.text.x=element_text(size=7, angle=90, margin = margin(t=0, r=0, b=0, l=0), hjust=0, vjust=0.5),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  ylab('somatic mutation calls (% of reads)') +
  scale_x_discrete(labels=spllabels)
cabsple2

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsbySamplePer.pdf'), cabsple2, width=80, height=70, units='mm')
}


# version with each sample's colour

# samples are ranked by % somatic, so cannot change that
# need to match the colour to each sample
splcols <- meta$colour[match(cabs$sample, meta$sample)]

cabsple3 <- ggplot(cabs2, aes(x=sample, y=somPer, fill=sample)) +
  geom_col() +
  scale_fill_manual(values=splcols) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(angle=90),
    legend.title=element_blank(),
    legend.position='none'
  ) +
  ylab('somatic mutation calls (% of reads)')
cabsple3

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsbySamplePer_col.pdf'), cabsple3, width=200, height=150, units='mm')
}


# Number of calls vs coverage ---------------------------------------------
# does number of calls for each sample correlate with coverage?
  # would not make the calls wrong (it would be expected we pick up more stuff with higher coverage)
  # but I think would mean we are not getting a full picture
callsvscoverage <- ggplot(cabs, aes(x=iedCov, y=nSom, colour=prion_disease)) +
  geom_point(size=1.5) +
  scale_colour_manual(values=cohcols) +
  theme_minimal() +
  theme(
    legend.title=element_blank()
  ) +
  theme(
    legend.position='none',
    #panel.grid.major.x=element_blank(),
    #panel.grid.minor.x=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  xlab('coverage') + ylab('somatic mutations calls') +
  scale_x_continuous(breaks=c(0, 10000, 20000, 30000), labels=c('0', '10,000', '20,000', '30,000')) +
  coord_cartesian(xlim=c(0, 30000))
callsvscoverage

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvscoverage.pdf'), callsvscoverage, width=60, height=60, units='mm')
}


cor(x=cabs$iedCov, y=cabs$nSom, method='pearson')
cor.test(x=cabs$iedCov, y=cabs$nSom, method='pearson')

cor(x=cabs$iedCov, y=cabs$nSom, method='spearman')
cor.test(x=cabs$iedCov, y=cabs$nSom, method='spearman')


# useful to have a version with labels
callsvscoverage2 <- ggplot(cabs, aes(x=iedCov, y=nSom, colour=prion_disease, label=sample)) +
  geom_point(size=2) +
  geom_text(size=1.8, hjust=1.0, vjust=2.0) +
  scale_colour_manual(values=cohcols) +
  theme_minimal() +
  theme(
    legend.title=element_blank()
  ) +
  xlab('coverage') + ylab('somatic mutations calls')
callsvscoverage2

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvscoverage_labels.pdf'), callsvscoverage2, width=120, height=100, units='mm')
}




# Number of calls vs flowcell ---------------------------------------------
  # anything weird about flowcell?
callsvsflowcell <- ggplot(cabs, aes(x=flowcell, y=nSom, colour=prion_disease)) +
  geom_quasirandom(width=0.1) +
  scale_colour_manual(values=cohcols) +
  scale_x_continuous(breaks=1:3) +
  theme_minimal() +
  theme(
    legend.title=element_blank()
  ) +
  ylab('somatic mutation calls')
callsvsflowcell

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsflowcell.pdf'), callsvsflowcell, width=80, height=100, units='mm')
}


# with percentages
callsvsflowcellPer <- ggplot(cabs, aes(x=flowcell, y=somPer, colour=prion_disease)) +
  geom_quasirandom(width=0.1) +
  scale_colour_manual(values=cohcols) +
  scale_x_continuous(breaks=1:3) +
  theme_minimal() +
  theme(
    legend.title=element_blank()
  ) +
  ylab('somatic mutation calls')
callsvsflowcellPer
ggsave(filename=here('needleSam', 'callsvsflowcellPercentage.pdf'), callsvsflowcellPer, width=80, height=100, units='mm')



# Number of calls vs tissue -----------------------------------------------
# first version with all the samples

callsvstissue <- ggplot(cabs, aes(x=gDNA_tissue, y=nSom, label=sample, colour=prion_disease)) +
  geom_quasirandom(width=0.2) +
  scale_colour_manual(values=cohcols, labels=c('inherited', 'sporadic', 'control')) +
  theme_minimal() +
  theme(
    legend.title=element_blank(),
    panel.grid.minor.x=element_blank(),
  ) +
  coord_cartesian(ylim=c(0,18)) +
  xlab('') + ylab('somatic calls')
callsvstissue

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvstissue.pdf'), callsvstissue, width=90, height=100, units='mm')
}


# with percentages
callsvstissuePer <- ggplot(cabs, aes(x=gDNA_tissue, y=somPer, label=sample, colour=prion_disease)) +
  geom_quasirandom(width=0.2) +
  scale_colour_manual(values=cohcols, labels=c('inherited', 'sporadic', 'control')) +
  theme_minimal() +
  theme(
    legend.title=element_blank(),
    panel.grid.minor.x=element_blank(),
  ) +
  xlab('') + ylab('somatic calls')
callsvstissuePer
ggsave(filename=here('needleSam', 'callsvstissuePercentage.pdf'), callsvstissuePer, width=90, height=100, units='mm')

# Chi2 test: is there a surplus/deficit of calls in brain or blood
cont <- table(ied$isSomatic, ied$gDNA_tissue)
chi <- chisq.test(cont, correct=FALSE)
chi
chi$observed
chi$expected # none below 5
# there are *more* calls from brain samples than we would have expected
# or less calls from blood samples than we would have expected

# this does not mean much
# ! we will see below OPR length predicts well number of somatic calls
# and brain samples are all long OPRs, so expect more somatic calls

# to demonstrate this, removing samples which have OPRD or reference OPR before Chi2 would be predicted to cancel out the effect
# = keep only OPRI samples
iedI <- subset(ied, substr(SV, 1, 4) == 'OPRI')
cont <- table(iedI$isSomatic, iedI$gDNA_tissue)
chi <- chisq.test(cont, correct=FALSE)
chi
chi$observed
chi$expected # none below 5
# correct, p-val is not significant anymore

# plot matched blood v brain samples as slope plot
# prepare data so one row for each individual
cabsma <- subset(cabs, cohort=='inherited_tissuecomparison') # cabs matched samples

# order them by individual then by tissue
cabsma <- cabsma[order(cabsma$individual, cabsma$gDNA_tissue),]

# and use that as factor level
cabsma$sample <- factor(cabsma$sample, levels=cabsma$sample)
# also fix factor level of gDNA tissue
cabsma$gDNA_tissue <- factor(cabsma$gDNA_tissue, levels=unique(cabsma$gDNA_tissue))

# just for the purpose of the plot, need to slightly move apart two dots to avoid overplotting but keep the connecting lines correct
tissueMatch <- ggplot(cabsma, aes(x=gDNA_tissue, y=nSom, fill=sample)) +
  geom_line(aes(group=individual), size=1.0, colour='#a7a7a7') +
  scale_fill_manual(values=cabsma2$colour) +
  geom_quasirandom(pch=21, size=5.0, colour='white', stroke=1.5, width=0.8, groupOnX=FALSE) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    legend.position='none'
  ) +
  coord_cartesian(ylim=c(0, 20)) +
  ylab('somatic mutations calls')
tissueMatch

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'tissueMatch.pdf'), tissueMatch, width=90, height=100, units='mm')
}

# with percentages
tissueMatchPer <- ggplot(cabsma, aes(x=gDNA_tissue, y=somPer, fill=sample)) +
  geom_line(aes(group=individual), size=1.0, colour='#a7a7a7') +
  scale_fill_manual(values=cabsma$colour) +
  coord_cartesian(ylim=c(0,0.3)) +
  geom_quasirandom(pch=21, size=5.0, colour='white', stroke=1.5, width=0.8, groupOnX=FALSE) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    legend.position='none',
    axis.text.x=element_text(size=9, margin = margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  ylab('somatic mutations calls (% of reads)')
tissueMatchPer

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'tissueMatchPercentage.pdf'), tissueMatchPer, width=60, height=68, units='mm')
}

# one-sample t-test

# make a simpler dataframe for this
tisco <- cabsma %>% # tissue comparison
  select(individual, gDNA_tissue, somPer) %>%
  pivot_wider(names_from=gDNA_tissue,
             values_from=somPer) %>%
  mutate(tissuediff=blood-brain) # compute delta blood - brain

# one-sample t-test: I am testing if the % decrease is significant
  # so null hypothesis = no decrease = 0% (hence one-sample)
  # I am specifically testing decrease, so one-sided
  # why specifically decrease: strand slippage model relies on replication,
    # hence should occur more often in tissues that divides a lot vs post-mitotic tissues
t.test(x=tisco$tissuediff, mu=0, alternative='greater')

# mean +- sd of difference
mean(tisco$tissuediff)
sd(tisco$tissuediff)


# Number of calls vs gender -----------------------------------------------

callsvsgender <- ggplot(cabs, aes(x=gender, y=nSom, colour=prion_disease)) +
  geom_quasirandom(width=0.2) +
  scale_colour_manual(values=cohcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.x=element_blank(),
    axis.title.x=element_blank()
  ) +
  scale_x_discrete(labels=c('female', 'male')) +
  ylab('somatic mutation calls')
callsvsgender

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsgender.pdf'), callsvsgender, width=90, height=100, units='mm')
}

# with percentages
callsvsgenderPer <- ggplot(cabs, aes(x=gender, y=somPer, colour=prion_disease)) +
  geom_quasirandom(width=0.2) +
  scale_colour_manual(values=cohcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.x=element_blank(),
    axis.title.x=element_blank()
  ) +
  scale_x_discrete(labels=c('female', 'male')) +
  ylab('somatic mutation calls (% of reads)')
callsvsgenderPer

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsgenderPercentage.pdf'), callsvsgenderPer, width=90, height=100, units='mm')
}



# Number of calls vs Sanger OPR -------------------------------------------

# by strand slippage model, would predict
  # more repeats, more opportunities for slippage, so more calls
  # less repeats, less opportunities for slippage, so fewer calls

# add a 'OPR difference' column
  # eg. OPRI1 = +1 / OPRD1 = -1
cabs$sanOPRdiff <- sapply(cabs$SV, function(o){ # Sanger OPR difference (eg. 1 OPRD = -1)
  
  if (is.na(o)) {
    return(0) # is NA, that means no SV, so OPRdiff = 0
  } else {
    diff <- as.numeric(substr(o, 5, 999)) # get the number; eg. OPRI8 >> 8
    
    if (substr(o, 4, 4)=='D') return(-diff) # if OPRD >> put a negative sign
    if (substr(o, 4, 4)=='I') return(diff) # if OPRI >> leave positive
  }
})
# Note; for PDG 46345 writing +5 for the 5 OPRI
# Note; could also do total OPR length (i.e. Sanger OPR difference + 5)
  # in a way it fits better the strand slippage model as what is determinant would be the total number of repeats
  # but I find the difference more intuitive as can still read the sample's genotype

# number of somatic calls vs Sanger OPR difference
callsvsSanOPRdiff <- ggplot(cabs, aes(x=sanOPRdiff, y=nSom, colour=prion_disease)) +
  geom_quasirandom(width=0.3) +
  scale_colour_manual(values=cohcols) +
  scale_x_continuous(breaks=-2:8) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.x=element_blank()
  ) +
  xlab('OPR change vs reference') + ylab('somatic mutation calls')
callsvsSanOPRdiff
cor(x=cabs$sanOPRdiff, y=cabs$nSom, method='pearson')
cor.test(x=cabs$sanOPRdiff, y=cabs$nSom, method='pearson')
cor(x=cabs$sanOPRdiff, y=cabs$nSom, method='spearman')
cor.test(x=cabs$sanOPRdiff, y=cabs$nSom, method='spearman')

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsSangerOPRdiff.pdf'), callsvsSanOPRdiff, width=120, height=100, units='mm')
}

# add linear regression
# ! lm (y ~ x)
# x = predictor // y = to be predicted
# here: x = OPR length // y = number of somatic calls
somfit <- lm(nSom ~ sanOPRdiff, data=cabs)
# formula is y = intercept + slope * x
# here: number of somatic calls = 3.1803 + 0.8938 * OPR difference

intercept <- somfit[[1]][1]
slope <- somfit[[1]][2]

# X start
# is shortest OPR
xstart <- min(cabs$sanOPRdiff)

# X stop
# is longest OPR
xstop <- max(cabs$sanOPRdiff)

# Y start
# >> y = intercept + slope * xstart
ystart = as.numeric(intercept + slope * xstart)

# Y stop
# >> y = intercept + slope * xstop
ystop = as.numeric(intercept + slope * xstop)

# Notes;
summary(somfit)

# version plot above with fitted line
callsvsSanOPRdiffFit <- ggplot(cabs, aes(x=sanOPRdiff, y=nSom, colour=prion_disease)) +
  geom_segment(aes(x=xstart, xend=xstop,
                   y=ystart, yend=ystop),
               linetype=2, colour='#ebebeb', size=1.0) +
  geom_quasirandom(width=0.3) +
  scale_colour_manual(values=cohcols) +
  scale_x_continuous(breaks=-2:8) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.x=element_blank()
  ) +
  xlab('OPR change vs reference') + ylab('somatic mutation calls')
callsvsSanOPRdiffFit

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsSangerOPRdiffFit.pdf'), callsvsSanOPRdiffFit, width=120, height=100, units='mm')
}

###

# same but with percentages
callsvsSanOPRdiffPer <- ggplot(cabs, aes(x=sanOPRdiff, y=somPer, colour=prion_disease)) +
  geom_quasirandom(width=0.3) +
  scale_colour_manual(values=cohcols) +
  scale_x_continuous(breaks=-2:8) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.x=element_blank()
  ) +
  xlab('OPR length vs reference') + ylab('somatic mutation calls (% of reads)')
callsvsSanOPRdiffPer
cor(x=cabs$sanOPRdiff, y=cabs$somPer, method='pearson')
cor.test(x=cabs$sanOPRdiff, y=cabs$somPer, method='pearson')
cor(x=cabs$sanOPRdiff, y=cabs$somPer, method='spearman')
cor.test(x=cabs$sanOPRdiff, y=cabs$somPer, method='spearman')

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsSangerOPRdiffPercentage.pdf'), callsvsSanOPRdiffPer, width=120, height=100, units='mm')
}
  
# add linear regression
# ! lm (y ~ x)
# x = predictor // y = to be predicted
# here: x = OPR length // y = number of somatic calls
somfit <- lm(somPer ~ sanOPRdiff, data=cabs)
# formula is y = intercept + slope * x
# here: number of somatic calls = 3.1803 + 0.8938 * OPR difference

intercept <- somfit[[1]][1]
slope <- somfit[[1]][2]

# X start
# is shortest OPR
xstart <- min(cabs$sanOPRdiff)

# X stop
# is longest OPR
xstop <- max(cabs$sanOPRdiff)

# Y start
# >> y = intercept + slope * xstart
ystart = as.numeric(intercept + slope * xstart)

# Y stop
# >> y = intercept + slope * xstop
ystop = as.numeric(intercept + slope * xstop)

# Notes;
summary(somfit)

# version plot above with fitted line
callsvsSanOPRdiffPerFit <- ggplot(cabs, aes(x=sanOPRdiff, y=somPer, colour=prion_disease, alpha=gDNA_tissue)) +
  geom_segment(aes(x=xstart, xend=xstop,
                   y=ystart, yend=ystop),
               linetype=2, colour='#ebebeb', size=1.0) +
  geom_quasirandom(width=0.3, size=1.2) +
  scale_colour_manual(values=cohcols) +
  scale_alpha_discrete(range=c(1.0, 0.5)) +
  scale_x_continuous(breaks=-2:8) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.x=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t=2, r=0, b=0, l=0)),
    axis.text.x=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  xlab('OPR length vs reference') + ylab('somatic mutation calls (% of reads)')
callsvsSanOPRdiffPerFit

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsSangerOPRdiffPercentageFit.pdf'), callsvsSanOPRdiffPerFit, width=80, height=70, units='mm')
}

# version with age labels
  # exploring more a possible relationship below
callsvsSanOPRdiffPerFitAgelabels <- ggplot(cabs, aes(x=sanOPRdiff, y=somPer, colour=prion_disease, label=age_atsample)) +
  geom_segment(aes(x=xstart, xend=xstop,
                   y=ystart, yend=ystop),
               linetype=2, colour='#ebebeb', size=1.0) +
  geom_quasirandom(width=0.3) +
  geom_text() +
  scale_colour_manual(values=cohcols) +
  scale_x_continuous(breaks=-2:8) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.x=element_blank()
  ) +
  xlab('OPR change vs reference') + ylab('somatic mutation calls (% of reads)')
callsvsSanOPRdiffPerFitAgelabels

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsSanOPRdiffAgelabelsPercentageFit.pdf'), callsvsSanOPRdiffPerFitAgelabels, width=120, height=100, units='mm')
}



# Number of calls vs Sanger OPR -- haplotype resolution -------------------

# same idea as above, but treating haplotypes as separate samples
# re-create cabs but including haplotype info
cabsh <- ied %>% # number of CAlls By Sample by Haplotype
  group_by(sample, haplotype) %>%
  tally(name='iedCov') %>% # coverage = simply number of rows per sample in ied
  mutate(sam_hap=paste(sample, haplotype, sep='_'), .before='sample') # add sample_haplotype ID column

# add all the metadata info
  # ! naturally, will be repeated for each haplotype
cabsh <- left_join(cabsh, meta, by='sample')

nsomps <- iedSom %>% # number of somatic mutation calls per sample per haplotype
  group_by(sample, haplotype) %>%
  tally(name='nSom') %>%
  mutate(sam_hap=paste(sample, haplotype, sep='_'), .before='sample') %>% # add sample_haplotype ID column
  ungroup() %>%
  select(-sample, -haplotype)

# join this to cabsh by matching sam_hap
cabsh <- left_join(x=cabsh, y=nsomps, by='sam_hap')

# nSom = NA actually means 0 somatic calls
cabsh[which(is.na(cabsh$nSom)) , 'nSom'] <- 0

# compute % of reads somatic
  # there are 2 options
  # 1- % of reads from that haplotype
  # 2- % of reads from that sample
# I think best to do 1- so haplotype proportions do not contribute
# i.e. number of somatic calls / coverage
cabsh$somPro <- cabsh$nSom / cabsh$iedCov
cabsh$somPer <- cabsh$somPro * 100

# make sure levels of prion_disease are in order of cohcols above
cabsh$prion_disease <- factor(cabsh$prion_disease, levels=c('inherited', 'sCJD', 'control'))

# order from sam_hap with largest number of somatic calls to lowest
cabsh <- cabsh[order(-cabsh$nSom),]
cabsh$sam_hap <- factor(cabsh$sam_hap, levels=cabsh$sam_hap)

# add sanger OPR difference column
cabsh$sanOPRdiff <- sapply(cabsh$SV, function(o){ # Sanger OPR difference (eg. 1 OPRD = -1)
  
  if (is.na(o)) {
    return(0) # is NA, that means no SV, so OPRdiff = 0
  } else {
    diff <- as.numeric(substr(o, 5, 999)) # get the number; eg. OPRI8 >> 8
    
    if (substr(o, 4, 4)=='D') return(-diff) # if OPRD >> put a negative sign
    if (substr(o, 4, 4)=='I') return(diff) # if OPRI >> leave positive
  }
})
# ! remember compound homozygous sample will show up as +5 (not -1)

# now do plot Sanger OPR vs number of calls again, but essentially treating separate haplotypes as separate samples
# number of somatic calls per haplotype vs Sanger OPR
cabsh$haplotype <- factor(cabsh$haplotype, levels=c(1, 2, 0))

callsvsSanOPRdiffHap <- ggplot(cabsh, aes(x=sanOPRdiff, y=nSom, colour=haplotype)) +
  geom_quasirandom(width=0.2) +
  scale_colour_manual(values=hapcols) +
  scale_x_continuous(breaks=-2:8) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.x=element_blank()
  ) +
  xlab('OPR change vs difference') + ylab('somatic mutation calls')
callsvsSanOPRdiffHap

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsSanOPRdiffHap.pdf'), callsvsSanOPRdiffHap, width=120, height=100, units='mm')
}

# with percentages
callsvsSanOPRdiffHapPer <- ggplot(cabsh, aes(x=sanOPRdiff, y=somPer, colour=haplotype)) +
  geom_quasirandom(width=0.2) +
  scale_colour_manual(values=hapcols) +
  scale_x_continuous(breaks=-2:8) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.x=element_blank()
  ) +
  xlab('OPR change vs difference') + ylab('somatic mutation calls')
callsvsSanOPRdiffHapPer
ggsave(filename=here('needleSam', 'callsvsSanOPRdiffHapPercentage.pdf'), callsvsSanOPRdiffHapPer, width=120, height=100, units='mm')

# Any striking if plotting sample per sample?
  # could be for example an OPRI sample where the reference haplotype is affected
# set OPRdiff column as factor for this, ggplot understands better what I am doing
cabsh$sanOPRdiff <- factor(cabsh$sanOPRdiff)

# list of plots, one per sample
hapod <- vector(mode='list', length=length(pdgs)) # one plot per sample, each will be histogram number of calls per haplotype

for (p in 1:length(pdgs)) {
  
  # get the Sanger genotype for that sample
  osan <- as.character(unlist(subset(meta, sample==pdgs[p], SV))) # OPR Sanger
  
  # if osan is NA >> it means reference
  if (is.na(osan)) {
    osan <- 'reference'
  }
  
  haplot <- ggplot( subset(cabsh, sample==pdgs[p]), aes(x=haplotype, y=nSom, fill=haplotype) ) +
    geom_col(position='dodge') +
    scale_fill_manual(values=hapcols) +
    theme_minimal() +
    ggtitle(osan) +
    coord_cartesian(ylim=c(0, 12)) +
    scale_y_continuous(breaks=seq(0, 12, 2)) +
    ylab('somatic mutation calls')
  
  # place the plot in the list
  hapod[[p]] <- haplot
  
}

hapodGrid <- ggarrange(plotlist=hapod, ncol=1, nrow=length(pdgs))

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsVsHap_bysample.pdf'), plot=hapodGrid, height=1000, width=100, units='mm')
}



# Is OPR mutated haplotype preferentially affected? -----------------------

# from needleSam.R: OPR mutated haplotype (if any) is set to 1
# Note; PDG 46345, 5 OPRI haplotype is considered as mutated
# ! here, should exclude samples where haplotype 1/2 distinction is not meaningful
  # i.e. the samples we could not haplotag (0 SNV) + the samples with reference OPR
pdgop <- as.character(unlist(subset(meta, !is.na(SV), sample))) # samples with a mutated OPR
pdgoptag <- intersect(pdgop, pdgs_tag) # samples with a mutated OPR + haplotagged

cabh <- iedSom %>% # calls by haplotype
  subset(sample %in% pdgoptag) %>% # only somatic calls from samples with mutated OPR + haplotagged
  group_by(haplotype) %>%
  tally(name='nSom')

cabh$haplotype <- factor(cabh$haplotype, levels=c(1, 2, 0))

# is there really a surplus of haplotype1 calls?

  # subset of ied, keeping only reads from samples pdgoptag
iedoptag <- ied %>%
  subset(sample %in% pdgoptag)

  # build contingency table
cont <- table(iedoptag$isSomatic, iedoptag$haplotype)

  # Chi2 test
chi <- chisq.test(cont, correct=FALSE)

chi$expected
chi$observed
# the expected count for haplotype0 is low so throws a Warning
# in this case it is usually recommended to use Fisher's exact test
fisher.test(cont)
# it gives very similar p-value (slightly lower)

# I do not actually care about the haplotype0 calls, so perhaps better solution to exclude haplotype0 calls and do Chi2 again
iedoptag2 <- iedoptag %>%
  subset(haplotype!=0)
iedoptag2$haplotype <- factor(iedoptag2$haplotype, levels=c(1,2))
cont <- table(iedoptag2$isSomatic, iedoptag2$haplotype)
chi <- chisq.test(cont, correct=FALSE) # p-value ~ 10x lower
chi
chi$expected
chi$observed

# all 3 approaches essentially give the same answer
# I think will use last (test only on haplotype1/haplotype2)

# get expected counts to add to plot
hap1exp <- chi$expected[ which(row.names(chi$expected)=='TRUE'), which(colnames(chi$expected)=='1') ]
hap2exp <- chi$expected[ which(row.names(chi$expected)=='TRUE'), which(colnames(chi$expected)=='2') ]

callsbyHaplotype <- ggplot(cabh, aes(x=haplotype, y=nSom, fill=haplotype)) +
  geom_col() +
  #geom_segment(aes(x=0.5, xend=1.5, y=hap1exp, yend=hap1exp), linetype=2, colour='#ebebeb', size=1) +
  #geom_segment(aes(x=1.5, xend=2.5, y=hap2exp, yend=hap2exp),linetype=2, colour='#ebebeb', size=1) +
  scale_fill_manual(values=hapcols) +
  scale_x_discrete(labels=c('OPR mutated', 'OPR reference', 'unassigned')) +
  theme_minimal() +
  theme(
    legend.position='none',
    #panel.grid.major.x=element_blank(),
    #panel.grid.minor.x=element_blank(),
    axis.title.x=element_blank(),
    # axis.text.x=element_text(angle=90, size=7, margin = margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  ylab('somatic mutation calls')
callsbyHaplotype

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsbyHaplotype.pdf'), callsbyHaplotype, width=40, height=60, units='mm')
}

# with percentages
  # (should look fairly similar as haplotype1/2 are roughly 50/50)
  # cannot use haplotype0 calls, very low count but of low total number of reads, so give higher percentage
  # definitely misleading, will simply exclude for the purpose of this plot

# add total count of each haplotype
iedoptag2$haplotype <- factor(iedoptag2$haplotype, levels=c(1, 2))

cabh2 <- subset(cabh, haplotype!=0)

cabh2 <- iedoptag2 %>%
  group_by(haplotype) %>%
  tally(name='hapReads') %>%
  left_join(cabh2, by='haplotype')

# add proportion of each haplotype of whole dataset
cabh2$hapPro <- cabh2$hapReads / nrow(iedoptag)

# compute somatic % of haplotype reads
cabh2$somPro <- (cabh2$nSom/cabh2$hapReads)
cabh2$somPer <- (cabh2$nSom/cabh2$hapReads) * 100

# compute expected counts & expected percentages
# expected count = proportion of haplotype * total number of somatic calls
cabh2 <- cabh2 %>%
  mutate(expCount=hapPro * sum(cabh2$nSom)) %>%
  mutate(expPer=(expCount/hapReads)*100)

cabh2$haplotype <- factor(cabh2$haplotype, levels=c(1, 2))

# plot as above
callsbyHaplotypePer <- ggplot(cabh2, aes(x=haplotype, y=somPer, fill=haplotype)) +
  geom_col() +
  geom_hline(yintercept=cabh2$expPer[1], linetype=2, colour='#ebebeb', size=2) +
  scale_fill_manual(values=hapcols) +
  scale_x_discrete(labels=c('OPR mutated', 'OPR reference', 'unassigned')) +
  theme_minimal() +
  theme(
    legend.position='none'
  ) +
  ylab('somatic mutation calls (% of reads)')
callsbyHaplotypePer

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsbyHaplotypePercentage.pdf'), callsbyHaplotypePer, width=100, height=100, units='mm')
}



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

callsvsnewOPR <- ggplot(sowp, aes(x=somOPR, y=nSom)) +
  geom_col(fill=NA, colour='#202222') +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90)
  ) +
  xlab('') + ylab('somatic mutation calls')
callsvsnewOPR

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsnewOPR.pdf'), callsvsnewOPR, width=150, height=100, units='mm')
}

# I do not think I can do in percentages here
# I do not see percentage of what I would do...
# actually I think it is not the most useful plot, as we know it is highly dependent on what the sample was originally

# can fill with haplotype as stacked barplot
sowph <- iedSom %>%
  group_by(somOPR, haplotype) %>%
  tally(name='nSom')
sowph$haplotype <- factor(sowph$haplotype, levels=c(1, 2, 0))

callsvsnewOPRHap <- ggplot(sowph, aes(x=somOPR, y=nSom, fill=haplotype)) +
  geom_col(position='stack') +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90),
    axis.title.x=element_blank(),
    legend.position='none'
  ) +
  scale_fill_manual(values=hapcols) +
  ylab('somatic mutation calls')
callsvsnewOPRHap

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvsnewOPRHap.pdf'), callsvsnewOPRHap, width=180, height=100, units='mm')
}




# Number of calls vs somatic OPR change -----------------------------------
# more interesting than above as includes 'what the sample was'
# eg. OPRI8 sample, somatic read is haplotype1 and says OPRI6 >>> -2
# essentially 'how far does it go from the genotype'

# first prepare column somOPRdiff
# this is somatic OPR difference vs reference (whatever the sample)
# eg. somatic calls is 8 OPRI >>> +8
iedSom$somOPRdiff <- as.integer(unlist(sapply(iedSom$somOPR, function(o) {
  if(substr(o, 1, 4) == 'OPRI') { # if somatic call is an insertion
    return(as.integer(substr(o, 5, 99))) # return the number of OPR
  }
  
  else if(substr(o, 1, 4) == 'OPRD') { # if somatic call is a deletion
    return( - as.integer(substr(o, 5, 99))) # return the number of OPR with a minus sign
  }
  
  else if (o == 'reference') { # if somatic call is reference (only possible for compound homozygous sample)
    return(0)
  }
  
})))

# second prepare exphapOPR column
# = Sanger OPR diff we would expect for this haplotype
  # eg. sample is 8 OPRI
  # haplotype1 read >>> expected +8
  # haplotype2 read >>> expected 0
tmp <- unlist(sapply(1:nrow(iedSom), function(r) { # for each somatic read
  
  # get that read's sample
  spl <- iedSom[r, 'sample']
  # get the SV of that read's sample
  sv <- iedSom[r, 'SV']
  # get that read's haplotype
  hap <- iedSom[r, 'haplotype']
  
  # if SV is NA, always return 0 (reference)
  # Note, it is also ok to return 0 for samples that we could not haplotag but do not have SV
  if (is.na(sv)) {
    return(0)
  } 
  
  # if haplotype is 0 and there is an SV, return NA (as we do not know what to expect)
  if (hap==0 & !is.na(sv)) {
    return(NA)
  }
  
  else if (spl=='pdg46345' & hap==1) { # if it is the compound homozygous sample and we are looking at haplotype1 (OPRI5)
    return(as.integer(substr(sv, 5, 99))) # return the number of OPR (+5)
  }
  
  else if (spl=='pdg46345' & hap==2) { # if it is the compound homozygous sample and we are looking at haplotype2 (OPRD1)
    return(-1) # return -1
  }
  
  else if (substr(sv, 1, 4) == 'OPRI' & hap==1) { # if SV is an insertion and we are looking at haplotype1 (mutated)
    return(as.integer(substr(sv, 5, 99))) # return the number of OPR
  }
  
  else if (substr(sv, 1, 4) == 'OPRI' & hap==2) { # if SV is an insertion but we are looking at haplotype2 (reference)
    return(0)
  }
  
  else if(substr(sv, 1, 4) == 'OPRD' & hap==1) { # if SV is a deletion and we are looking at haplotype1 (mutated)
    return( - as.integer(substr(sv, 5, 99))) # return minus the number of OPR
  }
  
  else if(substr(sv, 1, 4) == 'OPRD' & hap==2) { # if SV is a deletion but we are looking at haplotype2 (reference)
    return(0)
  }
  
  else{
    cat('\t \t \t \t >>> Something unexpected row', r, '\n')
  }
  
}))

length(tmp) == nrow(iedSom) # ok can add to iedSom
iedSom$exphapOPR <- tmp

# most NA should be from samples which have a mutated OPR but were not haplotagged (PDG 53747, 59060, 58648)
# as we do not know what to 'expect' from that read (was it from mutated OPR haplotype or not?)

# others NA are from unassigned reads (haplotype0) of samples which have a mutated OPR

# Note; samples which do not have a mutated OPR, we know both haplotypes are reference so what we expect the reads to be is reference = 0 OPRdiff

# third prepare column somOPRchange
  # = 'how far did the somatic call go from what is expected from the haplotype'
  # eg. 3 OPRI somatic call in a 4 OPRI sample >> -1
  # first get absolute difference between somatic OPR diff / Sanger OPR diff
  # then add the appropriate sign
tmp <- sapply(1:nrow(iedSom), function(re) { # for each read
  san <- iedSom[re, 'exphapOPR'] # Sanger OPR count for that haplotype
  som <- iedSom[re, 'somOPRdiff'] # somatic OPR count for that read
  
  # if exphapOPR is NA, return NA
  # (somOPRdiff is never NA)
  if (is.na(san)) {
    return(NA)
  }
  
  if (san < som) {
    # eg. Sanger = -2 / somatic = -1 >> should return +1
    # eg. Sanger = 0 / somatic = 4 >> should return +4
    return(abs(san - som))
  }
  
  else if (san > som) {
    # eg. Sanger = -1 / somatic = -2 >> should return -1
    # eg. Sanger = 4 / somatic = 0 >> should return -4
    return( - abs(san - som))
  }
})

iedSom$somOPRchange <- tmp # add the column



# prepare counts for plot
scoc <- iedSom %>% # somatic counts by somatic OPR change (eg. how many -1)
  group_by(somOPRchange) %>%
  tally(name='nSom')

# plot somatic OPR change vs number of calls
callsvssomaticOPRchange <- ggplot(scoc, aes(x=somOPRchange, y=nSom)) +
  geom_col(fill=NA, colour='#4D4D4D', width=0.8) +
  theme_minimal() +
  theme(
    panel.grid.minor.x=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t=2, r=0, b=0, l=0)),
    axis.text.x=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  scale_x_continuous(breaks=min(scoc$somOPRchange, na.rm=TRUE):max(scoc$somOPRchange, na.rm=TRUE)) +
  xlab('OPR change vs haplotype') + ylab('somatic mutation calls')
callsvssomaticOPRchange

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvssomaticOPRchange.pdf'), callsvssomaticOPRchange, width=74, height=70, units='mm')
}
  
# prepare counts, similar as above but adding haplotype resolution
scoch <- iedSom %>% # somatic counts by somatic OPR change (eg. how many -1) by haplotype
  group_by(somOPRchange, haplotype) %>%
  tally(name='nSom')

# make sure haplotype is in order I want
scoch$haplotype <- factor(scoch$haplotype, levels=c(1, 2, 0))

# would adding haplotype fill informative?
callsvssomaticOPRchangeHap <- ggplot(scoch, aes(x=somOPRchange, y=nSom , fill=haplotype)) +
  geom_col() +
  scale_fill_manual(values=hapcols) +
  theme_minimal() +
  theme(
    legend.position='none'
  ) +
  scale_x_continuous(breaks=min(scoc$somOPRchange, na.rm=TRUE):max(scoc$somOPRchange, na.rm=TRUE)) +
  xlab('OPR change vs haplotype') + ylab('somatic mutation calls')
callsvssomaticOPRchangeHap
# there are a few I am suspicious about
  # eg. -5 calls are haplotype1, reference call from PDG 46345 (5 OPRI/1 OPRD)
  # it may be more likely to be a +1 from 1 OPRD haplotype
  # but I do not see any objective reason to exclude those


# plot separately inherited (mostly mutated OPR) and sporadic/control (reference OPR)
scocc <- iedSom %>% # somatic counts by somatic OPR change (eg. how many -1) by cohort
  group_by(somOPRchange, prion_disease) %>%
  tally(name='nSom')

# add column inherited or not
scocc$inhTF <- scocc$prion_disease == 'inherited' # inherited TRUE or FALSE

# make sure inhTF is in order I want
scocc$inhTF <- factor(scocc$inhTF, levels=c(TRUE, FALSE))

# need to 'fill in' the dataframe so both inherited TRUE & FALSE have all OPRs
# otherwise, eg. -5 is all inherited, so makes a large bar which is inconsistent with the others
# or in other words, add the 0 counts (-5 non-inherited >> should be 0 nSom)

# make the template dataframe and left_join it to scocc
oprcgs <- min(scocc$somOPRchange, na.rm=TRUE) : max(scocc$somOPRchange, na.rm=TRUE) # OPR changes present in the data
tmpl <- as.data.frame(matrix(nrow=length(oprcgs)*2, ncol=2)) # nrow = all OPR changes for both inherited TRUE and FALSE, so x 2
colnames(tmpl) <- c('inhTF', 'somOPRchange')
tmpl$inhTF <- c(rep(TRUE, length(oprcgs)) , rep(FALSE, length(oprcgs)))
tmpl$inhTF <- factor(tmpl$inhTF, levels=c(TRUE, FALSE))
tmpl$somOPRchange <- c(oprcgs, oprcgs)

scocc <- left_join(tmpl, scocc, by=c('inhTF', 'somOPRchange'))

# replace nSom NA we added by 0
scocc[which(is.na(scocc$nSom)) , 'nSom'] <- 0

# add total number of reads from inherited / other
scocc <- ied %>%
  mutate(inhTF = factor((ied$prion_disease == 'inherited'))) %>%
  group_by(inhTF) %>%
  tally(name='totalInhTF') %>%
  left_join(scocc, by='inhTF')
# above changes the factor levels
scocc$inhTF <- factor(scocc$inhTF, levels=c(TRUE, FALSE))

# compute percentage of cohort (inherited / other)
scocc$somPer <- (scocc$nSom / scocc$totalInhTF) * 100

# colours for inherited TRUE / FALSE
inhTFcols <- c('#697a87', '#595E60')

# counts
callsvssomaticOPRchangeCoh <- ggplot(scocc, aes(x=somOPRchange, y=nSom , fill=inhTF)) +
  geom_col(position='dodge') +
  scale_fill_manual(values=inhTFcols) +
  theme_minimal() +
  theme(
    legend.position='none'
  ) +
  scale_x_continuous(breaks=min(scoc$somOPRchange, na.rm=TRUE):max(scoc$somOPRchange, na.rm=TRUE)) +
  xlab('OPR change vs haplotype') + ylab('somatic mutation calls')
callsvssomaticOPRchangeCoh

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvssomaticOPRchangeCoh.pdf'), callsvssomaticOPRchangeCoh, width=120, height=100, units='mm')
}
  
# percentages
callsvssomaticOPRchangeCohPer <- ggplot(scocc, aes(x=somOPRchange, y=somPer , fill=inhTF)) +
  geom_col(position='dodge') +
  scale_fill_manual(values=inhTFcols) +
  theme_minimal() +
  theme(
    legend.position='none'
  ) +
  scale_x_continuous(breaks=min(scoc$somOPRchange, na.rm=TRUE):max(scoc$somOPRchange, na.rm=TRUE)) +
  xlab('OPR change vs haplotype') + ylab('somatic mutation calls (% of reads)')
callsvssomaticOPRchangeCohPer

if (exportOrNo) {
  ggsave(filename=here('needleSam', 'callsvssomaticOPRchangeCohPercentage.pdf'), callsvssomaticOPRchangeCohPer, width=120, height=100, units='mm')
}




# number of somatic insertion/deletion by tissue --------------------------
# in https://pubmed.ncbi.nlm.nih.gov/16713248/
# there is an argument why expansions (insertions) may be more frequent in non-dividing cells

# add a column to iedSom whether OPR change vs expected from haplotype is a deletion or an insertion
  # i.e. whether iedSom$somOPRchange is positive or negative
iedSom$somOPRchange_type <- NA
iedSom$somOPRchange_type[which(iedSom$somOPRchange > 0)] <- 'insertion'
iedSom$somOPRchange_type[which(iedSom$somOPRchange < 0)] <- 'deletion'

# compare within tissue comparison cohort
idtc <- iedSom %>% # insertion/deletion tissue comparison
  subset(cohort=='inherited_tissuecomparison') %>%
  group_by(sample, individual, gDNA_tissue, somOPRchange_type) %>%
  tally(name='nSom') %>%
  pivot_wider(names_from=somOPRchange_type,
              values_from=nSom)
# one pair of samples could not be haplotagged, so cannot know if insertion or deletion
# can delete these samples
inds2del <- match(pdgs_notag, idtc$sample)
inds2del <- inds2del[!is.na(inds2del)]
idtc <- idtc[-inds2del,]

# >>> can see it will not be interesting. There are only 2 insertions, both in blood tissues


# Can chimeric reads/index hopping be an issue? ---------------------------
# stringent with filtering criteria in Guppy and before needleSam.R so I do not think likely
# but here let's take a few examples to see if likely or not

# somatic calls by new OPR by flowcell
oprflo <- iedSom %>%
  group_by(somOPR, flowcell) %>%
  tally(name='nSom')

# count total number of reads from each flowcell
# and the counts to oprflo
oprflo <- ied %>%
  group_by(flowcell) %>%
  tally(name='totreads') %>%
  left_join(x=oprflo, y=., by='flowcell')

# compute somatic calls in % of flowcell
oprflo$somPer <- (oprflo$nSom / oprflo$totreads) * 100

# convert flowcell to factor
oprflo$flowcell <- factor(oprflo$flowcell, levels=c(1, 2, 3))

# Example1/ flowcell3 has two samples OPRI8; others have 0
# >> are there more OPRI8 somatic calls from flowcell3 than flowcell 1 & 2?
subset(oprflo, somOPR=='OPRI8') # no 8 OPRI call from flowcell3

# Example2/ flowcell3 has four samples OPRI4; other have 0
nrow(subset(ied, flowcell==3 & opri=='OPRI4')) / nrow(subset(ied, flowcell==3)) # 20+ % of reads on flowcell3 are OPRI4
# >> are there more OPRI4 somatic calls from flowcell3 than flowcell 1 & 2?
subset(oprflo, somOPR=='OPRI4')
# no 4 OPRI call from flowcell3

# and how many somatic mutation from flowcell3 (only samples without a 4 OPRI)
nrow(subset(iedSom, sample %in% as.character(unlist(subset(meta, flowcell==3 & SV != 'OPRI4', sample)))))
# >>> I do not believe that chimeric read/index hopping can explain any somatic call





# write the somatic calls -------------------------------------------------

# export iedSom
write.xlsx(iedSom, file=here('needleSam', 'somaticCalls.xlsx'))

# export cabs
write.xlsx(cabs, file=here('needleSam', 'callsBySample.xlsx'))






# prepare a shortlist -----------------------------------------------------

# some of the criteria below cannot be applied 'equitably' to all samples
# only interested here in preparing a shortlist the best somatic mutation calls

# will do by read IDs, safer than row indices

# shortlist somatic calls
cat('\t \t \t \t >>>', nrow(iedSom), 'calls before shortlisting \n')

callsns <- c() # calls not shortlisted (read IDs)

# 1-- delete any haplotype0 calls from samples we could haplotag
# (the read is probably messy if Whatshap could not haplotag it)
tmp <- iedSom[which(iedSom$sample %in% pdgs_tag & iedSom$haplotype==0), 'read']
cat('\t \t \t \t >>>', lengths(tmp), 'calls are haplotype0 calls even though from a haplotagged sample \n')
callsns <- c(callsns, tmp)


# 2-- decrease allowed mismatches

  # original threshold: mean + sd of mismatch by bp of true reads
    # makes us keep the ~ 86% most accurate reads
length(which(ied$true_misbp < thrMis)) / length(which(!is.na(ied$true_misbp)))

# we saw somatic calls had a surplus of reads in the 5% mismatch range even after thresholding
  # >> more stringent threshold: just mean
    # makes us keep ~ 55% most accurate reads
thrMis2 <- mean(ied$true_misbp, na.rm=TRUE)
length(which(ied$true_misbp < thrMis2)) / length(which(!is.na(ied$true_misbp)))

tmp <- iedSom[which(iedSom$som_misbp > thrMis2) , 'read']

cat('\t \t \t \t >>>', length(tmp), 'calls are above stringent mismatch threshold \n')
callsns <- c(callsns, tmp)


# 3-- exclude any calls with low haplotype assignment score
tmp <- iedSom[which(iedSom$hapscore < 1.0), 'read'] # setting threshold at 0.5, basically 50% of the expected SNVs are not present
cat('\t \t \t \t >>>', length(tmp), 'calls have a low haplotype assignment score \n')
callsns <- c(callsns, tmp)


# remove these calls
callsns <- unique(callsns)
cat('\t \t \t \t >>> Excluding', length(callsns), 'calls \n')

iedSomsl <- iedSom[ - match(callsns, iedSom$read) ,] # sl for shortlist

cat('\t \t \t \t >>> Shortlisted to', nrow(iedSomsl), 'calls \n')

# in shortlist; any replicates?
  # i.e. more than 1 read same sample/haplotype/somOPR
iedSomsl %>%
  group_by(sample, haplotype, exphapOPR, somOPR, somOPRchange) %>%
  tally(name='nSom') %>%
  subset(nSom>1)




# reads for figure --------------------------------------------------------

# alignments in figure
# from top to bottom
figreads <- c('8f9756e2-1fb9-41d7-a6bf-a11463a55990',
              'a5b5b295-3e6d-4f48-8839-24aa3fd50718',
              'e1c59ac3-a547-4616-a120-4dd7d3463e24',
              '2f8a89cb-73d7-4591-af8b-f814ac373a24',
              'e67c028a-35e2-44b0-80b8-45ccaa100392',
              '9d387703-c39e-47df-bd96-49b0a1d62b4d',
              '86f5c651-7ce0-4778-a55f-bc01cacec79b',
              'e43b902e-cbaf-4feb-9170-51caecb1cf81',
              'aa058a24-7e53-4d55-ad91-f5dcdc54602a',
              '5a5e7247-0550-4f16-aaf2-e51ad2b5c9a5')

length(figreads)

figReds <- ied[match(figreads, ied$read),]

if (exportOrNo) {
  write.xlsx(figReds, file=here('needleSam', 'figureSomaticReads.xlsx'))
}



# write shortlist ---------------------------------------------------------

# write the shortlist

if (exportOrNo) {
  write.xlsx(iedSomsl, file=here('needleSam', 'somaticCallsShortlist.xlsx'))
}
