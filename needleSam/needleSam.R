# needle in a SAM

# finding a possible somatic mutation in the reads
# like a needle in a haystack!
# needle = somatic mutation / haystack = SAM files

oprlength <- 123
# threshold used in readSam() to exclude reads that start or terminate in the OPR by CIGAR lengths



# some common things we may need to change often --------------------------
hapcol <- 21 # nth field in the SAM file contains the HP tags (looks like HP:i:1; it was 20th then 21th in a newer version of minimap2)
hapcols <- c('#f1876b', '#4D4D4D', '#B3BDC4') # colours for the haplotypes, it goes 1, 2, 0



# packages and functions --------------------------------------------------


library(here)
library(stringr)
library(data.table)
library(openxlsx)
library(msa)
library(ggplot2)
library(tictoc)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggbeeswarm)

# substrEnding = last n characters of a string
# x = string (or vector of strings); n = last n characters
substrEnding <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# beforeFirstCharacter = given a string, takes everything before first instance of a character
# eg. beforeFirstCharacter('59060_file.bam', '_') >> '59060'
beforeFirstCharacter <- function(stri, cha) {
  as.character(sapply(stri, function(st){
    substr(st,
           1,
           which(strsplit(st, '')[[1]] == cha) [length(which(strsplit(st, '')[[1]] == cha))]-1)
  }))
}

###

# getReadlengthCigar() is a small helper function for readSam() below
  # the goal is to exlude reads which start or terminate in the OPR
  # these get categorised as reference because their CIGAR are for example 56M
# However, these reads do not span the whole OPR, so they should be excluded
  # but excluding by readlength does not work, for example a 103-bp read could be a OPRD1 or a reference read that terminated early in the OPR
# one solution is to compute the alignment length as told by the CIGAR (only Match/Insertion/Deletion)
  # this essentially corresponds to the window of hg38 that is covered by the read (as it will include deletions)

# get read length as told by CIGAR
getReadlengthCigar <- function(read) { # works for one row of a sam file
  cigar <- read[6]
  tmp <- str_extract_all(cigar, '[0-9]+[MDI]')[[1]] # all the Match/Deletion/Insertion components of the CIGAR
  rlci <- sum(as.numeric(unlist(strsplit(tmp, '[A-Z]')))) # rlci = readlength cigar
  return(rlci)
}

###

# readSam = reads SAM properly into R
readSam <- function(sampath) {
  
  # check the path actually looks like it leads to a SAM file
  if(substrEnding(sampath, 4) != '.sam') stop('\t \t \t \t >>> Error: that does not look like the path to a SAM file... \n')
  
  # if ok -- import the SAM file
  # numcols <- max(count.fields(sampath, sep='\t'), na.rm=TRUE) 
    # above is to count maximum number of columns to be used below as sprintf('col%i', numcols)
    # works well, but when importing more than one SAM, can have variable number of columns
    # alternative: set to a large number. It will create NA columns at the end, but at least we know all the SAMs look the same
    # (can delete the NA columns after rbind all the SAMs, if that is what we want to do)
  sam <- read.table(file=sampath, fill=TRUE, comment.char='@', sep='\t', na.strings=c('', 'NA'), 
                    col.names=sprintf('col%i', 1:30))
  
  # check what should be the haplotype tag column, can be messy for different reasons
  hpcol <- 21
  
  taghp <- substr(sam[,hpcol], 1, 2) # will take mostly NA or HP and we do not want anything else, hence below
  cat('\t \t \t \t \t >>> Found', length(which(taghp=='HP')), 'haplotags in the expected column \n')
  
  if (length(which(taghp!='HP')) > 0) {
    cat('\t \t \t \t \t >>> Removed', length(which(taghp!='HP')), 'unexpected tags from the HP column \n')
    sam[which(taghp!='HP'), hpcol] <- NA # the rows that are not NA or HP, replace them with NA
    # ! for these rows, it is possible that the HP tag was pushed to the next column
      # no need to look only at these rows, just look in the next column for any HP
    taghpnxt <- substr(sam[,hpcol+1], 1, 2)
    
    if(length(which(taghpnxt=='HP')) > 0) { # if found any HP tag to salvage from the next column
      cat('\t \t \t \t \t \t >>> for', length(which(taghpnxt=='HP')), 'of these reads, salvaged the HP tag found in the next column \n')
      sam[which(taghpnxt=='HP'), hpcol+1] <- sam[which(taghpnxt=='HP'), hpcol+1] # copy it to the HP column
      
    # ! Note this is far from an ideal solution, will copy the HP tag and it also means the columns are shifted for these rows
      # but at the moment I am not using any other tag than HP
    } 
  }
  
  # check there is no HP tag in other columns before or after HP column
  for (co in ((hpcol-4):(hpcol+4))[-5]) { # check all 4 columns before HP column, all 4 columns after HP column, except HP column itself
    tagco <- substr(sam[,co], 1, 2)
    if(length(which(tagco=='HP')) > 0) {
      cat('\t \t \t \t \t >>> found ', length(which(tagco=='HP')), ' HP tags in column', co ,'\n')
    }
  }
  
  # check there is no duplicate read ID
    # as I removed secondary alignments in filterBam.command it should be good
  if (sum(duplicated(sam[1])) != 0) {
    stop('\t \t \t \t \t >>> Error: found', sum(duplicated(sam[1])), 'duplicated read IDs \n')
  } else {
    cat('\t \t \t \t \t >>> All read IDs are unique \n')
  }
  
  # exclude reads that start or terminate in the OPR (see explanations for getReadlengthCigar() above)
  reads2del <- which(apply(sam, 1, getReadlengthCigar) < oprlength)
  
  if (length(reads2del) == 0) {
    cat('\t \t \t \t \t >>> 0 read that start or terminate in the OPR \n')
  } else {
    cat('\t \t \t \t \t >>> Deleting', length(reads2del),'reads that start or terminate in the OPR \n')
    sam <- sam[ - reads2del ,]
  }
  
  # exclude reads that include letters that are not ACGT
    # I think removing the reads that start/terminate reads above include this, but just to be safe
  totnotACGT <- sapply(1:nrow(sam), function(ri) { # for each read, count the total of letters that are not ACGT
    sum( ! strsplit(sam[ri,10], '')[[1]] %in% c('A', 'C', 'G', 'T'))
  })
  
  reads2del <- which(totnotACGT!=0)
  
  if (length(reads2del) != 0) {
    cat('\t \t \t \t \t >>> Deleting', length(reads2del),'reads with non-ACGT nucleotide \n')
    sam <- sam[ - reads2del ,]
  }
  
  # say final number of reads we are importing
  cat('\t \t \t \t \t >>> IMPORTING', nrow(sam),'READS \n')
  
  return(sam)
  
}


###

# from CIGAR of one read; get total insertion
getTotalInsertion <- function(read) {
  
  cigar <- read[6]
  
  if(length(grep(pattern='I', cigar))!=0) {
    tmp <- str_extract_all(cigar, '[0-9]+I')[[1]] # extracts the I components; eg. 59I129M32I >> 59I 32I
    totin <- sum(as.numeric(unlist(strsplit(tmp, '[A-Z]'))))
    # split at any letter (here will be only I), so example above: 59 32
    # convert to numeric & sum
  } else {
    totin <- 0
    return(totin)
  }
}

###

# from CIGAR of one read; get total deletion
getTotalDeletion <- function(read) {
  
  cigar <- read[6]
  
  if(length(grep(pattern='D', cigar))!=0) {
    tmp <- str_extract_all(cigar, '[0-9]+D')[[1]] # extracts the I components; eg. 59I129M32I >> 59I 32I
    totdel <- -sum(as.numeric(unlist(strsplit(tmp, '[A-Z]'))))
    # split at any letter (here will be only I), so example above: 59 32
    # convert to numeric & sum & turn to negative for deletion
    return(totdel)
  } else {
    totdel <- 0
    return(totdel)
  }
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

freqInBins <- function(df, datcol, grpcol, bins) {
  
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




# 1- Import what we need --------------------------------------------------


  # 1- Catalog of consensus OPR sequences
    # OPRConsensusCatalog.csv written by generateOPRCatalog.R
cato <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='OPRtemplates') # consensus sequence from OPRD4 up to OPRI24


  # 2- Meta information of samples
meta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='samples')
    # make the formatting match

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


  # 3- Import haplotype meta information
    # eg. how many SNVs available for haplotagging
hpmeta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='haplotypephasing')
    # sample: 52331 >> pdg52331, etc.
hpmeta$sample <- sprintf('pdg%s', hpmeta$sample)

  # and sequencing_summary (mainly for coverage)
seqmeta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='sequencing_summary')
# sample: 52331 >> pdg52331, etc.
seqmeta$sample <- sprintf('pdg%s', seqmeta$sample)


  # 4- Import the SAMs
    # preallocate list, each element = the SAM of one sample
    # 2 folders where to find the SAMs:
      # one for SAMs that could not be haplotagged
      # the other for SAMs that could be haplotagged
samfoldnotag <- here('haplotypephasing/bam_notag/OPRtrim_sam') # folder of no-tag SAMs
samfoldtag <- here('haplotypephasing/bam_haplotagged/OPRtrim_sam') # folder of tagged SAMs

samnms <- c(list.files(samfoldnotag), list.files(samfoldtag)) # names of all the SAMs
    # check expected number
length(samnms) == nrow(meta)

sams <- vector(mode='list', length=length(samnms)) # preallocate list that will store all the SAMs
pdgs <- sprintf('pdg%s', beforeFirstCharacter(samnms, '_')) # gives pdg58398 etc.
names(sams) <- pdgs # put the pdg58398 etc. as names of the elements in the list

    # get also pdgs of without-haplotag / with-haplotag, will be useful later
pdgs_notag <- sprintf('pdg%s', beforeFirstCharacter(list.files(samfoldnotag), '_'))
pdgs_tag <- sprintf('pdg%s', beforeFirstCharacter(list.files(samfoldtag), '_'))

    # gets the full paths where to find each SAM
sampaths <- c(list.files(samfoldnotag, full.names=TRUE),
              list.files(samfoldtag, full.names=TRUE))

    # fill in the list using function readSam from above
for (s in 1:length(sams)) {
  cat('\t \t \t \t >>>', s, '/', length(sams), 'Reading SAM of sample', names(sams)[s], '\n')
  sams[[s]] <- readSam(sampaths[s])
}




# 2- Calculate total insertion/deletion of each read ----------------------
# from the CIGAR; eg. 3I 2D 100M 3I 5D >> insertion = 6 bp // deletion = 5 bp

  # readSam() should have removed anything funky in HP column, but check before running below just to be sure
hpcol <- 21
sum(unlist(lapply(sams, function(sa) {
  which( substr(as.character(sa[,hpcol]), 1, 2) != 'HP') # any elements in all the HP columns that is not NA or does not start with HP
}))) == 0

  # preallocate list of indel dataframes
  # same length as sams
ies <- vector(mode='list', length=length(sams)) # ies for Insertion/dEletion dataframeS
names(ies) <- pdgs # put the pdg58398 etc. as names of the elements in the list

  # each element in the list will be a dataframe;
    # column1 = sample
    # column2 = haplotype (from whatshap haplotag)
    # column2 = read name
    # column3 = read sequence
    # column3 = total insertion (in bp, as per CIGAR)
    # column4 = total deletion (in bp, as per CIGAR)

for (i in 1:length(ies)) { # for each element in the list (each sample)
  
  cat('\t \t \t \t >>> Format & parsing CIGARs for SAM', i, '/', length(ies), '\n')
  
  # preallocate the dataframe
  colnms <- c('sample', 'haplotype', 'hapqual', 'read', 'samflag', 'seq', 'ins', 'del')
  tmp <- as.data.frame(matrix(nrow=nrow(sams[[i]]), ncol=length(colnms)))
  colnames(tmp) <- colnms
  
  # fill in sample name
  tmp$sample <- rep(names(ies)[i], nrow(sams[[i]]))
  
  # fill in haplotype
  tmp$haplotype <- as.numeric(substr(sams[[i]][,hpcol], 6, 6)) # simplifies the haplotype tag, from HP:i:1 / HP:i:2 to 1 / 2; will return NA if no present
  
  # fill in haplotype quality
  tmp$hapqual <- as.numeric(substr(sams[[i]][,hpcol+1], 6, 999)) # simplifies the haplotype quality, from eg. PC:i:150 to 150
  
  # fill in read names
  tmp$read <- sams[[i]][,1]
  
  # fill in sam flags
  tmp$samflag <- sams[[i]][,2]
  
  # fill in read sequences
  tmp$seq <- sams[[i]][,10]
  
  # fill in total insertions
  tmp$ins <- apply(sams[[i]], 1, getTotalInsertion)
  
  # fill in total deletions
  tmp$del <- apply(sams[[i]], 1, getTotalDeletion)
  
  # place that sample's dataframe in the list
  ies[[i]] <- tmp
}

# turn the ies list into one big dataframe
ied <- as.data.frame(rbindlist(ies)) # for Insertion/dEletion Dataframe

# make sure haplotype is stored as factor
ied$haplotype <- factor(ied$haplotype, levels=c(1, 2, 0))

# check no duplicated read ID after pooling all samples
sum(duplicated(ied$read))




# 4- Check haplotype tags -------------------------------------------------

# check a few things about haplotags
# sample that were not phased should haplotype all NA
subset(ied,
       sample %in% pdgs_notag &
         !is.na(haplotype))

# haplotypes can only be NA, 1, 2
subset(ied,
       (!is.na(haplotype) & ! haplotype %in% c(1,2))) # any haplotype not NA and not 1 or 2

# I think will be easier if we replace haplotype NA (i.e. unassigned read) by 0
ied$haplotype[which(is.na(ied$haplotype))] <- 0

# check no more NA in haplotype
sum(is.na(ied$haplotype))

# Are roughly 50% of the reads assigned to each haplotype?
  # for each sample, plot the number of reads / proportion of reads (barplot) from haplotype0 (unassigned), haplotype1, haplotype2

hco <- ied %>% # hco for haplotype counts
  group_by(sample, haplotype) %>%
  tally(name='nreads')

# add total number of reads in each sample
tmp <- ied %>%
  count(sample, name='totalsamplereads')

hco <- left_join(hco, tmp, by='sample')

# weirdly this seems to be the easiest way to order hco by pdgs
hco <- left_join(data.frame(sample=pdgs), hco, by='sample')

# normalise the read counts per haplotype
hco <- hco %>%
  mutate(haplotype_pro=nreads/totalsamplereads) %>% # haplotype proportion, i.e. proportion of reads assigned to haplotype 1, 2, or unassigned (0)
  mutate(sam_hap=paste(sample, haplotype, sep='_'), .after=haplotype) # add a unique sample_haplotype column, easier for plotting

# put levels of sample and sam_hap in right order for plotting
hco$sample <- factor(hco$sample, levels=unique(hco$sample))
hco$sam_hap <- factor(hco$sam_hap, levels=hco$sam_hap)
hco$haplotype <- factor(hco$haplotype, levels=c(1, 2, 0))

# haplotype assignments as stacked barplot, actual number of reads (so height = OPR coverage)
haploReads <- ggplot(hco, aes(x=sample, y=nreads, fill=haplotype)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=hapcols) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90),
    axis.title.x=element_blank()
  ) +
  ylab('number of reads')
haploReads

# haplotype assignments as stacked barplot, proportion of coverage (so height is 1 for all samples)

# I prefer samples without haplotag at the top

# add column tag or no tag
hco$htag <- TRUE
hco[which(hco$sample %in% pdgs_notag) , 'htag'] <- FALSE

hco <- hco %>%
  arrange(desc(htag))

hco$sample <- factor(hco$sample, levels=unique(hco$sample))

haploPro <- ggplot(hco, aes(x=sample, y=haplotype_pro, fill=haplotype)) +
  coord_flip() +
  geom_bar(stat='identity') +
  scale_fill_manual(values=hapcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    axis.title.x=element_blank(),
    axis.text.x=element_text(),
    axis.title.y=element_blank(),
    axis.text.y=element_text()
  ) +
  ylab('proportion of reads')
haploPro

ggsave(here('haplotypephasing', 'haploReads.pdf'), haploReads,width=200, height=100, units='mm')
ggsave(here('haplotypephasing', 'haploPro.pdf'), haploPro,width=200, height=100, units='mm')


# Why unassigned reads?
  # samples which did not have at least one heterozygous SNVs cannot be haplotagged, that is expected
  # but some samples have up to ~ 15% unassigned reads (i.e. haplotype0)
# worst case: PDG 55048 ~ 15% haplotype0
    # looking carefully at the reads in IGV, whatshap seems to do a good job, even considering possible small shifts in the alignment
    # there is only one usable SNV for PDG 55048, high unassigned % is mostly reads which do not cover the SNV
    # I think nothing unexpected
# does proportion of unassigned reads correlate strongly with how many SNVs available?

# add haplotype phasing meta info
  # keep all rows from hco (duplicate rows from hpmeta), so left_join
hco <- left_join(x=hco, y=hpmeta)

# add also the cohort of each sample
hco <- left_join(x=hco, y=select(meta, sample, prion_disease))


hcoun <- subset(hco, ! sample %in% pdgs_notag & haplotype==0) # keep only proportion of unassigned reads, from haplotagged samples

# ! any samples haplotagged but not there: their proportion of unassigned reads is 0, should add them to the plot
  # slightly clumsy, but will take one row for each (which will correspond to haplotype 1 or 2) and replace the few cells that need to be changed
rows2add <- hco[which(hco$sample == pdgs_tag[which(! pdgs_tag %in% hcoun$sample)]) , ]
rows2add$haplotype <- 0
rows2add$sam_hap <- paste(rows2add$sample, rows2add$haplotype, sep='_')
rows2add$nreads <- 0
rows2add$haplotype_pro <- 0
# rest is correct as dependent on the sample, not the haplotype
hcoun <- rbind(hcoun, rows2add)

hapUnassigned <- ggplot(hcoun, aes(x=genebody_heterozygousSNVs, y=haplotype_pro, colour=prion_disease)) +
  geom_quasirandom(width=0.3) +
  theme_minimal() +
  theme(
    panel.grid.minor.x=element_blank()
  ) +
  scale_x_continuous(breaks=1:16) +
  xlab('informative SNVs') + ylab('proportion unassigned reads')
hapUnassigned
ggsave(here('haplotypephasing/hapUnassigned_v4.pdf'), hapUnassigned,width=200, height=100, units='mm')
cor(hcoun$genebody_heterozygousSNVs, hcoun$haplotype_pro, method='pearson')
cor(hcoun$genebody_heterozygousSNVs, hcoun$haplotype_pro, method='spearman')




# 5- Plot insertions/deletions for each sample ----------------------------

# will then assign each read to a most likely OPR based on total insertion/deletion
# should leave some space for nanopore noise
  # eg. an OPRI often appears as a 23-bp insertion
# here: take the opportunity to look how much space we should leave

  # tentative thresholds
lownoise <- 6
highnoise <- 3
    # eg. for OPRI1: ok from 24 - lownoise up to 24 + highnoise
    # for example, if lownoise = 6 // highnoise = 3; considers potential OPRI1 anything from insertion of 18 bp to insertion of 27 bp
    # makes sense for lownoise to be larger, as there is always a drag towards lower lengths
    # (or in other words, nanopore makes more deletions errors than insertions)

# plot an insertion + a deletion histogram for each sample
# preallocate list to store the plots
# ! one list for all deletion plots and one list for all insertion plots
ips <- vector(mode='list', length=length(pdgs)) # Insertions PlotS
names(ips) <- pdgs

dps <- vector(mode='list', length=length(pdgs)) # Deletions PlotS
names(dps) <- pdgs

# make a table of target insertion/deletion lengths
  # eg. OPRD1 : deletion -24, etc.
# will do from OPRD4 up to OPRI24, so 4 OPRD + reference + 24 OPRI
tarl <- as.data.frame(matrix(nrow=29, ncol=2)) # target lengths
colnames(tarl) <- c('oprg', 'targetindel') # OPR genotype // target indel length
tarl$oprg <- c(sprintf('OPRD%i', 4:1) , 'reference' , sprintf('OPRI%i', 1:24))
tarl$targetindel <- c( seq(-24*4, -24*1, 24) , 0 , seq(24*1, 24*24, 24)) # (one OPR = 24 bp)

# first; DELETION plots
for (s in 1:length(dps)) {
  
  # PDG we are plotting
  pd <- names(dps)[s]
  
  # get its OPR variant (if any)
  sva <- as.character(subset(meta, sample==pd, SV))
  
  # and from it the target length
  target <- as.numeric(subset(tarl, oprg==sva, targetindel))
    # eg. PDG 1906 >> OPRD2 >> -48 bp
  
  # for that sample, count the reads with each deletion length
    # can also do it automatically with geom_histogram(), but I found it leave a bar when count = 0 if computed before, which I prefer
  tmp <- ied %>%
    subset(sample==pd) %>%
    group_by(del, haplotype) %>%
    tally(name='nreads') %>%
    mutate(freq=nreads/nrow(subset(ied, sample==pd))) 
      # above: compute frequency by dividing number of reads with each insertion length by total number of reads from that sample
  
  
  # if sample has a OPRD: version with marker lines
    # ! need to specify manually if PDG 46345 = compound homozygous 1 OPRD / 5 OPRI
  if (pd =='pdg46345') {
    
    target <- -24
    delp <- ggplot(tmp, aes(x=del, y=freq ,colour=haplotype)) +
      geom_col(fill='white', alpha=0.0, position='identity') +
      scale_color_manual(values=hapcols) +
      coord_cartesian(ylim=c(0,0.2), xlim=c(-250,0)) +
      geom_vline(xintercept=target) +
      geom_vline(xintercept=c(target-lownoise, target+lownoise), linetype=2) +
      theme_minimal() +
      theme(
        legend.position='none'
      ) +
      xlab('deletion') + ylab('frequency') +
      ggtitle(pd)

  } else if (!is.na(target) & target<0) { # if there is a target (i.e. sample has an OPRD) >> add the marker lines
    delp <- ggplot(tmp, aes(x=del, y=freq ,colour=haplotype)) +
      geom_col(fill='white', alpha=0.0, position='identity') +
      scale_color_manual(values=hapcols) +
      coord_cartesian(ylim=c(0,0.2), xlim=c(-250,0)) +
      geom_vline(xintercept=target) +
      geom_vline(xintercept=c(target-lownoise, target+lownoise), linetype=2) +
      theme_minimal() +
      theme(
        legend.position='none'
      ) +
      xlab('deletion') + ylab('frequency') +
      ggtitle(pd)
  } else { # if not (i.e. reference or OPRI), skip the marker lines
    delp <- ggplot(tmp, aes(x=del, y=freq ,colour=haplotype)) +
      geom_col(fill='white', alpha=0.0, position='identity') +
      scale_color_manual(values=hapcols) +
      coord_cartesian(ylim=c(0,0.2), xlim=c(-250,0)) +
      theme_minimal() +
      theme(
        legend.position='none'
      ) +
      xlab('deletion') + ylab('frequency') +
      ggtitle(pd)
  }
  
  # in any case, add deletion plot to the list
  dps[[s]] <- delp
  
}


# second; INSERTION plots
for (s in 1:length(ips)) {
  
  # PDG we are plotting
  pd <- names(ips)[s]
  
  # get its OPR variant (if any)
  sva <- as.character(subset(meta, sample==pd, SV))
  
  # and from it the target length
  target <- as.numeric(subset(tarl, oprg==sva, targetindel))
  # eg. PDG 1906 >> OPRD2 >> -48 bp
  
  # for that sample, count the reads with each insertion length
  # can also do it automatically with geom_histogram(), but I found it leave a bar when count = 0 if computed before, which I prefer
  tmp <- ied %>%
    subset(sample==pd) %>%
    group_by(ins, haplotype) %>%
    tally(name='nreads') %>%
    mutate(freq=nreads/nrow(subset(ied, sample==pd))) 
  # above: compute frequency by dividing number of reads with each insertion length by total number of reads from that sample
  
  # if sample has a OPRI: version with marker lines
  if (!is.na(target) & target>0) {
    insp <- ggplot(tmp, aes(x=ins, y=freq ,colour=haplotype)) +
      geom_col(fill='white', alpha=0.0, position='identity') +
      scale_color_manual(values=hapcols) +
      coord_cartesian(ylim=c(0,0.2), xlim=c(0,250)) +
      geom_vline(xintercept=target) +
      geom_vline(xintercept=c(target-lownoise, target+lownoise), linetype=2) +
      theme_minimal() +
      theme(
        legend.position='none'
      ) +
      xlab('insertion') + ylab('frequency') +
      ggtitle(pd)
  } else { # if not (i.e. reference or OPRI), skip the marker lines
    insp <- ggplot(tmp, aes(x=ins, y=freq ,colour=haplotype)) +
      geom_col(fill='white', alpha=0.0, position='identity') +
      scale_color_manual(values=hapcols) +
      coord_cartesian(ylim=c(0,0.2), xlim=c(0,250)) +
      theme_minimal() +
      theme(
        legend.position='none'
      ) +
      xlab('insertion') + ylab('frequency') +
      ggtitle(pd)
  }
  
  # in any case, add deletion plot to the list
  ips[[s]] <- insp
  
}

# arrange each list of plots in a grid of one column and export
  # grid of deletion plots

delgrid <- ggarrange(plotlist=dps, ncol=1, nrow=length(dps))
ggsave(filename=here('needleSam', 'deletionGrid.pdf'), plot=delgrid, height=1000, width=100, units='mm', useDingbats=FALSE)

insgrid <- ggarrange(plotlist=ips, ncol=1, nrow=length(ips))
ggsave(filename=here('needleSam', 'insertionGrid.pdf'), plot=insgrid, height=1000, width=100, units='mm', useDingbats=FALSE)

# Note there are 3 options for the histograms above:
  # position='stack' (default): each bar is a stacked barplot with proportion of each haplotype
    # in a sense it is the best way to plot all the information, but intuitively when I see the tip of the bars in orange (for example),
      # I imagine that there is another plot below it, which changes completely the interpretation
  # position='identity': histograms for each haplotypes are effectively separated, overlayed
  # position='dodge': same but the bars are a little bit shifted so there is less overlay
# I think 'identity' is the less likely to create confusion (even though it is the one with the most overlay)




# 6- Assign each read to a most likely OPR --------------------------------
# eg. 23-bp insertion >> could be a OPRI
# below also allows OPRI + OPRD in same read, although that is unlikely

# create the possible intervals for OPRDs
delInts <- vector(mode='list', length=4) # OPR deletion intervals
# minimum = 4 OPRD = - 24 * 4
# maximum = 1 OPRD = -24 * 1
delOpr <- seq(-24*4, -24*1, 24)
for (o in 1:length(delOpr)) {
  delInts[[o]] <- (delOpr[o] - lownoise) : (delOpr[o] + highnoise)
}
names(delInts) <- sprintf('OPRD%i', 4:1)

# create the possible intervals for OPRIs
inInts <- vector(mode='list', length=24) # OPR insertion intervals
# minimum = 1 OPRI = 24 * 1
# maximum = 24 OPRI = 24 * 24
inOpr <- seq(24*1, 24*24, 24)
for (o in 1:length(inOpr)) {
  inInts[[o]] <- (inOpr[o] - lownoise) : (inOpr[o] + highnoise)
}
names(inInts) <- sprintf('OPRI%i', 1:24)

# create the possible interval for reference
  # reference interval is a bit different as will overlap deletions/insertions
  # i.e. centered around 0, - lownoise / + highnoise
refInt <- (0-lownoise):(0+highnoise)

# given one total insertion (in bp), check if belongs to one OPRI interval
  # eg. +23 bp >> return OPRI1
checkOPRI <- function(totin) {
  
  check <- unlist(lapply(inInts, function(int){totin %in% int}))
  
  if (length(which(check==TRUE))==0){
    return(NA)
  } else(
    return(names(which(check) == TRUE))
  )
}

# given one total deletion (in bp), check if belongs to one OPRD interval
# eg. -23 bp >> return OPRD1
checkOPRD <- function(totdel) {
  
  check <- unlist(lapply(delInts, function(int){totdel %in% int}))
  
  if (length(which(check==TRUE))==0){
    return(NA)
  } else(
    return(names(which(check) == TRUE))
  )
}

# given both total deletion / total insertion, check if it should be categorised as reference
  # eg. -3 bp >> return reference
  # ! flipping around here: functions above are meant to be applied to all column; checkRef is meant to be apply to each row of ied
  # (because need both totdel & totins for each read)
checkRef <- function(read) { # read = row of ied
  
  totdel <- read['del']
  totins <- read['ins']
  
  # del and ins should be both in interval
    # eg. -5 del / +1 ins >> reference
    # eg. -12 del / +2 ins >> NA (in this case read will not have a category as deletion too small for an OPRD1)
  check <- (totdel %in% refInt) & (totdel %in% refInt) # check will be Logical, are both totdel / totins in Ref interval
  
  if (check) {
    read['del'] <- 'reference'
    read['ins'] <- 'reference'
    return(read)
  } else {
    return(read)
  }
}

# now assign reads to possible OPR mutations
ied$opri <- sapply(ied$ins, checkOPRI)
ied$oprd <- sapply(ied$del, checkOPRD)

# put reads that were not assigned to an OPRI and/or an OPRD in either 'reference' or 'unassigned'
  # 'reference' = total insertion & total deletion are in reference length interval (i.e. no or small insertion + no or small deletion)
  # 'unassigned' = none of total insertion or total deletion are in any interval (so anything falling between the cracks, eg. a 10 bp insertion)
# create a new column for these categories (prevents from plotting same reads in different plots below)
ied$othercategory <- NA
ied[which((ied$ins %in% refInt) & (ied$del %in% refInt)) , 'othercategory'] <- 'reference'
# explicitly call all the leftover reads 'unassigned'; leftover = 
  # do not have an OPRD
  # do not have an OPRI
  # are NA for othercategory (i.e. were not categorised as reference above)
ied[which(is.na(ied$opri) & is.na(ied$oprd) & is.na(ied$othercategory)) , 'othercategory'] <- 'unassigned'

# Note1; reads categorised as 'reference' or 'unassigned' will always be OPRI NA / OPRD NA
# Note2; it is a bit of a pain to keep 3 columns when really each read should be in a single category
  # (i.e. has an OPRI OR has OPRD OR is reference OR is unassigned)
  # the only advantage is to leave the opportunity for a single read to be both OPRD & OPRI
    # which are (hopefully) cases which will go away when aligning later

# Check
  # that means all reads should be in at least one category
  # i.e. all three columns cannot be NA at once
subset(ied, is.na(opri) & is.na(oprd) & is.na(othercategory))




# 7- Positive control for haplotype phasing -------------------------------
# eg. sample 8 OPRI: ~ 100% reads 8 OPRI should be one haplotype / ~ 100% reference should be the other haplotype
# similar plot as proportion of reads from each haplotype above, but one column for each sample / genotype

# make sure OPRI/OPRD/reference in order that makes sense
ied$opri <- factor(ied$opri, levels=sprintf('OPRI%i', 1:24))
ied$oprd <- factor(ied$oprd, levels=sprintf('OPRD%i', 4:1))
ied$othercategory <- factor(ied$othercategory, levels=c('reference', 'unassigned'))

# make haplotags easier to handle by making (for cases with a heterozygous OPR mutation):
  # mutated haplotype always haplotype1
  # reference haplotype always haplotype2
  # for PDG 46345: haplotype1 = 5 OPRI (as most important, probably causal mutation) / haplotype2 = 1 ORPD
# for samples that could not be haplotagged: all haplotype0 / for samples that are OPR reference/OPR reference (eg. control): leave it as it

# pdgsOPR = PDGs of sample with an OPR mutation (i.e. those for which we may need to swap haplotags)
pdgsOPR <- as.character(unlist(subset(meta, !is.na(SV), sample)))

  # ! do not pick PDGs we do not have here
  # (only applies if running the script on subset of the samples)
pdgsOPR <- intersect(pdgsOPR, pdgs)

# first set the whole haplotype column to integer (harder to change stuff is stays factor)
ied$haplotype <- as.integer(as.character(ied$haplotype))

for (p in 1:length(pdgsOPR)) {
  
  # check the PDG is not one we could not haplotag
  if (pdgsOPR[p] %in% pdgs_notag) {
    next # in which case, do nothing and try next PDG
  }
  
  osan <- as.character(subset(meta, sample==pdgsOPR[p], SV)) # take the Sanger genotype of that sample
  
  if (substr(osan, 1, 4) == 'OPRD') { # if Sanger genotype is a deletion
    
    # which haplotype matches the Sanger genotype?
    hap1 <- nrow(subset(ied, sample==pdgsOPR[p] & haplotype==1 & oprd==osan)) # how many reads from haplotype1 match the OPRD
    hap2 <- nrow(subset(ied, sample==pdgsOPR[p] & haplotype==2 & oprd==osan)) # how many reads from haplotype2 match the OPRD
    
    if (hap1 > hap2) { # if haplotype1 is mutated, no need to change anything
      cat('\t \t \t \t >>>', pdgsOPR[p], ': haplotype 1 is mutated one \n')
      next
    }
    
    else if (hap1 < hap2) {
      cat('\t \t \t \t >>>', pdgsOPR[p], ': haplotype 2 is mutated one \n')
      cat('\t \t \t \t \t swapping haplotags... \n')
      
      ied[which(ied$sample==pdgsOPR[p] & ied$haplotype==2), 'haplotype'] <- 99 # which() = these are the rows to change
      # ! do not 'cheat' and swap haplotags that are probably in the wrong haplotype
      # or in other words: do not look at the reads' OPR above, only at their tag
      # temporarily set these rows to 99
      
      # temporary set rows of the other haplotype to 9
      ied[which(ied$sample==pdgsOPR[p] & ied$haplotype==1), 'haplotype'] <- 9
      
      # now swap them
      ied[which(ied$sample==pdgsOPR[p] & ied$haplotype==99), 'haplotype'] <- 1 # haplotype 2 becomes haplotype 1
      ied[which(ied$sample==pdgsOPR[p] & ied$haplotype==9), 'haplotype'] <- 2 # and vice-versa
      
    } else stop('\t \t \t \t >>> Error: cannot find which haplotype matches the insertion for ', pdgsOPR[p], '\n')
  }
  
  ###
  
  if (substr(osan, 1, 4) == 'OPRI') { # if Sanger genotype is an insertion
    # same as above, just OPRI now
    
    # which haplotype matches best the Sanger genotype?
    hap1 <- nrow(subset(ied, sample==pdgsOPR[p] & haplotype==1 & opri==osan)) # how many reads from haplotype1 match the OPRI
    hap2 <- nrow(subset(ied, sample==pdgsOPR[p] & haplotype==2 & opri==osan)) # how many reads from haplotype2 match the OPRI
    
    if (hap1 > hap2) { # if haplotype1 is most likely the mutated one, no need to swap
      cat('\t \t \t \t >>>', pdgsOPR[p], ': haplotype 1 is mutated one \n')
      next
    }
    
    else if (hap1 < hap2) { # if haplotype2 is most likely the mutated one, no need to swap
      cat('\t \t \t \t >>>', pdgsOPR[p], ': haplotype 2 is mutated one \n')
      
      
      ied[which(ied$sample==pdgsOPR[p] & ied$haplotype==2), 'haplotype'] <- 99 # which() = these are the rows to change
      # ! do not 'cheat' and swap haplotags that are probably in the wrong haplotype
      # or in other words: do not look at the reads' OPR above, only at their tag
      # temporarily set these rows to 99
      
      # temporary set rows of the other haplotype to 9
      ied[which(ied$sample==pdgsOPR[p] & ied$haplotype==1), 'haplotype'] <- 9
      
      # now swap them
      cat('\t \t \t \t \t swapping haplotags... \n')
      ied[which(ied$sample==pdgsOPR[p] & ied$haplotype==99), 'haplotype'] <- 1 # haplotype 2 becomes haplotype 1
      ied[which(ied$sample==pdgsOPR[p] & ied$haplotype==9), 'haplotype'] <- 2 # and vice-versa
      
    } else stop('\t \t \t \t >>> Error: cannot find which haplotype matches the insertion for ', pdgsOPR[p], '\n')
  }
}

# turn haplotype column back into factor
ied$haplotype <- as.factor(ied$haplotype)
# in the order I want
ied$haplotype <- factor(ied$haplotype, levels=c(1, 2, 0))


# plot proportion of haplotype from each sample again
# between ### : same as above
###

hco <- ied %>% # hco for haplotype counts
  group_by(sample, haplotype) %>%
  tally(name='nreads')

# add total number of reads in each sample
tmp <- ied %>%
  count(sample, name='totalsamplereads')

hco <- left_join(hco, tmp, by='sample')

# order hco by pdgs
hco <- left_join(data.frame(sample=pdgs), hco, by='sample')

# normalise the read counts per haplotype
hco <- hco %>%
  mutate(haplotype_pro=nreads/totalsamplereads) %>% # haplotype proportion, i.e. proportion of reads assigned to haplotype 1, 2, or unassigned (0)
  mutate(sam_hap=paste(sample, haplotype, sep='_'), .after=haplotype) # add a unique sample_haplotype column, easier for plotting

# put levels of sample and sam_hap in right order for plotting
hco$sample <- factor(hco$sample, levels=unique(hco$sample))
hco$sam_hap <- factor(hco$sam_hap, levels=hco$sam_hap)
hco$haplotype <- factor(hco$haplotype, levels=c(1, 2, 0))

# haplotype assignments as stacked barplot, actual number of reads (so height = OPR coverage)
haploReads <- ggplot(hco, aes(x=sample, y=nreads, fill=haplotype)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=hapcols) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90),
    axis.title.x=element_blank()
  ) +
  ylab('number of reads')
haploReads

# haplotype assignments as stacked barplot, proportion of coverage (so height is 1 for all samples)

# I think best is same order as table all the way

# alternatively, samples without haplotag at the top
# add column tag or no tag
  # hco$htag <- TRUE
  # hco[which(hco$sample %in% pdgs_notag) , 'htag'] <- FALSE
  # hco <- hco %>%
  #   arrange(desc(htag))
  # hco$sample <- factor(hco$sample, levels=unique(hco$sample))

# same order as meta
hco <- left_join(data.frame(sample=meta$sample), hco, by='sample')
hco <- hco[seq(dim(hco)[1],1),]
hco$sample <- factor(hco$sample, levels=unique(hco$sample))

# sample names for y axis
splsyaxis <- paste0('#', substr(unique(hco$sample), 4, 99))

haploPro <- ggplot(hco, aes(x=sample, y=haplotype_pro, fill=haplotype)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_manual(values=hapcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    axis.title.x=element_text(size=9, margin = margin(t=2, r=0, b=0, l=0)),
    axis.text.x=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_blank(),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0), hjust=0)
  ) +
  ylab('proportion of reads') +
  scale_x_discrete(labels=splsyaxis)
haploPro

ggsave(here('haplotypephasing/haploReads2.pdf'), haploReads,width=200, height=100, units='mm')
ggsave(here('haplotypephasing/haploPro2.pdf'), haploPro, width=65, height=90, units='mm')

###

# plot OPRD / OPRI stacked barplots for each sample (fill = haplotype)
  # ! reads that have both an OPRI & an OPRD will be counted twice (i.e. once in OPRI plot + once in OPRD plot), but they are pretty rare

dir.create(here('needleSam', 'OPRHaploPlots'))

for (s in 1:length(pdgs)) {
  
  cat('\t \t \t \t >>> Plotting OPR by haplotype for sample ', s, '/', length(pdgs), '\n')
  
  # common parameters (easier for debugging)
  yaxislimit <- 1.0
  barwidth <- 1
  
  # get the sample we are plotting
  pd <- pdgs[s]
  
  # create a small list for the 4 plots (OPRD / reference / OPRI / unassigned)
  ohapPlots <- vector(mode='list', length=4)
  names(ohapPlots) <- c('OPRD', 'reference', 'OPRI', 'unassigned')
  
  # OPRD plot
    # first prepare the data
  tmp <- ied %>%
    subset(sample==pd) %>%
    group_by(oprd, haplotype) %>%
    tally(name='nreads') %>%
    mutate(freq=nreads/nrow(subset(ied, sample==pd))) # frequency = proportion to *total* number of reads of that sample
  
    # plot
  ohapPlots$OPRD <- ggplot(subset(tmp, !is.na(oprd)), aes(x=oprd, y=freq, fill=haplotype)) +
    geom_col(width=barwidth) +
    scale_fill_manual(values=hapcols) +
    theme_minimal() +
    theme(
      axis.title.x=element_blank()
    ) +
    coord_cartesian(ylim=c(0, yaxislimit))
  
  # reference plot
    # first prepare the data
  tmp <- ied %>%
    subset(sample==pd) %>%
    group_by(othercategory, haplotype) %>%
    tally(name='nreads') %>%
    mutate(freq=nreads/nrow(subset(ied, sample==pd))) # frequency = proportion to *total* number of reads of that sample
  
    # plot
  ohapPlots$reference <- ggplot(subset(tmp, othercategory=='reference'), aes(x=othercategory, y=freq, fill=haplotype)) +
    geom_col(width=barwidth) +
    scale_fill_manual(values=hapcols) +
    theme_minimal() +
    theme(
      axis.title.x=element_blank()
    ) +
    coord_cartesian(ylim=c(0, yaxislimit))
  
  # OPRI plot
  # first prepare the data
  tmp <- ied %>%
    subset(sample==pd) %>%
    group_by(opri, haplotype) %>%
    tally(name='nreads') %>%
    mutate(freq=nreads/nrow(subset(ied, sample==pd))) # frequency = proportion to *total* number of reads of that sample
  
  # plot
  ohapPlots$OPRI <- ggplot(subset(tmp, !is.na(opri)), aes(x=opri, y=freq, fill=haplotype)) +
    geom_col(width=barwidth) +
    scale_fill_manual(values=hapcols) +
    theme_minimal() +
    theme(
      axis.title.x=element_blank()
    ) +
    coord_cartesian(ylim=c(0, yaxislimit))
  
  # reference plot
    # first prepare the data
  tmp <- ied %>%
    subset(sample==pd) %>%
    group_by(othercategory, haplotype) %>%
    tally(name='nreads') %>%
    mutate(freq=nreads/nrow(subset(ied, sample==pd))) # frequency = proportion to *total* number of reads of that sample
  
    # plot
  ohapPlots$unassigned <- ggplot(subset(tmp, othercategory=='unassigned'), aes(x=othercategory, y=freq, fill=haplotype)) +
    geom_col(width=barwidth) +
    scale_fill_manual(values=hapcols) +
    theme_minimal() +
    theme(
      axis.title.x=element_blank()
    ) +
    coord_cartesian(ylim=c(0, yaxislimit))
  
  # prepare the grid (each plot alongside each other)
  ohapGrid <- ggarrange(plotlist=ohapPlots, ncol=length(ohapPlots), nrow=1)
  # save as pdf
  ggsave(filename=here('needleSam', 'OPRHaploPlots', paste0('OPRHap_', pd, '.pdf')), plot=ohapGrid, width=800, height=200, units='mm')
  
}

# these plots are not really good, but good enough for what I want to check now
  # I will improve these plots if we include them

# simplify plots above (make it a single plot) by
  # 1- plot only samples with an OPR & which have haplotags
  # 2- plot only that OPR + reference (or OPRD1 / OPRI5 for the compound homozygous sample)
# there are 14 samples with an OPR genotype; each 2 alleles, so 28 columns (42 if plot unassigned)
pdgoptag <- as.character(unlist(subset(meta, !is.na(SV) & ! sample %in% pdgs_notag, sample))) # these are the samples we want to plot

# ! do not pick PDGs we do not have here
  # (only applies if running the script on subset of the samples)
pdgoptag <- intersect(pdgoptag, pdgs)

# tally grossly the samples we need to plot, then will delete what we do not want to include in the plot
oprcount <- ied %>%
  subset(sample %in% pdgoptag) %>%
  group_by(sample, opri, oprd, othercategory, haplotype) %>%
  tally(name='nreads')

# for each sample, delete what does not correspond to its genotype
for (s in 1:length(pdgoptag)) {
  
  # get the sample PDG
  pd <- pdgoptag[s]
  
  # ! when reaching PDG 46345, consider the other OPR genotype and delete all the reference counts
  if (pd=='pdg46345') {
    
    ogen <- as.character(subset(meta, sample==pd, c(SV, SV2))) # take the 2 SVs
    
    # delete from oprcount everything that is not these 2 genotypes
    # or in other words, keep everything that matches
    oprcount <- oprcount[- which(oprcount$sample==pd &
                                   (!oprcount$opri %in% ogen) & (!oprcount$oprd %in% ogen)) , ]
    
    # there is sometimes a read or two that are eg. OPRI4 + OPRD1, and as OPRD1 is correct it passes filter above
      # remove also any read that has both an OPRI and an OPRD (if any)
    inds2del <- which(oprcount$sample==pd & !is.na(oprcount$opri) & !is.na(oprcount$oprd))
    if (length(inds2del)!=0) {
      oprcount <- oprcount[ - inds2del , ]
    }
    
  } else { # for the other samples
    # get its OPR genotype
    ogen <- as.character(subset(meta, sample==pd, SV))
    
    # delete from oprcount everything that is not that genotype or reference
    oprcount <- oprcount[- which(oprcount$sample==pd
                   & (oprcount$opri!=ogen | oprcount$oprd!=ogen | oprcount$othercategory!='reference')) , ]
    
    # same as above, remove any read that has both OPRD and OPRI (if any)
    inds2del <- which(oprcount$sample==pd & !is.na(oprcount$opri) & !is.na(oprcount$oprd))
    if (length(inds2del)!=0) {
      oprcount <- oprcount[ - inds2del , ]
    }
  }
  
}

# add proportion of each haplotype within genotype categories
oprcount <- oprcount %>%
  group_by(sample, opri, oprd, othercategory, .drop=FALSE) %>%
  mutate(totalreads_samplegenotype=sum(nreads)) %>%
  mutate(freq=nreads/totalreads_samplegenotype) %>%
  mutate(sam_ogen=paste(sample, opri, oprd, othercategory, sep='_')) # add an ID column

# plot
# will do in order of OPR, which should be in order of meta (minus the ones not present)
# same order as meta
oprcount <- left_join(data.frame(sample=meta$sample), oprcount, by='sample')
oprcount <- oprcount[seq(dim(oprcount)[1],1),] # reverse the rows
oprcount$sam_ogen <- factor(oprcount$sam_ogen, levels=unique(oprcount$sam_ogen))
# delete all the NA rows (from left_join above)
oprcount <- oprcount[-which(is.na(oprcount$sam_ogen)),]

oprHaplot <- ggplot(oprcount, aes(x=sam_ogen, y=freq, fill=haplotype)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=hapcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    axis.title.x=element_text(size=9, margin = margin(t=2, r=0, b=0, l=0)),
    axis.text.x=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
    axis.text.y=element_blank(),
    axis.title.y=element_blank()) +
  ylab('proportion of reads')
oprHaplot

ggsave(filename=here('needleSam', 'oprHapPlot.pdf'), plot=oprHaplot, width=54, height=90, units='mm')


# haplotags look ok but not perfect
# can be up to ~ 10% wrong haplotype, eg. OPRI8: 90% from haplotype1 / 10% from haplotype2
# could also be that error is from detecting OPRI8, but cannot believe that a ~ 192 bp insertion can be there 'by accident'
# additionally, errors are clearly symmetrical; in example above: reference for that sample will be 90% haplotype2 / 10% haplotype1
# so I think error is = wrong haplotype assignment, i.e. from whatshap haplotag

# there are 2 questions
  # 1- why the UNassigned reads (reads which get assigned haplotype 0)
  # 2- why the MISassigned reads (reads which get assigned what I think is the wrong haplotype from looking at the SV)

# Example
# PDG 47875 = 8 OPRI
# good example because
  # 14 SNVs available, so we do not expect whatshap to do many errors
  # from looking at phased VCF: all SNVs are from one haplotype
  # by looking at phased reads in IGV: OPRI8 haplotype is the one carrying the 14 SNVs

# I cannot see anything wrong with these 'haplotype-swapped' reads; they look real to me
# I see 3 options
  # 1- artefact on flowcell -- I do not think so, otherwise we would find tons of eg. OPRI6 on flowcell3 (as there are many such reads around)
  # 2- PCR artefact -- precisely template switching, i.e. the polymerase would be reading a hap2 molecule,
    # then switch just before the OPR to a hap1 molecule
  # 3- inter-chromosomal rearrangement
# I do not see a dry-lab way of differentiating the last two

# Looking deeper into swapped haplotypes

# proportion of swapped haplotypes per sample
  # for PDG 46345, we are interested in hap1 that are OPRD1 / hap2 that are OPRI5
swa <- subset(oprcount,
       (sample!='pdg46345' & haplotype==1 & othercategory=='reference') |
         (sample!='pdg46345' & haplotype==2 & !is.na(opri)) |
         (sample!='pdg46345' & haplotype==2 & !is.na(oprd)) |
         
         (sample=='pdg46345' & haplotype==1 & !is.na(oprd)) |
         (sample=='pdg46345' & haplotype==2 & !is.na(opri)) )

# is it true that it looks symmetrical?
swappedHap <- ggplot(swa, aes(x=sam_ogen, y=freq, fill=haplotype)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=hapcols) +
  theme_minimal() +
  theme(
    axis.title.y=element_blank()
  )
swappedHap
ggsave(filename=here('haplotypephasing', 'swappedHap.pdf'), plot=swappedHap, width=300, height=200, units='mm')
# yes, clear

# does it correlate with anything?

  # add every kind of sample info
swa <- left_join(swa, meta, by='sample')
swa <- left_join(swa, hpmeta, by='sample')

# add proportion left not haplotagged
proze <- select(hcoun, sample, haplotype_pro) # proportion haplotype0
colnames(proze) <- c('sample', 'prohap0')
swa <- left_join(swa, proze)

# just as sanity check
swappedHapVsnReads <- ggplot(swa, aes(x=nreads, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsnReads
# nreads is how frequency is calculated so obviously makes sense

# vs Coverage
swappedHapVsCov <- ggplot(swa, aes(x=totalreads_samplegenotype, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsCov
ggsave(filename=here('haplotypephasing', 'swappedHapVsCov.pdf'), plot=swappedHapVsCov, width=300, height=200, units='mm')

# vs tissue
swappedHapVsTissue <- ggplot(swa, aes(x=gDNA_tissue, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsTissue
ggsave(filename=here('haplotypephasing', 'swappedHapVsTissue.pdf'), plot=swappedHapVsTissue, width=300, height=200, units='mm')

# vs Cohort
  # cannot do as only inherited samples

# vs Flowcell
swappedHapVsFlow <- ggplot(swa, aes(x=flowcell, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols) +
  scale_x_continuous(breaks=c(1, 3))
swappedHapVsFlow
ggsave(filename=here('haplotypephasing', 'swappedHapVsFlow.pdf'), plot=swappedHapVsFlow, width=300, height=200, units='mm')

# vs OPR coverage
# genebody_coverage in in sequencing_summary sheet of AdditionalFile1.xlsx

swa <- seqmeta %>%
  select('sample', 'genebody_coverage') %>%
  left_join(swa, ., by='sample')

swappedHapVsOPRcov <- ggplot(swa, aes(x=genebody_coverage, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsOPRcov
ggsave(filename=here('haplotypephasing', 'swappedHapVsOPRcov.pdf'), plot=swappedHapVsOPRcov, width=300, height=200, units='mm')

# vs MRC score
swappedHapVsMRC <- ggplot(swa, aes(x=MRC_slope, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsMRC
ggsave(filename=here('haplotypephasing', 'swappedHapVsMRC.pdf'), plot=swappedHapVsMRC, width=300, height=200, units='mm')

# vs age
swappedHapVsAge <- ggplot(swa, aes(x=age_atsample, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsAge
ggsave(filename=here('haplotypephasing', 'swappedHapVsAge.pdf'), plot=swappedHapVsAge, width=300, height=200, units='mm')

# vs gender
swappedHapVsGender <- ggplot(swa, aes(x=gender, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsGender
ggsave(filename=here('haplotypephasing', 'swappedHapVsGender.pdf'), plot=swappedHapVsGender, width=300, height=200, units='mm')

# vs number of SNVs available for haplotagging
swappedHapVsHaplotagSNVs <- ggplot(swa, aes(x=genebody_heterozygousSNVs, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsHaplotagSNVs
ggsave(filename=here('haplotypephasing', 'swappedHapVsHaplotagSNVs.pdf'), plot=swappedHapVsHaplotagSNVs, width=300, height=200, units='mm')

# vs proportion not haplotagged
swappedHapVsProHap0 <- ggplot(swa, aes(x=prohap0, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsProHap0
ggsave(filename=here('haplotypephasing', 'swappedHapVsProHap0.pdf'), plot=swappedHapVsProHap0, width=300, height=200, units='mm')

# vs OPR
swa$SV <- factor(swa$SV, levels=c('OPRD2', 'OPRD1', 'OPRI1', 'OPRI2', 'OPRI4', 'OPRI5', 'OPRI6', 'OPRI8'))
swappedHapVsOPR <- ggplot(swa, aes(x=SV, y=freq, colour=haplotype)) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values=hapcols)
swappedHapVsOPR
ggsave(filename=here('haplotypephasing', 'swappedHapVsOPR.pdf'), plot=swappedHapVsOPR, width=300, height=200, units='mm')

# overall nothing striking


# do swapped reads have lower haplotag quality?

# first need to tag swapped reads

swais <- c() # here, will store indices of swapped reads (swapped indices)

for (p in 1:length(pdgoptag)) { # for each sample that was haplotagged + has an SV
  # if no haplotag or no SV: cannot tell if it is a swapped read
  
  # get the pdg we are working with
  pdg <- pdgoptag[p]
  
  # get its SV
  sv <- as.character(subset(meta, sample==pdg, SV))
  if (is.na(sv)) stop('\t \t \t \t Error: SV is not expected to be NA \n')
  
  # ! bit more complicated for PDG 46345 (compound homozygous)
  if (pdg == 'pdg46345') {
    
    # indices to change
      # any haplotype1 that says OPRD1
    inds <- which(ied$sample==pdg &
                    ied$haplotype==1 &
                    ied$oprd=='OPRD1')
    swais <- c(swais, inds)
    
      # any haplotype2 that says OPRI5
    inds <- which(ied$sample==pdg &
                    ied$haplotype==2 &
                    ied$oprd=='OPRI5')
    swais <- c(swais, inds)
    
  }
  
  else { # in the other cases
    
    # if SV is an OPRI
    if (substr(sv, 1, 4)=='OPRI') {
      
        # any haplotype1 that says reference
      inds <- which(ied$sample==pdg &
              ied$haplotype==1 &
              ied$othercategory=='reference')
      swais <- c(swais, inds)
      
        # any haplotype2 that says OPRI
      inds <- which(ied$sample==pdg &
              ied$haplotype==2 &
              ied$opri==sv)
      swais <- c(swais, inds)
      
    }
    
    # if SV is an OPRD
    if (substr(sv, 1, 4)=='OPRD') {
      
      # any haplotype1 that says reference
      inds <- which(ied$sample==pdg &
                      ied$haplotype==1 &
                      ied$othercategory=='reference')
      swais <- c(swais, inds)
      
      # any haplotype2 that says OPRD
      inds <- which(ied$sample==pdg &
                      ied$haplotype==2 &
                      ied$oprd==sv)
      swais <- c(swais, inds)
      
    }
  }
}

# few checks
  # how many swapped reads?
length(swais)
  # out of
length(pdgoptag) # samples
  # so around
length(swais) / length(pdgoptag) # per sample
  # sounds about right

  # there should not be duplicates
sum(duplicated(swais))

# add column swappedOPR
  # preallocate as all FALSE and just change the indices above to TRUE
swacol <- rep(FALSE, nrow(ied))
swacol[swais] <- TRUE
length(swacol) == nrow(ied)

ied$swappedOPR <- swacol


# now need to make haplotype quality more comparable between samples
  # currently: highly dependent on number of SNVs, usually number of SNVs correct * 30
# so if only 1 SNV: always 30 or 0
  # which cannot be compared to samples with many SNVs
  # eg. read with 14 SNVs, but only 2 match; I want a low haplotag score
  # but currently it would be 60, so higher than a read with 1 SNV that matches
# a possible solution: % of available SNVs that match (should be ~ % of maximum for that sample)
  # in examples above
    # 1 SNV: match = 1.0 / no-match = 0.0
    # 14 SNV: only 2 match = 0.14
besthq <- ied %>%
  group_by(sample) %>%
  summarise(besthapqual=max(hapqual, na.rm=TRUE))
# generates -Inf for the sample which do not have haplotags, turn into NA
besthq[which(besthq$besthapqual==-Inf), 'besthapqual'] <- NA

# add to ied
  # keep all rows of ied, so left_join
ied <- left_join(x=ied, y=besthq, by='sample')

# divide each hapqual by the best of its sample
ied$hapscore <- ied$hapqual/ied$besthapqual

# back to question: is haplotype quality lower for swapped reads?

# distribution of all scores
hapScoreDist <- ggplot(ied, aes(hapscore)) +
  geom_histogram(binwidth=0.1) +
  theme_minimal()
hapScoreDist
ggsave(filename=here('haplotypephasing', 'hapScoreDist.pdf'), plot=hapScoreDist, width=300, height=200, units='mm')

# same, but differentiated by swapped or not

hsbsw <- freqInBins(df=ied, datcol='hapscore', grpcol='swappedOPR', bins=seq(0, 1, by=0.1)) # haplotype scores by swapped or not

hapScoreSwapped <- ggplot(hsbsw, aes(x=bin, y=freq, fill=swappedOPR)) +
  geom_col(position='identity', alpha=0.5) +
  theme_minimal()
hapScoreSwapped
ggsave(filename=here('haplotypephasing', 'hapScoreSwapped.pdf'), plot=hapScoreSwapped, width=300, height=200, units='mm')

# yes swapped reads may have worse haplotype quality
t.test(ied$hapscore ~ ied$swappedOPR)




# Check strand directions -------------------------------------------------

# check nothing else than 0 (= forward) or 16 (= reverse complement) in samflag column
which(ied$samflag!=0 & ied$samflag!=16) # ok

# if yes -- replace all 0 to +, all 16 to -
ied[which(ied$samflag==0), 'samflag'] <- '+'
ied[which(ied$samflag==16), 'samflag'] <- '-'

# check all +/- now
which(ied$samflag!='+' & ied$samflag!='-')
# ok

# change column name to strand
colnames(ied)[which(colnames(ied)=='samflag')] <- 'strand'

# quick sanity checks with strand direction
strandcount <- ied %>%
  group_by(strand) %>%
  tally(name='count') %>%
  mutate(prop=count/nrow(ied)) 

strandPro <- ggplot(strandcount, aes(x=strand, y=prop)) +
  geom_col() +
  theme_minimal() +
  xlab('proportion of reads') + ylab('strand')
strandPro
ggsave(here('haplotypephasing', 'strand.pdf'), strandPro, width=100, height=100, units='mm')

# add a dummy column
strandcount$dummy <- factor(1)

# same but stack
strandProStack <- ggplot(strandcount, aes(x=dummy, y=prop, fill=strand)) +
  geom_col() +
  theme_minimal() +
  xlab('proportion of reads')
strandProStack
ggsave(here('haplotypephasing', 'strandProStack.pdf'), strandProStack, width=100, height=100, units='mm')


# vs haplotype
strandhap <- ied %>%
  group_by(strand, haplotype) %>%
  tally(name='count')

# counts
strandvsHapCount <- ggplot(strandhap, aes(x=haplotype, y=count, fill=strand)) +
  geom_col() +
  theme_minimal()
strandvsHapCount
ggsave(here('haplotypephasing', 'strandvsHapCount.pdf'), strandvsHapCount, width=100, height=100, units='mm')

# proportion
hapcount <- ied %>%
  group_by(haplotype) %>%
  tally(name='countHap')
strandhap <- left_join(x=strandhap, y=hapcount)
strandhap <- strandhap %>%
  mutate(propHap=count/countHap)

strandvsHapPro <- ggplot(strandhap, aes(x=haplotype, y=propHap, fill=strand)) +
  geom_col() +
  theme_minimal() +
  geom_hline(yintercept=0.5, linetype=2)
strandvsHapPro
ggsave(here('haplotypephasing', 'strandvsHapPro.pdf'), strandvsHapPro, width=100, height=100, units='mm')

# other way to look at it: assignment to haplotype for each strand direction
# add strand count
strandcount <- ied %>%
  group_by(strand) %>%
  tally(name='countStrand')
strandhap <- left_join(x=strandhap, y=strandcount)
strandhap <- strandhap %>%
  mutate(propStrand=count/countStrand)

hapvsStrandPro <- ggplot(strandhap, aes(x=strand, y=propStrand, fill=haplotype)) +
  geom_col() +
  theme_minimal() +
  scale_fill_manual(values=hapcols)
hapvsStrandPro
ggsave(here('haplotypephasing', 'hapvsStrandPro.pdf'), hapvsStrandPro, width=100, height=100, units='mm')


# vs flowcell
  # temporary copy of ied and add flowcell info to it
tmp <- ied
tmp <- left_join(x=tmp, y=meta)

flowstrand <- tmp %>%
  group_by(strand, flowcell) %>%
  tally(name='count')
  # add total number of reads on each flowcell
flowcount <- tmp %>%
  group_by(flowcell) %>%
  tally(name='countFlow')
flowstrand <- left_join(x=flowstrand, y=flowcount)
flowstrand <- flowstrand %>%
  mutate(propFlow=count/countFlow)
flowstrand$flowcell <- factor(flowstrand$flowcell, levels=c(1, 2, 3))

strandvsFlowPro <- ggplot(flowstrand, aes(x=flowcell, y=propFlow, fill=strand)) +
  geom_col() +
  theme_minimal() +
  geom_hline(yintercept=0.5, linetype=2) +
  scale_fill_manual(values=hapcols)
strandvsFlowPro
ggsave(here('haplotypephasing', 'strandvsFlowPro.pdf'), strandvsFlowPro, width=100, height=100, units='mm')


# vs cohort
  # can use same tmp copy
cohstrand <- tmp %>%
  group_by(strand, prion_disease) %>%
  tally(name='count')
  # add total number of reads in each cohort
cohcount <- tmp %>%
  group_by(prion_disease) %>%
  tally(name='countCoh')
cohstrand <- left_join(x=cohstrand, y=cohcount)
cohstrand <- cohstrand %>%
  mutate(propCoh=count/countCoh)

strandvsCohort <- ggplot(cohstrand, aes(x=prion_disease, y=propCoh, fill=strand)) +
  geom_col() +
  theme_minimal() +
  geom_hline(yintercept=0.5, linetype=2)
strandvsCohort
ggsave(here('haplotypephasing', 'strandvsCohort.pdf'), strandvsCohort, width=100, height=100, units='mm')

# vs sample
  # can use same tmp copy
samstrand <- tmp %>%
  group_by(strand, sample) %>%
  tally(name='count')
  # add total number of reads in each cohort
samcount <- tmp %>%
  group_by(sample) %>%
  tally(name='countSam')
samstrand <- left_join(x=samstrand, y=samcount)
samstrand <- samstrand %>%
  mutate(propSam=count/countSam)

strandvsSample <- ggplot(samstrand, aes(x=sample, y=propSam, fill=strand)) +
  geom_col() +
  geom_hline(yintercept=0.5, linetype=2) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90)
  )
strandvsSample
ggsave(here('haplotypephasing', 'strandvsSample.pdf'), strandvsSample, width=250, height=150, units='mm')

# related to above: would strand bias correlate with haplotype-swapped reads?
  # just keep positive strand proportion (negative is just 1 - that)
samplus <- samstrand %>%
  subset(strand=='+') %>%
  ungroup() %>%
  select(sample, propSam)
colnames(samplus) <- c('sample', 'propPlus')
  
swa <- left_join(x=swa, y=samplus)

strandvSwapped <- ggplot(swa, aes(x=propPlus, y=prohap0)) +
  geom_point() +
  theme_minimal()
strandvSwapped
ggsave(here('haplotypephasing', 'strandvSwapped.pdf'), strandvSwapped, width=100, height=100, units='mm')




# 8- Alignment to consensus OPR sequences ---------------------------------

# for example, say a read has an insertion of 92 bp >> it could be a OPRI4
# how can we measure how convincing it is?
# we will align to a consensus sequence of what an OPRI4 looks like
# this makes use of the consensus R sequence written by Simon, which uses flexible nucleotides (eg. Y = C or T)
# >> see OPRconsensus.xlsx / generateOPRCatalog.R / OPRConsensusCatalog.xlsx for more information

# previously: was trying alignment to published OPR sequences
  # but there is some logic to which nucleotides can change, so better to define a consensus sequence
  # advantage: it includes more possible OPRs (OPR sequences that are probably out there but were never detected or published)


# function msaConsensusFlexible takes an alignment and outputs a corrected consensus sequence
  # ! msa support alignment with flexible nucleotides, see https://github.com/UBod/msa/issues/12
  # but eg. G vs A returns R; which is correct but not what I want here
  # I want to leave flexibility only if the reference leaves that option open
  # where eg. Y vs T returns Y
# build look-up table for flexible nucleotides
fnu <- c('R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N') # flexible nucleotides
lut <- as.data.frame(matrix(nrow=4, ncol=length(fnu)))
colnames(lut) <- fnu
lut$R[1:2] <- c('A', 'G')
lut$Y[1:2] <- c('C', 'T')
lut$S[1:2] <- c('G', 'C')
lut$W[1:2] <- c('A', 'T')
lut$K[1:2] <- c('G', 'T')
lut$M[1:2] <- c('A', 'C')
lut$B[1:3] <- c('C', 'G', 'T')
lut$D[1:3] <- c('A', 'G', 'T')
lut$H[1:3] <- c('A', 'C', 'T')
lut$V[1:3] <- c('A', 'C', 'G')
lut$N[1:4] <- c('A', 'C', 'G', 'T')

# msaConsensusFlexible
  # use instead of msaConsensusSequence, will return the consensus sequence of an alignment, but
  # will replace ? in the consensus sequence which were thrown because of the presence of a flexible nucleotide, where the alignment was fine (eg. Y vs T)
  # ! Note; it assumes alignment one-to-one and only one of the sequence (typically the reference) allows flexible nucleotides
msaConsensusFlexible <- function (alignment){
  conm <- consensusMatrix(alignment)
  cons <- msaConsensusSequence(alignment)
  
  sapply(1:ncol(conm), function(pos) { # for each position of the alignment
    
    al <- conm[,pos] # get the alignment matrix at that position
    # count number of flexible nucleotides at that position
    nfnu <- sum(al[which(names(al) %in% fnu)])
    
    if (nfnu==0) { # if no flexible nucleotide at that position, do nothing
      return(NULL)
    }
    
    if (nfnu>0) { # if any flexible nucleotide at that position
      # get the flexible nucleotide
      f <- names(which(al[which(names(al) %in% fnu)] > 0))
      
      # get the standard nucleotide
      s <- names(which(al[which(!names(al) %in% fnu)] > 0))
      
      # is the alignment valid? Check in the look-up table; eg. T vs Y is ok
      if (!s %in% lut[,f]) { # if not -- do NOT modify the consensus sequence at that position (should leave ?)
        return(NULL)
      }
      
      else if (s %in% lut[,f]) { # if yes -- modify the consensus sequence at that position (should put the flexible nucleotide instead of the ?)
        substr(cons, pos, pos) <<- f
      }
    }
    
  })
  
  return(cons)
  
}


# function alignCon_NumMis:
  # takes a query DNA sequence (i.e. possible somatic mutation)
  # looks at the possible OPRD/OPRI for that read
  # takes consensus DNA sequence from the catalog of OPRs (by allele)
  # try aligning to the consensus sequence, gets number of mismatches
  # reports the number of mismatches vs consensus
alignCon_NumMis <- function(query, allele) {
  # eg. alignRef_NumMis('CGAAATGC', 'OPRD1'), i.e. try aligning CGAAATGC to the consensus OPRD1 sequence and tell me number of mismatches
  # can also be alignRef_NumMis('CGAAATGC', 'reference')
  
  if( is.na(allele) | allele=='unassigned') { # if allele is NA (eg. a 8OPRI read will usually have NA in oprd column; all unassigned reads are opri/oprd NA; etc.)
    return(NA)
  }
  
  if(! allele %in% cato$genotype) { # check the allele is in catalog (I think should always be ok as catalog covers all possible cases)
    stop('\t \t \t \t >>> Allele', allele, ' was not found in catalog \n')
  }
  
  
  # 1- take the consensus sequence (used as reference)
  ref <- as.character(subset(cato, genotype==allele, consensusseq))
  
  # 2- align query to that reference
  ali <- msa(c(ref, query), type='dna', method='ClustalW')
  
  # 3- correct the consensus sequence to allow flexible nucleotides
    # using msaConsensusFlexible()
  con <- msaConsensusFlexible(ali)
  
  # 4- count number of mismatches (?)
  nmis <- str_count(con, '\\?')
  
  # return that
  return(nmis)
  
}

# sanity checks
  # align reference OPR to the consensus sequence
alignCon_NumMis('CCTCAGGGCGGTGGTGGCTGGGGGCAGCCTCATGGTGGTGGCTGGGGGCAGCCTCATGGTGGTGGCTGGGGGCAGCCCCATGGTGGTGGCTGGGGACAGCCTCATGGTGGTGGCTGGGGTCAA',
                'reference')
# gives 0 mismatch

  # check eg. G vs A should NOT return R; only edit consensus if reference leaves the option open
alitest <- msa(c('CYGTK', 'CTAAT'), type='dna', method='ClustalW') # should return C Y ? ? K
msaConsensusFlexible(alitest) # ok


tic('\t \t \t \t Aligning reads took... \n')

# align all possible OPRI reads to their consensus sequence
opri_mis <- sapply(1:nrow(ied), function(ri){ # ri = read index; opri_mats = OPRI mismatches
  # Note; will be one element per row of ied, can add to ied as it
  cat('\t \t \t \t >>> Read', ri, 'out of', nrow(ied), '\n')
  
  readseq <- as.character(ied[ri, 'seq']) # read sequence
  oprgeno <- as.character(ied[ri, 'opri']) # possible OPR allele
  
  alignCon_NumMis(readseq, oprgeno)
})

# put the results in ied
ied$opri_nummis <- opri_mis



# align all possible OPRD reads to their consensus sequence
oprd_mis <- sapply(1:nrow(ied), function(ri){ # ri = read index; oprd_mis = OPRD mismatches
  cat('\t \t \t \t >>> Read', ri, 'out of', nrow(ied), '\n')
  
  readseq <- as.character(ied[ri, 'seq']) # read sequence
  oprgeno <- as.character(ied[ri, 'oprd']) # possible OPR genotype
  
  alignCon_NumMis(readseq, oprgeno)
})

# put the results in ied
ied$oprd_nummis <- oprd_mis



# align reference reads to the consensus sequence
ref_mis <- sapply(1:nrow(ied), function(ri){ # ri = read index; ref_mis = reference mismatches
  cat('\t \t \t \t >>> Read', ri, 'out of', nrow(ied), '\n')
  
  readseq <- as.character(ied[ri, 'seq']) # read sequence
  oprgeno <- as.character(ied[ri, 'othercategory']) # take what is othercategory column as allele
  # Note alignCon_NumMis takes care of 'unassigned' already, i.e. returns NA (as there is nothing to align to)
  
  alignCon_NumMis(readseq, oprgeno)
})

# put the results in ied
ied$reference_nummis <- ref_mis


toc()




# 9- Threshold for maximum number of mismatches? --------------------------
# what is a good threshold for maximum number of mismatches?
# eg. if 25% of the nucleotides in one read do not match the consensus OPR sequence, might not be so convincing
# >> look at the 'positive control' distribution of number of mismatches
  # eg. take sample genotyped as OPRD2 by Sanger
  # take all the reads that are consistent with an OPRD2
  # and look at their number of mismatches with OPRD2 consensus sequence

# add samples metadata
ied <- left_join(x=ied, y=meta, by='sample')
  # we do not need the colour column
ied$colour <- NULL

# plot mismatches plot
# list to store mismatches plots
# will basically do one list for each copy of the genome, i.e.
  # first list = all the SVs or reference for the samples w/o SV
  # second list = reference for the samples that have a SV / OPRD1 for the compound homozygous sample

# List 1
misPlots1 <- vector(mode='list', length=length(pdgs))
names(misPlots1) <- pdgs

for (s in 1:length(pdgs)) {
  
  # get the PDG we are plotting
  pdg <- pdgs[[s]]
  
  # get its Sanger genotype
  sang <- as.character(unique(subset(ied, sample==pdg, SV)))
  
  
  # if Sanger genotype is NA (example sCJD case)
  if (is.na(sang)) { # then plot the mismatches vs reference
    
    # then it means its genotype is reference
    sang <- 'reference'

    tmp <- subset(ied, sample==pdg & othercategory==sang) # subset of ied for plot
    misplot <- ggplot(tmp, aes(reference_nummis, after_stat(density))) + # mismatches plot
      geom_histogram(binwidth=1) +
      theme_minimal() +
      coord_cartesian(xlim=c(0,100), ylim=c(0, 0.3)) +
      labs(title=pdg,
           subtitle=sang)
    misPlots1[[s]] <- misplot # add it to the list
  
    
  # if Sanger genotype is an insertion...
  } else if (substr(sang, 1, 4) == 'OPRI') {
      
    tmp <- subset(ied, sample==pdg & opri==sang) # subset of ied for plot
    misplot <- ggplot(tmp, aes(opri_nummis, after_stat(density))) + # mismatches plot
      geom_histogram(binwidth=1) +
      theme_minimal() +
      coord_cartesian(xlim=c(0,100), ylim=c(0, 0.3)) +
      labs(title=pdg,
           subtitle=sang)
    misPlots1[[s]] <- misplot # add it to the list
    
    # if Sanger genotype is a deletion...
  } else if (substr(sang, 1, 4) == 'OPRD') {
    
    tmp <- subset(ied, sample==pdg & oprd==sang) # subset of ied for plot
    misplot <- ggplot(tmp, aes(oprd_nummis, after_stat(density))) + # mismatches plot
      geom_histogram(binwidth=1) +
      theme_minimal() +
      coord_cartesian(xlim=c(0,100), ylim=c(0, 0.3)) +
      labs(title=pdg,
           subtitle=sang)
    misPlots1[[s]] <- misplot # add it to the list
  }
}

# List 2
misPlots2 <- vector(mode='list', length=length(pdgs))
names(misPlots2) <- pdgs

for (s in 1:length(pdgs)) {
  
  # get the PDG we are plotting
  pdg <- pdgs[[s]]
  
  # get its first Sanger genotype
  sang <- as.character(unique(subset(ied, sample==pdg, SV)))
  
  if (is.na(sang)) { # if its first Sanger genotype is NA, we already plotted the reference reads above
    # so do not repeat the plot but instead plot an empty plot to fill the list
    misplot <- ggplot() + 
      theme_void() +
      labs(title=pdg)
    misPlots2[[s]] <- misplot # add it to the list
    
  # if first Sanger genotype is not NA, then we should plot the reference reads now
    # (or the 1 OPRD reads for the compound homozygous; will be in SV2 column)
  } else if(!is.na(sang)) {
    
    sang2 <- as.character(unique(subset(ied, sample==pdg, SV2))) # look at SV2 column, the only one not NA is the compound homozygous
    
    if(is.na(sang2)) {
      # then it means we should plot the reference reads
      sang2 <- 'reference'
      
      tmp <- subset(ied, sample==pdg & othercategory==sang2) # subset of ied for plot
      misplot <- ggplot(tmp, aes(reference_nummis, after_stat(density))) + # mismatches plot
        geom_histogram(binwidth=1) +
        theme_minimal() +
        coord_cartesian(xlim=c(0,100), ylim=c(0, 0.3)) +
        labs(title=pdg,
             subtitle=sang2)
      misPlots2[[s]] <- misplot # add it to the list
      
      # if not NA, i.e. we are dealing with the compound homozygous; so it will be 1 OPRD
    } else if (!is.na(sang2)) {
      if(sang2!='OPRD1') stop('\t \t \t \t >>> Something unexpected with the compound homozygous sample \n') # check that is correct
      
      tmp <- subset(ied, sample==pdg & oprd==sang2) # subset of ied for plot
      misplot <- ggplot(tmp, aes(oprd_nummis, after_stat(density))) + # mismatches plot
        geom_histogram(binwidth=1) +
        theme_minimal() +
        coord_cartesian(xlim=c(0,100), ylim=c(0, 0.3)) +
        labs(title=pdg,
             subtitle=sang2)
      misPlots2[[s]] <- misplot # add it to the list
    }
    
  }
}

# interleave the two lists
orderindex <- order(c(seq_along(misPlots1), seq_along(misPlots2))) # gives 1, 26, 2, 27, etc; so can concatenate both lists then order them
misPlots <- (c(misPlots1,misPlots2))[orderindex]

# arrange the interleaved list in a grid
misGrid <- ggarrange(plotlist=misPlots, nrow=length(pdgs), ncol=2) # should produce each row = one sample
ggsave(filename=here('needleSam', 'mismatches_vsconsensus.pdf'), misGrid, height=2000, width=100, units='mm', limitsize=FALSE)

# there is some effect of length on number of mismatches (expectedly); i.e. longer reads (eg. OPRI8) = more mismatches
# I think need to have a flexible threshold based on % error by nucleotide

# below: look at effect of length on number of mismatches

# for all 'expected genotype' reads, plot number of mismatches vs read length
  # add a read length column
ied$read_length <- nchar(ied$seq)

# subset we are interested in = all 'true' reads
  # = all reference reads (except for compound homozygous)
  # + all reads that match Sanger genotype (! 2 Sanger genotypes for compound homozygous)

# will add a 'true read' column and then subset based on that
ied$trueread <- NA

  # all reference reads, except for PDG 46345
ied[which(ied$sample != 'pdg46345' & ied$othercategory == 'reference'), 'trueread'] <- TRUE

  # all reads which have an OPRI that match their SV
ied[which(ied$opri == ied$SV), 'trueread'] <- TRUE

  # all reads which have an OPRI that match their SV2 (will only apply to PDG 46345)
ied[which(ied$opri == ied$SV2), 'trueread'] <- TRUE

  # all reads which have an OPRD that match their SV
ied[which(ied$oprd == ied$SV), 'trueread'] <- TRUE

  # all reads which have an OPRD that match their SV2 (will only apply to PDG 46345)
ied[which(ied$oprd == ied$SV2), 'trueread'] <- TRUE

  # total number of 'true reads'
sum(ied$trueread, na.rm=TRUE)
sum(ied$trueread, na.rm=TRUE) / nrow(ied)
  
  # will put FALSE for all unassigned reads
ied[which(ied$othercategory == 'unassigned'), 'trueread'] <- FALSE

  # how many NA left
sum((is.na(ied$trueread)))

# looks correct, essentially all reads left NA have a potential somatic mutation
  # mismatch threshold will say if convincing or not, though

  # put FALSE for these reads
ied[which(is.na(ied$trueread)), 'trueread'] <- FALSE

  # check no more NA in trueread
sum(is.na(ied$trueread))


# subset ied to keep only true reads
iedTrue <- subset(ied,
                  trueread)


# few checks
# there should not be any unassigned reads
subset(iedTrue, othercategory=='unassigned')
# there should not be any reads from pdg46345 that is reference
subset(iedTrue, sample=='pdg46345' & othercategory=='reference')
# there should not be any eg. 3 OPRI (as no sample has that as its Sanger genotype)
subset(iedTrue, opri=='OPRI3')
subset(iedTrue, opri=='OPRD3')

# looks good
  # ! number of mismatches are spread out between 3 columns
  # need to collapse both columns in one
  # first; check there is no case where both OPRD and OPRI are not NA
subset(iedTrue, !is.na(opri_nummis) & !is.na(oprd_nummis)) # check okay to delete... If only a few and one of opri_nummis or oprd_nummis is high
inds2del <- which(!is.na(iedTrue$opri_nummis) & !is.na(iedTrue$oprd_nummis))
iedTrue <- iedTrue[ - inds2del , ]


# quick sanity checks
subset(iedTrue, !is.na(opri) & is.na(opri_nummis))
subset(iedTrue, is.na(opri) & !is.na(opri_nummis))

subset(iedTrue, !is.na(oprd) & is.na(oprd_nummis))
subset(iedTrue, is.na(oprd) & !is.na(oprd_nummis))

identical(which(is.na(iedTrue$opri)) , which(is.na(iedTrue$opri_nummis)))
identical(which(is.na(iedTrue$oprd)) , which(is.na(iedTrue$oprd_nummis)))

# also, wherever read is reference, it is NA for both opri/oprd
subset(iedTrue, !is.na(reference_nummis) & (!is.na(opri) | !is.na(oprd)))

# or in summary, for each row there should only be one of OPRI/OPRD/othercategory that is not NA
checkNa <- apply(iedTrue, 1, function(re) {
  totna <- as.numeric(!is.na(re['opri_nummis'])) + as.numeric(!is.na(re['oprd_nummis'])) + as.numeric(!is.na(re['reference_nummis']))
  return(totna)
})
# should always be 1
which(checkNa != 1)
# ok so we can collapse into one column

# now to collapse both columns in one:
iedTrue$true_nummis <- NA # preallocate the column

# fill in the OPRI
iedTrue$true_nummis [which(!is.na(iedTrue$opri_nummis))] <- iedTrue$opri_nummis[which(!is.na(iedTrue$opri_nummis))]
iedTrue$true_nummis [which(!is.na(iedTrue$oprd_nummis))] <- iedTrue$oprd_nummis[which(!is.na(iedTrue$oprd_nummis))]
iedTrue$true_nummis [which(!is.na(iedTrue$reference_nummis))] <- iedTrue$reference_nummis[which(!is.na(iedTrue$reference_nummis))]
  # left which() gives the indices to fill in / right which() gets all the non-NA OPRI

# there should not be any NA left
sum(is.na(iedTrue$true_nummis)) # ok

# let's put this column back into ied, will make life easier
  # can simply left_join by read ID, it is safe as no duplicated read
  # check to be safe
sum(duplicated(ied$read))
tmp <- select(iedTrue, read, true_nummis) # keep only read/true_nummis columns
ied <- left_join(x=ied, y=tmp, by='read')

# easy check: all trueread TRUE should have a true_nummis
subset(ied,
       trueread & is.na(true_nummis))
  # few reads, which are the one which have both an OPRI and an OPRD (see above), it is fine we can leave them there

# or in other words, all FALSE trueread should be NA true_nummis
subset(ied,
       !trueread & !is.na(true_nummis))
# ok

# plot read length vs number of mismatches
  # can do from ied now; all not 'true read' will have NA for true_nummis and will not be plotted
lengthVsMis <- ggplot(ied, aes(x=read_length, y=true_nummis)) +
  geom_point() +
  theme_minimal() +
  coord_cartesian(xlim=c(0, 350))
lengthVsMis
ggsave(filename=here('needleSam', 'lengthVsMis.pdf'), lengthVsMis, width=300, height=300, units='mm')

# for any read length, most reads are with low number of mismatches which is good news
# however, there is a clear effect of read length as expected
  # it increases the variability of number of mismatches & the number of reads with many mismatches

# perhaps clearer if plotting average mismatch for each read length
meanMisbyLength <- ied %>%
  group_by(read_length) %>%
  summarise(meanMis=mean(true_nummis, na.rm=TRUE))

meanMisbyLengthPlot <- ggplot(meanMisbyLength, aes(x=read_length, y=meanMis)) +
  geom_point() +
  theme_minimal()
meanMisbyLengthPlot
ggsave(filename=here('needleSam', 'meanMisbyLength.pdf'), meanMisbyLengthPlot, width=300, height=300, units='mm')
# interestingly it makes like a chainsaw pattern
# I think it is simpl because reads tend to have lowest number of mismatches when their length match the target
  # eg. reference = 123 bp, reads that are exactly 123 bp are more likely to have low number of mismatches
  # a reference read that is eg. 118 bp will necessarily have at least 6 mismatches
# correlation read_length vs number of mismatches
cor(ied$read_length, ied$true_nummis, use='complete.obs')

# Normalise number of mismatches by doing % mismatch by nucleotide
  # compute % mismatches by nucleotide
    # simply need to divide mismatches by read_length
ied$true_misbp <- ied$true_nummis / ied$read_length # true read, mismatch by bp

# how does it look
  # mean / sd will give us an idea of the appropriate threshold
mismatchByNucleotide <- ggplot(ied, aes(true_misbp, after_stat(density))) +
  geom_histogram(binwidth=0.01) +
  geom_vline(xintercept=mean(ied$true_misbp, na.rm=TRUE), linetype=1) +
  geom_vline(xintercept=mean(ied$true_misbp, na.rm=TRUE) + sd(ied$true_misbp, na.rm=TRUE), linetype=2) +
  theme_minimal()
mismatchByNucleotide
ggsave(filename=here('needleSam', 'mismatchByNucleotide.pdf'), mismatchByNucleotide, width=100, height=100, units='mm')

# >> I am happy with mean + 1*sd as threshold 'mismatch by nucleotide'
thrMis <- mean(ied$true_misbp, na.rm=TRUE) + sd(ied$true_misbp, na.rm=TRUE)
cat('\t \t \t \t >>> Maximum mismatches allowed vs consensus sequence =', round(thrMis*100, 2), '% of the read length \n')
cat('\t \t \t \t \t or in other words,', round(100-thrMis*100, 2), '% of the read has to align correctly to its consensus sequence \n')

# plot how the threshold line for illustration
thrLine <- as.data.frame(cbind(min(ied$read_length):max(ied$read_length), min(ied$read_length):max(ied$read_length)))
colnames(thrLine) <- c('read_length', 'misauthor') # number of mismatches authorised
thrLine$misauthor <- thrLine$misauthor * thrMis
thrLinePlot <- ggplot(thrLine, aes(x=read_length, y=misauthor)) +
  geom_point()
thrLinePlot
ggsave(filename=here('needleSam', 'thrLinePlot.pdf'), thrLinePlot, width=100, height=100, units='mm')

# is there any correlation with read_length? there should not be any anymore
mismatchByNucleotide_vsReadLength <- ggplot(ied, aes(x=read_length, y=true_misbp)) +
  geom_point() +
  theme_minimal() +
  coord_cartesian(xlim=c(0, 350))
mismatchByNucleotide_vsReadLength
ggsave(filename=here('needleSam', 'mismatchByNucleotide_vsReadLength.pdf'), mismatchByNucleotide_vsReadLength, width=100, height=100, units='mm')
# I do not see any anymore

# do the same as above, mean and correlation
meanMisByNucbyLength <- ied %>%
  group_by(read_length) %>%
  summarise(mean_true_misbp=mean(true_misbp, na.rm=TRUE))

meanMisByNucbyLengthPlot <- ggplot(meanMisByNucbyLength, aes(x=read_length, y=mean_true_misbp)) +
  geom_point() +
  theme_minimal()
meanMisByNucbyLengthPlot
ggsave(filename=here('needleSam', 'meanMisByNucbyLength.pdf'), meanMisByNucbyLengthPlot, width=300, height=300, units='mm')

cor(ied$read_length, ied$true_misbp, use='complete.obs')




# 10- Filter by maximum number of mismatches -------------------------------

# thrMis is the maximum number of mismatches by nucleotide
# compute the maximum number of mismatches authorised for each read
  # i.e. read_length * thrMis
ied$authMis <- round(ied$read_length * thrMis) # number of mismatches authorised

# how many cases are there where one read is both OPRI + OPRD?
nrow(subset(ied, !is.na(opri) & !is.na(oprd)))

# now add a TRUE/FALSE to the reads that are below
  # because of above, will add 2 columns, opri_pass and oprd_pass
ied$opri_pass <- ied$opri_nummis < ied$authMis
ied$oprd_pass <- ied$oprd_nummis < ied$authMis

# ! for PDG 46345 (compound homozygous OPRD1/OPRI5), should consider reference reads as possible somatic mutations
# will simply add a column ref_pass for all; same as above
ied$ref_pass <- ied$reference_nummis < ied$authMis

# check it does what I expect to the PDG 46345 reads
subset(ied, sample=='pdg46345' & othercategory=='reference')


# can also compute % mismatch by nucleotide for each read
  # I am finding it less intuitive, but it has the advantage to be more flexible later as can be used with another threshold later
  # same as above, we need to do it in 3 columns
ied$opri_misbp <- ied$opri_nummis / ied$read_length # mismatches by bp for all OPRI reads
ied$oprd_misbp <- ied$oprd_nummis / ied$read_length # mismatches by bp for all OPRD reads
ied$ref_misbp <- ied$reference_nummis / ied$read_length # mismatches by bp for all reference reads




# write the database ------------------------------------------------------
# (as alignments take a while)
# and we will start exploring it in needleSam_exploration.R

# format a bit better ied

# remove metadata information, only keep read-level info
cols2remove <- c('individual', 'cohortID', 'cohort', 'gDNA_tissue', 'prion_disease',
                 'SV', 'SV2', 'codon129', 'otherSNVs', 'flowcell', 
                 'gender', 'birth_date', 'sample_date', 'death_date', 'death_age',
                 'onset_age', 'age_atsample', 'MRC_t0', 'MRC_intersect', 'MRC_slope')
cols2remove <- match(cols2remove, colnames(ied))

ied <- ied[, -cols2remove]

# arrange columns order
colsorder <- c('sample', 'read', 'seq', 'read_length', 'strand', 'haplotype', 'hapqual', 'besthapqual', 'hapscore',
               'swappedOPR', 'ins', 'del', 'opri', 'oprd', 'othercategory', 'opri_nummis', 'oprd_nummis',
               'reference_nummis', 'trueread', 'true_nummis', 'true_misbp', 'authMis',
               'opri_pass', 'oprd_pass', 'ref_pass', 'opri_misbp', 'oprd_misbp', 'ref_misbp')
length(colsorder) == ncol(ied)
colsorder <- match(colsorder, colnames(ied))
ied <- ied[, colsorder]

# write the database
write.csv(ied, file=here('needleSam', 'allOPRreads.csv'), row.names=FALSE)

# >>> move to needleSam_exploration.R