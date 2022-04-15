# work on SNVs called by nanopolish



library(openxlsx)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggbeeswarm)




# colours -----------------------------------------------------------------

# cohorts
  # inherited = yellow
  # sCJD = khaki
  # control = grey
cohcols <- c('#fcb505', '#78ac63', '#B3BDC4')




# import ------------------------------------------------------------------

# filtered SNVs
var <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='SNVs_filtered')

# unfiltered SNVs
varun <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='SNVs_unfiltered')

# samples metadata
meta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='samples')
meta$sample <- as.character(meta$sample)

# SNV allele frequency
af <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='SNVs_allelefrequency')




# how many total ----------------------------------------------------------

nrow(var)
# how many unique
length(unique(var$Existing_variation))




# how many outside of protein-coding --------------------------------------

# column Consequence include missense_variant or else
# but careful; if SNV is at third base of codon, could be in coding region but not missense

# all SNVs in coding region should have a number in column Protein_position
# (is not, it is '-', which will produce NA when converted to as.numeric, hence the Warning)
sum(is.na(as.numeric(var$Protein_position))) # number OUTSIDE coding region
sum(!is.na(as.numeric(var$Protein_position))) # number IN coding region

# check we get the same when searching for the '-'
sum(var$Protein_position=='-') # number OUTSIDE coding region
sum(var$Protein_position!='-') # number IN coding region

# check we are getting the same numbers using positions
# coding region = chr20:4,699,221--4,699,982
length(which(var$POS < 4699221 | var$POS > 4699982)) # number OUTSIDE coding region
length(which(var$POS >= 4699221 & var$POS <= 4699982)) # number IN coding region
# same

# add a tag to these SNVs for later
var$coding <- NA
var$coding[which(var$POS < 4699221 | var$POS > 4699982)] <- FALSE
var$coding[which(var$POS >= 4699221 & var$POS <= 4699982)] <- TRUE
# check no more NA
sum(is.na(var$coding))

# how many unique?

# all
# create composite chr_pos_ref_alt column
var$com <- paste(var$CHROM, var$POS, var$REF, var$ALT, sep='_')
length(unique(var$com))

# coding
nrow(unique(subset(var, coding, com)))

# non-coding
nrow(unique(subset(var, !coding, com)))




# how many in regulatory region -------------------------------------------

# regulatory region (according to Ensembl) = chr20:4,685,060–4,688,047
subset(var, POS>4685060 & POS<4688047)
nrow(subset(var, POS>4685060 & POS<4688047))
# how many unique?
nrow(unique(subset(var, POS>4685060 & POS<4688047, Existing_variation)))
# check same with com column
nrow(unique(subset(var, POS>4685060 & POS<4688047, com)))

# how many affect transcription binding site in regulatory region
nrow(subset(var, POS>4685060 & POS<4688047 & Consequence=='TF_binding_site_variant'))
# how many unique?
nrow(unique(subset(var, POS>4685060 & POS<4688047 & Consequence=='TF_binding_site_variant', Existing_variation)))
# check same with com column
nrow(unique(subset(var, POS>4685060 & POS<4688047 & Consequence=='TF_binding_site_variant', com)))



# how many non-coding per sample ------------------------------------------

nsnv_noncoding_persample <- var %>%
  subset(!coding) %>%
  group_by(SAMPLE) %>%
  summarise(nsnsv=n())




# versus Sanger -----------------------------------------------------------

# gather the Sanger SNV calls
snv129 <- meta[which(meta$codon129 == 'M129V'), c('individual', 'sample', 'codon129')]
colnames(snv129) <- c('invididual', 'sample', 'snv')
snv0 <- meta[which(meta$otherSNVs != '-'), c('individual', 'sample', 'otherSNVs')]
colnames(snv0) <- c('invididual', 'sample', 'snv')
san <- rbind(snv129, snv0)

# add the nanopolish SNV calls in coding
nan <- subset(var, coding) %>%
  select(SAMPLE, Protein_position, Amino_acids)
colnames(nan)[1] <- 'sample'
nan$sample <- as.character(nan$sample)

# ! left_join, keep all Sanger calls
sana <- left_join(x=san, y=nan, by='sample') # Sanger vs Nanopolish
# ok

# also do vice-versa so keep all nanopolish coding calls and add the Sanger
# (would pick up if anything extra from nanopolish, potentially false positive)
nasa <- left_join(x=nan, y=san, by='sample') # Nanopolish vs Sanger
# ok
# extra is PDG52331, codon129 was not tested by Sanger




# replicates when same invididual? ----------------------------------------

# i.e. matched blood vs brain samples

# add all of meta information to var (will include individual)
var$SAMPLE <- as.character(var$SAMPLE)
varm <- left_join(var, meta, by=c('SAMPLE'='sample')) # var with metadata

# check no individual is NA
varm[is.na(varm$individual),]

# individuals in tissue comparison cohort
idvs <- as.numeric(unlist(unique(subset(meta, cohort=='inherited_tissuecomparison', individual))))

# will (semi) manually check the SNVs calls for each
# tally by composite ID; if all have 2, ok
subset(varm, individual==idvs[1]) %>%
  group_by(com) %>%
  tally()

subset(varm, individual==idvs[2]) %>%
  group_by(com) %>%
  tally()

subset(varm, individual==idvs[3]) %>%
  group_by(com) %>%
  tally()

subset(varm, individual==idvs[4]) %>%
  group_by(com) %>%
  tally()

# SNV below detected only in blood sample
varm[which(varm$com=='chr20_4690801_T_C'),]
# below, it is actually present, but only in unfiltered
# as it is above SOR threshold
subset(varun, SAMPLE==53689 & POS==4690801)



# strand bias -------------------------------------------------------------

# from above, we know coding variants are all = Sanger
# (well, except PDG52331 which was not tested but will safe to assume it is M129V)

# what is max strand bias of true positive calls = StrandOddsRatio threshold
max(subset(varm, coding, SOR))

# we set threshold at 1
sorthr <- 1

# check we are obtaining same results
nrow(subset(varun, SOR<sorthr)) == nrow(var)

poscheck <- sort(subset(varun, SOR<sorthr)$POS)
identical(poscheck, sort(var$POS))


# number of individuals each unique SNV is found in -----------------------

# better than counting samples, as same individual can have multiple samples

# add a unique composite column: sample_chr_pos_ref_alt
varm$ucom <- paste(varm$SAMPLE, varm$CHROM, varm$POS, varm$REF, varm$ALT, sep='_')
sum(duplicated(varm$ucom)) # no duplicated

# about rs numbers

# some have multiple codes, just keep the rs
# it is the first one before the coma
varm$rs <- as.character(lapply(strsplit(varm$Existing_variation, ','), function(ex) {ex[1]}))

# do all SNVs here have an rs?
which(varm$Existing_variation == '-')
which(is.na(varm$Existing_variation))

# same number of unique rs as number of unique com ID?
length(unique(varm$rs)) == length(unique(varm$com))

# to do so, skip the calls from brain samples
varmi <- varm[-which(varm$DNA.from.tissue=='brain'), ] # single-nucleotide VARiants + Metadata + Individual

# now each individual is represented by a single sample
# group by SNV and count how many times called, should be = number of individuals who carry it
# (will also group by cohort to do stack, but same idea)
nindiv <- varmi %>%
  group_by(com, rs, prion_disease) %>%
  tally(name='nidvs') %>%
  arrange(desc(nidvs)) # in order of most individuals to least

# also need to group without prion_disease to get right descending order of SNVs
nin <- varmi %>%
  group_by(com, rs) %>%
  tally(name='nidvs') %>%
  arrange(desc(nidvs)) # in order of most individuals to least

# set the SNVs order as when grouped without looking at cohort (so most individuals >> least)
nindiv$rs <- factor(nindiv$rs, levels=nin$rs)

# cohort in right order
nindiv$prion_disease <- factor(nindiv$prion_disease, levels=c('inherited', 'sCJD', 'control'))

# y axis title
ytit <- paste0('number of individuals (of ', length(unique(meta$individual)), ')') # makes sure I do not do a mistake in the total number!

SNVvsnInd <- ggplot(nindiv, aes(x=rs, y=nidvs, fill=prion_disease)) +
  geom_col(position='stack', width=0.7) +
  scale_fill_manual(values=cohcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    axis.title.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text.x=element_text(size=7, angle=90, hjust=0, vjust=0.40, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  ylab(ytit) +
  scale_y_continuous(breaks=0:10)
SNVvsnInd
ggsave(here('SNVs', 'SNVvsnIndividuals.pdf'), plot=SNVvsnInd, width=120, height=80, units='mm')




# number of SNVs in each individual ---------------------------------------

# group by individual and count
# during grouping, keep all info that is individual-level
nsnvi <- varmi %>%
  group_by(individual, cohort, gDNA_tissue, prion_disease,
           SV, codon129, otherSNVs, gender) %>%
  tally(name='nsnvs')

# ! any individual in meta but not in nsnvi, means it had 0 SNV, should add them
unique(meta$individual[!meta$individual %in% nsnvi$individual])

toadd <- meta[!meta$individual %in% nsnvi$individual , ]
# skip any brain samples
toadd <- subset(toadd, gDNA_tissue=='blood')

# select the same column as above and add nsnvs = 0
toadd <- toadd %>%
  select(individual, cohort, gDNA_tissue, prion_disease,
         SV, codon129, otherSNVs, gender) %>%
  mutate(nsnvs=0)

# add to nsnvi
nsnvi <- rbind(nsnvi, toadd)


# now there should be same number of rows as number of individuals
length(unique(meta$individual)) == nrow(nsnvi)

# order individuals by number of SNVs
nsnvi <- nsnvi %>%
  arrange(desc(nsnvs))

nsnvi$individual <- factor(nsnvi$individual, levels=nsnvi$individual)

# now histogram SNV vs number of individuals who carry it
# (also colour by cohort)

# temporarily increase a little bit the individuals with count 0 so we see the colour
nsnvi2 <- nsnvi
nsnvi2[which(nsnvi2$nsnvs==0), 'nsnvs'] <- nsnvi2[which(nsnvi2$nsnvs==0) , 'nsnvs'] + 0.1

# order of cohort
nsnvi2$prion_disease <- factor(nsnvi2$prion_disease, levels=c('inherited', 'sCJD', 'control'))

individualvsnSNVs <- ggplot(nsnvi2, aes(x=individual, y=nsnvs, fill=prion_disease)) +
  geom_col(width=0.8) +
  scale_fill_manual(values=cohcols) +
  theme_minimal() +
  xlab('individual') + ylab('number of SNVs') +
  theme(
    panel.grid.minor.y=element_blank(),
    legend.position='none')
individualvsnSNVs
ggsave(here('SNVs', 'individualvsnSNVs.pdf'), plot=individualvsnSNVs, width=200, height=100, units='mm')




# any cohort with more SNVs? ----------------------------------------------

# use nsnvi abvove

nsnvi$prion_disease <- factor(nsnvi$prion_disease, levels=c('inherited', 'sCJD', 'control'))

# plot
nsnvbyCoh <- ggplot(nsnvi, aes(x=prion_disease, y=nsnvs, colour=prion_disease)) +
  geom_quasirandom(width=0.2, size=2) +
  scale_colour_manual(values=cohcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(angle=90, size=7, margin=margin(t=0, r=0, b=0, l=0), vjust=0.5),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  ylab('number of SNVs') +
  scale_x_discrete(labels=c('inherited', 'sporadic', 'control'))
nsnvbyCoh

ggsave(here('SNVs', 'SNVcountbyIndividualbyCohort.pdf'), plot=nsnvbyCoh, width=35, height=80, units='mm')

# test
summary(aov(nsnvs ~ prion_disease, data=nsnvi))




# same but only SNVs in regulatory region amplicon ------------------------

# i.e. before reverse primer of amplicon1
varmireg <- subset(varmi, POS<4688028)

# group by individual and count
# during grouping, keep all info that is individual-level
nsnvireg <- varmireg %>%
  group_by(individual, cohort, gDNA_tissue, prion_disease,
           SV, codon129, otherSNVs, gender) %>%
  tally(name='nsnvs')

# ! any individual in meta but not in nsnvireg, means it had 0 SNV in promoter, should add them
unique(meta$individual[!meta$individual %in% nsnvireg$individual])

toadd <- meta[!meta$individual %in% nsnvireg$individual , ]
# skip any brain samples
toadd <- subset(toadd, gDNA_tissue=='blood')

# select the same column as above and add nsnvs = 0
toadd <- toadd %>%
  select(individual, cohort, gDNA_tissue, prion_disease,
         SV, codon129, otherSNVs, gender) %>%
  mutate(nsnvs=0)

# add to nsnvi
nsnvireg <- rbind(nsnvireg, toadd)

# now there should be same number of rows as number of individuals
length(unique(meta$individual)) == nrow(nsnvireg)

# order individuals by number of SNVs
nsnvireg <- nsnvireg %>%
  arrange(desc(nsnvs))

nsnvireg$individual <- factor(nsnvireg$individual, levels=nsnvireg$individual)

# now histogram SNV vs number of individuals who carry it
# (also colour by cohort)

# temporarily increase a little bit the individuals with count 0 so we see the colour
nsnvireg2 <- nsnvireg
nsnvireg2[which(nsnvireg2$nsnvs==0), 'nsnvs'] <- nsnvireg2[which(nsnvireg2$nsnvs==0) , 'nsnvs'] + 0.1

# order of cohorot
nsnvireg2$prion_disease <- factor(nsnvireg2$prion_disease, levels=c('inherited', 'sCJD', 'control'))

individualvsnSNVsReg <- ggplot(nsnvireg2, aes(x=individual, y=nsnvs, fill=prion_disease)) +
  geom_col(width=0.8) +
  scale_fill_manual(values=cohcols) +
  theme_minimal() +
  xlab('individual') + ylab('number of SNVs') +
  theme(
    panel.grid.minor.y=element_blank(),
    legend.position='none')
individualvsnSNVsReg
ggsave(here('SNVs', 'individualvsnSNVsReg.pdf'), plot=individualvsnSNVsReg, width=200, height=100, units='mm')


# any cohort with more SNVs in regulatory region?
# use nsnvireg above
nsnvireg$prion_disease <- factor(nsnvireg$prion_disease, levels=c('inherited', 'sCJD', 'control'))

# plot
nsnvbyCohReg <- ggplot(nsnvireg, aes(x=prion_disease, y=nsnvs, colour=prion_disease)) +
  geom_quasirandom(width=0.2, size=2) +
  scale_colour_manual(values=cohcols) +
  theme_minimal() +
  theme(
    legend.position='none',
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  ylab('number of SNVs') +
  scale_x_discrete(labels=c('inherited', 'sporadic', 'control'))
nsnvbyCohReg

ggsave(here('SNVs', 'SNVcountbyIndividualbyCohortReg.pdf'), plot=nsnvbyCohReg, width=200, height=100, units='mm')

# test
summary(aov(nsnvs ~ prion_disease, data=nsnvireg))




# plot strand bias --------------------------------------------------------

# look at distribution
ggplot(varun, aes(SOR)) +
  geom_histogram(binwidth=0.1)

# add unique composite ID column
varun$ucom <- paste(varun$SAMPLE, varun$CHROM, varun$POS, varun$REF, varun$ALT, sep='_')
sum(duplicated(varun$uid)) # no duplicated

# add coding or not in unfiltered SNVs
varun$coding <- NA
varun$coding[which(varun$POS < 4699221 | varun$POS > 4699982)] <- FALSE
varun$coding[which(varun$POS >= 4699221 & varun$POS <= 4699982)] <- TRUE

varun <- varun[(order(varun$SOR)),] # order by SOR
varun$ucom <- factor(varun$ucom, levels=varun$ucom)

# fix rs numbers as above
# some have multiple codes, just keep the rs
# it is the first one before the coma
varun$rs <- as.character(lapply(strsplit(varun$Existing_variation, ','), function(ex) {ex[1]}))

# add a TRUE/FALSE column whether variation exists (has a RS number or not)
varun$variantHasRS <- sapply(varun$rs, function(v){
  if (v=='-') return(FALSE)
  else return(TRUE)
  })

# see below explanations for arrange()
# Update: after re-analysis, all have an rs number
# leaving code here, but it is not relevant anymore
var_SOR <- ggplot(varun %>%
                    arrange(desc(variantHasRS)), aes(x=ucom, y=SOR, colour=variantHasRS)) +
  geom_hline(yintercept=sorthr, linetype=2, size=1.0, colour='#b8b8b8') +
  geom_point(size=0.7) +
  scale_colour_manual(values=c('#4D4D4D', '#f1876b')) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.ticks.x=element_line(size=0.1),
    legend.position='none',
    axis.title.x=element_text(size=9, margin = margin(t=4, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  xlab('SNV call') + ylab('strand bias (SOR)') +
  scale_x_discrete(expand=c(0.01,0))
var_SOR
ggsave(here('SNVs', 'SOR.pdf'), plot=var_SOR, width=90, height=60, units='mm')


# same, but colour whether coding variant (so confirmed by Sanger, so true positive) or not
# add column 'confirmed by Sanger' or not
# basically, this is = coding SNVs from var

varun$truepos <- FALSE

# here are the SNVs we have in var
vfil <- unique(varm$ucom) # SNVs filtered

# SNV is coding + SNV is in var = true positives
varun[which(varun$ucom %in% vfil & varun$coding) , 'truepos'] <- TRUE

# ! need to plot the Sanger SNVs on top, otherwise will be hidden
# one way to do this is to put the coloured dots at the bottom of the dataframe, so gets drawn last
# this is what the arrange() below does
var_truepos <- ggplot(varun %>%
                        arrange(truepos), aes(x=ucom, y=SOR, colour=truepos)) +
  geom_hline(yintercept=sorthr, linetype=2, size=1.0, colour='#b8b8b8') +
  geom_point(size=0.7) +
  scale_colour_manual(values=c('#4D4D4D', '#f1876b')) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.ticks.x=element_line(size=0.1),
    legend.position='none',
    axis.title.x=element_text(size=9, margin = margin(t=4, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=2, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0))
  ) +
  xlab('SNV call') + ylab('strand bias (SOR)') +
  scale_x_discrete(expand=c(0.01,0))
var_truepos
ggsave(here('SNVs', 'SOR_truepositive.pdf'), plot=var_truepos, width=73, height=60, units='mm')




# number of individuals vs allele frequency -------------------------------

# use nin created above
# i.e. number of individuals who carry each unique SNV

# join allele frequency to nin built above
colnames(af)[1] <- 'rs'
nin <- left_join(af, nin, by='rs')

# check unique SNV
sum(duplicated(nin$rs))

# y axis title
ytit <- paste0('number of individuals (of ', length(unique(meta$individual)), ')') # makes sure I do not do a mistake in the total number!

nvsaf <- ggplot(nin, aes(x=allelefrequency, y=nidvs)) +
  geom_point(size=1.0) +
  theme_minimal() +
  scale_y_continuous(breaks=0:10) +
  theme(
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t=2, r=0, b=0, l=0)),
    axis.title.y=element_text(size=9, margin = margin(t=0, r=0, b=0, l=0)),
    axis.text.x=element_text(size=7, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=0, b=0, l=0)),
  ) +
  ylab(ytit) + xlab('allele frequency')
nvsaf

ggsave(filename=here('SNVs', 'SNVcountindividualsvsAF.pdf'), plot=nvsaf, width=73, height=60, units='mm')

cor(nin$nidvs, nin$allelefrequency, method='pearson')
cor.test(nin$nidvs, nin$allelefrequency, method='pearson')

# there is one vaguely interesting SNV

# low allele frequency but 8 individuals in our panel
intrs <- as.character(unlist(subset(nin, allelefrequency<0.15 & nidvs==8, rs))) # interesting rs
# who has it?
subset(varm, rs==intrs)



# Ns ----------------------------------------------------------------------

# number of individuals in each cohort
nsnvi %>%
  group_by(prion_disease) %>%
  tally()
