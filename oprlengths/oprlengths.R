# histogram of OPR lengths

library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)



# import metadata ---------------------------------------------------------

# will mainly use for the colours to use
meta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='samples')



# small functions ---------------------------------------------------------


beforeFirstChar <- function(stri, cha) {
  as.character(sapply(stri, function(st){
    substr(st,
           1,
           which(strsplit(st, '')[[1]] == cha) [length(which(strsplit(st, '')[[1]] == cha))] - 1)
  }))
}
# eg. beforeFirstChar('1906_counts.txt', '_') returns '1906'

# preallocate length vs count matrix --------------------------------------

# minimum length = 0 = 5 OPRD (deletion of the entire OPR region)
# maximum length = 123 + 100 * 3 * 8 = 2523 = 100 OPRI (can crop later)

minl <- 0
maxl <- 2523

# file paths we will loop through
fils <- c( list.files(here('haplotypephasing', 'bam_notag', 'readlengths'), full.names=TRUE) , 
           list.files(here('haplotypephasing', 'bam_haplotagged', 'readlengths'), full.names=TRUE))

# PDG # in the same order
pdgs <- c( list.files(here('haplotypephasing', 'bam_notag', 'readlengths')) , 
           list.files(here('haplotypephasing', 'bam_haplotagged', 'readlengths')))
pdgs <- beforeFirstChar(pdgs, '_')

# number of samples
numsamples <- length(fils)

# preallocate length vs count matrix
opr <- as.data.frame(matrix(nrow=length(minl:maxl), ncol=numsamples+1)) # +1 as first column will be read length

# add column names = sample's PDG#
colnames(opr) <- c('length', pdgs)

# fill in length column
opr$length <- minl:maxl


# fill in counts for each sample ------------------------------------------

# small dummy dataframe which will allow us to pad the counts table
pad <- as.data.frame(matrix(nrow=length(minl:maxl), ncol=1))
colnames(pad) <- 'length'
pad$length <- minl:maxl

for (s in 1:length(fils)) { # for sample s
  
  # read the table
  spl <- read.table(fils[s])
  colnames(spl) <- c('count', 'length')
  
  # pad it with NA, so have one row for every possible length
  splpad <- left_join(pad, spl, by='length')
  
  # now number of rows in sample counts should be same as OPR dataframe
  if (nrow(opr) != nrow(splpad)) stop('\t \t \t \t >>> Error')
  
  # if yes -- can add the counts to the OPR dataframe
    # match PDG# to find which column we need to fill
  opr[,match(pdgs[s], colnames(opr))] <- splpad$count
  
}

# now change the column names
# (numbers as column names creates issues)
colnames(opr)[2:ncol(opr)] <- paste0('pdg', colnames(opr)[2:ncol(opr)])


# trim OPR ----------------------------------------------------------------

# trim OPR to realistic lengths for plot
# minimum boundary = 123 bp (WT OPR) - 48 bp (2 OPRD) - 10 bp (to safely account for minION noise) = 65 bp
# maximum boundary = 123 bp (WT OPR) + 192 bp (8 OPRI) + 10 bp (minION noise) = 325 bp

lowbound <- 65
topbound <- 325

oprt <- opr[which(opr$length==lowbound) : which(opr$length==topbound) ,] # OPR trimmed



# convert to frequencies --------------------------------------------------

# loop through samples, compute total number of reads, divide all counts by that

count2freq <- function(sp) { # for one column of OPRt
  freq <- sp / sum(sp, na.rm=TRUE)
  return(freq)
}

oprtf <- as.data.frame(cbind(oprt$length, # OPR trimmed / frequencies
                             apply(oprt[,2:ncol(oprt)], 2, count2freq)))
colnames(oprtf)[1] <- 'length'



# compute possible lengths of OPR variants --------------------------------
# so can add some marks to plot

ref <- 9 * 1 * 3 + 8 * 4 * 3 # reference = 1 nonapeptide + 4 octapeptide
# OPR I / D are simply +/- full octapeptides, so +/- multiples of 8 * 3 = 24 bp
oct <- 8 * 3

possiblelengths <- c(ref-2*oct,
                     ref-1*oct,
                     ref,
                     ref+1*oct,
                     ref+2*oct,
                     ref+3*oct,
                     ref+4*oct,
                     ref+5*oct,
                     ref+6*oct,
                     ref+7*oct,
                     ref+8*oct,
                     ref+9*oct,
                     ref+10*oct)

# matching labels
oprlabels <- c('2OPRD',
               '1OPRD',
               'reference',
               '1OPRI',
               '2OPRI',
               '3OPRI',
               '4OPRI',
               '5OPRI',
               '6OPRI',
               '7OPRI',
               '8OPRI',
               '9OPRI',
               '10OPRI')


# plot for one sample -----------------------------------------------------

oprTrace <- ggplot(oprtf, aes(x=length, y=pdg1906)) +
  geom_col() +
  geom_vline(xintercept=possiblelengths, linetype=2) +
  theme_minimal() +
  theme(
    #panel.grid.major.x=element_blank(),
    #panel.grid.major.y=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    
  )
oprTrace



# plot for all samples ----------------------------------------------------

oprTracer <- function(samplenum) {
  
  # create subset of OPR to plot
  oprsub <- oprtf[, c(1,samplenum+1)]
  colnames(oprsub) <- c('length', 'freq')
  
  # find the colour we are supposed to use
  pdgcol <- meta$colour[which(meta$sample == pdgs[samplenum])]
  
  oprTrace <- ggplot(oprsub, aes(x=length, y=freq)) +
    geom_col(fill=pdgcol) +
    geom_vline(xintercept=possiblelengths, linetype=2) +
    theme_minimal() +
    theme(
      #panel.grid.major.x=element_blank(),
      #panel.grid.major.y=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank()
      
    ) +
    coord_cartesian(xlim=c(lowbound-20, topbound+20),
                    ylim=c(0, 0.1)) +
    ggtitle(pdgs[samplenum])
  
  return(oprTrace)
  
}

# preallocate list where to put all the plots
oprtraces <- vector(mode='list', length=numsamples)
names(oprtraces) <- pdgs

# fill it
for (s in 1:length(pdgs)) {
  oprtraces[[s]] <- oprTracer(s)
}

# list of plots is now in PDG order (from smallest PDG number to biggest)
# need to order according to order in Meta file
oprtraces <- oprtraces[as.character(meta$sample)]

# arrange in a grid & save
oprGrid <- ggarrange(plotlist=oprtraces, nrow=13, ncol=2)
ggsave(filename=here('oprlengths', 'oprtrace_grid.pdf'), plot=oprGrid, height=1000, width=240, units='mm')



# grid but actual number of reads -----------------------------------------


oprTracerReads <- function(samplenum) {
  
  # create subset of OPR to plot
  oprsub <- oprt[, c(1,samplenum+1)]
  colnames(oprsub) <- c('length', 'freq')
  
  # find the colour we are supposed to use
  pdgcol <- meta$colour[which(meta$sample == pdgs[samplenum])]
  
  oprTrace <- ggplot(oprsub, aes(x=length, y=freq)) +
    geom_col(fill=pdgcol) +
    geom_vline(xintercept=possiblelengths, linetype=2) +
    theme_minimal() +
    theme(
      #panel.grid.major.x=element_blank(),
      #panel.grid.major.y=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank()
      
    ) +
    coord_cartesian(xlim=c(lowbound-20, topbound+20),
                    ylim=c(0, 100)) +
    ggtitle(pdgs[samplenum])
  
  return(oprTrace)
  
}

# preallocate list where to put all the plots
oprreads <- vector(mode='list', length=numsamples)
names(oprreads) <- pdgs

# fill it
for (s in 1:length(pdgs)) {
  oprreads[[s]] <- oprTracerReads(s)
}

# list of plots is now in PDG order (from smallest PDG number to biggest)
# need to order according to order in Meta file
oprreads <- oprreads[as.character(meta$sample)]

# arrange in a grid & save
oprReads <- ggarrange(plotlist=oprreads, nrow=13, ncol=2)
ggsave(filename=here('oprlengths', 'oprtrace_reads.pdf'), plot=oprReads, height=1000, width=240, units='mm')





# prepare overlay plots ---------------------------------------------------

# melt the counts
oprtfm <- oprtf %>%
  pivot_longer(cols=-length,
               names_to='sample',
               values_to='freq')

# fix factors in oprtfm so same order of samples as in meta
sampleorder <- paste0('pdg', as.character(unlist(meta$sample)))

oprtfm$sample <- factor(oprtfm$sample,
                        levels=sampleorder)

# parameters common for all 3 plots
pdfwidth <- 170
pdfheight <- 42


# overlay plot – inherited part1 ------------------------------------------

samples2plot <- paste0('pdg', as.character(unlist(subset(meta, cohort=='inherited' & prion_disease!='control', sample))))
colours2use <- as.character(unlist(subset(meta, cohort=='inherited' & prion_disease!='control' , colour)))


opr_inherited <- ggplot(subset(oprtfm, sample %in% samples2plot),
                        aes(x=length, y=freq, fill=sample)) +
  geom_col(width=1, alpha=0.7, colour=NA) +
  #geom_vline(xintercept=possiblelengths, linetype=2) +
  #geom_vline(xintercept=c(lowbound-10, topbound+10)) +
  scale_fill_manual(values=colours2use) +
  scale_x_continuous(breaks=seq(100, 350, 50), labels=rep('', length(seq(100, 350, 50)))) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = 'white', colour=NA),
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    legend.position='none'
  ) +
  coord_cartesian(xlim=c(lowbound-10, topbound+10),
                  ylim=c(0, 0.2))

opr_inherited

ggsave(filename=here('oprlengths', 'opr_inherited.pdf'), plot=opr_inherited, height=pdfheight, width=pdfwidth, units='mm')



# overlay plot – sCJD -----------------------------------------------------

samples2plot <- paste0('pdg', as.character(unlist(subset(meta, cohort=='sporadic', sample))))
colours2use <- as.character(unlist(subset(meta, cohort=='sporadic', colour)))


opr_sporadic <- ggplot(subset(oprtfm, sample %in% samples2plot),
                        aes(x=length, y=freq, fill=sample)) +
  geom_col(width=1, alpha=0.7, colour=NA) +
  #geom_vline(xintercept=possiblelengths, linetype=2) +
  #geom_vline(xintercept=c(lowbound-10, topbound+10)) +
  scale_fill_manual(values=colours2use) +
  scale_x_continuous(breaks=seq(100, 350, 50), labels=rep('', length(seq(100, 350, 50)))) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = 'white', colour=NA),
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    legend.position='none'
  ) +
  coord_cartesian(xlim=c(lowbound-10, topbound+10),
                  ylim=c(0, 0.2))

opr_sporadic

ggsave(filename=here('oprlengths', 'opr_sporadic.pdf'), plot=opr_sporadic, height=pdfheight, width=pdfwidth, units='mm')



# overlay plot – inherited x tissue ---------------------------------------

samples2plot <- paste0('pdg', as.character(unlist(subset(meta, cohort=='inherited_tissuecomparison', sample))))
colours2use <- as.character(unlist(subset(meta, cohort=='inherited_tissuecomparison', colour)))

opr_tissue <- ggplot(subset(oprtfm, sample %in% samples2plot),
                   aes(x=length, y=freq, fill=sample)) +
  geom_col(width=1, alpha=0.7, colour=NA) +
  #geom_vline(xintercept=possiblelengths, linetype=2) +
  #geom_vline(xintercept=c(lowbound-10, topbound+10)) +
  scale_fill_manual(values=colours2use) +
  scale_x_continuous(breaks=seq(100, 350, 50), labels=rep('', length(seq(100, 350, 50)))) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = 'white', colour=NA),
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    legend.position='none'
  ) +
  coord_cartesian(xlim=c(lowbound-10, topbound+10),
                  ylim=c(0, 0.2))

opr_tissue

ggsave(filename=here('oprlengths', 'opr_tissue.pdf'), plot=opr_tissue, height=pdfheight, width=pdfwidth, units='mm')