# takes a sequence, plots GC%

library(here)
library(data.table)
library(ggplot2)

# GC percent --------------------------------------------------------------
# rolling sum in window

GCpercent <- function(sequence, window) {
  
  seqsplit <- strsplit(sequence, split='')[[1]]
  
  gc1 <- as.numeric(sapply(seqsplit, function(nu){
    if(nu=='A' | nu=='T') return(0)
    if(nu=='C' | nu=='G') return(1)
  }))
  
  gcroll <- frollsum(gc1, window,
                     fill=NA, align='right')
  
  gcpe <- gcroll/window
  
  return(gcpe)
  
}

# GC plot -----------------------------------------------------------------

GCplot <-  function(GCpercent) {
  
  gcdf <- cbind(1:length(GCpercent), as.data.frame(GCpercent))
  colnames(gcdf) <- c('pos', 'gcp')
  
  gcgg <- ggplot(gcdf, aes(x=pos, y=gcp)) +
    geom_line(colour='#5d5e5d', size=0.5) +
    coord_cartesian(ylim=c(0,1)) +
    theme_minimal() +
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      panel.grid.major.x=element_blank()
    )
    
  
  return(gcgg)
  
}


# GC plot – barplot -------------------------------------------------------

GCplotbar <-  function(GCpercent) {
  
  gcdf <- cbind(1:length(GCpercent), as.data.frame(GCpercent))
  colnames(gcdf) <- c('pos', 'gcp')
  
  gcgg <- ggplot(gcdf, aes(x=pos, y=gcp)) +
    geom_col(colour='#5d5e5d', size=0.5) +
    coord_cartesian(ylim=c(0,1)) +
    theme_minimal() +
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      panel.grid.major.x=element_blank()
    )
  
  
  return(gcgg)
  
}


# to use ------------------------------------------------------------------

# 1- import sequence in fasta
fasta <- here('GCplot', 'prnp_window1.fa')
seq <- as.character(read.table(fasta, skip=1))

# 2- run GCpercent & GCplot, with desired window
gcp1 <- GCplotbar(GCpercent(seq, 100))

ggsave((here('GCplot', 'gc_window1.pdf')), plot=gcp1, width=65, height=18, units='mm')

# second half
# 1- import sequence in fasta
fasta <- here('GCplot', 'prnp_window2.fa')
seq <- as.character(read.table(fasta, skip=1))

# 2- run GCpercent & GCplot, with desired window
gcp2 <- GCplotbar(GCpercent(seq, 100))

ggsave((here('GCplot', 'gc_window2.pdf')), plot=gcp2, width=52, height=18, units='mm')

# window3
# 1- import sequence in fasta
fasta <- here('GCplot', 'prnp_window3.fa')
seq <- as.character(read.table(fasta, skip=1))

# 2- run GCpercent & GCplot, with desired window
gcp3 <- GCplotbar(GCpercent(seq, 100))

ggsave(here('GCplot', 'gc_window3.pdf'), plot=gcp3, width=32, height=18, units='mm')
