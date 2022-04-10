# generate list of possible OPRs with R consensus sequence

  # we have a common R consensus making use of flexible nucleotide (eg. Y can be C or T)
  # here: generate possible OPRDs, OPRIs sequences making use of this consensus pattern
  # then these sequences will be used to align candidate somatic mutation reads in needleSam.R



# packages ----------------------------------------------------------------

library(ggplot2)
library(here)



# prepare the OPR catalog -------------------------------------------------

# 'Rx' consensus which works for R2, R3, R4, R2a, R2c, R3g and more that are likely possible but have not been observed (to my knowledge)
rxcon <- 'CCYCAKGGYGGYGGYTGGGRDCAR'
nchar(rxcon)

# consensus for R1 (and leaves 'same' nucleotides flexible; although R1 seems to never or rarely be altered)
r1con <- 'CCYCAKGGYGGYGGYGGYTGGGRDCAR'
nchar(r1con)

# essentially always start with R1, then gradually add more Rx

# genotypes will be
ges <- c(sprintf('OPRD%i', 4:1), # OPRD4 to ORD1
         'reference',
         sprintf('OPRI%i', 1:24)) # OPRI1 to OPRI24

# OPR catalog table
colnms_cato <- c('total_R', 'n_Rx', 'genotype', 'OPR_length', 'pattern', 'consensusseq')
cato <- as.data.frame(matrix(nrow=length(ges), ncol=length(colnms_cato))) # catalog
colnames(cato) <- colnms_cato

# fill in OPR catalog
cato$total_R <- 1:29 # total number of R (reference = 5 R)
cato$n_Rx <- cato$total_R - 1 # total number of Rx (= anything except R1) (reference = 4 Rx)
cato$genotype <- ges # genotypes: OPRD4 ...>>>... reference ...>>>... OPRI24

cato$pattern <- sapply(cato$n_Rx, function(nr) {
  paste0('R1/', paste(rep('Rx/', nr), sep='', collapse=''))
  })

cato$consensusseq <- sapply(cato$n_Rx, function(nr) {
  paste0(r1con, paste(rep(rxcon, nr), sep='', collapse=''))
  })

cato$OPR_length <- nchar(cato$consensusseq)

# as a sanity check
ggplot(cato, aes(x=n_Rx, y=OPR_length)) +
  geom_point() +
  geom_hline(yintercept=seq(27, 27+(4+24)*24, 24))
# ok



# write the OPR catalog ---------------------------------------------------

write.csv(cato, here('needleSam', 'OPRConsensusCatalog.csv'), row.names=FALSE)