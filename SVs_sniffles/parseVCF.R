library(here)
library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)




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




# import metadata ---------------------------------------------------------

# import metadata
meta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='samples')

# format meta
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




# import promoter VCFs ----------------------------------------------------

prosv <- as.data.frame(matrix(ncol=11, nrow=0))


# first import regulatory region amplicon calls
  # files finish by _promoterSnif.vcf
folder <- here('SVs_sniffles', 'vcf')

# find indices in list.files of the promoter (regulatory region) VCFs
proinds <- which(unlist(lapply(strsplit(list.files(folder), '_'), function(f){f[2]})) == 'promoterSnif.vcf') # promoter indices
profp <- list.files(folder, full.names=TRUE)[proinds] # promoter full paths

for (v in 1:length(profp)) {
  
  path <- profp[v]
  pdg <- strsplit(list.files(folder)[proinds][v], '_')[[1]][1]
  cat('\t \t \t \t >>> Sample #', pdg, '\n')
  
  ok <- try(read.table(path))
  
  if (!inherits(ok, 'try-error')) { # if no error, read the file; if error, skip to the next one
    tmp <- read.table((path))
    tmp <- cbind(rep(pdg, nrow(tmp)), tmp) # add sample where calls came from
    prosv <- rbind(prosv, tmp)
  }
  
}
# add source as first column
prosv <- cbind('promoter', prosv)




# import genebody VCFs ----------------------------------------------------

gbsv <- as.data.frame(matrix(ncol=11, nrow=0))


# second import genebody amplicon calls
  # files finish by _bodySnif.vcf

# find indices in list.files of the genebody VCFs
gbinds <- which(unlist(lapply(strsplit(list.files(folder), '_'), function(f){f[2]})) == 'bodySnif.vcf') # genebody indices
gbfp <- list.files(folder, full.names=TRUE)[gbinds] # genebody full paths

for (v in 1:length(gbfp)) {
  
  path <- gbfp[v]
  pdg <- strsplit(list.files(folder)[gbinds][v], '_')[[1]][1]
  cat('\t \t \t \t >>> Sample #', pdg, '\n')
  
  ok <- try(read.table(path))
  
  if (!inherits(ok, 'try-error')) { # if no error, read the file; if error, skip to the next one
    tmp <- read.table((path))
    tmp <- cbind(rep(pdg, nrow(tmp)), tmp) # add sample where calls came from
    gbsv <- rbind(gbsv, tmp)
  }
  
}
# add source as first column
gbsv <- cbind('genebody', gbsv)


# stick together both
colnames(prosv) <- sprintf('V%i', 1:ncol(prosv))
colnames(gbsv) <- sprintf('V%i', 1:ncol(gbsv))

sv <- rbind(gbsv, prosv)



# format ------------------------------------------------------------------

# column 10 is everything crammed; split in multiple columns
colsplitall <- as.data.frame(t(apply(sv, 1, function(row){
  tmp <- unlist(strsplit(row[10], ';'))
  # PRECISE/IMPRECISE is different
  # as.data.frame(t(c(tmp[1], as.character(lapply(sapply(tmp[2:length(tmp)], strsplit, '='), function(i) i[2])))))
  colsplit <- c(tmp[1], as.character(lapply(sapply(tmp[2:length(tmp)], strsplit, '='), function(i) i[2])))
  return(colsplit)
})))

sv <- cbind(sv,colsplitall)

colnames(sv) <- c('source',
                  'sample',
                  'chr',
                  'pos',
                  'id',
                  'ref',
                  'alt',
                  'qual',
                  'filter',
                  'info',
                  'format_nexttag',
                  'geno_readsRef_readsAlt',
                  'breakpoints_confidence',
                  'svmethod',
                  'chr2',
                  'end',
                  'std_quant_start',
                  'std_quant_stop',
                  'kurtosis_quant_start',
                  'kurtosis_quant_stop',
                  'svtype',
                  'subtype',
                  'svlen',
                  'strands',
                  'strands2',
                  're',
                  'ref_strand',
                  'strandbias_pval',
                  'af')

# now delete info column
sv$info <- NULL

# allele frequency column is stored as character, should be numeric
sv$af <- as.numeric(sv$af)
# same with strandbias_pval
sv$strandbias_pval <- as.numeric(sv$strandbias_pval)




# take only calls in sequenced window -------------------------------------

cat('\t \t \t \t >>> Total number of calls:', nrow(sv), 'calls \n')

# take only calls in PRNP genomic region
# earliest possible position
left <- 4685060
# latest possible position
right <- 4701756

svp <- subset(sv, pos>left & pos<right)
cat('\t \t \t \t >>> After filter = only PRNP genomic window:', nrow(svp), 'calls \n')




# look at calls in regulatory region amplicon -----------------------------
# i.e. 5' amplicon

prosv <- subset(svp, source=='promoter')

# filter them
cat('\t \t \t \t >>> Before filtering:', nrow(prosv), 'calls \n')

# filter based on allele frequency
# minimum = 10%
prosvf <- subset(prosv, af>0.10)

cat('\t \t \t \t >>> After filter = Allele Frequency > 0.10:', nrow(prosvf), 'calls \n')




# look at calls in genebody amplicon --------------------------------------

gbsv <- subset(svp, source=='genebody')

# filter them
cat('\t \t \t \t >>> Before filtering:', nrow(gbsv), 'calls \n')

# filter based on allele frequency
# minimum = 10%
gbsvf <- subset(gbsv, af>0.10)

cat('\t \t \t \t >>> After filter = Allele Frequency > 0.10:', nrow(gbsvf), 'calls \n')


# For each mutated OPR sample, is there at least one call?
pdgopr <- meta[!is.na(meta$SV), 'sample'] # PDG which have an OPR mutation
sum(pdgopr %in% unique(gbsvf$sample)) == length(pdgopr)


# Note; StrandBias flag is slightly problematic
# if look before filtering, many unique calls get called multiple times
# and as one read can only belong to one call, each call gets a low allele frequency
# it also creates imbalances in strand




# nice export -------------------------------------------------------------

# pool filtered calls
svf <- rbind(prosvf, gbsvf)

# add some sample information
meta$sample <- as.character(meta$sample)
sve <- meta %>% # SV Export
  select(sample, gDNA_tissue, prion_disease, SV) %>%
  right_join(., svf, by='sample')

sve$svlen <- as.numeric(sve$svlen)

# ! svlen does not always match the length of the inserted/deleted sequence
# see https://github.com/fritzsedlazeck/Sniffles/issues/267
# should look at svlen as this is calculated from multiple reads

# will add manually in Excel error vs expected SV length

# export and copy-paste to supplementary
write.xlsx(sve, file=here('SVs_sniffles', 'snifflescallsfilter.xlsx'))

# also export calls before filtering
write.xlsx(svp, file=here('SVs_sniffles', 'snifflescallsall.xlsx'))