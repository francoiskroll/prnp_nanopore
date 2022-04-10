# prepare VCFs for whatshap
  # from curated nanopolish SNVs calls

# packages ----------------------------------------------------------------

library(openxlsx)
library(here)
library(dplyr)
library(tidyr)




# import ------------------------------------------------------------------

# will make VCFs from filtered SNVs

# import filtered SNVs
snv <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='SNVs_filtered')
# will assume all these calls are true positives

# import unfiltered SNVs
snvunf <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='SNVs_unfiltered') # SNVs unfiltered

# import sample's metadata
meta <- read.xlsx(here('AdditionalFile1.xlsx'), sheet='samples')




# do we have enough SNVs for each sample? ---------------------------------
# for haplotype phasing, we need heterozygous variants in the gene-body reads

# remove all promoter SNVs from snv
snv <- snv[-which(snv$POS < 4687751),] # i.e. before end of genebody_Forward primer binding site

# are there any homozygous SNVs? These are not useful for haplotype phasing
subset(snv, Genotype=='1/1') # 0 homozygous SNP

# are all samples represented in SNVs? (i.e. have at least 1 heterozygous SNV that can be used for haplotype phasing)
pdgs <- unique(snv$SAMPLE)

meta[!meta$sample %in% pdgs,] # No; 5 samples have 0 SNV calls
nosv <- sort(meta$sample[!meta$sample %in% pdgs]) # PDG numbers of samples without useful SNVs

# for these samples, any usable SNVs thrown away durinig StrandBias filtering?

# delete SNVs in the promoter, for sure those are not useful here
snvunf <- snvunf[-which(snvunf$POS < 4687751),] # i.e. before end of genebody_Forward primer binding site

subset(snvunf, SAMPLE==nosv[1])
subset(snvunf, SAMPLE==nosv[2])
subset(snvunf, SAMPLE==nosv[3])
subset(snvunf, SAMPLE==nosv[4])
subset(snvunf, SAMPLE==nosv[5])
# all high strandbias
# >> nothing to salvage from the bin in terms of SNVs

# samples without SNVs are
# PDG 53747: 6 OPRI
# PDG 58648: 6 OPRI
# PDG 54890: no SV (sCJD)
# PDG 55050: no SV (sCJD)
# PDG 59060: 6 OPRI

# (59060/53747 are same individual)

# Note: I have tried doing haplotype phasing using a single SV (6 OPRI of PDG 58648), it runs but does not work well
  # whatshap haplotag seems too stringent about identifying the SV:
    # one haplotype gathers all the reads with a relatively 'clean' insertion
    # and all the other reads are considered other haplotype (reference), but I can see in IGV many contain the 6 OPRI
    # it may still be possible to do it more 'manually' in R by looking at the CIGARs, but probably a bit rough

# >> will only use SNV for haplotype phasing
# and keep OPR SVs as a way to check if the haplotype phasing is correct

# for the samples that present in snv:
bysmp <- snv %>%
  group_by(SAMPLE) %>%
  tally()



# prepare VCFs ------------------------------------------------------------

# create 2 folders
dir.create(here('haplotypephasing'))
dir.create(here('haplotypephasing/vcf_unphased')) # will store the vcf before phasing, i.e. the ones for which we have more than 1 SNV
dir.create(here('haplotypephasing/vcf_phased')) # will store the vcf after phasing
  # i.e. the ones for which we have a single SNV (we are writing it 'as if' it went through whatshap phase)
  # + whatshap phase will put the phased version of the VCF above in that folder as well

# loop through samples and prepare each VCF

for (p in 1:nrow(meta)) { # for each sample
  
  # get PDG number
  pdg <- meta$sample[p]
  
  # ! we do not have any SNV to do haplotype phasing for 5 samples (nosv above), so skip these
  if (pdg %in% nosv) {
    next
  } else if (pdg %in% unique(snv$SAMPLE)) { # for the samples that have enough SNV calls
    
    # get the filtered SNV calls for that sample
    vcf <- subset(snv, SAMPLE==pdg)
    
    # need to reconstruct the INFO column
    infocolumn <- as.character(apply(vcf, 1,
                                     function(ro){
                                       paste0('BaseCalledReadsWithVariant=', ro['BaseCalledReadsWithVariant'], ';',
                                              'BaseCalledFraction=', ro['BaseCalledFraction'], ';',
                                              'TotalReads=', ro['TotalReads'], ';',
                                              'AlleleCount=', ro['AlleleCount'], ';',
                                              'SupportFraction=', ro['SupportFraction'], ';',
                                              'SupportFractionByStrand=', ro['SupportFractionByStrand'], ';',
                                              'StrandFisherTest=', ro['StrandFisherTest'], ';',
                                              'SOR=', ro['SOR'], ';',
                                              'RefContext=', ro['RefContext'], ';')
                                     }))
    # will then add it below
    
    # format the VCF
    vcf <- vcf %>%
      select('CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'Genotype') %>% # take only the columns we need
      mutate(ID=row_number(), .before='REF') %>% # whatshap needs column ID, which was deleted somewhere; will just add the row number
      mutate(FORMAT='GT', .before='Genotype') %>% # whatshap needs column FORMAT, which was deleted somewhere; I think always GT in my case
      mutate(INFO=infocolumn, .before='FORMAT') %>% # whatshap needs column FORMAT, which was splitted into its components for readibility, will just add NA here
      rename('sample'='Genotype') %>% # column Genotype was originally called 'sample' I think, will change it back. It gives the zygosity
      rename('#CHROM'='CHROM') # add comment character so understands it is the header
    
    # ! whenever there is only SNV, need to edit the VCF to 'make it look like' it went through phasing
    # (whatshap phase does not work with a single SNV as the goal is to phase variants in relation to each other,
      # but whatshap haplotag works fine and that what we need to separate the reads from each haplotype)
    
    if (as.numeric(subset(bysmp, SAMPLE==pdg, n)) == 1) {
      
      # change FORMAT and sample to make it look like phased by whatshap phase
      vcf$FORMAT <- 'GT:PS'
      vcf$sample <- '0|1:4688888'
      
      # header is a bit different than below to write as if it went through whatshap phase
      header <- '##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##nanopolish_window=chr20:4680000-4705000
##INFO=<ID=TotalReads,Number=1,Type=Integer,Description="The number of event-space reads used to call the variant">
##INFO=<ID=SupportFraction,Number=1,Type=Float,Description="The fraction of event-space reads that support the variant">
##INFO=<ID=SupportFractionByStrand,Number=2,Type=Float,Description="Fraction of event-space reads that support the variant for each strand">
##INFO=<ID=BaseCalledReadsWithVariant,Number=1,Type=Integer,Description="The number of base-space reads that support the variant">
##INFO=<ID=BaseCalledFraction,Number=1,Type=Float,Description="The fraction of base-space reads that support the variant">
##INFO=<ID=AlleleCount,Number=1,Type=Integer,Description="The inferred number of copies of the allele">
##INFO=<ID=StrandFisherTest,Number=1,Type=Integer,Description="Strand bias fisher test">
##INFO=<ID=SOR,Number=1,Type=Float,Description="StrandOddsRatio test from GATK">
##INFO=<ID=RefContext,Number=1,Type=String,Description="The reference sequence context surrounding the variant call">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr20>
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">'
      
      # write the header
      vcfname <- paste0(pdg, '_whap.vcf') # notice p (below is _wha) to indicate already phased
      
      write.table(header, file=here('haplotypephasing/vcf_phased', vcfname),
                  quote=FALSE, row.names=FALSE, col.names=FALSE)
      
      # then append the VCF
      write.table(vcf, file=here('haplotypephasing/vcf_phased', vcfname),
                  quote=FALSE, row.names=FALSE, sep='\t', append=TRUE)
      
    } else if ((as.numeric(subset(bysmp, SAMPLE==pdg, n)) > 1)) {
      # in the other cases, i.e. when we have more than 1 SNP, will go through whatshap phase so write a more normal VCF
      
      
      # first write the header
      # ! deleted StrandSupport as there were formatting errors
      # SupportFractionByBase field also disappeared
      header <- '##fileformat=VCFv4.2
##nanopolish_window=chr20:4680000-4705000
##INFO=<ID=TotalReads,Number=1,Type=Integer,Description="The number of event-space reads used to call the variant">
##INFO=<ID=SupportFraction,Number=1,Type=Float,Description="The fraction of event-space reads that support the variant">
##INFO=<ID=SupportFractionByStrand,Number=2,Type=Float,Description="Fraction of event-space reads that support the variant for each strand">
##INFO=<ID=BaseCalledReadsWithVariant,Number=1,Type=Integer,Description="The number of base-space reads that support the variant">
##INFO=<ID=BaseCalledFraction,Number=1,Type=Float,Description="The fraction of base-space reads that support the variant">
##INFO=<ID=AlleleCount,Number=1,Type=Integer,Description="The inferred number of copies of the allele">
##INFO=<ID=StrandFisherTest,Number=1,Type=Integer,Description="Strand bias fisher test">
##INFO=<ID=SOR,Number=1,Type=Float,Description="StrandOddsRatio test from GATK">
##INFO=<ID=RefContext,Number=1,Type=String,Description="The reference sequence context surrounding the variant call">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
      
      vcfname <- paste0(pdg, '_wha.vcf')
      write.table(header, file=here('haplotypephasing/vcf_unphased', vcfname),
                  quote=FALSE, row.names=FALSE, col.names=FALSE)
      
      # then append the VCF
      write.table(vcf, file=here('haplotypephasing/vcf_unphased', vcfname),
                  quote=FALSE, row.names=FALSE, sep='\t', append=TRUE)
      
    }
  }
}