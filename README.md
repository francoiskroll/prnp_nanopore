# prnp_nanopore

Data and code for

François Kroll\*\, Athanasios Dimitriadis\*\, Tracy Campbell, Lee Darwent, John Collinge, Simon Mead, Emmanuelle Viré. 2022.

**Prion protein gene mutation detection using long-read Nanopore sequencing.**

_co-first authors*_

https://www.medrxiv.org/content/10.1101/2022.03.06.22271294v1


Below explains how to recreate the analyses.  

Please cite if you use some of the data or code.  

___

## about .command bash scripts

I use macOS.  
.command bash scripts are included in /utilities/.  

Add /utilities/ to PATH so scripts (and other scripts they depend on) are found.  

You probably need to add permissions for each .command script in /utilities/ with

    chmod u+x ~/.../utilities/XXX.command

___

Get in touch for questions

  * [![alt text][1.2]][1] [@francois_kroll](https://twitter.com/francois_kroll)

  * :email: francois@kroll.be

<!-- icons with padding -->
[1.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)

<!-- icons without padding -->
[1.2]: http://i.imgur.com/wWzX9uB.png (twitter icon without padding)

<!-- links to your social media accounts -->
[1]: https://twitter.com/francois_kroll

___

## GC% plot

Refers to GC% plot in Figure 1A.  

In directory /GCplot/, find `GCplotter.R`

Input =
* prnp_window1.fa
* prnp_window2.fa
* prnp_window3.fa

___

## Filtering the gene-body bam alignment files

Directory /genebodyBams/ contains .bam & .bam.bai files for each sample, gene-body amplicon.  
Please find /genebodyBams/ in Zenodo archive at XXX.  

In directory /utilities/, find `prnp_filterBam.command`.

Script loops through bam files in /genebodyBams/ and filters each with command  

    filterBam.command -i $bam -f 10000 -c 15000 -s 0.05 -p yes -o $out  

Read about filtering parameters in comments in `prnp_filterBam.command`.

`filterBam.command` is included in /utilities/.  

About `filterBam.command`:
* it calls a R script `readsToThrow.R`, included in /utilities/. Path is hard-coded.
* Note another hard-coded path to `picard.jar`  

In summary, to filter the gene-body bam files:

    cd ~/.../prnp_nanopore/
    prnp_filterBam.command

Which creates filtered version of each bam file in folder /bamfilt/.  
You can directly find this folder at Zenodo archive XXX.  

___

## Filtering the promoter bam alignment files

Same logic as for the gene-body amplicon above.  

Directory /promoterBams/ contains .bam & .bam.bai files for each sample, promoter amplicon.  
Please find /promoterBams/ in Zenodo archive at XXX.  

In directory /utilities/, find `promoter_filterBam.command`.

Script loops through bam files in /promoterBams/ and filters each with command  

    filterBam.command -i $bam -f 1000 -c 3500 -s 0.15 -p yes -o $out

Read about filtering parameters in comments in `promoter_filterBam.command`.

`filterBam.command` is included in /utilities/. It is the same script as used for the gene-body bam files above.  

In summary, to filter the promoter bam files:

    cd ~/.../prnp_nanopore/
    promoter_filterBam.command

Which creates filtered version of each bam file in folder /promoterbamfilt/.  
You can directly find this folder at Zenodo archive XXX.  

___

## Haplotype phasing

Haplotype phasing is done with `whatshap v1.1`.  

The first step is to phase the SNVs in each sample's VCF file.  

The SNVs were called by `nanopolish v0.13.2`. Find this step in athanadd's repo https://github.com/athanadd/prnp-nanopore-seq.

Typically, haplotype phasing would be performed directly on the VCF from the variant calling algorithm, here from nanopolish. Here, we re-wrote a VCF file for each sample containing only SNVs after strand bias filtering, using _AdditionalFile1.xlsx_ sheet _SNVs_filtered_ as input.   

In directory /utilities/ find `prepareVCFforWhatshap.R`.  
Input is _AdditionalFile1.xlsx_, especially sheet _SNVs_filtered_.

`prepareVCFforWhatshap.R` will write one VCF file for each sample, except for five samples which do not carry any SNV:
* #58648
* #54890
* #55050
* #59060
* #53747

As these samples do not carry any SNV, they cannot be haplotype-phased, hence there is no use creating a VCF file for them.  

Output folder for the new VCF files is /haplotypephasing/vcf_unphased/. Each VCF is called _SAMPLEID_wha.vcf_, wha is for whatshap.  

Note, whatshap cannot phase a single SNV, but a single heterozygous SNV is sufficient to haplotag the reads. Therefore, if the sample only carries a single SNV, `prepareVCFforWhatshap.R` will write the VCF 'as if' it already went through `whatshap`, i.e. it is written directly in /haplotypephasing/vcf_phased/ and is called _SAMPLEID_whap.vcf_, whap is for whatshap phased.  

Once the VCFs are ready, the actual phasing is done by script `prnp_haplotypePhasing.command`. Find it in /utilities/.  

The key steps `prnp_haplotypePhasing.command` are:  
* phase the SNVs in each VCF using `whatshap phase`. For each sample, it writes phased VCF in /haplotypephasing/vcf_phased/ called _SAMPLEID_whap.vcf_, whap is for whatshap phased.  
* using the phased VCF, haplotag the reads with `whatshap haplotag`. For each sample, it finds its filtered BAM in /bamfilt/.  

The above is called possibility 2 in `prnp_haplotypePhasing.command`. Possibility 1 is for the samples which carry a single SNV. For these, `prepareVCFforWhatshap.R` directly created a phased VCF, so we skip `whatshap phase` and directly haplotag the reads (the filtered BAM file) using `whatshap haplotag`.  

Note `prnp_haplotypePhasing.command` has hardcoded calls to the human reference genome hg38.fa. Cmd + F "hg38.fa" to find and modify them.  

For each sample which can be haplotype phased/haplotagged, the main output from `prnp_haplotypePhasing.command` should be a phased VCF file in /haplotypephasing/vcf_phased/ and a haplotagged BAM in /haplotypephasing/bam_haplotagged/ called _SAMPLEID_hp.bam_ (hp for haplotype).  

Five samples cannot be haplotype-phased as they do not carry any SNV. We want to manually move the BAM for these samples in a new folder /haplotypephasing/bam_notag/.  

    mkdir haplotypephasing/bam_notag/

    cp bamfilt/58648_f.bam haplotypephasing/bam_notag/58648_f.bam
    cp bamfilt/58648_f.bam.bai haplotypephasing/bam_notag/58648_f.bam.bai

    cp bamfilt/54890_f.bam haplotypephasing/bam_notag/54890_f.bam
    cp bamfilt/54890_f.bam.bai haplotypephasing/bam_notag/54890_f.bam.bai

    cp bamfilt/55050_f.bam haplotypephasing/bam_notag/55050_f.bam
    cp bamfilt/55050_f.bam.bai haplotypephasing/bam_notag/55050_f.bam.bai

    cp bamfilt/59060_f.bam haplotypephasing/bam_notag/59060_f.bam
    cp bamfilt/59060_f.bam.bai haplotypephasing/bam_notag/59060_f.bam.bai

    cp bamfilt/53747_f.bam haplotypephasing/bam_notag/53747_f.bam
    cp bamfilt/53747_f.bam.bai haplotypephasing/bam_notag/53747_f.bam.bai

You can find directly folders    
* /haplotypephasing/vcf_unphased/
* /haplotypephasing/vcf_phased/
* /haplotypephasing/bam_haplotagged/
* /haplotypephasing/bam_notag/

at the Zenodo archive XXX.  

___

## Trim the reads to keep only the octapeptide repeat region (OPR)

This is performed by `prnp_OPR.command`, found in /utilities/. It is ran once on all the BAM in /haplotypephasing/bam_haplotagged/ and once on all the /haplotypephasing/bam_notag/ (see below).  

For each sample (BAM file), the key steps are:  
* trim the reads in the BAM to keep only the OPR using `samtools ampliconclip`. Note, this step uses _OPRpos.bed_ in /needleSam/ folder, which are the positions to trim to keep only the OPR, i.e. the positions to clip/exclude. Beware, the path to _OPRpos.bed_ is hardcoded. The main output is the BAM containing only OPR reads, written in new folder /OPRtrim/ and named _SAMPLEID_opr.bam_. The step also creates a log file in new folder /clipLogs/, which we do not use further.  
* filter the trimmed BAM (only OPR reads) using `filterBam.command`. Read the comments for explanations about the filtering parameters. This step will partly repeat the first filtering performed on the BAM files.  
* count read lengths. This step will write a TXT file for each sample in new folder /readlengths/ named _SAMPLEID_lengths.txt_. The TXT file has two column: 1– number of reads of that length / 2– read length.  
* convert the trimmed BAM to SAM format, which we will use in `needleSam.R` (see below).    

First we run `prnp_OPR.command` on the haplotagged BAM files:  

    cd ~/Dropbox/nanopore/haplotypephasing/bam_haplotagged/
    prnp_OPR.command

Second we run `prnp_OPR.command` on the BAM files which could not be haplotagged:  

    cd ~/Dropbox/nanopore/haplotypephasing/bam_notag/
    prnp_OPR.command

Accordingly, folders
* /clipLogs/
* /OPRtrim/
* /OPRtrim_sam/
* /readlengths/

are found in the folder /bam_haplotagged/ or in the folder /bam_notag/ depending on the sample.  

As in the previous step, you can find directly folders    
* /haplotypephasing/vcf_unphased/
* /haplotypephasing/vcf_phased/
* /haplotypephasing/bam_haplotagged/
* /haplotypephasing/bam_notag/

at the Zenodo archive XXX.  

___

## Generate catalog of OPR templates

Finding candidate somatic mutations of the OPR will involve aligning OPR reads to template OPR sequences.  

Please look at _AdditionalFile1.xlsx_, sheet _OPRconsensus_ at this stage. The OPR templates are built from a consensus sequence of one R (R being a single repeat unit) which makes use of IPUPAC flexible nucleotide codes (see https://www.bioinformatics.org/sms/iupac.html). By consensus sequence, we mean that all of R2, R3, R4 are included in the consensus sequence. For example, the first codon of the repeat unit is always CCT or CCC. Accordingly, it is written as CCY in the consensus sequence, where Y can be either C or T. R1 has an extra codon compared to R2, R3, R4 so it can be written with exactly the same consensus sequence but it is also written with flexible nucleotides to allow for the same changes at Wobble positions.  

The catalog of OPR templates is generated by `generateOPRCatalog.R`, found in /needleSam/. Finding a somatic mutation of the OPR is like finding a needle in a SAM alignment file, hence 'needleSam'.  

`generateOPRCatalog.R` writes _OPRConsensusCatalog.csv_ in folder /needleSam/. It contains 29 OPR templates, from 4 OPRD (deletion of all R repeats except R1) up to 24 OPRI (insertion of 24 extra R repeats). The table was copied to _AdditionalFile1.xlsx_, sheet _OPRtemplates_.  


___

## Somatic mutation search

This is performed by `needleSam.R`, found in /needleSam/.  

Inputs to `needleSam.R` are some sheets from _AdditionalFile1.xlsx_ and, most importantly, the OPR trimmed SAM files from above, found in folders /haplotypephasing/bam_haplotagged//OPRtrim_sam/ and /haplotypephasing/bam_notag/OPRtrim_sam/.  

Follow the code/comments in `needleSam.R`. The key steps are:  
* Import the SAM files  
* Calculate the total insertion/deletion of each read by parsing the CIGAR  
* Assign each read to a most likely OPR genotype based on its total insertion/deletion  
* Align each read to its most likely OPR template sequence  
* Calculate the maximum number of mismatches allowed based on the set of 'true reads', i.e. the reads which confirm the Sanger genotype of their sample.  
* Identify candidate somatic OPR reads. A read is a candidate somatic OPR read if it is unexpected from its sample (typically it is different than reference or the known mutated OPR) and passes the mismatch threshold.  

Note, the last step is actually performed by `needleSam_exploration.R`, see below.  

The main output of `needleSam.R` is _allOPRreads.csv_, found in /needleSam/. It stores all the OPR reads from the various SAM files, with various information about each read, e.g. which is its most likely OPR and how many mismatches did it have with that OPR template.

This the main analysis which I hope could be useful to someone else, so do not hesitate to get in touch for questions!

francois@kroll.be or twitter @francois_kroll.

___

## Exploring candidate somatic mutations

This is performed by `needleSam_exploration.R`, found in /needleSam/.  

The main input is the database _allOPRreads.csv_ created by `needleSam.R`.  

`needleSam_exploration.R` first identifies the candidate somatic mutations (see above), then does various analyses and plots.  It also writes _somaticCalls.xlsx_, which was copied in _AdditionalFile1.xlsx_, sheet _somaticmutationcalls_, minus some unnecessary columns. Read about the columns in _AdditionalFile1.xlsx_, sheet _TableofContents_.  

Note, some personal information including precise ages were removed from _AdditionalFile1.xlsx_, sheet _samples_ to protect anonymity. Accordingly, some sub-analyses (not included in the publication) of `needleSam_exploration.R` will throw errors.  

___

## Control PCR

___

## Literature survey of OPR genotypes  

Before having the idea of the OPR consensus sequence, we created a survey of all mutated OPRs reported in literature. The approach using the OPR consensus sequence is better as it generalises to more possible cases, but we are including this survey if it is any useful. It may not be exhaustive, but it should be close.  

Find _OPRLitCatalog.xlsx_ in /OPR_litsurvey/.

Publications typically give the R1/R2/... pattern and not the full nucleotide sequence. Script `seqFromOPRpattern.R` generates the full OPR sequence from the R1/R2/... pattern. It uses `eachOPR.xlsx` for the sequence of each repeat.  

___
