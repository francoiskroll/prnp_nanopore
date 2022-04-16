# prnp_nanopore

Data and code for

François Kroll\*\, Athanasios Dimitriadis\*\, Tracy Campbell, Lee Darwent, John Collinge, Simon Mead, Emmanuelle Viré. 2022.  

_co-first authors*_

**Prion protein gene mutation detection using long-read Nanopore sequencing.**

https://doi.org/10.1101/2022.03.06.22271294

Below explains how to recreate the analyses.  

Please cite if you use some of the data or code.  

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

## about fast5/fastq/bam sequencing data

fast5 files were too heavy to upload even for Zenodo. Please get in touch if you would like to access them.  

fastq files (i.e. after basecalling) are available at the Zenodo archive https://doi.org/10.5281/zenodo.6427185, folder /fastq/.  

bam files (i.e. after alignment to human reference genome hg38) are available at the Zenodo archive https://doi.org/10.5281/zenodo.6427185, folders /genebodyBams/ and /promoterBams/.  

___

## about .command bash scripts

I used macOS.  
.command bash scripts are included in /utilities/.  

Add /utilities/ to PATH so scripts (and other scripts they depend on) are found.  

You probably need to add permissions for each .command script in /utilities/ with

    chmod u+x ~/.../utilities/XXX.command

___

## about .Rproj and package here()

In R scripts, I tried to avoid hard-coded paths using the package `here()`.  

`here()` starts the path wherever the .Rproj file of the R project opened is.  

So, for this to work, make sure to get file _prnp_nanopore.Rproj_ included and open it first in RStudio before opening a script. You can then run `here()` to check where it is starting. It should start at the folder containing the folders from this repositories, i.e. containing /utilities/, /GCplot/,  /SNVs/, etc.   

___

## GC% plot

Refers to GC% plot in Figure 1A.  

In directory /GCplot/, find `GCplotter.R`

Input =
* _prnp_window1.fa_
* _prnp_window2.fa_
* _prnp_window3.fa_

___

## Pilot experiment sequencing protein-coding amplicon

This refers to Figure S1.  

#### basecalling

Basecalling of the fast5 files was done with `guppy` with command:

    guppy_basecaller --input_path ~/.../fast5_pooled \
      --save_path ~/.../pilotexp/ \
      --flowcell FLO-MIN106 \
      --kit SQK-LSK109 \
      --calib_detect \
      --verbose_logs \
      --records_per_fastq 0 \
      --compress_fastq

* I used a calibration strand, so added flag `--calib_detect` so it detects it and puts it a separate output.
* `--verbose_logs` to output a log file.
* By default, puts 4000 reads per fastq file. By setting to 0, output all reads in one file.
* `--compress_fastq`: gzip compression of output fastq files (~ 50% reduced file size).

Default settings include
* High Accuracy basecall model
* Q-score filtering OFF

Generates three fastq files:
* _fastq_runid_d3373c24347bdd069710849e408f6251f1da2c57_0_0.fastq.gz_ (1060 reads)
* _fastq_runid_dbccf13003c5dd981e72427f54717c842ee74709_0_0.fastq.gz_ (33 reads)
* in /calibration_strands/_fastq_runid_d3373c24347bdd069710849e408f6251f1da2c57_0_0.fastq.gz_ (22 reads)

To count number of reads in each fastq.gz:

    numberReadsFastqGz.command fastq_runid_d3373c24347bdd069710849e408f6251f1da2c57_0_0.fastq.gz
    numberReadsFastqGz.command fastq_runid_dbccf13003c5dd981e72427f54717c842ee74709_0_0.fastq.gz
    numberReadsFastqGz.command calibration_strands/fastq_runid_d3373c24347bdd069710849e408f6251f1da2c57_0_0.fastq.gz

`numberReadsFastqGz.command` is included in /utilities/.  

Can merge the first two fastq.gz. The smaller file (33 reads) contains reads generated during mux scan, which are fine to use.  

To merge the first two fastq.gz:  

    cd ~/.../pilotexp/
    cat *.fastq.gz > 180421_pool.fastq.gz

I did this analysis on 18/04/2021, hence the filename.  

And check the merged file:

    numberReadsFastqGz.command 180421_pool.fastq.gz  

Merged file _180421_pool.fastq.gz_ now has 1093 reads, so correct.  

#### alignment

Alignment to hg38 was done using `minimap2` with command:

    minimap2 -ax map-ont ~/.../hg38.fa \
      180421_pool.fastq.gz > 180421.sam

Then to index/sort:

    samtools view 180421.sam -o 180421.bam
    samtools sort 180421.bam > 180421s.bam
    samtools index 180421s.bam

#### consensus sequence

Consensus sequence was computed by `canu`.  

Prior to this, I removed all reads above 2000 bp using `nanofilt` (https://github.com/wdecoster/nanofilt):

    gunzip -c 180421_pool.fastq.gz | NanoFilt --maxlength 2000 | gzip > 180421_max2000.fastq.gz  

Could do this with script `filterBam.command` included in /utilities/ (see below). I had not written it at the time.  

    canu \
      -p con -d canuconsensus \
      genomeSize=1.5k \
      -nanopore 180421_max2000.fastq.gz

This writes _con.contigs.fasta_ which contains the consensus sequence.  

#### polish the consensus sequence

Using `nanopolish`.  

First to align the fastq reads to their consensus sequence:

    minimap2 -ax map-ont canucon/con.contigs.fasta \
      180421_max2000.fastq.gz > 180421con.sam

    samtools view 180421con.sam -o 180421con.bam
    samtools sort 180421con.bam > 180421cons.bam
    samtools index 180421cons.bam

Then to polish the consensus sequence:

    nanopolish index -d fast5_pooled/ \
      180421_max2000.fastq.gz

    nanopolish variants --consensus -o canuconpolish.vcf \
        -r 180421_max2000.fastq.gz \
        -b 180421cons.bam \
        -g con.contigs.fasta

This writes _canuconpolish.vcf_ which lists the edits that `nanopolish` would make to the consensus sequence.  

To make these edits:

    nanopolish vcf2fasta --skip-checks -g con.contigs.fasta canuconpolish.vcf > canuconpolish.fa

I then trimmed the adapters manually in Benchling and aligned the consensus to the Sanger sequence to create Figure S1b.  

### call single-nucleotide variants

Sample was M129V, called by Sanger sequencing.  

Indexing fastq to fast5 files:

    nanopolish index -d ~/.../fast5_pooled 180421_pool.fastq.gz

Then calling variants with:

    nanopolish variants \
      --reads 180421_pool.fastq.gz \
      --bam 180421s.bam \
      --genome ~/.../hg38.fa \
      --outfile 180421_nanopolish.vcf \
      --ploidy 2 \
      --window "chr20:4680000:4705000" \
      --threads 4 \
      --min-candidate-frequency 0.1 \
      --min-candidate-depth 10 \
      --calculate-all-support

Which correctly called M129V chr20:4699605A>G.  

fast5 files and _180421_pool.fastq.gz_ are included at the Zenodo archive https://doi.org/10.5281/zenodo.6427185 in folder /pilotexp/.  

___

## All panel: basecalling/demultiplexing/alignment/SNV calling

All panel = the 25 samples listed in Table 1.  

Please find these steps in co-first author Athanasios Dimitriadis' repository at https://github.com/athanadd/prnp-nanopore-seq.  

The below starts after these steps, namely it uses SNV calls from `nanopolish` and bam files after alignment with `minimap2`.  

___

## Single-nucleotide variants (SNVs)

SNVs were called on all samples of the panel, see co-first author Athanasios Dimitriadis' repository at https://github.com/athanadd/prnp-nanopore-seq.  

Find script `SNVscalls.R` in folder /SNVs/. Main input is _AdditionalFile1.xlsx_. It plots Figure 2a,b; panels of Figure S2; and extra plots not included in the manuscript.  

Question from reviewer (paraphrased):  
> M129V individuals tend to have a higher number of non-coding SNVs (Table 1).  
Are these SNVs occuring on the same haplotype?  

Wrote script `SNVs_M129V.R` to answer.  
It expects phased vcf files in /haplotypephasing/vcf_phased/ as input, see below how these were generated.   

___

## Filtering the gene-body bam alignment files

Directory /genebodyBams/ contains bam & bam.bai files for each sample, gene-body amplicon.  
Please find /genebodyBams/ in Zenodo archive at https://doi.org/10.5281/zenodo.6427185.  

In directory /utilities/, find `prnp_filterBam.command`.

Script loops through bam files in /genebodyBams/ and filters each with command  

    filterBam.command -i $bam -f 10000 -c 15000 -s 0.05 -p yes -o $out  

Read about filtering parameters in comments in `prnp_filterBam.command`.

`filterBam.command` is included in /utilities/.  

About `filterBam.command`:
* it calls a R script `readsToThrow.R`, included in /utilities/. Path is hard-coded.
* note another hard-coded path to `picard.jar`  

In summary, to filter the gene-body bam files:

    cd ~/.../prnp_nanopore/
    prnp_filterBam.command

Which creates filtered version of each bam file in folder /bamfilt/.  
You can directly find this folder at Zenodo archive https://doi.org/10.5281/zenodo.6427185.  

___

## Filtering the promoter bam alignment files

> Note, we have not always been consistent with regards to how this amplicon (the 2988-bp amplicon in Figure 1a) is called. It is sometimes referred to as promoter amplicon or regulatory region amplicon.

Same logic as for the gene-body amplicon above.  

Directory /promoterBams/ contains .bam & .bam.bai files for each sample, promoter amplicon.  
Please find /promoterBams/ in Zenodo archive at https://doi.org/10.5281/zenodo.6427185.  

In directory /utilities/, find `promoter_filterBam.command`.

Script loops through bam files in /promoterBams/ and filters each with command  

    filterBam.command -i $bam -f 1000 -c 3500 -s 0.15 -p yes -o $out

Read about filtering parameters in comments in `promoter_filterBam.command`.

`filterBam.command` is included in /utilities/. It is the same script as used for the gene-body bam files above.  

In summary, to filter the promoter bam files:

    cd ~/.../prnp_nanopore/
    promoter_filterBam.command

Which creates filtered version of each bam file in folder /promoterbamfilt/.  
You can directly find this folder at Zenodo archive https://doi.org/10.5281/zenodo.6427185.  

___

## Calling structural variants (SVs)

Including (inherited) OPR mutations insertions/deletions.  

This is done by tool `sniffles`, ran on filtered bam alignment files from above.  

Find script `runSnifflesloop.command` in /utilities/.  

Main steps are:
* Loop through bam files of promoter/regulatory regions in folder /promoterbamfilt/, writing a vcf file _SAMPLEID_promoterSnif.vcf_ for each sample  
* Loop through bam files of gene-body in folder /bamfilt/, writing a vcf file _SAMPLEID_bodySnif.vcf_ for each sample  

For gene-body amplicon, and hence to call the OPR mutations, the final `sniffles` command was:

    sniffles --mapped_reads "$BAM" \
        --vcf "$VCF" \
        --genotype \
        --min_support 300 \
        --min_length 20 \
        --max_distance 50

`min_support 300`: minimum 300 reads to support the variant. It corresponds to ~ 10% of lowest coverage (see _AdditionalFile1.xlsx_, sheet _sequencing_summary_, column _genebody_coverage_, minimum is 3165x for sample #52331).  

`min_length 20`: minimum insertion/deletion length 20 bp. 1 OPRD or 1 OPRI is 24 bp, plus allow ~ 4 bp error.  

`max_distance 50`: controls how far SVs from different reads can be from each other to be called together (I think). Preprint https://www.biorxiv.org/content/10.1101/2021.05.27.445886v1.full.pdf benchmarked this parameter on their data and decided on 50. We followed their recommendation.  

You can read further about this parameter in thread: https://github.com/fritzsedlazeck/Sniffles/issues/267  

Calls were filtered further with script `parseVCF.R`, in folder `SVs_sniffles`.  

Filters:  
* only calls within _PRNP_ genomic region, precisely chr20:4685060–4701756
* only calls with allele frequency > 0.1  

The filtered calls were then added to _AdditionalFile1.xlsx_, sheet _SVs_.  

##### some comments about OPR SV calling with sniffles

I found StrandBias flag implemented by sniffles difficult to use/trust in our case, see comment in `parseVCF.R`. In the final calls, three still had a StrandBias flag raised (see _AdditionalFile1.xlsx_, sheet _SVs_, column _filter_) while we knew by gel/Sanger sequencing that they were true positives.

Sniffles also tended to call unique mutations multiple times, even after benchmarking the `max_distance` parameter. This created artificially low allele frequencies, as reads carrying a unique mutation are splitted between two calls. In the final calls, this was the case for 4 OPRI sample #57265 (see _AdditionalFile1.xlsx_, sheet _SVs_). Each call has an allele frequency of 0.10 and 0.13 while they represent the same mutation and the call should thus have an allele frequency ~ 0.23.  

In summary, if you are reading this specifically to call OPR mutations in similar data, I would probably recommend using the parameters we benchmarked, but also to inspect carefully both the calls before/after filtering and the reads in IGV. Overall, varying parameters usually made sniffles calls multiple times the same SV, rather than missing SVs entirely.

You are welcome to use our data if you want to test this further/try another SV calling algorithm. Feel free to get in touch if you need help getting started.  

___

## Haplotype phasing

Haplotype phasing is done with `whatshap`.  

The first step is to phase the SNVs in each sample's VCF file.  

The SNVs were called by `nanopolish`. Please find this step in co-first author Athanasios Dimitriadis' repository at https://github.com/athanadd/prnp-nanopore-seq.  

Typically, haplotype phasing would be performed directly on the VCF from the variant calling algorithm, here from `nanopolish`. Here, we re-wrote a VCF file for each sample containing only SNVs after strand bias filtering, using _AdditionalFile1.xlsx_ sheet _SNVs_filtered_ as input.   

In directory /utilities/ find `prepareVCFforWhatshap.R`.  
Input is _AdditionalFile1.xlsx_, especially sheet _SNVs_filtered_.

`prepareVCFforWhatshap.R` will write one vcf file for each sample, except for five samples which do not carry any SNV:
* #58648
* #54890
* #55050
* #59060
* #53747

As these samples did not carry any SNV, they cannot be haplotype-phased, hence there is no use creating a vcf file for them.  

Output folder for the new vcf files is /haplotypephasing/vcf_unphased/. Each vcf is called _SAMPLEID_wha.vcf_, _wha_ is for whatshap.  

Note, whatshap cannot phase a single SNV, but a single heterozygous SNV is sufficient to haplotag the reads. Therefore, if the sample only carries a single SNV, `prepareVCFforWhatshap.R` will write the vcf "as if" it already went through `whatshap`, i.e. it is written directly in /haplotypephasing/vcf_phased/ and is called _SAMPLEID_whap.vcf_, _whap_ is for whatshap phased.  

Once the vcf files are ready, the actual phasing is done by script `prnp_haplotypePhasing.command`. Find it in /utilities/.  

The key steps `prnp_haplotypePhasing.command` are:  
* phase the SNVs in each vcf using `whatshap phase`. For each sample, it writes phased vcf in /haplotypephasing/vcf_phased/ called _SAMPLEID_whap.vcf_, _whap_ is for whatshap phased.  
* using the phased vcf, haplotag the reads with `whatshap haplotag`. For each sample, it finds its filtered bam in /bamfilt/.  

The above is called possibility #2 in `prnp_haplotypePhasing.command`. Possibility #1 is for the samples which carry a single SNV. For these, `prepareVCFforWhatshap.R` directly created a phased vcf, so we skip `whatshap phase` and directly haplotag the reads (the filtered bam file) using `whatshap haplotag`.  

Note `prnp_haplotypePhasing.command` has hard-coded calls to the human reference genome _hg38.fa_. Cmd + F "hg38.fa" to find and modify them.  

For each sample which can be haplotype phased/haplotagged, the main output from `prnp_haplotypePhasing.command` should be a phased vcf file in /haplotypephasing/vcf_phased/ and a haplotagged BAM in /haplotypephasing/bam_haplotagged/ called _SAMPLEID_hp.bam_, _hp_ is for haplotype.  

Five samples cannot be haplotype-phased as they do not carry any SNV. We manually move the bam for these samples in a new folder /haplotypephasing/bam_notag/.  

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

at the Zenodo archive https://doi.org/10.5281/zenodo.6427185.  

___

## Trim the reads to keep only the octapeptide repeat region (OPR)

This is performed by `prnp_OPR.command`, found in /utilities/. It is ran once on all the BAM in /haplotypephasing/bam_haplotagged/ and once on all the /haplotypephasing/bam_notag/ (see below).  

For each sample (bam file), the key steps are:  
* trim the reads in the bam to keep only the OPR using `samtools ampliconclip`. Note, this step uses _OPRpos.bed_ in /needleSam/ folder, which are the positions to trim to keep only the OPR, i.e. the positions to clip/exclude. Note, the path to _OPRpos.bed_ is hard-coded. The main output is the bam containing only OPR reads, written in new folder /OPRtrim/ and named _SAMPLEID_opr.bam_. The step also creates a log file in new folder /clipLogs/, which we do not use further.  
* filter the trimmed bam (only OPR reads) using `filterBam.command`. Read the comments for explanations about the filtering parameters. This step will partly repeat the first filtering performed on the bam files.  
* count read lengths. This step will write a txt file for each sample in new folder /readlengths/ named _SAMPLEID_lengths.txt_. The txt file has two column: 1: number of reads of that length // 2: read length.  
* convert the trimmed BAM to SAM format, which we will use in `needleSam.R` (see below).    

First we run `prnp_OPR.command` on the haplotagged bam files:  

    cd ~/Dropbox/nanopore/haplotypephasing/bam_haplotagged/
    prnp_OPR.command

Second we run `prnp_OPR.command` on the bam files which could not be haplotagged:  

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

at the Zenodo archive https://doi.org/10.5281/zenodo.6427185.  

___

## OPR lengths

Find script _oprlengths.R_ in folder /oprlengths/.  

It loops through the _SAMPLEID_lengths.txt_ files in the /readlengths/ folders created in the previous section, i.e. /haplotypephasing/bam_haplotagged/readlengths/ and /haplotypephasing/bam_notag/readlengths/, depending on whether the sample could be haplotagged or not.  

For each sample, it creates a histogram frequency of reads vs read length. These are then overlayed by cohort to create Figure 2c,d.  
___

## Generate catalog of OPR templates

Finding candidate somatic mutations of the OPR will involve aligning OPR reads to template OPR sequences.  

Please look at _AdditionalFile1.xlsx_, sheet _OPRconsensus_ at this stage. The OPR templates are built from a consensus sequence of one R (R being a single repeat unit) which makes use of IPUPAC flexible nucleotide codes (see https://www.bioinformatics.org/sms/iupac.html). By consensus sequence, we mean that all of R2, R3, R4 are included in the consensus sequence we called Ri. For example, the first codon of any R is always CCT or CCC. Accordingly, it is written as CCY in the consensus sequence, where Y can be either C or T. R1 has an extra codon compared to R2, R3, R4 so it cannot be written with exactly the same consensus sequence but it is also written with flexible nucleotides to allow for the same changes at wobble positions.  

The catalog of OPR templates is generated by `generateOPRCatalog.R`, found in /needleSam/.

> Finding a somatic mutation of the OPR is like finding a needle in a sam alignment file, hence "needleSam" !  

`generateOPRCatalog.R` writes _OPRConsensusCatalog.csv_ in folder /needleSam/. It contains 29 OPR templates, from 4 OPRD (deletion of all R repeats except R1) up to 24 OPRI (insertion of 24 extra R repeats). The table was copied to _AdditionalFile1.xlsx_, sheet _OPRtemplates_.  


___

## Somatic mutation search

This is performed by `needleSam.R`, found in /needleSam/.  

Inputs to `needleSam.R` are some sheets from _AdditionalFile1.xlsx_ and, most importantly, the OPR trimmed sam files from above, found in folders /haplotypephasing/bam_haplotagged//OPRtrim_sam/ and /haplotypephasing/bam_notag/OPRtrim_sam/.  

Follow the code/comments in `needleSam.R`. The key steps are:  
* Import the sam files  
* Calculate the total insertion/deletion of each read by parsing the CIGAR  
* Assign each read to a most likely OPR genotype based on its total insertion/deletion  
* Align each read to its most likely OPR template sequence  
* Calculate the maximum number of mismatches allowed based on the set of 'true reads', i.e. the reads which confirm the Sanger genotype of their sample.  
* Identify candidate somatic OPR reads. A read is a candidate somatic OPR read if it is unexpected from its sample (typically it is different than reference or the known mutated OPR) and passes the mismatch threshold.  

Note, the last step is actually performed by `needleSam_exploration.R`, see below.  

The main output of `needleSam.R` is _allOPRreads.csv_. Note, it is too heavy for GitHub so please find in Zenodo archive https://doi.org/10.5281/zenodo.6427185, folder /needleSam/. It stores all the OPR reads from the various SAM files, with various information about each read, e.g. which is its most likely OPR and how many mismatches did it have with that OPR template.

This the main analysis which I hope could be useful to someone else, so do not hesitate to get in touch for questions!

francois@kroll.be or twitter @francois_kroll.

___

## Exploring candidate somatic mutations

This is performed by `needleSam_exploration.R`, found in /needleSam/.  

The main input is the database _allOPRreads.csv_ created by `needleSam.R`.  

`needleSam_exploration.R` first identifies the candidate somatic mutations (see above), then does various analyses and plots.  It also writes _somaticCalls.xlsx_, which was copied in _AdditionalFile1.xlsx_, sheet _somaticmutationcalls_, minus some unnecessary columns. Read about the columns in _AdditionalFile1.xlsx_, sheet _TableofContents_.  

___

## Control PCR

These are the data obtained from amplification and sequencing of a reference OPR to test for PCR-introduced errors.  

See Methods, _Control PCR to test for PCR-introduced errors_ for more details about the experiment.  

The analysis follows closely the approach above, only meaningful difference is the input data.  

Reads aligned to hg38 are _controlPCR_56635_sorted_mdfix.bam_.  

It was too heavy to be opened in IGV on my laptop, so subsampled following https://bioinformatics.stackexchange.com/questions/402/how-can-i-downsample-a-bam-file-while-keeping-both-reads-in-pairs/5648  

    cd ~/.../prnp_nanopore/controlPCR
    samtools view -bs 42.1 controlPCR_56635_sorted_mdfix.bam > controlPCR_sample.bam
    samtools index controlPCR_sample.bam

Files  
* controlPCR_56635_sorted_mdfix.bam  
* controlPCR_56635_sorted_mdfix.bam.bai  
* controlPCR_sample.bam  
* controlPCR_sample.bam.bai  

are included in the Zenodo archive https://doi.org/10.5281/zenodo.6427185, folder /controlPCR/.  

bam alignment file is filtered with

    filterBam.command -i controlPCR_56635_sorted_mdfix.bam -f 900 -c 1300 -s 0.20 -p yes -o 56635conPCR_f.bam

_filterBam.command_ is found in /utilities/.  

Note, primers were forward #106 and reverse #46  
Primer #106 positions: chr20:4699127–4699147  
Primer #46 positions: chr20:4700121–4700141  
So amplicon length is 4700141 − 4699127 = 1014 bp, hence `-f 900 -c 1300` above.  

Next, we move _56635conPCR_f.bam_ and _56635conPCR_f.bam.bai_ in a new folder /filt/.  

Then we apply `prnp_OPR.command` on _56635conPCR_f.bam_. See above section _Trim the reads [...]_ for explanations about the script.  

    cd ~/Dropbox/nanopore/controlPCR/filt
    prnp_OPR.command

Outputs from `prnp_OPR.command` are
* _/controlPCR/clipLogs/56635conPCR_clipLog.txt_, i.e. a log file produced by `samtools ampliconclip`, not used further
* _/controlPCR/OPRtrim/56635conPCR_opr.bam_, i.e. bam alignment file with reads trimmed to keep only the OPR  
* _/controlPCR/OPRtrim/56635conPCR_opr.bam.bai_, i.e. index file of above  
* _/controlPCR/OPRtrim/OPRtrim_sam/OPRtrim_sam/56635conPCR_opr.sam_, i.e. sam version of bam file above, will serve as input to `conPCR_needleSam.R` (see below)    
* _/controlPCR/OPRtrim/OPRtrim_sam/56635conPCR_lengths_, i.e. read length vs number of reads, will not use here  

You can find all of the above directly at Zenodo archive https://doi.org/10.5281/zenodo.6427185.  

We then go through `conPCR_needleSam.R`, which takes _56635conPCR_opr.sam_ as main input.  
Main output is `conPCR_allOPRreads.csv`, as previously.  

Note, _conPCR_allOPRreads.csv_ is too heavy for GitHub and can be found instead at the Zenodo archive https://doi.org/10.5281/zenodo.6427185, folder /controlPCR/.  

Next, we go through `conPCR_needleSam_exploration.R`, which takes _conPCR_allOPRreads.csv_ as input. Script is a short version of `needleSam_exploration.R` (see above) mainly to create Figure 3e.  

_controlPCR_somaticcalls.xlsx_: somatic mutation calls from the control PCR.  

___

## Literature survey of OPR genotypes  

Before having the idea of the OPR consensus sequence, we created a survey of all mutated OPRs reported in literature. The approach using the OPR consensus sequence is better as it generalises to more possible cases, but we are including this survey if it is any useful. It may not be exhaustive, but it should be close.  

Find _OPRLitCatalog.xlsx_ in /OPR_litsurvey/.

Publications typically give the R1/R2/... pattern and not the full nucleotide sequence. Script `seqFromOPRpattern.R` generates the full OPR sequence from the R1/R2/... pattern. It uses `eachOPR.xlsx` for the sequence of each repeat.  

___
