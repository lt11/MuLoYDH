We are currently upgrading the pipeline. If you want to try it and you find any issue please contact us writing an email to "nicolo.tellini.2@gmail.com" or "lorenzo.tattini@gmail.com".

"The mule [mulo in Italian] always appears to me a most surprising animal. That a hybrid should possess more reason, memory, obstinacy, social affection, powers of muscular endurance, and length of life, than either of its parents, seems to indicate that art has here outdone nature."

The Voyage of the Beagle, Charles Darwin

# README #

This README describes the steps required to configure and run the MuLoYDH pipeline.

### What is this repository for? ###

* Quick summary

The MuLoYDH pipeline performs the analysis of paired-end short-read sequencing experiments of clonal samples from yeast diploid hybrids.

MuLoYDH requires two inputs: (1) a dataset of paired-end short-read sequencing experiments of clonal samples from yeast diploid hybrids and (2) the two parental genomes which were used to produce the hybrids in fasta format as well as the corresponding annotations in the "general feature format" (gff). 

A collection of high quality fasta assemblies and the corresponding annotations are provided in the "Rep" folder of the package. For more details regarding the assemblies please refer to Yue, Jia-Xing, et al. "Contrasting evolutionary genome dynamics between domesticated and wild yeasts." Nature genetics 49.6 (2017): 913-924 or visit https://yjx1217.github.io/Yeast_PacBio_2016/welcome/.

Given the two aforementioned inputs, MuLoYDH performs calling and annotation of:

1. de novo small variants (single-nucleotide variants and insertions/deletions smaller than 50 bp)
2. copy-number variants
3. loss-of-heterozygosity (LOH) regions

All the tools embedded in MuLoYDH are run automatically by the wrapper "RunMuLoYDH.sh".

The package contains three folders, the aforementioned "Rep" folder, the "Scr" folder (with the the wrapper and the scripts embedded in it), and the "Lib" folder (with addition scripts, templates and databases).

* Version

This README refers to version 5.

### How do I get set up? ###

* Summary of set up

Once you have downloaded the package you will need to: (1) check for requirements and dependencies and eventually install them (if you do not have root privilegies you will need the help of your system administrator or you can install locally at least the dependencies), (2) prepare input data and folders, (3) configure the wrapper, (4) run the wrapper.

* Requirements

The following linux (openSUSE, https://www.opensuse.org/) packages are required: vi, wget, unzip, gzip, tar, bzip2, git, R, time, curl, java-1_8_0-openjdk-devel, gcc, gcc-c++, perl, make, python, autoconf, automake, zlib-devel, libbz2-devel, xz-devel, libcurl-devel, libopenssl-devel, ncurses-devel, tcsh, gcc-fortran, libxml2, libxml2-devel.

* Dependencies

Command-line tools: SAMtools, freebayes, FastQC, bedtools, MUMmer, BWA, htslib, tabix, bcftools, GEM, VCFtools, SnpEff, Control-FREEC

R libraries: vcfR, ggplot2, IRanges, scales, ape, rtracklayer

A step-by-step guide to install all the dependecies on an openSUSE linux machine is provided in the "Supplementary information" section "Dependencies" (Supplementary file 1.pdf) of the paper and can be downlaoded at http://ircan.unice.fr/~ltattini/SuppData.zip.

* Database configuration

SnpEff database is pre-configured. Instructions to configure it for novel organisms are reported in BuildSnpEff.txt.

* Input data

Input fastq file name must be set in the format: "exp_id.R1.fastq.gz", "exp_id.R2.fastq.gz". The prefix "exp_id" may contain only alphanumeric characters and it must not contain the string "ParHyb".

Input fasta files must be named "strain_id.genome.fa" properly sorted, e.g. chrI, chrII, chrIII, chrIV, chrV, chrVI. The prefix "strain_id" may contain only alphanumeric characters (e.g. S288C.genome.fa).

* Wrapper configuration

In order to run the pipeline 10 variables must be defined in the "user settings" section of the wrapper. A description of the variables defined in the header section of the wrapper is reported here:

assemblies labels (must match prefix of assemblies stored in "Rep" folder)

Ref1Label="SK1"

Ref2Label="YPS128"

short labels (used to name file and plotting)

Ref1="SK1"

Ref2="NA"

prefix of optional parental short read data

ParentalZero="-"

prefix of optional control short read data or "-" to filter variants shared by more than 95% of the samples (the threshold value can be modified through the variable shareThreshold of the FiltSharedVar.R script) or "skip" to bypass this subtraction step

HybridZero="A887R79 A887R80"

label of optional chromosome bearing selection marker (will be filtered before calculation of LOH statistics); if no selection marker has been used the variable must be set to "-"

MarkerChrom="chrII"

"collinear" or "rearranged" parental genomes

GenomesStructure="collinear"

reciprocal overlap for LOH filtering

OverlapThreshold=0.5

read length of the experiments

ReadLength="150"

Using a wrapper configured as above means that:

1. The "Rep/Ann"" folder contains the following files: SK1.all_feature.gff, SK1.subtel.txt, SK1.centromere.txt, YPS128.all_feature.gff, YPS128.subtel.txt, YPS128.centromere.txt
2. The "Rep/Asm" contains the following files: SK1.genome.fa, YPS128.genome.fa
3. The plot produced will bear the labels "SK1" and "NA" for the two parental genomes
4. No parental paired-end short-read data is provided (default settings)
5. The "Exp" folder contain 2 control samples (4 files: A887R79.R1.fastq.gz, A887R79.R2.fastq.gz, A887R80.R1.fastq.gz, A887R80.R2.fastq.gz), and a number of test samples with other IDs
6. Chromosome II was used to select hybrids thus LOH statistics will be calculated filtering out events from chromosome II
7. Parental genomes are known a priori to be collinear
8. One LOH in a test sample will be considered the same event as one LOH occurring in a control sample (and filtered out) if the reciprocal overlap is larger than 0.5

* How to run tests

Experiments must be placed in the "Exp" folder under the main folder (i.e. at the same level of "Scr" and "Rep"). Once the experiments have been placed in the "Exp" folder you can change directory to folder "Scr" and run a test using:

$ bash MuLoYDH.V*.sh &

### Modules Description ###

* MuLoYDH general description

The MuLoYDH pipeline requires as input: (1) a dataset of short-read sequencing experiments from yeast diploid hybrids and (2) the two parental genomes which were used to produce the hybrids in FASTA format as well as the corresponding annotations in the "general feature format" (GFF). Reads from hybrid data are mapped against the assemblies of the two parental genomes separately (standard mappings) and against the union of the two aforementioned assemblies (namely a multi-FASTA obtained concatenating the two original assemblies) to produce the competitive mappings. In the latter case, reads from parent 1 are expected to map to the assembly of parent 1 on the basis of the presence of single-nucleotide markers. Conversely, reads from parent 2 are expected to map to the assembly of parent 2. Standard mappings are used to determine the presence of LOH regions and CNVs. The latter are also exploited to discriminate LOHs due to recombination from those resulting by deletion of one parental allele. The markers between the parental assemblies are determined by the NUCmer algorithm and are exploited to map LOH segments. Markers are genotyped from standard mappings. De novo small variants are determined from both competitive and standard mappings. Competitive mapping allows for direct variant phasing in heterozygous regions. Variant calling from competitive mapping is performed setting ploidy = 1 in heterozygous regions and ploidy = 2 in LOH blocks. Regions characterized by reads with low mapping quality (MAPQ < 5 in the competitive mapping) are assessed from standard mapping using arbitrarily the assembly from parent 1.

* Determination of single-nucleotide marker positions

Single-nucleotide marker positions are determined through the NUCmer algorithm. In order to obtain reliable marker positions and take advantage of the "seed and extend" strategy of the algorithm, markers are calculated in both direct (assembly 1 vs assembly 2) and reverse (assembly 2 vs assembly 1) ways. The intersection of the two sets is retained for LOH detection and to calculate statistics. In the collinear mode markers are determined chromosome-by-chromosome, aligning a chromosome of parent 1 against the corresponding homologous from parent 2, whereas with the rearranged option they are calculated through a single whole-genome alignment of parental assemblies. Running MuLo in collinear mode (a) provides more uniform marker distribution in core chromosomal regions but it can be used only when collinearity between the two parental subgenomes holds.

* Classification of single-nucleotide markers

Markers are classified as lying in collinear or rearranged regions as determined by MUMmer and custom R scripts. The fraction of markers lying within rearranged regions is calculated taking into account inter- and intra-chromosome inversions and translocations.

* Markers genotyping, small variants calling, annotation and filtering

Markers calling and genotyping is performed using SAMtools (mpileup) and BCFtools (call) from standard mappings. Markers are quality filtered removing those with quality <(μ−σ), where μ is the sample marker mean quality value and σ is the corresponding standard deviation as well as those lying in telomeric and subtelomeric regions. The strategy implemented in MuLoYDH for calling de novo small variants relies on a stringent procedure to limit the number of false positives and keep the number of false negatives as low as possible. De novo SNVs and indels are called with (1) SAMtools and BCFtools and (2) FreeBayes. Only variants called by both are retained. Both callers are exploited using competitive and standard mappings as described above. Regions characterized by reads with MAPQ <5 in competitive mappings are determined by custom R scripts, bash scripts and BEDTools. Parental and control hybrid variation is subtracted from hybrids data using custom bash scripts, VCFtools and tabix The resulting variants are quality filtered masking those characterized by quality <(μ−σ), where μ is the sample markers mean quality value and σ is the corresponding standard deviation. Variants bearing marker alleles are filtered out, while those lying within (sub)telomeric regions are masked. Small variants are annotated by means of SnpEff.

* Copy-number variants calling and annotation

Copy-number variants are estimated by means of Control-FREEC with no matched normal samples, using standard mappings against both parental genomes. Read-count data are normalized by GC-content and mappability. Mappability is calculated with GEM-mappability. Results are annotated with p-values calculated with both Kolmogorov-Smirnov and Wilcoxon Rank-Sum tests.

* B-allele frequency calculation

B-allele frequency (BAF) values are calculated from standard mapping as Na/(Nr+Na), where: Na is the number of read bearing the most abundant alternative (non-reference) allele and Nr is the number of reads bearing the reference allele at each marker position.

* Loss-of-heterozygosity detection and annotation

Loss-of-heterozygosity regions are determined and annotated using custom R scripts. Considering standard mappings of each hybrid against both parental assemblies, marker positions characterized by non-matching genotype or alternate allele are filtered out, as well as multiallelic sites and those lying in telomeric or subtelomeric regions. Markers involved in large deletions, as predicted by Control-FREEC, are masked. Finally, stretches of consecutive marker positions are grouped in LOH regions. Genomic coordinates of each LOH event are determined using both the "first/last coordinates" and the "start/end coordinates". First/last coordinates are determined using the coordinates of the first and the last markers of the event. Start/end coordinates are calculated using the average coordinate of the first (last) marker and the last (first) marker of the adjacent event. LOH regions are annotated as terminal/interstitial as well as with genomic features embedded and those potentially involved in breakpoints. Annotation is performed based on the genomic features downloaded from the "Population-level Yeast Reference Genomes" website. Interstitial LOHs are defined as homozygous segments that are flanked on both sides by heterozygous markers. Terminal LOHs are defined as homozygous regions extended to the end of the chromosmal arm.

* Calculation of low-marker-density-regions

Regions characterized by less than one marker in 300 bp are calculated using custom R scripts which are embedded in the MuLoYDH pipeline.

### Results ###

The pipeline produces a number of folders with your results. In the following sections you will find:

* Type of result

Folder/: a brief description of the results reported in the "Folder".

* Reference assemblies

Ref/: fasta sequences used in the analysis.

* Assemblies alignment

MUMmer/: markers detected with MUMmer. Those which are used to call LOHs are reported in the intersect.snps file (e.g., SK1_S288C.intersect.snps). The regions characterized by low marker density are reported in the files within the folder "LowDensityRegions".

* Mappings

MapDdp/: bam files.

* Coverage

DepthDdp/: coverage statistics calculated from de-duplicated bam files.
DepthRaw/: coverage statistics calculated from raw (with duplicates) bam files.

* Markers

Markers/: vcf files of samples and controls.

* CNVs

CNV/: configuration files for all the sample and controls (folder "Config"), mappability data (folder "Mappability") and results for both the assemblies (folder "Results").

* LOHs

LOH/: sample-specific folders (e.g. folder "A887R79") contain the graphical representation of the LOH regions called. Plots of markers distributions are reported in folder "Markers" along with plots reporting markers distances, quality statistics (e.g. file "Stat.Markers.SK1-BY.txt"), and a list of marker classified as lying in structural-rearranged regions. Heatmaps of the LOHs detected in the dataset are reported in folder "AllSegments", e.g. file "HeatMap.SK1.pdf" and "HeatMap.SK1.jpg". The subfolder "Bin25_HistoInterTerminal" contains plots with the length distribution of terminal and interstitial events.

LOH/Summary/: this subfolder reports the main results. Folder "Stat" reports the raw mean statistics per sample for LOH regions detect in the dataset. Under "Plot" you can find plots of BAF data and plots of the mean density of interstitial events (e.g. file "MDIE.SK1.pdf"). Folder "Tables" contains the file "LOH.features.txt" that reports all the LOH segments annotated with upstream/downstream/overlapping genomic features. File "LOH.filtered.variants.txt" reports the LOH segments detected in the dataset. The columns reported in the file refer to:

sample -> sample ID

ref -> reference genome

chr -> chromosome of the segment

start -> start coordinate of the segment calculated as the average position from the first marker of the segment and the last one of the upstream segment

first -> coordinate of the first marker of the segment

last -> coordinate of the last marker of the segment

end -> end coordinate of the segment calculated as the average position from the last marker of the segment and the first one of the downstream segment

status -> status of the segment, 0 = homozygous reference, 1 = heterozygous, 0* events which were merged into a larger LOH

len -> number of markers lying within the segment

denES -> average number of bp between markers of the segment calculated using "start" and "end" coordinates

evrES -> average base-pairs separating markers of the segment calculated using "start" and "end" coordinates

distLF -> average number of bp between markers of the segment calculated using "last" and "first" coordinates

evrLF -> average base-pairs separating markers of the segment calculated using "last" and "first" coordinates

SNV -> number of single-nucleotide variants detected within the segment

denSNV -> density of single-nucleotide variants detected within the segment, calculated using "last" and "first" coordinates

InDel -> number of indels detected within the segment 

denInDel -> density of indels detected within the segment, calculated using "last" and "first" coordinates

TI -> annotation for interstitial LOH segments (int) and terminal segments (ter); events detected also in the control (and thus subtracted) are annotated as "sub"

### Who do I talk to? ###

* Repo owner or admin

For any question or for reporting bugs please write to "nicolo.tellini.2@gmail.com" and "lorenzo.tattini@gmail.com".

### License ###

MuLoYDH is distributed under the MIT license. A number of dependencies are under more restrictive licenses, for which commerical use of the software needs to be discussed with the corresponding developers.

"Ah, pardon! Tarapio tapioco come se fosse antani, la supercazzola prematurata con dominus vobiscum blinda??"
[Earl Conte Raffaello 'Lello' Mascetti]
