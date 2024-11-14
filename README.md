<img width="2309" alt="image" src="https://github.com/user-attachments/assets/5458ee28-a07f-4a5a-b6d6-aa5c896d837e">

# METAFORA
**Metafora** (**Met**hylation **a**nalysis **f**or **o**utlier **r**egion **a**nnotation) is a workflow to call methylation outlier regions--segments of the genome where methylation levels in an individual deviate significantly from the normal levels seen in the general population--from long-read whole genome sequencing data. Metafora has been tested for both PacBio and ONT methylation calls and can further correct for technology-specific biases to allow joint analysis of across different sequencing modalities.

Methylation outliers identified by Metafora can be used to pinpoint regions of the epignome that are being dysregulated. For rare disease genomes, this can help find aberrantly methylated imprinting regions that may lead to imprinting disorders as well as other diagnostic epivariants. In addition, methylation outliers can help interpret the functional impacts of nearby rare variants. 

Benefits of Metafora include its ability to model uncertainty of methylation proportions from variable read depth and its ability to correct for technical covariates and hidden factors when detecting methylation outliers. Metafora uses an unbiased segmentation approach that learns methylation outlier regions from the data rather than relying on annotations of CpG islands or promoters. It also has broad applicability across many sequencing modalities and tools.


## Methods
<img width="1662" alt="image" src="https://github.com/user-attachments/assets/60f7ba2e-14b1-4243-a859-ef14a1f92881">


To call outliers, Metafora uses a four step approach detailed below:
1. First, a mean population methylation reference is generated for a given tissues by aggregating beta values (methylation proportions) over each CpG across all samples sequenced for that tissue. This gives us an expectation for the "normal" level of methylation that we'd predict at any given CpG.
2. Next, the observed methylation profile of a given sample is compared to the population mean reference by conducting a hypothesis test at each CpG to determine if sample-level methylation equals the population mean. We account for uncertainty in the sample-level beta estimate due to variable read depth by modeling these methylation proportion with a beta distribution. P-values from these hypothesis tests are converted to signed z-scores which represents a population deviance score. Values of this score near zero correspond to no difference, lower negative values suggest stronger hypomethylation, and higher positive values suggest stronger hypermethylation.
3. After generating individual CpG population deviance scores, the whole deviance score profile is segmented to find contiguous blocks of CpG that are all consistently deviant. Methylation deviance scores are summarized over these regions. This segmentation gives us candidate outlier regions, blocks of CpGs where a sample deviates significanlty from the population mean.
4. Finally, once we have these candidate outlier regions, methylation proportions are summarized over each region across all samples. These methylation betas are then converted to M-values and then corrected by regression out the effect of known and hidden covariates. Corrected M-values are then scaled and centered to generate methylation z-scores, which can then be thresholded to determine the methylation outliers that exist at the tails of these methylation distributions.

### Automated Global Outlier Detection

Finding meaningful methylation outlier regions assumes these methylation outliers are rare in any given genome and that outside these regions, global methylation profiles are highly concordant. We usually think of methylation outliers being driven in *cis* by nearby rare genetic variants. If a sample's global methylation profile doesnt match the population mean across a substantial portion of the genome--due to technical factors, mismatching tissue, extreme deviance in expected cell-type proprotions, etc.--calling individual methylation outlier regions is not as useful. Further, including these global outliers could bias the results of other samples within the cohort by skewing our estimates of mean methylation. To this end, one important step of Metafora is to identify and exclude global methylation outliers, samples who have global--potentially *trans*-acting--differences in their methylation profiles. 

This process is automated in Metafora by calculating pair-wise correlations of sample methylation profiles. As a stable epigenetic mark, we have found methylation to be highly conserved and correlated across samples. We find on average any two samples should be have pearson correlation values near ~0.96 even when sequenced from different sequencing modalities, though this might vary based on tissue heterogeneity and other factors. To identify global outliers, Metafora calculates a correlation matrix across methylation profiles of all samples, and summarizes each sample by its mean correlation with all other samples. This mean should be centered near 1 for most concordant samples, while global outliers will have much lower outlying values. To call global outliers then, we simply take the distribution of this mean correlation and use tukey's method to call individuals with outlying low levels. These samples will be removed from hidden factors calculation and methylation outlier regions will not be called for them. 

As an example here is a Metafora run of a cohort of 43 PBMC samples, into which 3 "global outliers" were spiked in that were seqeunced from different tissues and incorrectly labeled as PBMC: 2 from Cerebellum tissue, and 1 from a cultured fibroblast sample. As methylation should be tissue-specific, methylation profiles from different tissues should display global methylation deviation and these samples should be significantly less correlated. Metafora's pairwise corrleation analysis detects all three of these samples as global outliers and will automatically exclude them from the analysis. 

<img width="1056" alt="image" src="https://github.com/user-attachments/assets/ec5265df-0853-47b6-aff1-6b86aea62315">

In addition to mismatching tissues, this automated outlier process can also detect samples which have low sequencing or DNA quality, or have global dysregulation of methylation due to mutations in chromatin or methylation regulation genes. Further analysis to understand why a sample is a global outlier by a user is warranted! 

In future releases of Metafora, automated outlier detection will be able to be skipped to force methylation outlier calling on all samples (though not recommended as some global outlier samples can have 10k+ outliers which makes interpretation nearly impossible, and will significantly slow down runtimes). Also it will be possible to skip the automated detection and simply provide a list of outlier samples to exclude based on manual inspection of the data (or these samples can simply be removed from the input table). 

### Calling Methylation Outliers on Sex Chromosomes

To call methylation outliers on the sex chromosomes, we take a sex-stratified approach to bypass the complexities of dealing with varying sex chromosome copy number. So, X-chromosome outliers are called separately in XX and in XY individuals, and Y-chromosome outliers will only be called in XY individuals. Sex chromosome copy number is automatically estimated based on read depth over both chromosome X and chromosome Y. If a sample does not cluster well with either XX or XY individuals, because they, for instance, have a rare sex chromosome combination (X, XXX, XYY, etc.), they will be flagged and excluded from analyses.

<img width="937" alt="image" src="https://github.com/user-attachments/assets/db392990-fbef-4d61-8c5c-e6eeab1f0167">

## How large a sample size is needed to run Metafora?
As with all outlier analyses, the larger the reference population, the more power you have to detect robust outliers. We recommend at least **30 samples** per tissue to start using Metafora to call outliers, though the more samples the better. Sample size could be increased by using publically available sequencing data of a matched tissue (and adjusting for sequencing modality as a `Batch` covariate if public data are coming from a different sequencing type). 

By decreasing the default MAX_ABS_ZSCORE threshold, it could be possible to detect outliers with fewer samples, but these outliers will be less reliable, and one might also find more false positives. 

For sex-stratified outliers on the X and Y chromosome, it is recommended to have at least 20-30 samples per sex as well. 

## Running Metafora

Metafora is written as a snakemake workflow wrapper around a series of R scripts. To run, clone this github repository, enter the directory, and then update a config file to specify an input table of samples / methylation filepaths, reference genome, runtime parameters, and optional covariates file. Then the workflow can be run using the following snakemake command
```
snakemake -pr --snakefile Metafora.snakefile --configfile config.yml --use-conda --profile <cluster_profile> --jobs 100
```
Metafora was designed to be run parrallelized across many jobs on a slurm controller or other HPC job scheduler. This can be configured in snakemake by setting up a cluster profile. 



### Input
Input paths are supplied to the workflow in a tsv file provided in the `sample_table:` directive of the config file. An example input table is included: `example_input.sample_table.txt`

Input table is expected to have the following 4 columns: 
* `Sample_name` (Unique sample ID)
* `Tissue` (Name of tissue this sample was sequenced from i.e. WholeBlood, PBMC, LCL, Fibroblast, Brain, etc.)
* `Technology` (Long read methylation input type. Currently can be one of 3 values: [**PacBio**, **ONT**, **Metafora**]).
* `Input_file` (Path to input file corresponding to input type)

```
Sample_name  Tissue  Technology  Input_file
SampleA  WholeBlood  ONT  /path/to/SampleA.ONT.bam
SampleB WholeBlood PacBio /path/to/SampleB.PacBio.bam
SampleC Brain Metafora /path/to/SampleC.methylation.METAFORA_formatted.bed.gz
```

#### Input Types
For **PacBio**, expected input is a methylated bam file. As part of pre-processing, Metafora will run cpg-tools to get a methylation bed file and properly format this file for downstream analyses. For PacBio preprocessing, a path to a static binary of cpg-tools and a downloaded CpG methylation model must be downloaded and included in the config file under the `pb_cpg_tools:` directive. These can be downloaded from pb_cpg_tools github page. 

For **ONT**, expected input is a methylated bam file. As part of pre-processing, Metafora will run modBam2Bed on the methylated file to generate a bed file of methylation proportions. **ONT** is specified for nanopore genomes where methylation has been called with guppy or dorado and methylation is stored in the Bam file as MM/ML tags. For nanopolish called ONT methylation you must manually format to be consistent with the **Metafora** input type.

If using other tools or technologies to call methylation, it is possible to manually format methylation calls to be consistent with Metafora's expected formatting and input them with the input type **Metafora**. Metafora formatting requires a 5-column, tsv-separated bed file that is sorted by position, bgzipped, and tabix indexed. 5 columns should be `chromosome`, `start`, `end`, `depth`, and `beta` as seen in the example below. Using this formatting it would hypothetically be possible to call outliers from whole-genome bisulfite sequencing cohorts, though this has not been tested. 

```
chromosome      start   end     depth   beta
chr1    10468   10468   57      0.75439649122807
chr1    10470   10470   66      0.757554545454545
chr1    10483   10483   65      0.861544615384615
chr1    10488   10488   66      0.924259090909091
chr1    10492   10492   57      0.982456140350877
chr1    10496   10496   66      0.984854545454546
chr1    10524   10524   74      0.93245
chr1    10541   10541   70      0.942882857142857
chr1    10562   10562   70      0.899971428571429
```
### Reference Genome Specification
In addition to input bams, the user must also specificy a reference genome that was used for aligning bams / calling methylation. Reference is specified as a fasta in the `reference_fasta:` directive of the config yaml file. In addition to this reference genome, the user may also need to point to using the `reference_seqnames_table:` directive a table of chromosome seqnames mapped to whether they are an **autosome**, **sex_chromosome_X**, or **sex_chromosome_Y**. Metafora will only run on chromosomes listed in this file, so make sure seqnames exactly match between reference fasta and seqnames table. By default, Metafora is set up to run on GRCh38 where autosomes are named `chr1`, `chr2`, ... , `chr22` and sex chromosomes are `chrX` and `chrY`. For any other reference genome or version make sure this table is updated!! **sex_chormosome_X** and **sex_chromosome_Y** must both be specified or metafora will automatically run with `SKIP_SEX_CHROMSOME_ESTIMATION` turned on and will not call outliers on sex chromosomes. 

reference_seqnames_table example:
```
chr1 autosome
chr2 autosome
 . . .
chr20 autosome
chr21 autosome
chrX sex_chromosome_X
chrY sex_chromosome_Y
```

**Disclaimer**: Metafora has only been tested on build hg38, but scripts *should* be general enough to work on other genome builds as well. If you run into troubles, please reach out :) 

### Specifying Additional Covariates to Correct (Optional)
By default, Metafora controls for technical covariates, biological factors and batch effects that could confound outlier analysis. It acheives this by computing hidden factors using PCA of the methylation proprotion matrix. Hidden factors will capture the largest sources of variability between samples, usually coming from technical factors such as sequencing modality, sequencing depth, library quality, etc. Because the PCs are capturing this source of variation, we do not believe it is necesarry to explicity correct for them in the model, though we give users the options to explicity provide a matrix of covariates to correct for and regress out for methylation z-score estimation. 

To explicitly correct for covariates, one simply needs to specify a matrix of covariates in the `covariates:` directive of the config file. Covariates matrix should be all numeric with the exception of an optional `Batch` column that contains a factor representing a "Batch-like" variable. Batch need not be sequencing batch per se, as we have witnessed that, as a stable epigenetic mark, methylation is much less prone to batch effects resulting from different sequencing times/sites. `Batch`, instead, can be any group that could explain a major variation in the data set: ie. sequencing technology (PacBio vs ONT), flowcell chemistry (ONT_R9 vs ONT_R10), methylation/base calling model used (dorado sup 5mc vs dorado hac 5mc_5hmc), or some combination there of. The covariates matrix must also have a `Sample_name` column that matches exactly sample IDs present in the input data table. The rest of the columns should be **linearly independent** numeric columns of additional variables that the user wants to be corrected for when calling outliers, for instance, Age, cell type proprotions, culture time for cell lines, post-mortem interval, DNA Integrity (DIN), etc. An example covariates file is present in the github as `example_input.covariates.txt`. 

Sex should not be explicitly included in the covariates matrix, as Metafora automatically estimates sex chromosome copy number by measuring depth across sex chromosomes and then controls for this estimated sex during outlier calling, but the user can explicitly include sex if running with `SKIP_SEX_CHROMOSOME_ESTIMATION` turned on. 

### Output
Metafora generates two final output files per sample to represent the methylation outlier in the sample. These output files can be found in sample-level directories generated in the `output_dir` specified in the config file. 

`{Sample_name}.GRCh38.tissue_{Tissue}.METAFORA.outlier_regions.bed`: is a tsv file containing outlier region coordinates and statistics per line. 
```
seqnames start end width strand ID num.mark seg.mean startRow endRow seg_id pop_median delta zscore
chr2 1133646 1134450 21 * SampleA 21 4.69754522381825 1564 1584 chr2_SampleA_39 0.181064229024943 0.458492567091562 3.07107383009539
chr2 188033870 188034821 27 * SampleA 27 -3.48184554362126 27 53 chr2_SampleA_6361 0.833792189044445 -0.334283018312737 -3.39982706663906
chr2 241732357 241733011 21 * SampleA 21 4.55240254572716 7827 7847 chr2_SampleA_8192 0.477200448180219 0.451191637246917 4.63241419813897
chr5 177116845 177117161 31 * SampleA 31 4.03587755330864 2955 2985 chr5_SampleA_6638 0.0531266709887945 0.359516369485046 4.68047407558207
chr12 27709846 27710286 21 * SampleA 21 -4.1008044265898 612 632 chr12_SampleA_969 0.589461168278312 -0.367239640500534 -3.34552920325453
chr20 31546859 31548400 49 * SampleA 49 5.74149646908452 1727 1775 chr20_SampleA_1018 0.54845890348952 0.439128638276589 7.47994322352052
```

* `seqnames`,`start`,and `end` represent coordinates of outlier region
* `ID` marks sample ID this outlier was found in
* `num.mark` counts number of CpG sites present in this outlier region
* `seg.mean` represents mean of the population deviation score across the outlier region
* `seg_id` is a unique ID for each candidate outlier that was identified
* `pop_median` represents the median methylation level across all samples for this region
* `delta` is the effect size difference of the methylation outlier in this sample. How much lower/higher the outlier is compared to the `pop_median`
* `zscore` is the signed methylation zscore representing how significantly different methylation levels in this sample are compared to rest of population

`{Sample_name}.tissue_{Tissue}.METAFORA.outlier_regions.zscore.mat`: is a matrix containing z-scores for all samples across each outlier region. Columns represent samples while row represent outlier regions that can be cross-linked to the outlier_regions.bed file. Each cell represents methylation z-score for that samples in that methylation region.

In addition to these summary files, we also plot a pdf image of each methylation outlier in the `outlier_plots` directory. These figures show 1. the population mean methylation with confidence intervals showing standard error of that estimate, as well as the observed sample-level methlyation. 2. the population deviance score across the same region with red line highlighting the segmented candidate outlier region. 3. the population methylation z-score distribution after correction and normalization with the outlier sample highlighted in red. An example plot for an outlier is included below: 

<img width="300" alt="image" src="https://github.com/user-attachments/assets/383de49a-5755-49d2-af5b-2af2585de0f0">

In addition to methylation outlier regions per sample, Metafora also outputs a tissue methylation reference, the estimate of mean methylation proportions at each cpg as well as a matrix of beta and depth values per chrom across the whole cohort in the `Population_Methylation.tissue_{tissue}/` directory. Also data tables of estimated hidden factors and sex chromosome copy number / sex assignments can be found as well as plots for correlation outlier, global autosomal PCs, and sex chromosome copy numbers in the `Global_Methylation_PCA_tissue_{tissue}` directory.

### Dependencies
To run Metafora, one must install snakemake and (ideally) mamba through conda. Conda environments for each step are then specified for each snakemake rule and will be automatically installed and activated when running snakemake with the `--use-conda` parameter. 

```
conda create -n snakemake -c conda-forge snakemake mamba
conda activate snakemake

#run metafora with following command after updating config.yml file
snakemake -pr --snakefile Metafora.snakefile --configfile config.yml --use-conda --profile <cluster_profile> --jobs 100
```

## Additional Workflow Parameters
The following parameters can be set in the `params:` directive of the config file, giving the user additional control on how to run Metafora
* `MAX_DEPTH`: Maxmimum read depth that methylation estimates are capped at (Default: **30**). Metafora caps read depth at this value to ensure one sample is unduly influencing mean methylation profile because of higher sequencing depth. Capping read depth also ensures p-values from CpG methylation proportion hypothesis testing are not being inflated.
* `MIN_SEG_SIZE`: Minimum number of cpgs within a segment during the segmentation algorithm (Default: **20**). This value applies both to the mean profile segmentation to calculate hidden factors, as well as the deviance score segmentation to find candidate methylation outlier regions. By default, we set MIN_SEG_SIZE to 20 because we believe larger outlier blocks of CpGs tend to be easier to interpret, but this value can be lowered to identify more shorter-sized outliers.
* `MIN_ABS_ZSCORE`: Threshold for absolute methylation z-score above which an outlier is called (Default: **3**). Increasing threshold will lead to more stringent outliers, while decreasing threshold could be helpful with small sample sizes where there is more limited power. User can always start with a low-threshold and manually filter outliers based on reported z-score in the output tables.
* `MIN_ABS_DELTA`: Threshold of absolute effect difference **delta** between outlier sample and population median methylation (Default: **0.25**). Outliers with absolute deltas lower than this value will be filtered out and not reported in the output tables. When looking for methlyation outliers, we are usually interested in large effect differences corresponding often to allele-specific effects. While this value is arbitrary, we have found it effective at filtering out the significant outliers that are not substantively that different in terms of methylation levels differences.
* `SKIP_SEX_CHROM_ESTIMATION`: Whether to run sex chromosome estimation and call methylation outliers on sex chromosomes (Default: **FALSE**). Set to TRUE if user wants to skip this estimation step and only call methylation outliers on autosomes. This will automatically be set to TRUE if the user did not specify both a **sex_chromosome_X** and **sex_chromosome_Y** in the `reference_seqnames_table`. 

# Disclaimer

**Warning**: Metafora is still actively being developed, if you are testing it out and run into any problems or have specific feature requests, please reach out to me or submit issue! tannerj@stanford.edu :) 
