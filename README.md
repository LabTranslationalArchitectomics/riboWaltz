![riboWaltz](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/ribowaltz.png)

# riboWaltz

R package for calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data
------------------------------------------------------------------------

## Overview

Ribosome profiling is a powerful technique used to study translation at the genome-wide level, generating unique information concerning ribosome positions along RNAs. Optimal localization of ribosomes requires the proper identification of the ribosome P-site in each ribosome protected fragment, a crucial step to determine trinucleotide periodicity of translating ribosomes and draw correct conclusions concerning where ribosomes are located. To determine the P-site within ribosome footprints at nucleotide resolution, the precise estimation of its offset with respect to the protected fragment is necessary. 

__riboWaltz__ is an R package for calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data. Taking advantage of a two-step algorithm where offset information is passed through populations of reads with different length in order to maximize offset coherence, __riboWaltz__ computes with high precision the P-site offset. __riboWaltz__ also provides a variety of graphical representations, laying the foundations for further accurate RiboSeq analyses and improved interpretation of positional information.

------------------------------------------------------------------------

## License

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg?maxAge=2592000?style=flat)](https://opensource.org/licenses/MIT)

------------------------------------------------------------------------

## Dependencies

__riboWaltz__ requires R version >= 3.3.0 and the following packages to work:

* CRAN
    + data.table (>= 1.10.4.3)
    + ggplot2 (>= 2.2.1)
    + ggrepel (>= 0.6.5)
* Bioconductor
    + Biostrings (>= 2.46.0)
    + GenomicAlignments (>= 1.14.1)
    + GenomicFeatures (>= 1.24.5)
    + GenomicRanges (>= 1.24.3)
    + IRanges (>= 2.12.0)

All the dependencies are automatically installed running the code in the next section.

------------------------------------------------------------------------

## Installation

To install __riboWaltz__ directly from GitHub the _devtools_ package is required. If it is not installed, run
    
    install.packages("devtools")
	
Otherwise, load _devtools_ and install __riboWaltz__ typing
	
	library(devtools)
    install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
  
Please note: to install __riboWaltz__ generating the vignette replace the last command with:
  
    install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE, build_vignettes = TRUE)

------------------------------------------------------------------------

## Loading

To load __riboWaltz__ type

	library(riboWaltz)
	
------------------------------------------------------------------------

## Getting help

The following section explains how to make use of riboWalz by introducing all the functions included in the package and reporting most of the data structures and graphical outputs obtained with the default options. For additional examples and further details about the usage of each parameter in the functions please refer to their documenation by running
 
	?function_name

or

	help(package = riboWaltz)
 
in the R console.
 
A complete reference manual can be found [here](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/ReferenceManual.pdf).   

------------------------------------------------------------------------

## Usage

### Acquiring input files

  One or more BAM files can be read and converted into a list of data tables or into a GRangesList object by running the `bamtolist` function. To run `bamtolist`, only the path to the BAM file(s) and an annotation file (see below for additional information) are required. For convenience, it is suggested to to rename the BAM files before their acquisition in order to maintain the same nomenclature of the samples through the whole analysis. However, it is also possible to assign the desired name to the samples thanks to the *list_name* option. Pay attention to the order in which their are provided: the first string is assigned to the first file, the second string to the second one and so on. If the original BAM files come from an alignment on transcripts (intended as an alignment based on a reference FASTA of all the transcript sequences), no reads associated to the negative strand should be present and they are therefore removed. Moreover, multiple options for treating the read lengths are available (see next chapter for more information). The syntax of the `bamtolist` function is as follow:
  
	reads_list <- bamtolist(bamfolder = path_to_bam, annotation = annotation_file)
	
  The resulting data structures contain for each read the name of the reference transcript on which it aligns, the leftmost and rightmost position of the read and its length. Two additional columns are also attached, reporting the leftmost and rightmost position of the CDS of the reference sequence with respect to its 1st nuclotide. An example of the final output of the `bedtolist` function is provided by the *reads_list* dataset included in the package, that contains the data for a sample called *Samp1* (a subset of the original dataset mainly composed by reads aligning on the translation initiation site is provided. Please contact the authors for more information). Here the first rows:
  
  |  transcript  |  end5  |  end3  |  length  |  start_pos  |  stop_pos  |
  |:------:|:-----:|:------:|:------:|:------:|:------:|
  |  ENSMUST00000000001.4  |  92  |  119  |  28  |  142  |  1206  |
  |  ENSMUST00000000001.4  |  94  |  122  |  29  |  142  |  1206  |
  |  ENSMUST00000000001.4  |  131  |  159  |  29  |  142  |  1206  |
  |  ENSMUST00000000001.4  |  132  |  158  |  27  |  142  |  1206  |
  |  ENSMUST00000000001.4  |  140  |  167  |  28  |  142  |  1206  |
  |  ENSMUST00000000001.4  |  142  |  170  |  29  |  142  |  1206  |
  
  Alternatively, the BAM file can be first converted into BED files and then into a list of data tables or into a GRangesList object through the functions `bamtobed` and `bedtolist`. The `bamtobed` function calls the bamtobed utility of the BEDTools suite (for the installation follow the instructions at http://bedtools.readthedocs.io/en/latest/content/installation.html). The BEDTools suite has been developed for command line environments, so `bamtobed` can be only run on UNIX, LINUX and Apple OS X operating systems. Once the bed files are produced, it is possible to switch to any other machine without further restrictions. To run `bamtobed`, only the path to the BAM file(s) is required, possibly coupled with the location of the directory where the BED files should be saved. The command for running `bamtobed` appears as follow

    bamtobed(bamfolder=path_to_bam, bedfolder = path_to_bed)

   Then, the `bedtolist` function loads and reads the BED files merging them into a list. `bedtolist` only requires the path to the BED file(s) and an annotation file. The syntax of this function is the same as for `bamtolist`:

    reads_list <- bedtolist(bedfolder = path_to_bed, annotation = annotation_file)

#### Selection of read lengths

  Different lengths of ribosome protected fragments may derive from alternative ribosome conformations. Therefore, the researcher should be free to modify the tolerance for the selection of the read length according to the aim of the experiment. For this reason, __riboWaltz__ has multiple options for treating read lengths specified by the parameter *length_filter_mode* included in both `bamtolist` (used in the examples below) and `bedtolist`:
  
1. all read lengths are included in the analysis (all-inclusive mode, default)
    
	reads_list <- bamtolist(bamfolder = path_to_bam, annotation = annotation_file,
							length_filter_mode = "none")

2. only read lengths specified by the user are included (manual mode)

	reads_list <- bamtolist(bamfolder = path_to_bam, annotation = annotation_file,
							length_filter_mode = "custom", length_filter_vector = 27:30)

3. only read lengths satisfying a periodicity threshold are included in the analysis (periodicity threshold mode). The user can change the desired threshold (the default is 50%). This mode enables the removal of all the reads without periodicity.

	reads_list <- bamtolist(bamfolder = path_to_bam, annotation = annotation_file,
							length_filter_mode = "periodicity", periodicity_threshold = 50)

For additional details please referes to the documentation provided by ?bamtolist or ?bedtolist.
	
### Annotation data table
  
  A reference annotation file is required to attach to the data tables two additional columns containing the position of the start and the stop codons with respect to the beginning of the transcript, two crucial information for localizing the reads within the three region of the transcrips (5' UTR, the CDS and the 3' UTR) and computing the P-site offsets. To do this, the annotation file must contain at least five columns reporting the name of the transcripts and the length of the whole transcript and of the annotated 5' UTR, the CDS and the 3' UTR. The annotation file can be either provided by the user or generated starting from a GTF file by using the `create_annotation` function. In the latter case, the name of the transript in the annotation data table are composed by the ENST ID and version, dot separated. During the generation of the annotation file, the input GTF is converted in a TxDb object and than in a data table. Therefore, a TxDb object can also be directly used as input of `create_annotation`. Here an example of the output:
  
  |  transcript  |  l_tr  |  l_utr5  |  l_cds  |  l_utr3  |
  |:------:|:---------:|:------:|:---------:|:------:|
  |  ENSMUST00000000001.4  |  3262  |  141  |  1065  |  2056  |
  |  ENSMUST00000000003.11  |  902  |  140  |  525  |  237  |
  |  ENSMUST00000000010.8  |  2574  |  85  |  753  |  1736  |
  |  ENSMUST00000000028.11  |  2143  |  313  |  1701  |  129  |
  |  ENSMUST00000000033.9  |  3708  |  115  |  543  |  3050  |
  |  ENSMUST00000000049.5  |  1190  |  51  |  1038  |  101  |
  
### Overview of the data

  Two graphical outputs can be produced before the identification of the P-site offset, in order to have an overview of the whole read sets. The first plot shows the distribution of the length of the reads for a specified sample and can be exploited to identify one or more populaton of reads i.e. one or more conformation of the ribosomes bound to the mRNAs. This plot is provided by `rlength_distr`. This function, as all the other contained in __riboWaltz__ producing a graphical output, returns a list containing both the data to generate the plot and the plot itself. For more details about the data tables for generating the plot and for an example of their structure please refer to the vignette of the package and to the documentation of its functions
  
  Note that a wide range of read lengths can make hard to read the plot. This issue can be easily solved specifying by the *cl* option a confidence level that restricts the distribution to a more narrow range of lengths

    example_length_dist_zoom <- rlength_distr(reads_list, sample = "Samp1", cl = 99)
    example_length_dist_zoom[["plot"]]
<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_length_dist_zoom.png" width="300" />
</p>

  The second plot consists of 4 metaheatmaps that show the abundance of the 3' the and 5' end of the reads mapping around the start and the stop codons, stratified by their length. This plot, generated by `rends_heat`, is particularly useful for understanding which extremity of the reads is the best choice for the computation of the P-site offset. Even if __riboWaltz__ is able to automatically recognize the best read end to use for the P-site identification, in some cases it may be necessary to provide this information. Also in this case it is possible to restrict the output to a subset of read lengths specified by a confidence level *cl*.
  
  Note that as for all the metaprofiles generated by __riboWaltz__, the data table associated to the plot contains three main columns: the first one indicating the distance (in nucleotides) from either the start or the stop codon, the second one reporting the value of the plot corresponding to that position and the third one specifying if the line of the data table refers either to the inital or to the final region of the coding sequence. If the metaprofile is stratified for the length of the reads, an additional column with this information is present. Only for `rends_heat' the data table includes an additional column that specifies the extremity of the reads involved in the plot.

    example_ends_heatmap <- rends_heat(reads_list, mm81cdna, sample = "Samp1", cl = 85,
                                       utr5l = 25, cdsl = 40, utr3l = 25)
    example_ends_heatmap[["plot"]]
![example_ends_heatmap](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_ends_heatmap.png)

  In our example, looking at the reads aligning around the translation initiation site (TIS) it is clearly visible a different trend of the signal coming from the 5' and the 3' extremity. In fact, the distance between the 5' end and the TIS varies depending on the read length (shorter the reads, closer are the 5' ends to TIS), while the 3' end often alignes on a specific nucleotide. This may suggest that the more stable extremity (i.e. the best option for the identification of the P-site offset) is the latter one and this information can be passed to the function `psite`. Nevertheless, in our example we are going to employ the automatic selection of the extremity to show how __riboWaltz__ works without any restriction (see below).
  
### P-site offset

  Starting from the list produced by `bedtolist` the P-site offsets can be computed using the function `psite`. This function compute the P-site identification starting from the reads that align on the start codon of any transcript, exploiting the knowledge that their associated P-sites corresponds to the translation start site. However, it is possible to remove reads which extremities are too close to the start codon by specifying the *flanking* option.
  
  `psite` processes one sample at a time, printing on the screen the extremity chosen for the computation of the P-site offsets and the value upon which the correction is based.

    psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")

  The result is a data table containing for all the read lengths of each sample the percentage of reads in the whole dataset and the percentage of reads aligning on the start codon (if any; used for the computation of the P-site offsets as described above). The data table also reports the distance of the P-site from the two extremities of the reads before and after the correction step. An additional column contains the name of the sample.

  |  length  |  total_percentage  |  start_percentage  |  around_start  |  offset_from_5  |  offset_from_3  |  adj_offset_from_5  |  adj_offset_from_3  |  sample  |
  |:------:|:---------:|:------:|:---------:|:------:|:------:|:------:|:------:|:------:|
  |  19  |  0.715  |  0.110  |  T  |  12  |  6  |  8  |  10  |  Samp1  |
  |  20  |  0.811  |  0.194  |  T  |  12  |  7  |  7  |  12  |  Samp1  |
  |  21  |  0.992  |  0.223  |  T  |  15  |  5  |  7  |  13  |  Samp1  |
  |  22  |  0.921  |  0.364  |  T  |  6  |  15  |  10  |  11  |  Samp1  |
  |  23  |  1.112  |  0.974  |  T  |  6  |  16  |  12  |  10  |  Samp1  |
  |  24  |  1.892  |  1.722  |  T  |  7  |  16  |  7  |  16  |  Samp1  |
    
  For every read length of each sample a plot of the ribosome occupancy profile for the 5' and the 3’ extremity around the start codon is produced. The optimal offsets (dotted black line) as well as the inferred offsets before (dashed vertical lines) after the correction (continuous vertical lines) are reported. The regions used for their computation (depending on the *flanking* option) are shaded. Here two examples for reads of 28 and 31 nucleotides:

<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/meta_psite_length28.png" width="750" />
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/meta_psite_length31.png" width="750" />
</p>

  The initial dataset must be updated with new information resulting from the identification of the P-site offset. The function `psite_info` to attaches to the exsisting data tables the localization of the P-site along the transcript and its position with respect to the start and stop codons. The associated region of the transcript (5' UTR, CDS, 3' UTR) and, optionally, the sequence of the triplet covered by the P-site are also added. These information are required for facilitating further analyses and in particular to compute the number of in-frame Psites, the percentage of P-sites falling in the three transcript regions, verify the trinucleotide periodicity of the reads along the coding sequence and generating metplots, as discussed below.

    reads_psite_list <- psite_info(reads_list, psite_offset)
	
  Here an example of the resulting data table for one sample:
  
  |  transcript  |  end5  |  psite  |  end3 |  length  |  start_pos  |  stop_pos  |  psite_from_start  |  psite_from_stop  |  psite_region  |
  |:------:|:-----:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|
  |  ENSMUST00000000001.4  |  92  |  103  |  119  |  28  |  142  |  1206  |  -39  |  -1103  |  5utr  |
  |  ENSMUST00000000001.4  |  94  |  106  |  122  |  29  |  142  |  1206  |  -36  |  -1100  |  5utr  |
  |  ENSMUST00000000001.4  |  131  |  143  |  159  |  29  |  142  |  1206  |  1  |  -1063  |  cds  |
  |  ENSMUST00000000001.4  |  132  |  142  |  158  |  27  |  142  |  1206  |  0  |  -1064  |  cds  |
  |  ENSMUST00000000001.4  |  140  |  151  |  167  |  28  |  142  |  1206  |  9  |  -1055  |  cds  |
  |  ENSMUST00000000001.4  |  142  |  154  |  170  |  29  |  142  |  1206  |  12  |  -1052  |  cds  |

#### Sequence data

  Optionally, a file containing transcript or genome sequence information in FASTA format can be provided as input for `psite_info` to perform P-site specific codon sequence analysis. By specifying a genome build, the corresponding BSGenome object in R will be used for sequence retrieval. Genome sequences and genome builds must be coupled with a TxDb object or a GTF file providing genomic annotation. For further datails please refer to the documentation of `psite_info` running ?psite_info. 
  
### Codon coverage

  The updated dataset can be used to compute the coverage of the triplets. To this aim the function `codon_coverage` can be run, using the *psite* option to specify if the codon coverage should be based either on  number of read footprints per triplet or on the number of P-sites per triplet.
  
	coverage_dt <- codon_coverage(reads_psite_list, mm81cdna, psite = FALSE)

  The result is a data table as the following one
  
  |  transcript  |  start  |  end  |  start_dist |  stop_dist  |  region  |  Samp1  |
  |:------:|:-----:|:------:|:------:|:------:|:------:|:------:|
  |  ENSMUST00000000001.4  |  90  |  93  |  -17  |  -371  |  5utr  |  1  |
  |  ENSMUST00000000001.4  |  93  |  96  |  -16  |  -370  |  5utr  |  2  |
  |  ENSMUST00000000001.4  |  96  |  99  |  -15  |  -369  |  5utr  |  2  |
  |  ENSMUST00000000001.4  |  99  |  102  |  -14  |  -368  |  5utr  |  2  |
  |  ENSMUST00000000001.4  |  102  |  105  |  -13  |  -367  |  5utr  |  2  |
  |  ENSMUST00000000001.4  |  105  |  108  |  -12  |  -366  |  5utr  |  2  |
  
### P-sites per region

  The dataset containing the position of the identified P-sites and the associated information can also be used to compute the percentage of P-sites falling in the three annotated regions of the transcripts (5' UTR, CDS and 3'UTR) expoiting the `region_psite` function that generates a barplot of the resulting values. Moreover, the function calculates and plots the percentage of region length (reported in column "RNAs"). 

	example_psite_region <- region_psite(reads_psite_list, mm81cdna, sample = "Samp1")
	example_psite_region[["plot"]]
<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_psite_per_region.png" width="300" />
</p>

### In-frame P-sites

  For generating a data tables containing, for each transcript, the number of ribosome protected fragments with in-frame P-site mapping on the CDS the function `psite_per_cds` can be exploited. The resulting data table can be emploied to estimate transcript-specific translation levels and perform differential analysis comparing multiple conditions.
	
	psite_cds <- psite_per_cds(reads_psite_list, mm81cdna)
	
  The result is a data table as the following one
  
  |  transcript  |  length  |  Samp1  |
  |:------:|:------:|:------:|
  |  ENSMUST00000000001.4  |  1065  |  31  |
  |  ENSMUST00000000003.11  |  525  |  0  |
  |  ENSMUST00000000010.8  |  753  |  0  |
  |  ENSMUST00000000028.11  |  1701  |  2  |
  |  ENSMUST00000000033.9  |  543  |  26  |
  |  ENSMUST00000000049.5  |  1038  |  0  |

### 3-nucleotide periodicity

  To visualize the percentage of identified P-sites (i.e. ribosomes) in the three reading frames along the coding sequence, the functions `frame_psite_length` and `frame_psite` can be exploited. Both of them compute how many ribosomes are in the three frames for the 5' UTR, the CDS and the 3' UTR, with the following difference: the first one divides the results depending on the read length, while the second one works handling the reads all together.

    example_frames_stratified <- frame_psite_length(reads_psite_list, sample = "Samp1",
                                                    region = "all", cl = 90)
    example_frames_stratified[["plot"]]
<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_frames_stratified.png" width="550" />
</p>

    example_frames <- frame_psite(reads_psite_list, sample = "Samp1", region = "all")
    example_frames[["plot"]]
<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_frames.png" width="550" />
</p>

  In both cases it is possible to specified a set of samples of interest, a range of read lengths and to display the results for either all the three regions of the transcript or just one of them. Depending on this choices the box containing the plots are differently arranged, to optimise the organization and the visualization of the data.

### Metaplots

  The `metaprofile_psite` function generates metaprofiles (sum of single, transcript-specific profiles) based on the P-sites previously identified. This plots are useful to verify the so-called 3-nt periodicity of ribosomes along transcripts at genome-wide scale. The contribution from many replicates can be combined in a single plot, taking into account possible scale factors coming from any normalization of the data chosen by the user. It is possible to use the whole transcriptome (as in the example below), restrict the analysis to a subset of transcripts and even look at single RNAs.

    example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
                                             utr5l = 20, cdsl = 40, utr3l = 20,
											 plot_title = "auto")
    example_metaprofile[["plot"]]
![example_metaprofile](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_metaprofile.png)

  The `metaprofile_psite` utility also provides the *length_range* option, which allows to plot metaprofiles employing specified sub-populations of reads, depending on their length. Here an example using reads of 28 nucleotides.

    example_metaprofile_28 <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
                                                length_range = 28, utr5l = 20, cdsl = 40,
                                                utr3l = 20, plot_title = "auto")
    example_metaprofile_28[["plot"]]
![example_metaprofile_28](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_metaprofile_28.png)

  Another way to visualize the previous data is a metahaetmap where the intensity of the colour depends on the abundance of P-sites in a specific position (corresponding to the height of the line for the metaprofiles). This representation is the optimal choice in order to compare multiple samples or look at the behavior of the data if different populations of reads are considered. The `metaheatmap_psite` function provides a collection of metaheatmaps generated using the same set of transcripts while the reads are associated to either a variety of biological conditions or subsets of the same dataset.

  To show how it works, let's suppose we want to check how the reads of 28 nucleotdes (the most frequent in our example) behave with respect to the whole dataset. In principle it would be sufficient to compare the two metaprofiles displayed above, but in order to have a better view of the data we will exploit the `metaheatmap_psite` function. The first step is to create a list of data tables containing the data of interest:

    comparison_dt <- list()
    comparison_dt[["subsample_28nt"]] <- reads_psite_list[["Samp1"]][length == 28]
    comparison_dt[["whole_sample"]] <- reads_psite_list[["Samp1"]]

  Then the list with the names of the data tables of interest can be defined. Pay attention: here we don't have any replicate, thus each element of *names_list* is a vector with just one string.

    names_list <- list("Only_28" = c("subsample_28nt"),
                       "All" = c("whole_sample"))

  At this point it is sufficient to run `metaheatmap_psite` using the data tables and the list of names previously defined.

    example_metaheatmap <- metaheatmap_psite(comparison_dt, mm81cdna, sample = names_list,
                                             utr5l = 20, cdsl = 40, utr3l = 20, log = F)
    example_metaheatmap[["plot"]]
![example_metaheatmap](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_metaheatmap.png)

### Codon usage

  To understand what codons display higher or lower ribosome density, __riboWaltz__ provides the user with the analysis of the empirical codon usage, i.e. the frequency of in-frame P-sites along the coding sequence codon by codon, normalized for the frequency in sequences of each codon. The empirical condon usage is provided by the `codon_usage_psite` function which also returns a bar plot reporting the computed values, highlighting the start and the stop codon and labeling each bar with the corresponding amino acid. To this aim, the path to the fasta file employed during the the alignment step must be provided by the user through the *fastapath* option. The syntax for generating the bar plot is the following (note that due to its dimension, the fasta file used by the author is not included among the example data of the package).
  
   	codon_usage_barplot <- codon_usage_psite(reads_psite_list, mm81cdna, sample = "Samp1",
						   fastapath = path_to_fasta)
![codon_usage_barplot](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/codon_usage_barplot.png)
  
  To unravel possible defects in ribosome elongation at specific codons or aa-tRNAs use is it possible to exploit `codon_usage_psite` to compare empirical usage from two conditions or organisms or to asses potenial differences between the empirical codon usage and the more diffused codon usage bias based on codon frequencies. To this aim a set of 64 values can be passed to the function by the user through the option *codon_usage*. The structure of the required data table (here called *cub_mouse* and reporting the codon usage bias in mouse downloaded from http://www.kazusa.or.jp/codon) is as follow:

  |  codon  |  usage_index  |
  |:------:|:------:|
  |  UUU  |  17.2  |
  |  UCU  |  16.2  |
  |  UAU  |  12.2  |
  |  UGU  |  11.4  |
  |  UUC  |  21.8  |
  |  UCC  |  18.1  |
  
  If such a data table is provided, `codon_usage_psite` returns a second graphical output: a scatter plot where each codon is represented by a dot. The following image shows the comparison between the empirical codon usage reported in the previous figure and the codon usage bias in mouse contained in the data table *cub_mouse* (as for the fasta file, cub_mouse is not included among the example data of the package):
  
	codon_usage_scatter <- codon_usage_psite(reads_psite_list, mm81cdna, sample = "Samp1",
						   fastapath = path_to_fasta, codon_usage = cub_mouse)
<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/codon_usage_scatter.png" width="300" />
</p>

------------------------------------------------------------------------

### Contact

fabio.lauria@unitn.it

t.tebaldi@unitn.it

gabriella.viero@cnr.it
