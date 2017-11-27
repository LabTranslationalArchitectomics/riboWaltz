![riboWaltz](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/ribowaltz.png)

# riboWaltz

R package for calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data
------------------------------------------------------------------------

### Overview

Ribosome profiling is a powerful technique used to study translation at the genome-wide level, generating unique information concerning ribosome positions along RNAs. Optimal localization of ribosomes requires the proper identification of the ribosome P-site in each ribosome protected fragment, a crucial step to determine trinucleotide periodicity of translating ribosomes and draw correct conclusions concerning where ribosomes are located. To determine the P-site within ribosome footprints at nucleotide resolution, the precise estimation of its offset with respect to the protected fragment is necessary. 

__riboWaltz__ is an R package for calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data. Taking advantage of a two-step algorithm where offset information is passed through populations of reads with different length in order to maximize offset coherence, __riboWaltz__ computes with high precision the P-site offset. __riboWaltz__ also provides a variety of graphical representations, laying the foundations for further accurate RiboSeq analyses and better interpretation of positional information.

------------------------------------------------------------------------

### License

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg?maxAge=2592000?style=flat)](https://opensource.org/licenses/MIT)

------------------------------------------------------------------------

### Installation

To install __riboWaltz__, directly from GitHub: 
    
    install.packages("devtools")
    library("devtools")
    install_github("LabTranslationalArchitectomics/riboWaltz")
  
To install __riboWaltz__ and generate the vignette, substitute the latter command with:
  
    install_github("LabTranslationalArchitectomics/riboWaltz", build_vignettes=TRUE)


------------------------------------------------------------------------

### Usage

#### Aquiring input files

  One or more BAM files can be read and converted into BED files by using the `bamtobed` function which in turn calls the bamtobed utility of the BEDTools suite (for the installation follow the instructions at http://bedtools.readthedocs.io/en/latest/content/installation.html). The BEDTools suite has been developed for command line environments, so this step of __riboWaltz__ can be only run on UNIX, LINUX and Apple OS X operating systems. Once the bed files are produced, it is possible to switch to any other machine without further restrictions.

To run `bamtobed`, only the path to the BAM file(s) is required, possibly coupled with the location of the output directory. For convenience, it is suggested to to rename the BAM files before their acquisition, in order to maintain the same nomenclature of the samples through the whole process (however, it is possible to modify them while running `bedtolist`, as explained below).

    bamtobed(bamfolder=path_to_bam, bedfolder=path_to_bed)

  This worked example is base on BAM files aligned to a GENCODE transcriptome reference. The resulting files will contain, for each read, the name of the aligned transcript, the start and end coordinate, its length and the associated strand. The BED files will appear like the following (the header of the columns are here included for clarity):

  |  transcript  |  start  |  end  |  length  |  strand  |
  |:------:|:-----:|:------:|:------:|:------:|
  |  ENSMUST00000000001.4  |  92  |  119  |  28  |  +  |
  |  ENSMUST00000000001.4  |  94  |  122  |  29  |  +  |
  |  ENSMUST00000000001.4  |  131  |  159  |  29  |  +  |
  |  ENSMUST00000000001.4  |  132  |  158  |  27  |  +  |
  |  ENSMUST00000000001.4  |  140  |  167  |  28  |  +  |
  |  ENSMUST00000000001.4  |  142  |  170  |  29  |  +  |

   The next step loads and reads the BED files, merging them in a list. The `bedtolist` function only requires the path to the BED files and an annotation file. Moreover, the *list_name* option allows to assign the desired name to the samples. Pay attention to the order in which their are provided: the first string will be assigned to the first file, the second string to the second one and so on.

    reads_list <- bedtolist(bedfolder=path_to_bed, annotation=annotation_file)
	
  In case the original BAM files come from an alignment on transcripts, the reads associated to the negative strand should be present in a low percentage, and they will be removed. An example of the final output of the `bedtolist` function is provided by the *reads_list* dataset included in the package, that contains the data for a single sample called *Samp1* (a subset of the original data is here provided. Please contact the authors for the whole datset).

#### Annotation data frame
  
  A reference annotation file is required to attach to the data frames two additional columns containing the position of the start and the stop codons with respect to the beginning of the transcript. To do this, the annotation file must contain at least five columns reporting the name of the transcripts and the length of the whole transcript and of the annotated 5' UTR, the CDS and the 3' UTR. Here an example:
  
  |  transcript  |  l_tr  |  l_utr5  |  l_cds  |  l_utr3  |
  |:------:|:---------:|:------:|:---------:|:------:|
  |  ENSMUST00000000001.4  |  3262  |  141  |  1065  |  2056  |
  |  ENSMUST00000000003.11  |  902  |  140  |  525  |  237  |
  |  ENSMUST00000000010.8  |  2574  |  85  |  753  |  1736  |
  |  ENSMUST00000000028.11  |  2143  |  313  |  1701  |  129  |
  |  ENSMUST00000000033.9  |  3708  |  115  |  543  |  3050  |
  |  ENSMUST00000000049.5  |  1190  |  51  |  1038  |  101  |

  The annotation file can be either provided by the user, or generated starting from a GTF file by using the `create_annotation` function.
  
  Optionally, a third file containing transcript sequence information in FASTA format can be provided as input to perform P-site specific codon sequence analysis

#### Overview of the data

  Two graphical outputs can be produced before the identification of the P-site offset, in order to have an overview of the whole read sets. The first plot shows the distribution of the length of the reads for a specified sample, provided by `rlength_distr`. This function, as all the other contained in __riboWaltz__ producing a graphical output, returns a list containing both the data to generate the plot and the plot itself.
  
  Note that a wide range of read lengths can make hard to read the plot. This issue can be easily solved specifying by the *cl* option a confidence level that restricts the distribution to a more narrow range of lengths.

    example_length_dist_zoom <- rlength_distr(reads_list, sample="Samp1", cl=99)
    example_length_dist_zoom[["plot"]]
<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_length_dist_zoom.png" width="300" />
</p>

  The second plot consists of 4 metaheatmaps that show the abundance of the 3' the and 5' end of the reads mapping around the start and the stop codons, stratified by their length. Also in this case it is possible to restrict the output to a subset of read lengths specified by a confidence level *cl*.

    example_ends_heatmap <- rends_heat(reads_list, mm81cdna, sample="Samp1", cl=85,
                                       utr5l = 25, cdsl = 40, utr3l = 25)
    example_ends_heatmap[["plot"]]
![example_ends_heatmap](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_ends_heatmap.png)

  The latter plot is particularly useful for understanding which extremity of the reads is the best choice for the computation of the P-site offset. Even if __riboWaltz__ is able to automatically recognize the best read end to use for the P-site identification, in some cases it might be useful to provide directly this information. In our example, looking at the reads aligning around the translation initiation site (TIS), a different trend of the signal coming from the 5' and the 3' extremity is clearly visible. In fact, the distance between the 5' end and the TIS varies depending on the read length (shorter the reads, closer are the 5' ends to the TIS), while the 3' end often alignes on a specific nucleotide. This may suggest that the more stable extremity (i.e. the best option for the identification of the P-site offset) is the latter one and this information can be passed to the function `psite`. Nevertheless, in our example we are going to employ the automatic selection of the extremity to show how __riboWaltz__ works by default (see below).
  
#### P-site offset

  Starting from the list produced by `bedtolist` the P-site offsets can be computed using the function `psite`. This function computes the P-site starting from the reads that align on the start codon of any transcript, exploiting the knowledge that their associated P-sites correspond to the translation start site. However, it is possible to remove reads whose extremities are too close to the start codon by specifying the *flanking* option.
  
  `psite` processes one sample at a time, printing on the screen the extremity chosen for the computation of the P-site offsets and the value upon which the correction is based.

    psite_offset <- psite(reads_list, flanking = 6, extremity="auto")

  The result is a data frame containing for all the read lengths of each sample the percentage of reads in the whole dataset and the percentage of reads aligning on the start codon (if any). The data frame also reports the distance of the P-site from the two extremities of the reads before and after the correction step. An additional column contains the name of the sample.

  For every read length of each sample a plot of the ribosome occupancy profile for the 5' and the 3’ extremity around the start codon is produced. The optimal offsets (dotted black line) as well as the inferred offsets before (dashed vertical lines) and after the correction (continuous vertical lines) are reported. The regions used for their computation (depending on the *flanking* option) are shaded. Here two examples for reads of 28 and 31 nucleotides:

<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/meta_psite_length28.png" width="750" />
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/meta_psite_length31.png" width="750" />
</p>

  The initial dataset must be updated with new information resulting from the identification of the P-site offset. The function `psite_info` attaches to the exsisting data frames the localization of the P-site along the transcript and its position with respect to the start and stop codons. The associated region of the transcript (5' UTR, CDS, 3' UTR) and, optionally, the sequence of the triplet covered by the P-site are also added for facilitating further analyses.

    reads_psite_list <- psite_info(reads_list, psite_offset)
	
  Here an example of the resulting data frame for one sample:
  
  |  transcript  |  end5  |  psite  |  end3 |  length  |  start_pos  |  stop_pos  |  psite_from_start  |  psite_from_stop  |  psite_region  |
  |:------:|:-----:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|
  |  ENSMUST00000000001.4  |  92  |  103  |  119  |  28  |  142  |  1206  |  -39  |  -1103  |  5utr  |
  |  ENSMUST00000000001.4  |  94  |  106  |  122  |  29  |  142  |  1206  |  -36  |  -1100  |  5utr  |
  |  ENSMUST00000000001.4  |  131  |  143  |  159  |  29  |  142  |  1206  |  1  |  -1063  |  cds  |
  |  ENSMUST00000000001.4  |  132  |  142  |  158  |  27  |  142  |  1206  |  0  |  -1064  |  cds  |
  |  ENSMUST00000000001.4  |  140  |  151  |  167  |  28  |  142  |  1206  |  9  |  -1055  |  cds  |
  |  ENSMUST00000000001.4  |  142  |  154  |  170  |  29  |  142  |  1206  |  12  |  -1052  |  cds  |

  The updated dataset can be used as input to `psite_per_cds` for generating a list of data frames containing, for each transcript, the number of ribosome protected fragments with in-frame P-site mapping on the CDS. This data frame can be used to estimate transcript-specific translation levels and perform differential analysis comparing multiple conditions.
	
	psite_cds_list <- psite_per_cds(reads_psite_list, mm81cdna)
	
  The result is a data frame as the following:
  
  |  transcript  |  l_region  |  psite_count  |
  |:------:|:------:|:------:|
  |  ENSMUST00000000001.4  |  3262  |  31  |
  |  ENSMUST00000000003.11  |  525  |  0  |
  |  ENSMUST00000000010.8  |  753  |  0  |
  |  ENSMUST00000000028.11  |  1701  |  2  |
  |  ENSMUST00000000033.9  |  543  |  26  |
  |  ENSMUST00000000049.5  |  1038  |  0  |

#### 3-nucleotide periodicity

  To verify if the identified P-sites (i.e. ribosomes) are in the correct frame along the coding sequence, the functions `frame_psite_length` and `frame_psite` can be exploited. Both of them compute how many ribosomes are in the three frames for the 5' UTR, the CDS and the 3' UTR, with the following difference: the first one divides the results depending on the read length, while the second one works handling the reads all together.

    example_frames_stratified <- frame_psite_length(reads_psite_list, sample="Samp1",
                                                    region="all", cl=90)
    example_frames_stratified[["plot"]]
<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_frames_stratified.png" width="550" />
</p>

    example_frames <- frame_psite(reads_psite_list, sample="Samp1", region="all")
    example_frames[["plot"]]
<p align="center">
<img src="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_frames.png" width="550" />
</p>

  In both cases it is possible to specify a set of samples of interest, a range of read lengths and to display the results for either all the three regions of the transcript or just one of them. Depending on this choices the box containing the plots are differently arranged, to optimise the organization and the visualization of the data.

#### Metaplots

  The `metaprofile_psite` function generates metaprofiles (sum of single, transcript-specific profiles) based on the P-sites previously identified. These plots are useful to verify the so-called 3-nt periodicity of ribosomes along transcripts at genome-wide scale. The contribution from many replicates can be combined in a single plot, taking into account possible scale factors coming from any normalization of the data chosen by the user. It is possible to use the whole transcriptome (as in the example below), restrict the analysis to a subset of transcripts and even look at single RNAs.

    example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
                                             utr5l = 20, cdsl = 40, utr3l = 20)
    example_metaprofile[["plot"]]
![example_metaprofile](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_metaprofile.png)

  The `metaprofile_psite` utility also provides the *length_range* option, which allows to plot metaprofiles employing specified sub-populations of reads, depending on their length. Here an example using reads of 28 nucleotides.

    example_metaprofile_28 <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
                                                length_range = 28, utr5l = 20, cdsl = 40,
                                                utr3l = 20)
    example_metaprofile_28[["plot"]]
![example_metaprofile_28](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_metaprofile_28.png)

  Another way to visualize the previous data is a heatmap where the intensity of the colour depends on the abundance of P-sites in a specific position (corresponding to the height of the line for the meta-profiles). This representation is the optimal choice in order to compare multiple samples or look at the behavior of the data if different populations of reads are considered. The `metaheatmap_psite` function provides a collection of heatmaps generated using the same set of transcripts while the reads are associated to either a variety of biological conditions or subsets of the same dataset.

  To show how the function works, let's suppose we want to check how the reads of 28 nucleotdes (the most frequent length in our example) behave with respect to the whole dataset. In principle it would be sufficient to compare the two metaprofiles displayed above, but in order to have a better view of the data we will exploit the `metaheatmap_psite` function. The first step is to create a list of data frames containing the data of interest:

    comparison_df <- list()
    comparison_df[["subsample_28nt"]] <- subset(reads_psite_list[["Samp1"]], length == 28)
    comparison_df[["whole_sample"]] <- reads_psite_list[["Samp1"]]

Then the list with the names of the data frames of interest can be defined. Pay attention: here we don't have any replicate, thus each element of *names_list* is a vector with just one string:

    names_list <- list("Only_28" = c("subsample_28nt"),
                       "All" = c("whole_sample"))

At this point it is sufficient to run `metaheatmap_psite` using the data frames and the list of names previously defined.

    example_metaheatmap <- metaheatmap_psite(comparison_df, mm81cdna, sample = names_list,
                                             utr5l = 20, cdsl = 40, utr3l = 20, log=F)
    example_metaheatmap[["plot"]]
![example_metaheatmap](https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/vignettes/example_metaheatmap.png)

------------------------------------------------------------------------

### Contact

fabio.lauria@unitn.it

t.tebaldi@unitn.it

gabriella.viero@cnr.it
