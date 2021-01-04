## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")

## ---- eval = FALSE------------------------------------------------------------
#  library(devtools)
#  install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE,
#                 build_opts = c("--no-resave-data", "--no-manual"))

## ---- eval = TRUE, warning = FALSE--------------------------------------------
library(riboWaltz)

## ---- eval = FALSE------------------------------------------------------------
#  ?function_name

## ---- eval = FALSE------------------------------------------------------------
#  help(package = riboWaltz)

## -----------------------------------------------------------------------------
head(reads_list[["Samp1"]])

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom",
#                             length_filter_vector = 27:30)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  filtered_list <- length_filter(data = reads_list, length_filter_mode = "periodicity",
#                             periodicity_threshold = 70)

## ---- eval = TRUE, echo = FALSE, message = FALSE, warning = FALSE-------------
data(reads_list)
data(mm81cdna)
library(data.table)

## -----------------------------------------------------------------------------
head(mm81cdna)

## ---- echo = TRUE-------------------------------------------------------------
psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")

## ---- echo = TRUE-------------------------------------------------------------
head(psite_offset, 10)

## ---- out.width = '690px', fig.retina = NULL, echo = FALSE--------------------
knitr::include_graphics("meta_psite_length28.png")
knitr::include_graphics("meta_psite_length31.png")

## ---- echo = TRUE, eval = TRUE------------------------------------------------
reads_psite_list <- psite_info(reads_list, psite_offset)
head(reads_psite_list[["Samp1"]])

## ---- echo = FALSE, eval = TRUE-----------------------------------------------
head(reads_psite_list[["Samp1"]])

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  codon_coverage_example <- codon_coverage(reads_psite_list, mm81cdna, psite = FALSE)
#  head(codon_coverage_example[Samp1 > 0])

## ---- echo = FALSE, eval = TRUE-----------------------------------------------
codon_coverage_example <- data.table(transcript = rep("ENSMUST00000000001.4", 6),
                                     start = c(90, 93, 96, 99, 102, 105),
                                     end = c(93, 96, 99, 102, 105, 108),
                                     from_cds_start = -17:-12,
                                     from_cds_stop = -371:-366,
                                     region = rep("5utr", 6),
                                     Samp1 = c(1, 2, 2, 2, 2, 2))
head(codon_coverage_example)

## ---- echo = FALSE, eval = TRUE, include = FALSE------------------------------
cds_coverage_example <- cds_coverage(reads_psite_list, mm81cdna)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  cds_coverage_example <- cds_coverage(reads_psite_list, mm81cdna)
#  head(cds_coverage_example)

## ---- echo = FALSE, eval = TRUE-----------------------------------------------
head(cds_coverage_example)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_length_dist <- rlength_distr(reads_list, sample = "Samp1")
#  example_length_dist[["plot_Samp1"]]

## ---- out.width = '275px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_length_dist_zoom <- rlength_distr(reads_list, sample = "Samp1", cl = 99)
#  example_length_dist_zoom[["plot_Samp1"]]

## ---- out.width = '275px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_zoom.png")

## ---- echo = FALSE------------------------------------------------------------
    set.seed(10)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  reads_list[["Samp2"]] <- reads_list[["Samp1"]][sample(.N, 1000)]

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_length_dist_rep <-  rlength_distr(reads_list,
#                                            sample = list("Samp_avg" = c("Samp1", "Samp2")),
#                                            cl = 99, multisamples = "average",
#                                            colour = "gray70")
#  example_length_dist_rep[["plot"]]

## ---- out.width = '275px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_rep.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  comparison_list <- list()
#  comparison_list[["start_codon"]] <- reads_list[["Samp1"]][end5 <= cds_start & end3 >= cds_start]
#  comparison_list[["whole_sample"]] <- reads_list[["Samp1"]]

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  sample_list <- list("Only_start" = c("start_codon"),
#                     "All" = c("whole_sample"))

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_length_dist_split <-  rlength_distr(comparison_list,
#                                              sample = sample_list,
#                                              cl = 99, multisamples = "average",
#                                              plot_style = "split",
#                                              colour = c("dodgerblue", "gray70"))
#  example_length_dist_split[["plot"]]

## ---- out.width = '500px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_split.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_length_dist_dodged <-  rlength_distr(comparison_list,
#                                               sample = sample_list,
#                                               cl = 99, multisamples = "average",
#                                               plot_style = "dodged",
#                                               colour = c("dodgerblue", "gray70"))
#  example_length_dist_dodged[["plot"]]

## ---- out.width = '450px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_dodged.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_length_dist_mirrored <-  rlength_distr(comparison_list,
#                                                 sample = sample_list,
#                                                 cl = 99, multisamples = "average",
#                                                 plot_style = "mirrored",
#                                                 colour = c("dodgerblue", "gray70"))
#  example_length_dist_mirrored[["plot"]]

## ---- out.width = '450px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_mirrored.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  comparison_list <- list()
#  comparison_list[["start_codon"]] <- reads_list[["Samp1"]][end5 <= cds_start & end3 >= cds_start]
#  comparison_list[["whole_sample1"]] <- reads_list[["Samp1"]]
#  comparison_list[["whole_sample2"]] <- reads_list[["Samp2"]]
#  
#  sample_list <- list("Only_start" = c("start_codon"),
#                     "All" = c("whole_sample1", "whole_sample2"))

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_length_dist_split_rep <-  rlength_distr(comparison_list,
#                                                  sample = sample_list,
#                                                  cl = 99, multisamples = "average",
#                                                  plot_style = "split",
#                                                  colour = c("dodgerblue", "gray70"))
#  example_length_dist_split_rep[["plot"]]

## ---- out.width = '500px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_split_rep.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_ends_heatmap <- rends_heat(reads_list, mm81cdna, sample = "Samp1", cl = 85,
#                                     utr5l = 25, cdsl = 40, utr3l = 25)
#  example_ends_heatmap[["plot"]]

## ---- out.width = '700px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_ends_heatmap.png")

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  example_psite_region <- region_psite(reads_psite_list, mm81cdna, sample = "Samp1")
#  example_psite_region[["plot"]]

## ---- out.width = '260px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_psite_per_region.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_frames_stratified <- frame_psite_length(reads_psite_list, sample = "Samp1",
#                                                  region = "all", cl = 90)
#  example_frames_stratified[["plot"]]

## ---- out.width = '450px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_frames_stratified.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_frames <- frame_psite(reads_psite_list, sample = "Samp1", region = "all")
#  example_frames[["plot"]]

## ---- out.width = '450px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_frames.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
#                                           utr5l = 20, cdsl = 40, utr3l = 20,
#                                           plot_title = "sample.transcript")
#  example_metaprofile[["plot_Samp1"]]

## ---- out.width = '690px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaprofile.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_metaprofile_28 <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
#                                              length_range = 28, utr5l = 20, cdsl = 40,
#                                              utr3l = 20, colour = "aquamarine4",
#                                              plot_title = "sample.transcript.length_range")
#  example_metaprofile_28[["plot_Samp1"]]

## ---- out.width = '690px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaprofile_28.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  comparison_list <- list()
#  comparison_list[["subsample_28nt"]] <- reads_psite_list[["Samp1"]][length == 28]
#  comparison_list[["whole_sample"]] <- reads_psite_list[["Samp1"]]

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  sample_list <- list("Only_28" = c("subsample_28nt"),
#                     "All" = c("whole_sample"))

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_metaprofile_split <- metaprofile_psite(comparison_list, mm81cdna, sample = sample_list,
#                                                 multisamples = "average", plot_style = "split",
#                                                 utr5l = 20, cdsl = 40, utr3l = 20,
#                                                 frequency = TRUE, plot_title = "transcript",
#                                                 colour = c("aquamarine4", "gray70"))
#  example_metaprofile_split[["plot"]]

## ---- out.width = '690px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaprofile_split.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_metaprofile_overlaid <- metaprofile_psite(comparison_list, mm81cdna, sample = sample_list,
#                                                    multisamples = "average", plot_style = "overlaid",
#                                                    utr5l = 20, cdsl = 40, utr3l = 20,
#                                                    frequency = TRUE, plot_title = "transcript",
#                                                    colour = c("aquamarine4", "gray70"))
#  example_metaprofile_overlaid[["plot"]]

## ---- out.width = '690px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaprofile_overlaid.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_metaprofile_mirrored <- metaprofile_psite(comparison_list, mm81cdna, sample = sample_list,
#                                                    multisamples = "average", plot_style = "mirrored",
#                                                    utr5l = 20, cdsl = 40, utr3l = 20,
#                                                    frequency = TRUE, plot_title = "transcript",
#                                                    colour = c("aquamarine4", "gray70"))
#  example_metaprofile_mirrored[["plot"]]

## ---- out.width = '690px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaprofile_mirrored.png")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  example_metaheatmap <- metaheatmap_psite(comparison_list, mm81cdna, sample = sample_list,
#                                           utr5l = 20, cdsl = 40, utr3l = 20, log_colour = F,
#                                           plot_title = "Comparison metaheatmap")
#  example_metaheatmap[["plot"]]

## ---- out.width = '700px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaheatmap.png")

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  example_cu_barplot <- codon_usage_psite(reads_psite_list, mm81cdna, sample = "Samp1",
#                                          fastapath = "path/to/transcriptome/FASTA/file",
#                                          fasta_genome = FALSE,
#                                          frequency_normalization = FALSE)
#  example_cu_barplot[["plot"]]

## ---- out.width = '690px', fig.retina = NULL, echo = FALSE--------------------
knitr::include_graphics("example_cu_barplot.png")

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  example_cu_scatter_2samples <- codon_usage_psite(comparison_list, mm81cdna, sample = c("All", "Only_28"),
#                                                   fastapath = "path/to/transcriptome/FASTA/file",
#                                                   fasta_genome = FALSE,
#                                                   frequency_normalization = FALSE)
#  example_cu_scatter_2samples[["plot_comparison"]]

## ---- out.width = '300px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_cu_scatter_2samples.png")

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
cub_mouse <- data.table(codon = c("UUU", "UCU", "UAU", "UGU", "UUC", "UCC", "UAC", "UGC", "UUA", "UCA", "UAA", "UGA", "UUG", "UCG", "UAG", "UGG", "CUU", "CCU", "CAU", "CGU", "CUC", "CCC", "CAC", "CGC", "CUA", "CCA", "CAA", "CGA", "CUG", "CCG", "CAG", "CGG", "AUU", "ACU", "AAU", "AGU", "AUC", "ACC", "AAC", "AGC", "AUA", "ACA", "AAA", "AGA", "AUG", "ACG", "AAG", "AGG", "GUU", "GCU", "GAU", "GGU", "GUC", "GCC", "GAC", "GGC", "GUA", "GCA", "GAA", "GGA", "GUG", "GCG", "GAG", "GGG"),
                        value = c(17.2, 16.2, 12.2, 11.4, 21.8, 18.1, 16.1, 12.3, 6.7, 11.8, 1.0, 1.6, 13.4, 4.2, 0.8, 12.5, 13.4, 18.4, 10.6, 4.7, 20.2, 18.2, 15.3, 9.4, 8.1, 17.3, 12.0, 6.6, 39.5, 6.2, 34.1, 10.2, 15.4, 13.7, 15.6, 12.7, 22.5, 19.0, 20.3, 19.7, 7.4, 16.0, 21.9, 12.1, 22.8, 5.6, 33.6, 12.2, 10.7, 20.0, 21.0, 11.4, 15.4, 26.0, 26.0, 21.2, 7.4, 15.8, 27.0, 16.8, 28.4, 6.4, 39.4, 15.2))

## -----------------------------------------------------------------------------
head(cub_mouse)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  example_cu_scatter_cub <- codon_usage_psite(reads_psite_list, mm81cdna, sample = "Samp1",
#                                              fastapath = "path/to/transcriptome/FASTA/file",
#                                              fasta_genome = FALSE, codon_values = cub_mouse,
#                                              frequency_normalization = FALSE,
#                                              label_scatter = TRUE, label_number = 5)
#  example_cu_scatter_cub[["plot_comparison"]]

## ---- out.width = '300px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_cu_scatter_cub.png")

