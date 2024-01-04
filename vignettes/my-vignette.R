## ----eval = FALSE-------------------------------------------------------------
#  install.packages("devtools")

## ----eval = FALSE-------------------------------------------------------------
#  library(devtools)
#  install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE,
#                 build_vignettes = TRUE)

## ----eval = TRUE, warning = FALSE---------------------------------------------
library(riboWaltz)

## ----eval = FALSE-------------------------------------------------------------
#  ?function_name

## ----eval = FALSE-------------------------------------------------------------
#  help(package = riboWaltz)

## ----eval = FALSE-------------------------------------------------------------
#    reads_list <- bamtolist(bamfolder = "path/to/BAM/files", annotation = annotation_dt)

## -----------------------------------------------------------------------------
head(reads_list[["Samp1"]])

## ----eval = FALSE-------------------------------------------------------------
#    bamtobed(bamfolder = "path/to/BAM/files", bedfolder = "path/to/output/directory")

## ----eval = FALSE-------------------------------------------------------------
#    reads_list <- bedtolist(bedfolder = "path/to/BED/files", annotation = annotation_dt)

## ----eval = TRUE, echo = FALSE------------------------------------------------
library(data.table)
example_reads_list <- list()
example_reads_list[["Samp_example"]] <- data.table(row_ID = c(1, 2, 3, 4, 5, 6),
                                                   transcript = rep("ENSMUST00000000001.4", 6),
                                                   end5 = c(92, 92, 92, 94, 94, 95),
                                                   end3 = c(119, 119, 122, 122, 123, 123))
example_reads_list[["Samp_example"]][, length := end3 - end5 + 1
                                     ][, cds_start := 142
                                       ][, cds_stop := 1206]

## ----eval = TRUE, echo = TRUE-------------------------------------------------
print(example_reads_list[["Samp_example"]], row.names=FALSE)

## ----eval = TRUE, echo = TRUE, results = "hide"-------------------------------
filtered_list <- duplicates_filter(data = example_reads_list,
                                   extremity = "both")
filtered_list[["Samp_example"]]

## ----echo = FALSE, eval = TRUE------------------------------------------------
print(filtered_list[["Samp_example"]], row.names=FALSE)

## ----eval = TRUE, echo = TRUE, results = "hide"-------------------------------
filtered_list <- duplicates_filter(data = example_reads_list,
                                   extremity = "5end",
                                   keep = "shortest")
filtered_list[["Samp_example"]]

## ----echo = FALSE, eval = TRUE------------------------------------------------
print(filtered_list[["Samp_example"]], row.names=FALSE)

## ----eval = TRUE, echo = TRUE, results = "hide"-------------------------------
filtered_list <- duplicates_filter(data = example_reads_list,
                                   extremity = "3end",
                                   keep = "longest")
filtered_list[["Samp_example"]]

## ----echo = FALSE, eval = TRUE------------------------------------------------
print(filtered_list[["Samp_example"]], row.names=FALSE)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  filtered_list <- length_filter(data = example_reads_list,
#                                 length_filter_mode = "periodicity",
#                                 periodicity_threshold = 70)

## ----eval = TRUE, echo = TRUE, results = "hide"-------------------------------
filtered_list <- length_filter(data = example_reads_list,
                               length_filter_mode = "custom",
                               length_range = 29:30)
filtered_list[["Samp_example"]]

## ----echo = FALSE, eval = TRUE------------------------------------------------
print(filtered_list[["Samp_example"]], row.names=FALSE)

## ----eval = TRUE, echo = FALSE, message = FALSE, warning = FALSE--------------
data(reads_list)
data(mm81cdna)

## -----------------------------------------------------------------------------
head(mm81cdna)

## ----echo = TRUE--------------------------------------------------------------
psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")

## ----echo = TRUE--------------------------------------------------------------
head(psite_offset, 10)

## ----out.width = '690px', fig.retina = NULL, echo = FALSE---------------------
knitr::include_graphics("meta_psite_length28.png")
knitr::include_graphics("meta_psite_length31.png")

## ----echo = TRUE, eval = TRUE, results = "hide"-------------------------------
reads_psite_list <- psite_info(reads_list, psite_offset)
head(reads_psite_list[["Samp1"]])

## ----echo = FALSE, eval = TRUE------------------------------------------------
head(reads_psite_list[["Samp1"]])

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  codon_coverage_example <- codon_coverage(reads_psite_list, mm81cdna, psite = FALSE)
#  head(codon_coverage_example[Samp1 > 0])

## ----echo = FALSE, eval = TRUE------------------------------------------------
codon_coverage_example <- data.table(transcript = rep("ENSMUST00000000001.4", 6),
                                     start = c(90, 93, 96, 99, 102, 105),
                                     end = c(93, 96, 99, 102, 105, 108),
                                     from_cds_start = -17:-12,
                                     from_cds_stop = -371:-366,
                                     region = rep("5utr", 6),
                                     Samp1 = c(1, 2, 2, 2, 2, 2))
head(codon_coverage_example)

## ----echo = FALSE, eval = TRUE, include = FALSE-------------------------------
cds_coverage_example <- cds_coverage(reads_psite_list, mm81cdna)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  cds_coverage_example <- cds_coverage(reads_psite_list, mm81cdna)
#  head(cds_coverage_example)

## ----echo = FALSE, eval = TRUE------------------------------------------------
head(cds_coverage_example)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  for(i in 2:6){
#      samp_name <- paste0("Samp", i)
#  	set.seed(i)
#  	reads_list[[samp_name]] <- reads_list[["Samp1"]][sample(.N, 5000)]
#  }

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  psite_offset <- psite(reads_list, flanking = 6, extremity = "3end")
#  reads_psite_list <- psite_info(reads_list, psite_offset)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  example_length_dist <- rlength_distr(reads_list, sample = "Samp1",
#                                       colour = "#333f50")
#  example_length_dist[["plot_Samp1"]]

## ----out.width = '300px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_basic.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  example_length_dist <- rlength_distr(reads_list, sample = "Samp1",
#                                       cl = 99, colour = "#333f50")
#  example_length_dist[["plot_Samp1"]]

## ----out.width = '300px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_basic_zoom.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  input_samples <- c("Samp1", "Samp2")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  example_length_dist <- rlength_distr(reads_list,
#                                       sample = input_samples,
#                                       multisamples = "independent",
#                                       plot_style = "split",
#                                       cl = 99, colour = c("#333f50", "#39827c"))
#  example_length_dist[["plot_Samp1"]]
#  example_length_dist[["plot_Samp2"]]

## ----out.width = '300px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_1_1.png")
knitr::include_graphics("example_length_dist_1_2.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  example_length_dist <- rlength_distr(reads_list,
#                                       sample = input_samples,
#                                       multisamples = "independent",
#                                       plot_style = "facet",
#                                       cl = 99, colour = c("#333f50", "#39827c"))
#  example_length_dist[["plot"]]

## ----out.width = '450px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_2.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  example_length_dist <- rlength_distr(reads_list,
#                                       sample = input_samples,
#                                       multisamples = "independent",
#                                       plot_style = "dodge",
#                                       cl = 99, colour = c("#333f50", "#39827c"))
#  example_length_dist[["plot"]]

## ----out.width = '400px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_3.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  example_length_dist <- rlength_distr(reads_list,
#                                       sample = input_samples,
#                                       multisamples = "independent",
#                                       plot_style = "mirror",
#                                       cl = 99, colour = c("#333f50", "#39827c"))
#  example_length_dist[["plot"]]

## ----out.width = '330px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_4.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  input_samples <- c("Samp1", "Samp2")
#  example_length_dist <- rlength_distr(reads_list,
#                                       sample = input_samples,
#                                       multisamples = "average",
#                                       cl = 99, colour = "#333f50")
#  example_length_dist[["plot_Average"]]

## ----out.width = '300px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_rep1_1.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"))
#  example_length_dist <- rlength_distr(reads_list,
#                                       sample = input_samples,
#                                       multisamples = "average",
#                                       cl = 99, colour = "#333f50")
#  example_length_dist[["plot_S1"]]

## ----out.width = '300px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_rep1_2.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                     "S2" = c("Samp3", "Samp4", "Samp5"))
#  example_length_dist <- rlength_distr(reads_list,
#                                       sample = input_samples,
#                                       multisamples = "average",
#                                       plot_style = "split",
#                                       cl = 99, colour = c("#333f50", "#39827c"))
#  example_length_dist[["plot_S1"]]
#  example_length_dist[["plot_S2"]]

## ----out.width = '300px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_rep2_1_1.png")
knitr::include_graphics("example_length_dist_rep2_1_2.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#  								"S2" = c("Samp3", "Samp4", "Samp5"),
#  						    	"S3" = c("Samp6"))
#  example_length_dist <- rlength_distr(reads_list,
#                                       sample = input_samples,
#                                       multisamples = "average",
#                                       plot_style = "dodge",
#                                       cl = 99,
#                                       colour = c("#333f50", "#39827c", "gray70"))
#  example_length_dist[["plot"]]

## ----out.width = '400px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_rep2_2.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#  						     "S3" = c("Samp6"))
#  example_length_dist <- rlength_distr(reads_list,
#                                       sample = input_samples,
#                                       multisamples = "average",
#                                       plot_style = "mirror",
#                                       cl = 99, colour = c("#333f50", "gray70"))
#  example_length_dist[["plot"]]

## ----out.width = '330px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_length_dist_rep2_3.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  example_ends_heatmap <- rends_heat(reads_list, mm81cdna,
#                                     sample = "Samp1",
#                                     cl = 85, utr5l = 25, cdsl = 40, utr3l = 25)
#  example_ends_heatmap[["plot_Samp1"]]

## ----out.width = '700px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_ends_heatmap.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                        "S2" = c("Samp3", "Samp4", "Samp5"),
#                        "S3" = c("Samp6"))
#  
#  example_psite_per_region <- region_psite(reads_psite_list, mm81cdna,
#                                           sample = input_samples,
#                                           multisamples = "average",
#                                           plot_style = "stack",
#                                           cl = 85,
#                                           colour = c("#333f50", "gray70", "#39827c"))
#  example_psite_per_region[["plot"]]

## ----out.width = '350px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_psite_per_region1.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                        "S2" = c("Samp3", "Samp4", "Samp5"),
#                        "S3" = c("Samp6"))
#  
#  example_psite_per_region <- region_psite(reads_psite_list, mm81cdna,
#                                           sample = input_samples,
#                                           multisamples = "average",
#                                           plot_style = "dodge",
#                                           cl = 85,
#                                           colour = c("#333f50", "gray70", "#39827c"))
#  example_psite_per_region[["plot"]]

## ----out.width = '375px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_psite_per_region2.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                        "S2" = c("Samp3", "Samp4", "Samp5"),
#                        "S3" = c("Samp6"))
#  
#  example_frames_stratified <- frame_psite_length(reads_psite_list, mm81cdna,
#                                                  sample = input_samples,
#                                                  multisamples = "average",
#                                                  plot_style = "facet",
#                                                  region = "all",
#                                                  cl = 85, colour = "#333f50")
#  example_frames_stratified[["plot"]]

## ----out.width = '375px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_frames_stratified.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                        "S2" = c("Samp3", "Samp4", "Samp5"))
#  
#  example_frames <- frame_psite(reads_psite_list, mm81cdna,
#                             sample = input_samples,
#                             multisamples = "average",
#                             plot_style = "facet",
#                             region = "cds",
#                             colour = c("#333f50", "#39827c"))
#  example_frames[["plot"]]

## ----out.width = '375px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_frames1.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                        "S2" = c("Samp3", "Samp4", "Samp5"))
#  
#  example_frames <- frame_psite(reads_psite_list, mm81cdna,
#                                sample = input_samples,
#                                multisamples = "average",
#                                plot_style = "mirror",
#                                region = "all",
#                                colour = c("#333f50", "#39827c"))
#  example_frames[["plot"]]

## ----out.width = '575px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_frames2.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2")))
#  
#  example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna,
#                                           sample = input_samples,
#                                           multisamples = "average",
#                                           utr5l = 20, cdsl = 40, utr3l = 20,
#                                           colour = "#333f50")
#  example_metaprofile[["plot_S1"]]

## ----out.width = '690px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaprofile1.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                        "S2" = c("Samp3", "Samp4", "Samp5"))

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna,
#                                           sample = input_samples,
#                                           multisamples = "average",
#                                           plot_style = "facet",
#                                           utr5l = 20, cdsl = 40, utr3l = 20,
#                                           colour = c("#333f50", "#39827c"))
#  example_metaprofile[["plot"]]

## ----out.width = '650px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaprofile2.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna,
#                                           sample = input_samples,
#                                           multisamples = "average",
#                                           plot_style = "overlap",
#                                           utr5l = 20, cdsl = 40, utr3l = 20,
#                                           colour = c("#333f50", "#39827c"))
#  example_metaprofile[["plot"]]

## ----out.width = '650px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaprofile3.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna,
#                                           sample = input_samples,
#                                           multisamples = "average",
#                                           plot_style = "mirror",
#                                           utr5l = 20, cdsl = 40, utr3l = 20,
#                                           colour = c("#333f50", "#39827c"))
#  example_metaprofile[["plot"]]

## ----out.width = '650px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaprofile4.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                        "S2" = c("Samp3", "Samp4", "Samp5"),
#                        "S3" = c("Samp6"))
#  
#  example_metaheatmap <- metaheatmap_psite(reads_psite_list, mm81cdna,
#                                           sample = input_samples,
#                                           multisamples = "average",
#                                           utr5l = 20, cdsl = 40, utr3l = 20,
#                                           colour = "#333f50"))
#  example_metaheatmap[["plot"]]

## ----out.width = '650px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_metaheatmap.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"))
#  
#  example_cu_barplot <- codon_usage_psite(reads_psite_list, mm81cdna,
#  		                                        sample = input_samples,
#  		                                        multisamples = "average",
#  		                                        plot_style = "facet",
#  		                                        fastapath = "path/to/transcriptome/FASTA/file",
#  		                                        fasta_genome = FALSE,
#  		                                        frequency_normalization = FALSE)
#  example_cu_barplot[["plot_S1"]]

## ----out.width = '650px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_cu_barplot1.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                        "S2" = c("Samp3", "Samp4", "Samp5"))
#  
#  example_cu_barplot <- codon_usage_psite(reads_psite_list, mm81cdna,
#  		                                        sample = input_samples,
#  		                                        multisamples = "average",
#  		                                        plot_style = "facet",
#  		                                        fastapath = "path/to/transcriptome/FASTA/file",
#  		                                        fasta_genome = FALSE,
#  		                                        frequency_normalization = FALSE)
#  example_cu_barplot[["plot"]]

## ----out.width = '650px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_cu_barplot2.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"))
#  
#  example_cu_barplot <- codon_usage_psite(reads_psite_list, mm81cdna,
#  		                                        sample = input_samples,
#  		                                        multisamples = "average",
#  		                                        plot_style = "facet",
#  		                                        fastapath = "path/to/transcriptome/FASTA/file",
#  		                                        fasta_genome = FALSE,
#  		                                        frequency_normalization = TRUE,
#  		                                        include_stop_codons = FALSE)
#  example_cu_barplot[["plot_S1"]]

## ----out.width = '650px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_cu_barplot3.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  input_samples <- list("S1" = c("Samp1", "Samp2"),
#                        "S2" = c("Samp3", "Samp4", "Samp5"))

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  example_cu_barplot <- codon_usage_psite(reads_psite_list, mm81cdna,
#                                          sample = input_samples,
#                                          contrast_sample = c("S1", "S2"),
#                                          fastapath = "path/to/transcriptome/FASTA/file",
#                                          fasta_genome = FALSE,
#                                          frequency_normalization = FALSE)
#  example_cu_barplot[["plot"]]

## ----out.width = '350px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_cu_scatterplot1.png")

## ----eval=TRUE, echo=FALSE----------------------------------------------------
cub_mouse <- data.table(codon = c("UUU", "UCU", "UAU", "UGU", "UUC", "UCC", "UAC", "UGC", "UUA", "UCA", "UAA", "UGA", "UUG", "UCG", "UAG", "UGG", "CUU", "CCU", "CAU", "CGU", "CUC", "CCC", "CAC", "CGC", "CUA", "CCA", "CAA", "CGA", "CUG", "CCG", "CAG", "CGG", "AUU", "ACU", "AAU", "AGU", "AUC", "ACC", "AAC", "AGC", "AUA", "ACA", "AAA", "AGA", "AUG", "ACG", "AAG", "AGG", "GUU", "GCU", "GAU", "GGU", "GUC", "GCC", "GAC", "GGC", "GUA", "GCA", "GAA", "GGA", "GUG", "GCG", "GAG", "GGG"),
                        value = c(17.2, 16.2, 12.2, 11.4, 21.8, 18.1, 16.1, 12.3, 6.7, 11.8, 1.0, 1.6, 13.4, 4.2, 0.8, 12.5, 13.4, 18.4, 10.6, 4.7, 20.2, 18.2, 15.3, 9.4, 8.1, 17.3, 12.0, 6.6, 39.5, 6.2, 34.1, 10.2, 15.4, 13.7, 15.6, 12.7, 22.5, 19.0, 20.3, 19.7, 7.4, 16.0, 21.9, 12.1, 22.8, 5.6, 33.6, 12.2, 10.7, 20.0, 21.0, 11.4, 15.4, 26.0, 26.0, 21.2, 7.4, 15.8, 27.0, 16.8, 28.4, 6.4, 39.4, 15.2))

## -----------------------------------------------------------------------------
head(cub_mouse)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  example_cu_barplot <- codon_usage_psite(reads_psite_list, mm81cdna,
#                                          sample = input_samples,
#                                          contrast_sample = "S1",
#                                          codon_values = cub_mouse,
#                                          fastapath = "path/to/transcriptome/FASTA/file",
#                                          fasta_genome = FALSE,
#                                          frequency_normalization = FALSE,
#                                          label_scatter = TRUE, label_number = 5)
#  example_cu_barplot[["plot"]]

## ----out.width = '350px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("example_cu_scatterplot2.png")

