## ---- eval = FALSE-------------------------------------------------------
#  library(devtools)
#  install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)

## ---- eval = TRUE, warning = FALSE---------------------------------------
library(riboWaltz)

## ------------------------------------------------------------------------
head(reads_list[["Samp1"]])

## ---- eval = TRUE, echo = FALSE------------------------------------------
data(reads_list)
data(mm81cdna)
library(data.table)

## ------------------------------------------------------------------------
head(mm81cdna)

## ---- echo = TRUE, fig.show = 'hold', fig.width = 4, fig.asp = 1, fig.align = 'center', out.width = '250px'----
example_length_dist <- rlength_distr(reads_list, sample = "Samp1")
head(example_length_dist[["dt"]])
example_length_dist[["plot"]]

## ---- echo = TRUE, fig.show = 'hold', fig.width = 4, fig.asp = 1, fig.align = 'center', out.width = '250px', warning = FALSE----
example_length_dist_zoom <- rlength_distr(reads_list, sample = "Samp1", cl = 99)
example_length_dist_zoom[["plot"]]

## ---- echo = TRUE, fig.show = 'hold', fig.width = 15, fig.asp = 1/3, fig.align = 'center', out.width = '700px', message = FALSE, warning = FALSE----
example_ends_heatmap <- rends_heat(reads_list, mm81cdna, sample = "Samp1", cl = 85,
                                      utr5l = 25, cdsl = 40, utr3l = 25)
head(example_ends_heatmap[["dt"]])
example_ends_heatmap[["plot"]]

## ---- echo = TRUE--------------------------------------------------------
psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")

## ---- echo = TRUE--------------------------------------------------------
head(psite_offset, 10)

## ---- out.width = '690px', fig.retina = NULL, echo = FALSE---------------
knitr::include_graphics("meta_psite_length28.png")
knitr::include_graphics("meta_psite_length31.png")

## ---- echo = TRUE--------------------------------------------------------
reads_psite_list <- psite_info(reads_list, psite_offset)
head(reads_psite_list[["Samp1"]])

## ---- echo = TRUE--------------------------------------------------------
psite_cds_list <- psite_per_cds(reads_psite_list, mm81cdna)

## ---- echo = TRUE--------------------------------------------------------
head(psite_cds_list[["Samp1"]])

## ---- echo = TRUE, fig.show = 'hold', fig.width = 10, fig.asp = 1/2.1, fig.align = 'center', out.width = '450px', message = FALSE, warning = FALSE----
example_frames_stratified <- frame_psite_length(reads_psite_list, sample = "Samp1",
                                                   region = "all", cl = 90)
head(example_frames_stratified[["dt"]])
example_frames_stratified[["plot"]]

## ---- echo = TRUE, fig.show = 'hold', fig.width = 10, fig.asp = 1/2.3, fig.align = 'center', out.width = '450px', message = FALSE, warning = FALSE----
example_frames <- frame_psite(reads_psite_list, sample = "Samp1", region = "all")
head(example_frames[["dt"]])
example_frames[["plot"]]

## ---- echo = TRUE, fig.show = 'hold', fig.width = 15, fig.asp = 1/3, fig.align = 'center', out.width = '690px', message = FALSE, warning = FALSE----
example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
                                            utr5l = 20, cdsl = 40, utr3l = 20)
example_metaprofile[["plot"]]

## ---- echo = TRUE, fig.show = 'hold', fig.width = 15, fig.asp = 1/3, fig.align = 'center', out.width = '690px', message = FALSE, warning = FALSE----
example_metaprofile_28 <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
                                               length_range = 28, utr5l = 20, cdsl = 40,
                                               utr3l = 20)
example_metaprofile_28[["plot"]]

## ---- echo = TRUE--------------------------------------------------------
comparison_dt <- list()
comparison_dt[["subsample_28nt"]] <- reads_psite_list[["Samp1"]][length == 28]
comparison_dt[["whole_sample"]] <- reads_psite_list[["Samp1"]]

## ---- echo = TRUE--------------------------------------------------------
names_list <- list("Only_28" = c("subsample_28nt"),
                   "All" = c("whole_sample"))

## ---- echo = TRUE, fig.show = 'hold', fig.width = 15, fig.asp = 1/2.5, fig.align = 'center', out.width = '700px', message = FALSE, warning = FALSE----
example_metaheatmap <- metaheatmap_psite(comparison_dt, mm81cdna, sample = names_list,
                                         utr5l = 20, cdsl = 40, utr3l = 20, log = F)
example_metaheatmap[["plot"]]

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  codon_usage_barplot <- codon_usage_psite(reads_psite_list, mm81cdna, sample = "Samp1",
#  										  	 fastapath = path_to_fasta)

## ---- out.width = '690px', fig.retina = NULL, echo = FALSE---------------
knitr::include_graphics("codon_usage_barplot.png")

## ---- eval=TRUE, echo=FALSE----------------------------------------------
cub_mouse <- data.table(codon = c("UUU", "UCU", "UAU", "UGU", "UUC", "UCC", "UAC", "UGC", "UUA", "UCA", "UAA", "UGA", "UUG", "UCG", "UAG", "UGG", "CUU", "CCU", "CAU", "CGU", "CUC", "CCC", "CAC", "CGC", "CUA", "CCA", "CAA", "CGA", "CUG", "CCG", "CAG", "CGG", "AUU", "ACU", "AAU", "AGU", "AUC", "ACC", "AAC", "AGC", "AUA", "ACA", "AAA", "AGA", "AUG", "ACG", "AAG", "AGG", "GUU", "GCU", "GAU", "GGU", "GUC", "GCC", "GAC", "GGC", "GUA", "GCA", "GAA", "GGA", "GUG", "GCG", "GAG", "GGG"),
                        value = c(17.2, 16.2, 12.2, 11.4, 21.8, 18.1, 16.1, 12.3, 6.7, 11.8, 1.0, 1.6, 13.4, 4.2, 0.8, 12.5, 13.4, 18.4, 10.6, 4.7, 20.2, 18.2, 15.3, 9.4, 8.1, 17.3, 12.0, 6.6, 39.5, 6.2, 34.1, 10.2, 15.4, 13.7, 15.6, 12.7, 22.5, 19.0, 20.3, 19.7, 7.4, 16.0, 21.9, 12.1, 22.8, 5.6, 33.6, 12.2, 10.7, 20.0, 21.0, 11.4, 15.4, 26.0, 26.0, 21.2, 7.4, 15.8, 27.0, 16.8, 28.4, 6.4, 39.4, 15.2))

## ------------------------------------------------------------------------
head(cub_mouse)

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  codon_usage_scatter <- codon_usage_psite(reads_psite_list, mm81cdna, sample = "Samp1",
#  											   fastapath = path_to_fasta, codon_usage = cub_mouse)

## ---- out.width = '260px', fig.retina = NULL, echo = FALSE, fig.align = "center"----
knitr::include_graphics("codon_usage_scatter.png")

