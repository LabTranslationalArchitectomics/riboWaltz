## ---- eval=TRUE, warning=FALSE-------------------------------------------
library(riboWaltz)

## ------------------------------------------------------------------------
head(reads_list[["Samp1"]])

## ---- eval=TRUE, echo=FALSE----------------------------------------------
data(reads_list)
data(mm81cdna)

## ------------------------------------------------------------------------
head(mm81cdna)

## ---- echo=TRUE, fig.show='hold', fig.width=4, fig.asp=1, fig.align='center', out.width='250px'----
example_length_dist <- rlength_distr(reads_list, sample="Samp1")
head(example_length_dist[["df"]])
example_length_dist[["plot"]]

## ---- echo=TRUE, fig.show='hold', fig.width=4, fig.asp=1, fig.align='center', out.width='250px'----
example_length_dist_zoom <- rlength_distr(reads_list, sample="Samp1", cl=99)
example_length_dist_zoom[["plot"]]

## ---- echo=TRUE, fig.show='hold', fig.width=15, fig.asp=1/3, fig.align='center', out.width='700px', message=FALSE, warning=FALSE----
example_ends_heatmap <- rends_heat(reads_list, mm81cdna, sample="Samp1", cl=85,
                                      utr5l = 25, cdsl = 40, utr3l = 25)
example_ends_heatmap[["plot"]]

## ---- echo=TRUE----------------------------------------------------------
psite_offset <- psite(reads_list, flanking = 6, extremity="auto")

## ---- echo=TRUE----------------------------------------------------------
head(psite_offset, 10)

## ---- out.width = '700px', fig.retina = NULL, echo =FALSE----------------
knitr::include_graphics("meta_psite_length28.png")
knitr::include_graphics("meta_psite_length31.png")

## ---- echo=TRUE----------------------------------------------------------
reads_psite_list <- psite_info(reads_list, psite_offset)
head(reads_psite_list[["Samp1"]])

## ---- echo=TRUE, fig.show='hold', fig.width=10, fig.asp=1/2.1, fig.align='center', out.width='450px', message=FALSE, warning=FALSE----
example_frames_stratified <- frame_psite_length(reads_psite_list, sample="Samp1",
                                                   region="all", cl=90)
example_frames_stratified[["plot"]]

## ---- echo=TRUE, fig.show='hold', fig.width=10, fig.asp=1/2.3, fig.align='center', out.width='450px', message=FALSE, warning=FALSE----
example_frames <- frame_psite(reads_psite_list, sample="Samp1", region="all")
example_frames[["plot"]]

## ---- echo=TRUE, fig.show='hold', fig.width=15, fig.asp=1/3, fig.align='center', out.width='700px', message=FALSE, warning=FALSE----
example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
                                            utr5l = 20, cdsl = 40, utr3l = 20)
example_metaprofile[["plot"]]

## ---- echo=TRUE, fig.show='hold', fig.width=15, fig.asp=1/3, fig.align='center', out.width='700px', message=FALSE, warning=FALSE----
example_metaprofile_28 <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1",
                                               length_range = 28, utr5l = 20, cdsl = 40,
                                               utr3l = 20)
example_metaprofile_28[["plot"]]

## ---- echo=TRUE----------------------------------------------------------
comparison_df <- list()
comparison_df[["subsample_28nt"]] <- subset(reads_psite_list[["Samp1"]], length == 28)
comparison_df[["whole_sample"]] <- reads_psite_list[["Samp1"]]

## ---- echo=TRUE----------------------------------------------------------
names_list <- list("Only_28" = c("subsample_28nt"),
                   "All" = c("whole_sample"))

## ---- echo=TRUE, fig.show='hold', fig.width=15, fig.asp=1/2.5, fig.align='center', out.width='700px', message=FALSE, warning=FALSE----
example_metaheatmap <- metaheatmap_psite(comparison_df, mm81cdna, sample = names_list,
                                         utr5l = 20, cdsl = 40, utr3l = 20, log=F)
example_metaheatmap[["plot"]]

