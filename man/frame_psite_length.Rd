% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frames.R
\name{frame_psite_length}
\alias{frame_psite_length}
\title{Percentage of P-sites per reading frame stratified by read length.}
\usage{
frame_psite_length(
  data,
  annotation,
  sample,
  multisamples = "average",
  plot_style = "split",
  transcripts = NULL,
  region = "all",
  length_range = NULL,
  cl = 100,
  colour = "#172969"
)
}
\arguments{
\item{data}{Either list of data tables or GRangesList object from
\code{\link{psite_info}}.}

\item{annotation}{Data table as generated by \code{\link{create_annotation}}.}

\item{sample}{Either character string, character string vector or named list
of character string(s)/character string vector(s) specifying the name of
the sample(s) and replicate(s) of interest. If a list is provided, each
element of the list is considered as an independent sample associated with
one ore multiple replicates. Multiple samples and replicates are handled
and visualised according to \code{multisample} and \code{plot_style}.}

\item{multisamples}{Either "average" or "independent". It specifies how to
handle multiple samples and replicates stored in \code{sample}:
\itemize{
\item if \code{sample} is a character string vector and \code{multisample} is
set to "average" the elements of the vector are considered as replicates
of one sample and a single heatmap is returned.
\item if \code{sample} is a character string vector and \code{multisample} is
set to "independent", each element of the vector is analysed independently
of the others. The number of plots returned and their organization is
specified by \code{plot_style}.
\item if \code{sample} is list, \code{multisample} must be set to "average".
Each element of the list is analysed independently of the others, its
replicates averaged and its name reported in the plot. The number of plots
returned and their organization is specified by \code{plot_style}.
Note: when this parameter is set to "average" the heatmap associated with
each sample displays the length- and frame- specific mean signal
computed across the replicates. Default is "average".
}}

\item{plot_style}{Either "split" or "facet". It specifies how to organize and
display multiple heatmaps:
\itemize{
\item "split": heatmap(s) for each sample are returned as an independent
ggplot object.
\item "facet": the heatmap(s) for each sample are placed one below the other
(when \code{region} is "all") or one next to the other (when \code{region}
is set to either "5utr", "cds", "3utr") in independent boxes. Default is
"split".
}}

\item{transcripts}{Character string vector listing the name of transcripts to
be included in the analysis. Default is NULL, i.e. all transcripts are used.}

\item{region}{Character string specifying the region(s) of the transcripts to
be analysed. It can be either "5utr", "cds", "3utr", respectively for 5'
UTRs, CDSs and 3' UTRs, or "all", i.e. all regions are considered. Note:
according to this parameter the heatmaps are differently arranged to
optimise the organization and the visualization of the data.}

\item{length_range}{Integer or integer vector for restricting the plot to
a chosen range of read lengths. Default is NULL, i.e. all read lengths are
used.}

\item{cl}{Integer value in \link{1,100} specifying a confidence level for
restricting the plot to an automatically-defined range of read lengths. The
new range is computed according to the most frequent read lengths, which
accounts for the cl\% of the sample and is defined by discarding the
(100-cl)\% of read lengths falling in the tails of the read lengths
distribution. If multiple samples are analysed, a single range of read
lengths is computed such that at least the cl\% of all sample are
represented. Default is 100.}

\item{colour}{Character string specifying the colour of the plot. The colour
scheme is as follow: tiles corresponding to the lowest signal are always
white, tiles corresponding to the highest signal are of the specified
colour and the progression between these two colours follows a linear
gradient. Default is dark blue.}
}
\value{
List containing: one or more ggplot object(s) and the data table with
the corresponding x- and y-axis values and the z-values, defining the color
of the tiles ("plot_dt"); an additional data table with raw and scaled
number of read extremities mapping around the start and the stop codon, per
length, for each sample ("count_dt").
}
\description{
This function computes, for each transcript region (5' UTR, coding sequence
and 3' UTR), the percentage of P-sites falling in the three possible
translation reading frames, stratified by read length. It takes into account
only transcripts with annotated regions of length > 0 and the resulting
values are visualized as heatmaps. Multiple samples and replicates can be
handled in several ways.
}
\examples{
## data(reads_list)
## data(mm81cdna)
##
## ## Generate fake samples and replicates
## for(i in 2:6){
##   samp_name <- paste0("Samp", i)
##   set.seed(i)
##   reads_list[[samp_name]] <- reads_list[["Samp1"]][sample(.N, 5000)]
## }
##
## ## Compute and add p-site details
## psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")
## reads_psite_list <- psite_info(reads_list, psite_offset)
##
## ## Define the list of samples and replicate to use as input
## input_samples <- list("S1" = c("Samp1", "Samp2"),
##                       "S2" = c("Samp3", "Samp4", "Samp5"),
##                       "S3" = c("Samp6"))
##
## Generate heatmaps
## example_frames_stratified <- frame_psite_length(reads_psite_list, mm81cdna,
##                                                 sample = input_samples,
##                                                 multisamples = "average",
##                                                 plot_style = "split",
##                                                 region = "all",
##                                                 colour = "#333f50")
}
\seealso{
\code{\link{frame_psite}} for a similar plot where read lengths are
collapsed.
}
