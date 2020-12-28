#' Read length distributions.
#'
#' This function generates read length distributions.
#'
#' @param data List of data tables from \code{\link{bamtolist}},
#'   \code{\link{bedtolist}}, \code{\link{length_filter}} or
#'   \code{\link{psite_info}}.
#' @param sample Character string specifying the name of the sample of interest.
#' @param transcripts Character string vector listing the name of transcripts to
#'   be included in the analysis. Default is NULL i.e. all transcripts are used.
#' @param cl Integer value in [1,100] specifying a confidence level for
#'   restricting the plot to a sub-range of read lengths i.e. to the cl% of
#'   read lengths associated to the highest signals. Default is 100.
#' @return List containing a ggplot2 object ("plot") and the data table with the
#'   associated data ("dt").
#' @examples
#' data(reads_list)
#'
#' ## Generate the length distribution for all read lengths:
#' lendist_whole <- rlength_distr(reads_list, sample = "Samp1", cl = 100)
#' lendist_whole[["plot"]]
#'
#' ## Generate the length distribution for a sub-range of read lengths:
#' lendist_sub95 <- rlength_distr(reads_list, sample = "Samp1", cl = 95)
#' lendist_sub95[["plot"]]
#' @import data.table
#' @import ggplot2
#' @export
rlength_distr <- function(data, sample, transcripts = NULL, cl = 100) {

  if(length(transcripts) == 0) {
    dt <- data[[sample]]
  } else {
    dt <- data[[sample]][transcript %in% transcripts]
  }

  xmin <- quantile(dt$length, (1 - cl/100)/2)
  xmax <- quantile(dt$length, 1 - (1 - cl/100)/2)

  setkey(dt, length)
  dist <- dt[CJ(unique(dt$length)), list(count = .N), by = .EACHI
             ][, percentage := (count / sum(count)) * 100]

  p <- ggplot(dist, aes(as.numeric(length), percentage)) +
    geom_bar(stat = "identity", fill = "gray80") +
    labs(title = sample, x = "Read length", y = "Count (%)") +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(hjust = 0.5), panel.grid.minor.x = element_blank()) +
    scale_x_continuous(limits = c(xmin - 0.5, xmax + 0.5), breaks = seq(xmin + ((xmin) %% 2), xmax, by = max(c(1, floor((xmax - xmin)/7)))))

  output <- list()
  output[["plot"]] <- p
  output[["dt"]] <- dist
  return(output)
}
