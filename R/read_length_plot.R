#' Plot read length distributions.
#'
#' For a specified sample this function plots the read length distribution. It
#' is possible to visualise the distribution for all the read lengths or to
#' restrict the graphical output to a sub-range of read lengths.
#'
#' @param data A list of data tables from either \code{\link{bamtolist}} or
#'   \code{\link{bedtolist}}.
#' @param sample A character string specifying the name of the sample of
#'   interest.
#' @param cl An integer value in \emph{[1,100]} specifying the confidence level
#'   for restricting the plot to a sub-range of read lengths. By default it is
#'   set to 100, meaning that the whole distribution is displayed.
#' @return A list containing a ggplot2 object, and a data table with the
#'   associated data.
#' @examples
#' data(reads_list)
#'
#' ## Visualise distribution for all the read lengths
#' lendist_whole <- rlength_distr(reads_list, sample = "Samp1", cl = 100)
#' lendist_whole[["plot"]]
#'
#' ## Visualise the metaheatmaps for a sub-range of read lengths (the middle 95%)
#' lendist_sub95 <- rlength_distr(reads_list, sample = "Samp1", cl = 95)
#' lendist_sub95[["plot"]]
#' @import data.table
#' @import ggplot2
#' @export
rlength_distr <- function(data, sample, cl = 100) {
  dt <- data[[sample]]
  
  xmin <- quantile(dt$length, (1 - cl/100)/2)
  xmax <- quantile(dt$length, 1 - (1 - cl/100)/2)
  
  setkey(dt, length)
  dist <- dt[CJ(min(dt$length) : max(dt$length)), list(count = .N), by = length
             ][, count := (count / sum(count)) * 100]

  p <- ggplot(dist, aes(as.numeric(length), count)) +
    geom_bar(stat = "identity", fill = "gray80") +
    labs(title = sample, x = "Read length", y = "Count (%)") +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(xmin - 0.5, xmax + 0.5), breaks = seq(xmin + ((xmin) %% 2), xmax, by = max(c(1, floor((xmax - xmin)/7)))))
  
  output<-list()
  output[["plot"]] <- p
  output[["dt"]] <- dist
  return(output)
}
