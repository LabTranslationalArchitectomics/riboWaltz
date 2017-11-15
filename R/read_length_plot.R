#' Plot read length distributions.
#'
#' Plots the read length distribution for a specified sample of the input list. 
#' It is possible to visualise the whole length distribution or to restrict it
#' to a chosen range of read lengths.
#'
#' @param data A list of data frames from \code{\link{bamtolist}}.
#' @param sample A character string specifying the name of the sample of
#'   interest.
#' @param cl An integer with value in \emph{[1,100]} specifying the read length 
#'   confidence level for restricting the distribution to a chosen range of 
#'   lengths. By default it is set to 100, i.e. the whole distribution will be
#'   displayed.
#' @return A list containing a ggplot2 object, and a data frame with the
#'   associated data.
#' @examples
#' data(reads_list)
#'
#' ## Visualise the whole read length distribution
#' lendist_whole <- rlength_distr(reads_list, sample = "Samp1", cl = 100)
#' lendist_whole[["plot"]]
#'
#' ## Visualise the middle 95% of the read length distribution
#' lendist_sub95 <- rlength_distr(reads_list, sample = "Samp1", cl = 95)
#' lendist_sub95[["plot"]]
#' @import ggplot2
#' @export
rlength_distr <- function(data, sample, cl = 100) {
  df <- data[[sample]]
  dist <- table(factor(df$length, levels = c(min(df$length):max(df$length))))
  dist <- data.frame(length = names(dist), count = as.vector((dist/sum(dist)) * 100))

  xmin <- quantile(df$length, (1 - cl/100)/2)
  xmax <- quantile(df$length, 1 - (1 - cl/100)/2)

  p <- ggplot(dist, aes(as.numeric(as.character(length)), count)) +
    geom_bar(stat = "identity", fill = "gray80") +
    labs(title = sample, x = "Read length", y = "Count (%)") +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(xmin-0.5, xmax+0.5), breaks = seq(xmin + ((xmin) %% 2), xmax, by=floor((xmax-xmin)/7)))

  output<-list()
  output[["plot"]]<-p
  output[["df"]]<-dist
  return(output)
}