#' Plot metaheatmaps of the read ends.
#'
#' Plots 4 metaheatmaps that show the abundance of the 5' the and 3' end of the
#' reads mapping around the start and the stop codons, stratified by their
#' length. It is possible to visualise the heatmaps for all the reads or to
#' restrict the graphical output to a sub-range of read lengths.
#'
#' @param data A list of data frames from \code{\link{bamtolist}}.
#' @param annotation A data frame with a reference annotation of the transripts.
#'   It must contain at least four columns named \emph{transcript},
#'   \emph{l_utr5}, \emph{l_cds}, \emph{l_utr3} containing the name of the
#'   transcripts (the same as in the reference transcriptome), the position of
#'   the first nucleotide of the \emph{5' UTR}, the \emph{CDS} and the  \emph{3'
#'   UTR}, respectively. No specific order is required.
#' @param sample A character string specifying the name of the sample of
#'   interest.
#' @param cl An integer with value in \emph{[1,100]} specifying the read length
#'   confidence level for restricting the distribution to a chosen range of
#'   lengths. By default it is set to 99.
#' @param utr5l A positive integer specifying the length (in nucleotides) of the
#'   5' UTR portion that will flank the start codon in the plot. The default
#'   value is 50.
#' @param cdsl A positive integer specifying the length (in nucleotides) of the
#'   coding sequence portion that will flank both the start and the stop codon
#'   in the plot. The default value is 50.
#' @param utr3l A positive integer specifying the length (in nucleotides) of the
#'   3' UTR portion that will flank the stop codon in the plot. The default
#'   value is 50.
#' @param log A logical value whether or not to use a logarithmic scale colour
#'   (it is suggested for data with a strong difference between the lowest and
#'   the highest signal. Default is FALSE.
#' @param colour A character string specifying the colour to be used for the
#'   plot.
#' @return A list containing a ggplot2 plot object, and a data frame with the
#'   associated data.
#' @examples
#' data(reads_list)
#'
#' ## Visualise the heatmap for the whole range of read lengths
#' heatend_whole <- rends_heat(reads_list, sample = "Samp1", cl = 100)
#' heatend_whole[["plot"]]
#'
#' ## Visualise the heatmap for the middle 95% of the range of read lengths and
#' reducing the regions flanking the start and the stop codons
#' heatend_sub95 <- rends_heat(reads_list, sample = "Samp1", cl = 95,
#' utr5l = 30, cdsl = 40, utr3l = 30)
#' heatend_sub95[["plot"]]
#' @import ggplot2
#' @export
rends_heat <- function(data, annotation, sample, cl = 99, utr5l = 50, cdsl = 50, utr3l = 50,
                       log = F, colour = "black") {
  df <- data[[sample]]
  df$start.dist.end5 <- df$end5 - df$start_pos
  df$stop.dist.end5 <- df$end5 - df$stop_pos
  df$start.dist.end3 <- df$end3 - df$start_pos
  df$stop.dist.end3 <- df$end3 - df$stop_pos
  minlen <- ceiling(quantile(df$length, (1 - cl/100)/2))
  maxlen <- ceiling(quantile(df$length, 1 - (1 - cl/100)/2))
  
  rownames(annotation) <- as.character(annotation$transcript)
  l.transcripts <- rownames(annotation)[which(annotation$l_utr5 >= utr5l &
                                                annotation$l_cds >= 2 * (cdsl + 1) &
                                                annotation$l_utr3 >= utr3l)]
  c.transcripts <- intersect(unique(df$transcript), l.transcripts)
  
  # 5' end
  start.sub <- df[which(df$transcript %in% c.transcripts & df$start.dist.end5 %in% seq(-utr5l, cdsl)), ]
  stop.sub <- df[which(df$transcript %in% c.transcripts & df$stop.dist.end5 %in% seq(-cdsl, utr3l)), ]
  start.tab <- aggregate(rep(1, nrow(start.sub)), by = list(x = start.sub$length, y = start.sub$start.dist.end5), sum, drop = F)
  stop.tab <- aggregate(rep(1, nrow(stop.sub)), by = list(x = stop.sub$length, y = stop.sub$stop.dist.end5), sum, drop = F)
  colnames(start.tab) <- colnames(stop.tab) <- c("length", "dist", "count")
  start.tab$region <- "start"
  stop.tab$region <- "stop"
  final.tab5 <- rbind(start.tab, stop.tab)
  final.tab5$end <- "5end"
  
  # 3' end
  start.sub <- df[which(df$transcript %in% c.transcripts & df$start.dist.end3 %in% seq(-utr5l, cdsl)), ]
  stop.sub <- df[which(df$transcript %in% c.transcripts & df$stop.dist.end3 %in% seq(-cdsl, utr3l)), ]
  start.tab <- aggregate(rep(1, nrow(start.sub)), by = list(x = start.sub$length, y = start.sub$start.dist.end3), sum, drop = F)
  stop.tab <- aggregate(rep(1, nrow(stop.sub)), by = list(x = stop.sub$length, y = stop.sub$stop.dist.end3), sum, drop = F)
  colnames(start.tab) <- colnames(stop.tab) <- c("length", "dist", "count")
  start.tab$region <- "start"
  stop.tab$region <- "stop"
  final.tab3 <- rbind(start.tab, stop.tab)
  final.tab3$end <- "3end"
  
  final.tab <- rbind(final.tab5, final.tab3)
  final.tab$region <- factor(final.tab$region, levels = c("start", "stop"), labels = c("Distance from start (nt)", "Distance from stop (nt)"))
  final.tab$end <- factor(final.tab$end, levels = c("5end", "3end"), labels = c("5' end", "3' end"))
  
  max<-max(final.tab$count)
  p <- ggplot(final.tab, aes(dist, length)) +
    geom_tile(aes(fill = count)) +
    labs(title = paste(sample, "5' / 3' read end metaheatmaps", sep = " - "), y = "Read length") +
    theme_bw(base_size = 20) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.title.x = element_blank()) +
    facet_grid(end ~ region, scales = "free", switch = "x") +
    theme(strip.background = element_blank()) +
    scale_y_continuous(limits = c(minlen-0.5, maxlen+0.5), breaks = seq(minlen + ((minlen) %% 2), maxlen, by=max(2,floor((maxlen-minlen)/7)))) +
    geom_vline(xintercept = 0, linetype = 2, color = "red")
  
  if (log == F) {
    p <- p +
      scale_fill_gradient("Number\nof read\nextremities\n", low = "white", high = colour, limits = c(0.1, max), breaks = c(0.1, max/2, max), labels = c("0", floor(max/2), floor(max)), na.value = "white")
  } else {
    p <- p +
      scale_fill_gradient("Number\nof read\nextremities\n", low = "white", high = colour, limits = c(0.1, max), breaks = c(0.1, 10^(log10(max)/2 - 0.5), floor(max)), labels = c("0", floor(10^(log10(max)/2 - 0.5)), floor(max)), trans = "log", na.value = "transparent")
  }
  
  output<-list()
  output[["plot"]]<-p
  output[["df"]]<-final.tab
  return(output)
  return(p)
}
