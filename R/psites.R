#' Identify the ribosome P-site position along the reads.
#' 
#' Identifies, within each read, the position of the ribosome P-site, determined
#' by the localisation of its first nucleotide. The function returns the 
#' position of the P-site specifically inferred for the many read lengths by 
#' processing the samples separately. It also allows to plot a collection of
#' length-specific occupancy metaprofiles of the read ends aligning around the
#' start codon, displaying the identified P-sites offsets.
#'
#' @param data A list of data frames from \code{\link{bamtolist}}.
#' @param flanking An integer that specifies, for all the reads aligning on the
#'   start codon, the minimum number of nucleotides that must flank the
#'   translation initiation site in both directions. Default is 6.
#' @param extremity A character string specifing which extremity of the reads
#'   should be used for the correction step of the P-site identification. It can
#'   be either "5end" or "3end" for the 5' and the 3' extremity, respectively.
#'   Default is "auto", meaning that the best extremity is automatically chosen.
#' @param plot A logical value whether or not to plot the occupancy metaprofiles
#'   for the computation of the offsets. Default is FALSE.
#' @param plotdir A character string specifying the (existing or not) folder
#'   where the occupancy metaprofiles shuold be saved. This parameter is
#'   considered only if \code{plot} is TRUE. By default this argument is NULL,
#'   which implies the folder is set as a subfolder of the working directory,
#'   called \emph{offset_plot}.
#' @param plotformat Either "png" (the default) or "pdf". This parameter
#'   specifies the file format in which the occupancy metaprofiles shuold be
#'   saved. It is considered only if \code{plot} is TRUE.
#' @param cl An integer with value in \emph{[1,100]} specifying the read length
#'   confidence level for restricting the plot of the occupancy metaprofiles to
#'   a range of lengths. By default it is set to 99. This parameter is
#'   considered only if \code{plot} is TRUE.
#' @details This function compute the P-site identification starting from the 
#'   reads that align on the start codon of any annotated coding sequences,
#'   exploiting the knowledge that their associated P-sites corresponds to the
#'   triplet AUG. The P-site identification is then divided in two steps: i)
#'   computation of the offsets between the extremities of the reads and the
#'   start codons based on the alignment of 5' and the 3' end around the
#'   translation initiation site ii) correction of some offsets based on the
#'   global results of the previous step.
#' @return A data frame.
#' @examples
#' data(reads_list)
#'
#' ## Compute the P-site offset automatically not plotting the metaprofiles
#' based on the alignment of the read ends
#' psite(reads_list, flanking = 6, extremity="auto")
#'
#' ## Compute the P-site offset specifying the extremity used for the correction
#' step, plotting the metaprofiles based on the alignment of the read ends, only
#' for the middle 95% of the read length. The plots will be placed in the
#' current working directory.
#' psite_offset <- psite(reads_list, flanking = 6, extremity="3end",plot = TRUE, cl = 95)
#' psite_offset
#' @import ggplot2
#' @import cowplot
#' @export
psite <- function(data, flanking = 6, extremity="auto", plot = FALSE,
                  plotdir = NULL, plotformat="png", cl = 99) {
  names <- names(data)
  offset <- NULL
  for (n in names) {
    cat(sprintf("processing %s\n", n))
    df <- data[[n]]
    lev <- sort(unique(df$length))
    df$start.dist.end5 <- df$end5 - df$start_pos
    df$stop.dist.end5 <- df$end5 - df$stop_pos
    df$start.dist.end3 <- df$end3 - df$start_pos
    df$stop.dist.end3 <- df$end3 - df$stop_pos
    start.sub <- df[which(df$start.dist.end5 <= -flanking & df$start.dist.end3 >= flanking - 1), ]
    minlen <- min(start.sub$length)
    maxlen <- max(start.sub$length)
    t <- table(factor(start.sub$length, levels = lev))

    # offset
    offset.temp <- data.frame(length = as.numeric(as.character(names(t))), percentage = (as.vector(t)/sum(as.vector(t))) * 100)
    rownames(offset.temp) <- offset.temp$length
    offset.temp$around_start <- ifelse(offset.temp$percentage == 0, "F", "T")

    offset.temp$offset_from_5 <- as.numeric(as.vector(by(start.sub$start.dist.end5, factor(start.sub$length, levels = lev), function(x) names(which.max(table(x))))))
    offset.temp$offset_from_3 <- as.numeric(as.vector(by(start.sub$start.dist.end3, factor(start.sub$length, levels = lev), function(x) names(which.max(table(x))))))

    # adjusted offset
    best.offset.from3.tab <- by(offset.temp$percentage, offset.temp$offset_from_3, function(x) sum(x))
    best.offset.from5.tab <- by(offset.temp$percentage, offset.temp$offset_from_5, function(x) sum(x))
    best.offset.from3 <- names(which.max(best.offset.from3.tab))
    best.offset.from5 <- names(which.max(best.offset.from5.tab))
    if((extremity == "auto" &
         best.offset.from3.tab[best.offset.from3] > best.offset.from5.tab[best.offset.from5] &
         as.numeric(best.offset.from3) <= minlen - 2) |
        (extremity == "auto" &
         best.offset.from3.tab[best.offset.from3] <= best.offset.from5.tab[best.offset.from5] &
         as.numeric(best.offset.from5) > minlen - 1) |
         extremity == "3end") {
      best.offset <- as.numeric(best.offset.from3)
      line_plot <- "from3"
      cat(sprintf("best offset: %i nts from the 3' end\n", best.offset))
      adj_offset_from_3 <- as.numeric(as.character(do.call(rbind, list(by(start.sub, factor(start.sub$length, levels = lev), function(x) {
        t <- table(factor(x$start.dist.end3, levels = seq(min(x$start.dist.end3) - 2, max(x$start.dist.end3))))
        t[1:2]<-t[3]+1
        locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) == -2)]))) + 1
        adjoff <- locmax[which.min(abs(locmax - best.offset))]
        ifelse(length(adjoff) != 0, as.numeric(as.character(adjoff)), best.offset)
      })))))
      adj_offset_from_3[is.na(adj_offset_from_3)] <- best.offset
      offset.temp$adj_offset_from_5 <- - adj_offset_from_3 + offset.temp$length - 1
      offset.temp$adj_offset_from_3 <- adj_offset_from_3
    } else {
      if((extremity == "auto" &
            best.offset.from3.tab[best.offset.from3] <= best.offset.from5.tab[best.offset.from5] &
            as.numeric(best.offset.from5) <= minlen - 1) |
           (extremity == "auto" &
            best.offset.from3.tab[best.offset.from3] > best.offset.from5.tab[best.offset.from5] &
            as.numeric(best.offset.from3) > minlen - 2) |
           extremity == "5end") {
        best.offset <- as.numeric(best.offset.from5)
        line_plot <- "from5"
        cat(sprintf("best offset: %i nts from the 5' end\n", -best.offset))
        adj_offset_from_5 <- as.numeric(as.character(do.call(rbind, list(by(start.sub, factor(start.sub$length, levels = lev), function(x) {
          t <- table(factor(x$start.dist.end5, levels = seq(min(x$start.dist.end5) - 2, max(x$start.dist.end5) + 1 )))
          t[1:2]<-t[3]+1
          locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) == -2)]))) + 1
          adjoff <- locmax[which.min(abs(locmax - best.offset))]
          ifelse(length(adjoff) != 0, as.numeric(as.character(adjoff)), best.offset)
        })))))
        adj_offset_from_5[is.na(adj_offset_from_5)] <- best.offset
        offset.temp$adj_offset_from_5 <- abs(adj_offset_from_5)
        offset.temp$adj_offset_from_3 <- abs(offset.temp$adj_offset_from_5 - offset.temp$length + 1)
      }
    }

    t <- table(factor(df$length, levels = lev))
    offset.temp$offset_from_5 <- - offset.temp$offset_from_5
    offset.temp$total_percentage <- as.numeric(format(round((as.vector(t)/sum(as.vector(t))) * 100, 3), nsmall=4))
    offset.temp$start_percentage <- as.numeric(format(round(offset.temp$percentage, 3), nsmall=4))
    offset.temp$sample <- n
    offset.temp <- offset.temp[, c("length", "total_percentage", "start_percentage", "around_start", "offset_from_5", "offset_from_3", "adj_offset_from_5", "adj_offset_from_3", "sample")]

    # plot
    if (plot == T || plot == TRUE) {
      options(warn=-1)
      if (length(plotdir) == 0) {
        dir <- getwd()
        plotdir <- paste(dir, "/offset_plot", sep = "")
      }
      if (!dir.exists(plotdir)) {
        dir.create(plotdir)
      }
      minlen <- ceiling(quantile(start.sub$length, (1 - cl/100)/2))
      maxlen <- ceiling(quantile(start.sub$length, 1 - (1 - cl/100)/2))
      for (len in minlen:maxlen) {
        progress <- ceiling(((len + 1 - minlen)/(maxlen - minlen + 1)) * 100)
        cat(sprintf("plotting   %s\r", paste("|", paste(rep("-", progress/2), collapse = ""), ">  ", as.character(progress),
                                             "%", sep = "")))
        start.sub <- df[which(df$start.dist.end5 %in% seq(-len + 1, 0) & df$length == len), ]
        start.tab5 <- as.data.frame(table(factor(start.sub$start.dist.end5, levels = (-len + 1):(len))))
        start.sub <- df[which(df$start.dist.end3 %in% seq(0, len - 2) & df$length == len), ]
        start.tab3 <- as.data.frame(table(factor(start.sub$start.dist.end3, levels = -len:(len - 2))))
        colnames(start.tab5) <- colnames(start.tab3) <- c("distance", "reads")
        start.tab5$distance <- as.numeric(as.character(start.tab5$distance))
        start.tab3$distance <- as.numeric(as.character(start.tab3$distance))
        start.tab5$extremity = "5'end"
        start.tab3$extremity = "3'end"

        final.tab <- rbind(start.tab5[start.tab5$distance <= 0, ], start.tab3[start.tab3$distance >= 0, ])
        final.tab$extremity <- factor(final.tab$extremity, levels = c("5'end", "3'end"))

        p <- ggplot(final.tab, aes(distance, reads, color = extremity)) +
          geom_line() +
          geom_vline(xintercept = seq(-round(len/3) * 3, round(len/3) * 3, 3), linetype = 2, color = "gray90") +
          geom_vline(xintercept = 0, color = "gray50") +
          geom_vline(xintercept = - offset.temp[as.character(len), "offset_from_5"], color = "#D55E00", linetype = 2, size = 1.1) +
          geom_vline(xintercept = offset.temp[as.character(len), "offset_from_3"], color = "#56B4E9", linetype = 2, size = 1.1) +
          geom_vline(xintercept = - offset.temp[as.character(len), "adj_offset_from_5"], color = "#D55E00", size = 1.1) +
          geom_vline(xintercept = offset.temp[as.character(len), "adj_offset_from_3"], color = "#56B4E9", size = 1.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - len, xmax = -flanking , fill = "#D55E00", alpha = 0.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - 1 , xmax = len - flanking - 1, fill = "#56B4E9", alpha = 0.1) +
          labs(x = "Distance from start (nt)", y = "Number of read extremities", title = paste(n, " - length=", len, " nts", sep = ""), color="Extremity") +
          theme_bw(base_size = 20) +
          scale_fill_discrete("") +
          theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_x_continuous(breaks = seq(-floor(len/5) * 5, floor(len/5) * 5, 5))

        if(line_plot == "from3"){
          p <- p + geom_vline(xintercept = best.offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best.offset - len + 1, color = "black", linetype = 3, size = 1.1)
        } else {
          p <- p + geom_vline(xintercept = best.offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best.offset + len -1, color = "black", linetype = 3, size = 1.1)
        }

        subplotdir <- paste(plotdir, n, sep = "/")
        dir.create(subplotdir)
        save_plot(paste(subplotdir, "/", len, ".", plotformat, sep = ""), p, base_height = 5, base_aspect_ratio = 3)
      }
      cat(sprintf("plotting   %s\n", paste("|", paste(rep("-", progress/2), collapse = ""), ">  ", as.character(progress),
                                           "%", sep = "")))
      options(warn=0)
    }

    offset <- rbind(offset, offset.temp)
  }
  return(offset)
}

#' Updates reads information adding features associated to the inferred P-sites.
#' 
#' Updates the data frames of the input list with the P-site position identfied 
#' by \code{\link{psite}}. It attaches to the data frames 4 columns containing 
#' the P-site position with respect to the beginning of the transcript and to 
#' both its start and stop codon and the region of the transcript (5' UTR, CDS, 
#' 3' UTR) that includes the P-site. Please note: if a transcript doesn't 
#' present any annotated CDS then the positions of the P-site from both the 
#' start and the stop codon will be set to NA. If the FASTA file with the
#' transcriptome sequences is provided, an additional column containing the
#' three nucleotides covered by the P-site will be attached.
#'
#' @param data A list of data frames from \code{\link{bamtolist}}
#' @param offset A data frame from \code{\link{psite}}.
#' @param fastapath A character string specifying the path to the FASTA file
#'   containing the nucleotide sequence of the transcripts. For each mRNA the 
#'   the record description line must contain the transcript name as in the 
#'   reference transcriptome and the sequence must derive from the same relase 
#'   of the genome. Use this argument to attach to the data frames an additional
#'   column containing the three nucletotides corresponding to the identified
#'   P-sites. Default is NULL.
#' @return A list of data frames.
#' @examples
#' data(reads_list)
#' data(psite_offset)
#' data(mm81cdna)
#'
#' reads_psite_list <- psite_info(reads_list, psite_offset)
#' @export
psite_info <- function(data, offset, fastapath = NULL) {
  names <- names(data)
  for (n in names) {
    cat(sprintf("processing %s\n", n))
    df <- data[[n]]
    suboff <- offset[which(offset$sample == n), ]
    cat("adding p-site position\n")
    df <- dplyr::left_join(df, suboff[, c("length", "adj_offset_from_3")], by = "length")
    colnames(df)[colnames(df) == "adj_offset_from_3"] <- "psite"
    df$psite <- df$end3 - df$psite
    df <- df[, c("transcript", "end5", "psite", "end3", "length", "start_pos", "stop_pos")]
    df$psite_from_start <- ifelse(df$stop_pos == 0, 0, df$psite - df$start_pos)
    df$psite_from_stop <- ifelse(df$stop_pos == 0, 0, df$psite - df$stop_pos)
    cat("adding region\n")
    df$psite_region <- ifelse(df$stop_pos == 0,
                              NA,
                              ifelse(df$psite_from_start >= 0
                                     & df$psite_from_stop <= 0,
                                     "cds",
                                     ifelse(df$psite_from_start < 0
                                            & df$psite_from_stop < 0,
                                            "5utr",
                                            "3utr")))
    if(length(fastapath) != 0) {
      cat("adding codon\n\n")
      sequences_biost <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
      df$psite_codon <- as.character(subseq(sequences_biost[as.character(df$transcript)],
                                            start = df$psite,
                                            end = df$psite + 2))
    }

    data[[n]] <- df
  }
  return(data)
}
