#' Identify the ribosome P-site position within the reads.
#'
#' This function identifies within each read the position of the ribosome
#' P-site, determined by the localisation of its first nucleotide. The function
#' processes the samples separately starting from the reads aligning on the
#' reference codon (selected by the user between the start codon and the second
#' to last codon) of any annotated coding sequence. It then returns the position
#' of the P-site specifically inferred for all the read lengths. It also allows
#' to plot a collection of read length-specific occupancy metaprofiles
#' showing the P-sites offsets computed throughout the two steps of the
#' algorithm.
#'
#' @param data A list of data frames from either \code{\link{bamtolist}} or
#'   \code{\link{bedtolist}}.
#' @param flanking An integer that specifies how many nucleotides, at least, of
#'   the reads mapping on the reference codon must flank the reference codon in
#'   both directions. Default is 6.
#' @param start A logical value whether ot not to compute the P-site offsets
#'   starting from the reads aligning on the translation initiation site. FALSE
#'   implies that the reads mapping on the last triplet before the stop codon
#'   are used instead. Default is TRUE.
#' @param extremity A character string specifing which extremity of the reads
#'   should be used in the correction step of the algorithm. It can be either
#'   "5end" or "3end" for the 5' and the 3' extremity, respectively. Default is
#'   "auto", meaning that the best extremity is automatically selected.
#' @param plot A logical value whether or not to plot the occupancy metaprofiles
#'   showing the P-sites offsets computed throughout the two steps of the
#'   algorithm. Default is FALSE.
#' @param plotdir A character string specifying the (existing or not) location
#'   of the directory where the occupancy metaprofiles shuold be stored. This
#'   parameter is considered only if \code{plot} is TRUE. By default this
#'   argument is NULL, which implies it is set as a subfolder of the working
#'   directory, called \emph{offset_plot}.
#' @param plotformat Either "png" (the default) or "pdf", this parameter
#'   specifies the file format of the generated metaprofiles. It is considered
#'   only if \code{plot} is TRUE.
#' @param cl An integer value in \emph{[1,100]} specifying the confidence level
#'   for restricting the generation of the occupancy metaprofiles to a sub-range
#'   of read lengths. By default it is set to 99. This parameter is considered
#'   only if \code{plot} is TRUE.
#' @return A data frame.
#' @examples
#' data(reads_list)
#'
#' ## Compute the P-site offset automatically selecting the otimal read
#' extremity for the correction step and not plotting any metaprofile
#' psite(reads_list, flanking = 6, extremity="auto")
#'
#' ## Compute the P-site offset specifying the extremity used in the correction
#' step and plotting the metaprofiles only for a sub-range of read lengths (the
#' middle 95%). The plots will be placed in the current working directory.
#' psite_offset <- psite(reads_list, flanking = 6, extremity="3end", plot = TRUE, cl = 95)
#' psite_offset
#' @import ggplot2
#' @import cowplot
#' @export
psite <- function(data, flanking = 6, start = TRUE, extremity="auto", plot = FALSE,
                  plotdir = NULL, plotformat="png", cl = 99) {
  names <- names(data)
  offset <- NULL
  for (n in names) { 
    n = "Samp1"
    cat(sprintf("processing %s\n", n))
    df <- data[[n]]
    lev <- sort(unique(df$length))
    
    if(start == T | start == TRUE){
      base <- 0
      df$site.dist.end5 <- df$end5 - df$start_pos
      df$site.dist.end3 <- df$end3 - df$start_pos
    } else {
      base <- -5
      df$site.dist.end5 <- df$end5 - df$stop_pos - base
      df$site.dist.end3 <- df$end3 - df$stop_pos - base
    }
    site.sub <- df[which(df$site.dist.end5 <= -flanking & df$site.dist.end3 >= flanking - 1), ]
    minlen <- min(site.sub$length)
    maxlen <- max(site.sub$length)
    t <- table(factor(site.sub$length, levels = lev))
    
    # offset
    offset.temp <- data.frame(length = as.numeric(as.character(names(t))), percentage = (as.vector(t)/sum(as.vector(t))) * 100)
    rownames(offset.temp) <- offset.temp$length
    offset.temp$around_site <- ifelse(offset.temp$percentage == 0, "F", "T")
    
    offset.temp$offset_from_5 <- as.numeric(as.vector(by(site.sub$site.dist.end5, factor(site.sub$length, levels = lev), function(x) names(which.max(table(x))))))
    offset.temp$offset_from_3 <- as.numeric(as.vector(by(site.sub$site.dist.end3, factor(site.sub$length, levels = lev), function(x) names(which.max(table(x))))))
    
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
      adj_offset_from_3 <- as.numeric(as.character(do.call(rbind, list(by(site.sub, factor(site.sub$length, levels = lev), function(x) {
        t <- table(factor(x$site.dist.end3, levels = seq(min(x$site.dist.end3) - 2, max(x$site.dist.end3))))
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
        adj_offset_from_5 <- as.numeric(as.character(do.call(rbind, list(by(site.sub, factor(site.sub$length, levels = lev), function(x) {
          t <- table(factor(x$site.dist.end5, levels = seq(min(x$site.dist.end5) - 2, max(x$site.dist.end5) + 1 )))
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
    offset.temp$site_percentage <- as.numeric(format(round(offset.temp$percentage, 3), nsmall=4))
    offset.temp$sample <- n
    offset.temp <- offset.temp[, c("length", "total_percentage", "site_percentage", "around_site", "offset_from_5", "offset_from_3", "adj_offset_from_5", "adj_offset_from_3", "sample")]
    if(start == TRUE){
      colnames(offset.temp) <- c("length", "total_percentage", "start_percentage", "around_start", "offset_from_5", "offset_from_3", "adj_offset_from_5", "adj_offset_from_3", "sample")
    } else {
      colnames(offset.temp) <- c("length", "total_percentage", "stop_percentage", "around_stop", "offset_from_5", "offset_from_3", "adj_offset_from_5", "adj_offset_from_3", "sample")
    }
    
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
      minlen <- ceiling(quantile(site.sub$length, (1 - cl/100)/2))
      maxlen <- ceiling(quantile(site.sub$length, 1 - (1 - cl/100)/2))
      for (len in minlen:maxlen) {
        progress <- ceiling(((len + 1 - minlen)/(maxlen - minlen + 1)) * 25)
        cat(sprintf("\rplotting   %s\r", paste(paste(rep(c(" ", "<<", "-"), 
                                                         c(25 - progress, 1, progress)), collapse = ""), " ", as.character(progress*4),
                                               "% ", paste(rep(c("-", ">>", " "), c(progress, 1, 25 - progress)), collapse = ""), sep = "")))
        site.temp <- df[which(df$site.dist.end5 %in% seq(-len + 1, 0) & df$length == len), ]
        site.tab5 <- as.data.frame(table(factor(site.temp$site.dist.end5, levels = (-len + 1):(len))))
        site.temp <- df[which(df$site.dist.end3 %in% seq(0, len - 2) & df$length == len), ]
        site.tab3 <- as.data.frame(table(factor(site.temp$site.dist.end3, levels = (-len):(len - 2))))
        colnames(site.tab5) <- colnames(site.tab3) <- c("distance", "reads")
        site.tab5$distance <- as.numeric(as.character(site.tab5$distance))
        site.tab3$distance <- as.numeric(as.character(site.tab3$distance))
        site.tab5$extremity = "5'end"
        site.tab3$extremity = "3'end"
        
        final.tab <- rbind(site.tab5[site.tab5$distance <= 0, ], site.tab3[site.tab3$distance >= 0, ])
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
          theme(plot.title = element_text(hjust = 0.5))
        
        if(line_plot == "from3"){
          p <- p + geom_vline(xintercept = best.offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best.offset - len + 1, color = "black", linetype = 3, size = 1.1)
        } else {
          p <- p + geom_vline(xintercept = best.offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best.offset + len -1, color = "black", linetype = 3, size = 1.1)
        }
        
        p <- p + 
          scale_x_continuous(breaks = seq(-floor(len/5) * 5, floor(len/5) * 5, 5), labels = as.character(seq(-floor(len/5) * 5, floor(len/5) * 5, 5) + base))
        
        subplotdir <- paste(plotdir, n, sep = "/")
        dir.create(subplotdir)
        save_plot(paste(subplotdir, "/", len, ".", plotformat, sep = ""), p, base_height = 5, base_aspect_ratio = 3)
      }
      cat(sprintf("\rplotting   %s\n",
                  paste(paste(rep(c(" ", "<<", "-"), c(25 - progress, 1, progress)), collapse = ""), " ", 
                        as.character(progress*4), "% ", 
                        paste(rep(c("-", ">>", " "), c(progress, 1, 25 - progress)), collapse = ""), sep = "")))
      options(warn=0)
    }
    
    offset <- rbind(offset, offset.temp)
  }
  return(offset)
}

#' Update reads information according to the inferred P-sites.
#' 
#' Starting ftom the P-site position identfied by \code{\link{psite}}, this
#' function updates the data frames that contains information about the reads.
#' It attaches to the data frames 4 columns reporting the P-site position with
#' respect to the 1st nucleotide of the transcript, the start and the stop codon
#' of the annotated coding sequence (if any) and the region of the transcript
#' (5' UTR, CDS, 3' UTR) that includes the P-site. Please note: if a transcript
#' is not associated to any annotated CDS then the positions of the P-site from
#' both the start and the stop codon is set to NA. If either a FASTA file or a
#' BSgenome data package with the nucleotide  sequences is provided, an
#' additional column reporting the three nucleotides covered by the P-site is
#' attached.
#'
#' @param data A list of data frames from either \code{\link{bamtolist}} or
#'   \code{\link{bedtolist}}.
#' @param offset A data frame from \code{\link{psite}}.
#' @param fastapath A character string specifying the path to the FASTA file
#'   containing the reference nucleotide sequences. Please make sure that the
#'   sequences and their names derive from the same release of and are in
#'   agreement with the annotation file used in the
#'   \code{\link{create_annotation}} function. Either \code{fastapath} or
#'   \code{bsgenome_dp} coulped with \code{txdb} must be specified to attach an
#'   additional column reporting the three nucletotides covered by the
#'   identified P-sites. Default is NULL.
#' @param bsgenome_dp A character string specifying the name of the BSgenome
#'   data package to be loaded. If the specified data package is not already
#'   present in your system, it is installed through the biocLite.R script.
#'   Please check the data packages available in the Bioconductor repositories
#'   for your version of R/Bioconductor using the
#'   \code{\link[BSgenome]{available.genomes}} function from the BSgenome
#'   package. Either \code{fastapath} or \code{bsgenome_dp} coulped with
#'   \code{txdb} must be specified to attach an additional column reporting the
#'   three nucletotides covered by the identified P-sites. Default is NULL.
#' @param txdb A TxDb object. This parameter is considered only if
#'   \code{bsgenome_dp} is specified. Default is NULL.
#' @param granges A logical value whether or not to return a GRangesList object.
#'   Default is FALSE, meaning that a list of data frames (the required input
#'   for the downstream analyses and graphical outputs provided by riboWaltz) is
#'   returned instead.
#' @return A list of data frames or a GRangesList object.
#' @examples
#' data(reads_list)
#' data(psite_offset)
#' data(mm81cdna)
#' 
#' reads_psite_list <- psite_info(reads_list, psite_offset)
#' @import GenomicFeatures
#' @import BSgenome
#' @export
psite_info <- function(data, offset, fastapath = NULL, bsgenome_dp = NULL,
                       txdb = NULL,  granges = FALSE) {
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
    
    if(length(fastapath) != 0 & length(bsgenome_dp) != 0){
      warning("fastapath and bsgenome_dp are both specified. Only fastapath will be considered\n")
      bsgenome_dp = NULL
    }
    
    if(length(bsgenome_dp) != 0 & length(txdb) == 0){
      cat("\n")
      stop("\nERROR: txdb is not specified \n\n")
    }
    
    if(length(fastapath) != 0 | length(bsgenome_dp) != 0){
      if(length(fastapath) != 0) {
        cat("adding codon\n\n")
        sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
      } else {
        if(length(bsgenome_dp) != 0){
          if(bsgenome_dp %in% installed.genomes()){
            library(bsgenome_dp, character.only = TRUE)
          } else {
            source("http://www.bioconductor.org/biocLite.R")
            biocLite(bsgenome_dp, suppressUpdates = TRUE)
            library(bsgenome_dp, character.only = TRUE)
          }
        }
        sequences <- extractTranscriptSeqs(get(bsgenome_dp), txdb, use.names=T)
      }
      df$psite_codon <- as.character(subseq(sequences[as.character(df$transcript)],
                                            start = df$psite,
                                            end = df$psite + 2))
    }
    
    if (granges == T || granges == TRUE) {
      df <- GenomicRanges::makeGRangesFromDataFrame(df,
                                                    keep.extra.columns=TRUE,
                                                    ignore.strand=TRUE,
                                                    seqnames.field=c("transcript"),
                                                    start.field="end5",
                                                    end.field="end3",
                                                    strand.field="strand",
                                                    starts.in.df.are.0based=FALSE)
      strand(df) <- "+"
    }
    
    data[[n]] <- df
  }
  
  if (granges == T || granges == TRUE) {
    data <- GenomicRanges::GRangesList(data)
  }
  
  return(data)
}
