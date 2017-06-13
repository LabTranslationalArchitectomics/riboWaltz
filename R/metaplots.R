#' Plot ribosome occupancy metaprofiles at single-nucleotide resolution.
#'
#' Builds a metaprofile based on the identified P-site of the reads. It sums up
#' the number of P-sites per nucleotide for the chosen region, starting from one
#' ore more replicates of a sample.
#'
#' @param data A list of data frames from \code{\link{psite_info}}.
#' @param annotation A data frame with a reference annotation of the transripts.
#'   It must contain at least four columns named \emph{transcript},
#'   \emph{l_utr5}, \emph{l_cds}, \emph{l_utr3} containing the name of the
#'   transcripts (the same as in the reference transcriptome), the position of
#'   the first nucleotide of the \emph{5' UTR}, the \emph{CDS} and the  \emph{3'
#'   UTR}, respectively. No specific order is required.
#' @param sample A character string vector specifying the name of the sample (or
#'   of its replicates) of interest.
#' @param scale_factors A numeric vector of scale factors to be used for merging
#'   the replicates of the specified sample (if any). The vector must contain at
#'   least a set of values named after the strings listed in \code{sample}. No
#'   specific order is required.
#' @param length_range Either "all", an integer or an integer vector. In the
#'   first case all the reads will be used to generate the metaprofile.
#'   Otherwise, only the reads matching the specified length(s) will be
#'   employed.
#' @param transcripts A character string vector specifying the name of the
#'   transcripts to be included in the metaprofile. By default this argument is
#'   NULL, which implies all the transcripts in \code{data} will be used. Note
#'   that if either the 5' UTR, the coding sequence or the 3' UTR of a
#'   transcript is shorther than what is specified by \code{utr5l},
#'   \eqn{2*}\code{cdsl} and \code{utr3l} respectively, the transcript will not
#'   be cosidered.
#' @param utr5l A positive integer specifying the length (in nucleotides) of the
#'   5' UTR portion that will flank the start codon in the plot. The default
#'   value is 25.
#' @param cdsl A positive integer specifying the length (in nucleotides) of the
#'   coding sequence portion that will flank both the start and the stop codon
#'   in the plot. The default value is 50.
#' @param utr3l A positive integer specifying the length (in nucleotides) of the
#'   3' UTR portion that will flank the stop codon in the plot. The default
#'   value is 25.
#' @param plot_title Either "auto" (the default), NULL or any character string.
#'   When \emph{auto}, the title of the plot is the name of the sample followed
#'   by the number of transcript and the read lengths employed for generating
#'   the metaprofile.
#' @return A list containing a ggplot2 plot object, a data frame with the
#'   associated data and the transcripts employed for generating the plot.
#' @examples
#' data(reads_psite_list)
#' data(mm81cdna)
#'
#' ## Generate the metaprofile employing the whole dataset
#' metaprof_whole <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1")
#' metaprof_whole[["plot"]]
#'
#' ## Generate the metaprofile employing reads of 27, 28 and 29 nucleotides and
#' a specified set of transcripts (for example those with at least one P-site on
#' the first nucleotide of the coding sequence)
#' sample_name <- "Samp1"
#' sub_reads_psite_list <- subset(reads_psite_list[[sample_name]], psite_from_start == 0)
#' transcript_names <- as.character(sub_reads_psite_list$transcript)
#' metaprof_sub <- metaprofile_psite(reads_psite_list, mm81cdna, sample = sample_name,
#' length_range = 27:29, transcripts = transcript_names)
#' metaprof_sub[["plot"]]
#' @import ggplot2
#' @export
metaprofile_psite <- function(data, annotation, sample, scale_factors = NULL,
                              length_range = "all", transcripts = NULL,
                              utr5l = 25, cdsl = 50, utr3l = 25,
                              plot_title = "auto") {
  rownames(annotation) <- as.character(annotation$transcript)
  l.transcripts <- rownames(annotation)[which(annotation$l_utr5 >= utr5l &
                                                annotation$l_cds >= 2 * (cdsl + 1) &
                                                annotation$l_utr3 >= utr3l)]
  if (length(transcripts) == 0) {
    c.transcripts <- l.transcripts
    ntr <- length(c.transcripts)
  } else {
    c.transcripts <- intersect(l.transcripts, transcripts)
    ntr <- length(transcripts)
  }

  if(!identical(length_range, "all") & !inherits(length_range, "numeric") & !inherits(length_range, "integer")){
    warning("length_range is invalid. Set to default \"all\"\n")
    length_range="all"
  }

  for (samp in sample) {
    df <- data[[samp]][which(data[[samp]]$transcript %in% c.transcripts), ]
    if (identical(length_range, "all")) {
      start.sub <- df[which(df$psite_from_start %in% seq(-utr5l, cdsl)), ]
      stop.sub <- df[which(df$psite_from_stop %in% seq(-cdsl, utr3l)), ]
    } else {
      start.sub <- df[which(df$psite_from_start %in% seq(-utr5l, cdsl) & df$length%in%length_range), ]
      stop.sub <- df[which(df$psite_from_stop %in% seq(-cdsl, utr3l) & df$length%in%length_range), ]
    }
    start.tab <- as.data.frame(table(factor(start.sub$psite_from_start, levels = -utr5l:cdsl)))
    start.tab$reg <- "start"
    stop.tab <- as.data.frame(table(factor(stop.sub$psite_from_stop, levels = -cdsl:utr3l)))
    stop.tab$reg <- "stop"
    colnames(start.tab) <- colnames(stop.tab) <- c("distance", "reads", "reg")
    samp.tab <- rbind(start.tab, stop.tab)

    if (length(scale_factors) != 0 & length(scale_factors) == length(sample)) {
      samp.tab$reads <- samp.tab$reads * scale_factors[samp]
    }

    if (exists("final.tab.psitemetaprofile")) {
      final.tab.psitemetaprofile$reads <- final.tab.psitemetaprofile$reads + samp.tab$reads
    } else {
      final.tab.psitemetaprofile <- samp.tab
    }
  }
  final.tab.psitemetaprofile$reg <- factor(final.tab.psitemetaprofile$reg, levels = c("start", "stop"), labels = c("Distance from start (nt)", "Distance from stop (nt)"))

  linestart <- data.frame(reg = rep(c("Distance from start (nt)", "Distance from stop (nt)"), times = c(length(c(rev(seq(-3, -utr5l, -3)), seq(3, cdsl, 3))), length(c(rev(seq(-2, -cdsl, -3)), seq(1, utr3l, 3))))), line = c(rev(seq(-3, -utr5l, -3)), seq(3, cdsl, 3), rev(seq(-2, -cdsl, -3)), seq(1, utr3l, 3)))
  linered <- data.frame(reg = c("Distance from start (nt)", "Distance from stop (nt)"), line =c(0, 1))
  
  if(length(plot_title) != 0 && plot_title == "auto"){
    if (identical(length_range, "all")) {
      plot_title <- paste(paste(sample, collapse = "+"), " (", ntr, " tr)", sep = "")
    } else {
      if(min(length_range)==max(length_range)) {
        plot_title <- paste(paste(sample, collapse = "+"), " (", ntr, " tr) - Read length: ", min(length_range), " nts", sep = "")
      } else {
        if(identical(length_range, min(length_range):max(length_range)) | identical(length_range, seq(min(length_range),max(length_range),1))){
          plot_title <- paste(paste(sample, collapse = "+"), " (", ntr, " tr) - Read lengths: ", min(length_range), "-", max(length_range), " nts",sep = "")
        } else {
          plot_title <- paste(paste(sample, collapse = "+"), " (", ntr, " tr) - Read lengths: ", paste(length_range, collapse=","), " nts",sep = "")
        }
      }
    }
  }

  plot <- ggplot(final.tab.psitemetaprofile, aes(as.numeric(as.character(distance)), reads)) +
    geom_line() +
    geom_vline(data = linered, aes(xintercept = line), linetype = 1, color = "red") +
    labs(x = "", y = "P-site", title = plot_title) +
    theme_bw(base_size = 20) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    facet_grid(. ~ reg, scales = "free", switch = "x") +
    theme(strip.background = element_blank()) +
    geom_vline(data = linestart, aes(xintercept = line), linetype = 3, color = "gray60")
  plot

  output<-list()
  output[["plot"]]<-plot
  output[["df"]]<-final.tab.psitemetaprofile
  output[["transcripts"]] <- c.transcripts
  return(output)
}

#' Plot ribosome occupancy metaheatmap at single-nucleotide resolution.
#'
#' Plots a heatmap-like metaprofile based on the identified P-site of the reads.
#' It is similar to \code{\link{metaprofile_psite}}, but here the intensity of
#' the signal along the chosen region is represented by a continuous color scale
#' rather than by the height of a line chart. This graphical output is the
#' optimal choice in order to compare multiple samples or look at the behaviour
#' of the data if different populations of reads are considered.
#'
#' @param data A list of data frames from \code{\link{psite_info}}.
#' @param annotation A data frame with a reference annotation of the transripts.
#'   It must contain at least four columns named \emph{transcript},
#'   \emph{l_utr5}, \emph{l_cds}, \emph{l_utr3} containing the name of the
#'   transcripts (the same as in the reference transcriptome), the position of
#'   the first nucleotide of the \emph{5' UTR}, the \emph{CDS} and the  \emph{3'
#'   UTR}, respectively. No specific order is required.
#' @param sample A list of character string vectors specifying the name of the
#'   samples of interest. The samples contained in each elements of the list
#'   will be merged toghether, using the scale factors specified by
#'   \code{scale_factors} (if any). The name of the list elements will be used
#'   as labels for the raws of the heatmaps.
#' @param scale_factors A numeric vector of scale factors to be used for merging
#'   the replicates of the specified sample (if any). The vector must contain at
#'   least a set of values named after the strings listed in \code{sample}. No
#'   specific order is required.
#' @param length_range Either "all", an integer or an integer vector. In the
#'   first case all the reads will be used to generate the metaheatmap.
#'   Otherwise, only the reads matching the specified length(s) will be
#'   employed.
#' @param transcripts A character string vector specifying the name of the
#'   transcripts to be included in the metaprofile. By default this argument is
#'   NULL, which implies all the transcripts in \code{data} will be used. Note
#'   that if either the 5' UTR, the coding sequence or the 3' UTR of a
#'   transcript is shorther than what is specified by \code{utr5l},
#'   \eqn{2*}\code{cdsl} and \code{utr3l} respectively, the transcript will not
#'   be cosidered.
#' @param utr5l A positive integer specifying the length (in nucleotides) of the
#'   5' UTR portion that will flank the start codon in the plot. The default
#'   value is 25.
#' @param cdsl A positive integer specifying the length (in nucleotides) of the
#'   coding sequence portion that will flank both the start and the stop codon
#'   in the plot. The default value is 50.
#' @param utr3l A positive integer specifying the length (in nucleotides) of the
#'   3' UTR portion that will flank the stop codon in the plot. The default
#'   value is 25.
#' @param log A logical value whether or not to use a logarithmic scale colour
#'   (it is suggested for data with a strong difference between the lowest and
#'   the highest signal. Default is FALSE.
#' @param colour A character string specifying the colour to be used for the
#'   plot.
#' @param plot_title Either NULL (the default), "auto" or any character string.
#'   When \emph{auto}, the title of the plot reports the number of transcript
#'   and the read lengths employed for generating the heatmap.
#' @return A list containing a ggplot2 plot object, a data frame with the
#'   associated data and the transcripts employed for generating the plot.
#' @examples
#' data(reads_psite_list)
#'
#' ## Generate the metaheatmap employing the whole dataset
#' metaheat_whole <- metaheatmap_psite(reads_psite_list, mm81cdna, sample = list("Whole"=c("Samp1")))
#' metaheat_whole[["plot"]]
#'
#' ## Generate the metaheatmap employing reads of 27, 28 and 29 nucleotides and
#' a specified set of transcripts (for example those with at least one P-site on
#' the first nucleotide of the coding sequence)
#' sample_name <- "Samp1"
#' sub_reads_psite_list <- subset(reads_psite_list[[sample_name]], psite_from_start == 0)
#' transcript_names <- as.character(sub_reads_psite_list$transcript)
#' metaheat_sub <- metaheatmap_psite(reads_psite_list, mm81cdna, sample = list("sub"=sample_name),
#' length_range = 27:29, transcripts = transcript_names, plot_title = "auto")
#' metaheat_sub[["plot"]]
#'
#' ## Generate two metaheatmaps in the same plot. In this exampe one data frame
#' contains the whole dataset while in the other only reads of 28 nucleotides
#' are present
#' sample_name <- "Samp1"
#' metaheat_df <- list()
#' metaheat_df[["subsample_28nt"]] <- subset(reads_psite_info[[sample_name]], length == 28)
#' metaheat_df[["whole_sample"]] <- reads_psite_info[[sample_name]]
#' names_list <- list("Only_28" = c("subsample_28nt"), "All" = c("whole_sample"))
#' metaheat_comparison <- metaheatmap_psite(metaheat_df, mm81cdna, sample = names_list)
#' metaheat_comparison[["plot"]]
#' @import ggplot2
#' @export
metaheatmap_psite <- function(data, annotation, sample, scale_factors = NULL,
                              length_range = "all", transcripts = NULL,
                              utr5l = 25, cdsl = 50, utr3l = 25, log = F,
                              colour = "black", plot_title = NULL) {
  rownames(annotation) <- as.character(annotation$transcript)
  l.transcripts <- rownames(annotation)[which(annotation$l_utr5 >= utr5l &
                                                annotation$l_cds >= 2 * (cdsl + 1) &
                                                annotation$l_utr3 >= utr3l)]
  if (length(transcripts) == 0) {
    c.transcripts <- l.transcripts
    ntr <- length(c.transcripts)
  } else {
    c.transcripts <- intersect(l.transcripts, transcripts)
    ntr <- length(transcripts)
  }

  if(!identical(length_range, "all") & !inherits(length_range, "numeric") & !inherits(length_range, "integer")){
    warning("length_range is invalid. Set to default \"all\"\n")
    length_range="all"
  }

  if(inherits(length_range, "numeric") | inherits(length_range, "integer")){
    for(sampgroup in names(sample)){
      len_check_vec <- numeric()
      for(samp in sample[[sampgroup]]){
        len_check <- unique(data[[samp]]$length)
        if(sum(length_range %in% len_check) == 0) {
          warning(sprintf("data frame \"%s\" does not contain reads of the specified lengths", samp))
        }
        len_check_vec <- c(len_check_vec, len_check)
      }
      if(sum(length_range %in% len_check_vec) == 0) {
        sample[[sampgroup]] <- NULL
        warning(sprintf("none of the data frames in \"%s\" contain reads of the specified lengths: sample \"%s\" will not be considered ", sampgroup, sampgroup))
      }
    }
  }

  for(sampgroup in names(sample)){
    for (samp in sample[[sampgroup]]) {
      df <- data[[samp]][which(data[[samp]]$transcript %in% c.transcripts), ]
      if (identical(length_range, "all")) {
        start.sub <- df[which(df$psite_from_start %in% seq(-utr5l, cdsl)), ]
        stop.sub <- df[which(df$psite_from_stop %in% seq(-cdsl, utr3l)), ]
      } else {
        start.sub <- df[which(df$psite_from_start %in% seq(-utr5l, cdsl) & df$length%in%length_range), ]
        stop.sub <- df[which(df$psite_from_stop %in% seq(-cdsl, utr3l) & df$length%in%length_range), ]
      }
      start.tab <- as.data.frame(table(factor(start.sub$psite_from_start, levels = -utr5l:cdsl)))
      start.tab$reg <- "start"
      stop.tab <- as.data.frame(table(factor(stop.sub$psite_from_stop, levels = -cdsl:utr3l)))
      stop.tab$reg <- "stop"
      colnames(start.tab) <- colnames(stop.tab) <- c("distance", "reads", "reg")
      samp.tab <- rbind(start.tab, stop.tab)

      if (length(scale_factors) != 0 & length(scale_factors) == length(sample)) {
        samp.tab$reads <- samp.tab$reads * scale_factors[samp]
      }

      if (exists("temp.tab.heatmap")) {
        temp.tab.heatmap$reads <- temp.tab.heatmap$reads + samp.tab$reads
      } else {
        temp.tab.heatmap <- samp.tab
      }
      temp.tab.heatmap$sample <- sampgroup
    }

    if (exists("final.tab.heatmap")) {
      final.tab.heatmap <- rbind(final.tab.heatmap, temp.tab.heatmap)
    } else {
      final.tab.heatmap <- temp.tab.heatmap
    }
    rm(temp.tab.heatmap)
  }

  final.tab.heatmap$reg <- factor(final.tab.heatmap$reg, levels = c("start", "stop"), labels = c("Distance from start (nt)", "Distance from stop (nt)"))
  final.tab.heatmap$sample <- factor(final.tab.heatmap$sample, levels = rev(unique(final.tab.heatmap$sample)))

  linestart <- data.frame(reg = rep(c("Distance from start (nt)", "Distance from stop (nt)"), times = c(length(c(rev(seq(-3, -utr5l, -3)), seq(3, cdsl, 3))), length(c(rev(seq(-2, -cdsl, -3)), seq(1, utr3l, 3))))), line = c(rev(seq(-3, -utr5l, -3)), seq(3, cdsl, 3), rev(seq(-2, -cdsl, -3)), seq(1, utr3l, 3)))
  linered <- data.frame(reg = c("Distance from start (nt)", "Distance from stop (nt)"), line =c(0, 1))

  if(length(plot_title) != 0 && plot_title == "auto"){
    if (identical(length_range, "all")) {
      plot_title <- paste(ntr, " tr", sep = "")
    } else {
      if(min(length_range)==max(length_range)) {
        plot_title <- paste(ntr, " tr - Read length: ", min(length_range), " nts", sep = "")
      } else {
        if(identical(length_range, min(length_range):max(length_range)) | identical(length_range, seq(min(length_range),max(length_range),1))){
          plot_title <- paste(ntr, " tr - Read lengths: ", min(length_range), "-", max(length_range), " nts",sep = "")
        } else {
          plot_title <- paste(ntr, " tr - Read lengths: ", paste(length_range, collapse=","), " nts",sep = "")
        }
      }
    }
  }

  max <- max(final.tab.heatmap$reads)

  plot <- ggplot(final.tab.heatmap, aes(as.numeric(as.character(distance)), sample)) +
    geom_vline(data = linestart, aes(xintercept = line), linetype = 3, color = "gray60") +
    geom_vline(data = linered, aes(xintercept = line), linetype = 1, color = "red") +
    geom_tile(aes(fill = reads), height = 0.7) +
    labs(x = "", y = "", title = plot_title) +
    theme_bw(base_size = 20) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    facet_grid(. ~ reg, scales = "free", switch = "x") +
    theme(strip.background = element_blank())

  if (log == F) {
    plot <- plot +
      scale_fill_gradient("P-site signal\n", low = "white", high = colour, limits = c(0.1, max), breaks = c(0.1, max/2, max), labels = c("0", floor(max/2), floor(max)), na.value = "white")
  } else {
    plot <- plot +
      scale_fill_gradient("P-site signal\n", low = "white", high = colour, limits = c(0.1, max), breaks = c(0.1, 10^(log10(max)/2 - 0.5), floor(max)), labels = c("0", floor(10^(log10(max)/2 - 0.5)), floor(max)), trans = "log", na.value = "transparent")
  }

  output<-list()
  output[["plot"]] <- plot
  output[["df"]] <- final.tab.heatmap
  output[["transcripts"]] <- c.transcripts
  return(output)
}
