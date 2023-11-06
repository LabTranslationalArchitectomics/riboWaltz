#' Read length distributions.
#'
#' This function generates read length distributions, displayed as bar plots.
#' Multiple samples can be handled in several ways.
#'
#' @param data Either list of data tables or GRangesList object from
#'   \code{\link{bamtolist}}, \code{\link{bedtolist}},
#' \code{\link{length_filter}} or \code{\link{psite_info}}.
#' @param sample Either character string or character string vector specifying
#'   the name of the sample(s) of interest. A named list of one or more
#'   character strings and/or character string vectors can be provided. In this
#'   case i) each list element should include the name of the replicate(s)
#'   related to the sample of interest and ii) the name assigned to the elements
#'   of the list are displayed in the plot. Multiple replicates are handled
#'   according to \code{multisample}.
#' @param multisamples Either "separated" or "average". It specifies how to handle
#'   multiple samples and replicates. If "saparated", one bar plot for each
#'   sample included in \code{sample} is returned as an independent ggplot
#'   object. If \code{sample} is a list, it is unlisted, coerced to character
#'   string and handled accordingly. If "average" i) one barplot is returned if
#'   \code{sample} is a character string vector or ii) one bar plot is built for
#'   each element of \code{sample} when it is a list. If "average", the bar plot
#'   displays for each length the mean signal and the corresponding standard
#'   error computed among the replicates. In this case a single ggplot object is
#'   returned, where multiple bar plots are organized and displayed according to
#'   \code{plot_style}. Default is "separated".
#' @param plot_style Either "split", "dodged" or "mirrored". It specifies how to
#'   organize and display multiple bar plots. If "split", the bar plots are
#'   placed one next to the other, in independent boxes. If "dodged", all bar
#'   plots are included in one box, with each bar side by side with the others.
#'   If "mirrored", \code{sample} must be a list of exactly two elements and the
#'   two bar plots are mirrored along the x axis. Default is "split".
#' @param transcripts Character string vector listing the name of transcripts to
#'   be included in the analysis. Default is NULL i.e. all transcripts are used.
#' @param cl Integer value in [1,100] specifying a confidence level for
#'   restricting the plot to a sub-range of read lengths. The new range is
#'   associated to the most abundant populations of reads accounting for the cl%
#'   of the sample. If multiple sample names are provided one range of read
#'   lengths is computed, such that at least the cl% of all sample are
#'   represented. Default is 100.
#' @param colour Character string or character string vector specifying the
#'   colour of the bar plot(s). If \code{sample} is a list of multiple
#'   elements and \code{multisamples} is set to "average", a colour for
#'   each element of the list is required. If this parameter is not specified
#'   the R default palette is employed. Default is NULL.
#' @return List containing one or more ggplot object(s) and the data table with
#'   the associated data ("dt").
#' @examples
#' data(reads_list)
#'
#' ## Generate the length distribution for all read lengths:
#' lendist_whole <- rlength_distr(reads_list, sample = "Samp1", cl = 100)
#' lendist_whole[["plot_Samp1"]]
#'
#' ## Generate the length distribution for a sub-range of read lengths:
#' lendist_sub95 <- rlength_distr(reads_list, sample = "Samp1", cl = 95)
#' lendist_sub95[["plot_Samp1"]]
#' @import data.table
#' @import ggplot2
#' @export
rlength_distr <- function(data, sample, multisamples = "separated",
                          plot_style = "split", transcripts = NULL, cl = 100,
                          colour = NULL) {
  
  if(class(data[[1]])[1] == "GRanges"){
    data_tmp <- list()
    for(i in names(data)){
      data_tmp[[i]] <- as.data.table(data[[i]])[, c("width", "strand") := NULL
                                                ][, seqnames := as.character(seqnames)]
      setnames(data_tmp[[i]], c("seqnames", "start", "end"), c("transcript", "end5", "end3"))
    }
    data <- data_tmp
  }
  
  check_sample <- setdiff(unlist(sample), names(data))
  if(length(check_sample) != 0){
    cat("\n")
    stop(sprintf("incorrect sample name(s): \"%s\" not found\n\n",
                 paste(check_sample, collapse = ", ")))
  }
  
  if(length(sample) == 0){
    cat("\n")
    stop("at least one sample name must be spcified\n\n")
  }
  
  if(multisamples == "separated" & is.list(sample)) {
    cat("\n")
    warning("parameter multisamples is set to \"separated\" but a list of samples is provided:\nparameter sample will be unlisted and coerced to character string\n", call. = FALSE)
    sample <- as.character(unlist(sample))
  }
  
  if(is.list(sample) & length(sample) > 2 & plot_style == "mirrored") {
    cat("\n")
    warning("parameter sample is a list of dimension > 2.\nparameter plot_style set to default \"split\"\n", call. = FALSE)
    plot_style <- "split"
  }
  
  if(!(multisamples %in% c("average", "separated"))){
    cat("\n")
    warning("parameter multisamples must be either \"separated\" or \"average\"\nset to default \"separated\"\n", call. = FALSE)
    multisamples <- "separated"
  }
  
  if(multisamples == "average" & 
     !(plot_style %in% c("mirrored", "dodged", "split"))){
    cat("\n")
    warning("parameter plot_style must be either \"split\", \"mirrored\" or \"dodged\"\nset to default \"split\"\n", call. = FALSE)
    plot_style <- "split"
  }
  
  if(multisamples == "average" & ((is.list(sample) & length(sample) == 1) |
                                  (is.character(sample) & length(sample) > 1))){
    plot_style <- "dodged"
    if(is.character(sample) & length(sample) > 1){
      warning("consider to use a list of one named element to provide a title to the plot\n", call. = FALSE)
    }
  } 
  
  if(multisamples == "average" & ((is.character(sample) & length(sample) == 1)  |
                                  (is.list(sample) & length(as.character(unlist(sample))) == 1))){
    cat("\n")
    warning("parameter multisamples is set to \"average\" but only one sample is provided\nset to default \"separated\"", call. = FALSE)
    multisamples <- "separated"
  }
  
  if(is.character(sample) | (is.list(sample) & length(sample) == 1)){
    if(is.list(sample) & length(sample) == 1){
      samp_name <- names(sample)
      sample <- as.character(unlist(sample))
      sample_l <- list()
      sample_l[[samp_name]] <- sample
    } else {
      sample_l <- list("sample_name" = sample)
    }
  } else {
    sample_l <- sample
  }
  
  if(!is.null(colour)){
    if(length(sample_l) == 1){
      colour <- colour[1]
    } else {
      if(length(sample_l) == length(colour)) {
        names(colour) <- names(sample_l)
      } else {
        if(length(sample_l) > length(colour)){
          warning(sprintf("at least %s colours must be specified\nsystem default colours will be used",
                          length(sample_l)), call. = FALSE)
          colour = NULL
        }
      }
    }
  }
  
  #define length range
  for(samp in as.character(unlist(sample_l))){
    if(length(transcripts) == 0) {
      dt <- data[[samp]]
    } else {
      dt <- data[[samp]][transcript %in% transcripts]
    }
    
    if(!exists("length_range")){
      length_range <- seq(quantile(dt$length, (1 - cl/100)/2),
                          quantile(dt$length, 1 - (1 - cl/100)/2))
    } else {
      xmin <- min(min(length_range), quantile(dt$length, (1 - cl/100)/2))
      xmax <- max(max(length_range), quantile(dt$length, 1 - (1 - cl/100)/2))
      length_range <- seq(xmin, xmax)
    }
  }
  
  xmin = min(length_range)
  xmax = max(length_range)
  
  # compute count and percentage of reads of defined lengths
  for(samp in as.character(unlist(sample_l))){
    if(length(transcripts) == 0) {
      dt <- data[[samp]]
    } else {
      dt <- data[[samp]][transcript %in% transcripts]
    }
    
    setkey(dt, length)
    dist_dt <- dt[CJ(length_range), .(count = .N), by = .EACHI
                  ][, percentage := (count / sum(count)) * 100]
    
    if(!exists("final_dt")){
      final_dt <- dist_dt
      setnames(final_dt, c("count", "percentage"), c(paste0(samp, "_count"), paste0(samp, "_percentage")))
    } else {
      final_dt <- final_dt[, (paste0(samp, "_count")) := dist_dt[, count]
                           ][, (paste0(samp, "_percentage")) := dist_dt[, percentage]]
    }
  }
  
  col_plot_sel <- list()
  col_plot_se_sel <- list()
  
  if(length(sample_l) == 1){
    samp_group <- names(sample_l)
    samp <- sample_l[[samp_group]]
    if(length(samp) != 1 & multisamples == "average"){
      samp <- paste0(samp, "_percentage")
      final_dt[, (paste0("mean_", samp_group, "_percentage"))
               := apply(.SD, 1, mean), .SDcols = samp
               ][, (paste0("se_", samp_group, "_percentage"))
                 := apply(.SD, 1, sd)/sqrt(length(samp)), .SDcols = samp]
      col_plot_sel[[samp_group]] <- paste0("mean_", samp_group, "_percentage")
      col_plot_se_sel[[samp_group]] <- paste0("se_", samp_group, "_percentage")
    } else {
      col_plot_sel[[samp_group]] <- samp
    }
  } else {
    for(samp_group in names(sample_l)){
      samp <- paste0(sample_l[[samp_group]], "_percentage")
      if(length(samp) != 1 & multisamples == "average"){
        final_dt[, (paste0("mean_", samp_group, "_percentage")) :=
                   apply(.SD, 1, mean), .SDcols = samp
                 ][, (paste0("se_", samp_group, "_percentage")) :=
                     apply(.SD, 1, sd)/sqrt(length(samp)), .SDcols = samp]
        col_plot_sel[[samp_group]] <- paste0("mean_", samp_group, "_percentage")
        col_plot_se_sel[[samp_group]] <- paste0("se_", samp_group, "_percentage")
      } else {
        col_plot_sel[[samp_group]] <- samp
      }
    }
  }
  
  output <- list()
  output[["dt"]] <- final_dt
  
  if(multisamples != "separated"){
    plot_dt <- copy(final_dt)
    if(plot_style == "mirrored"){
      plot_dt[, (col_plot_sel[[2]]) := - get(col_plot_sel[[2]])]
    }
    col_sel <- c("length", as.character(unlist(col_plot_sel)))
    melt_plot_dt <- melt.data.table(plot_dt[, ..col_sel],
                                    id.vars = c("length"),
                                    variable.name = "sample",
                                    value.name = "mean"
    )[, sample := factor(sample,
                         levels = as.character(unlist(col_plot_sel)),
                         labels = names(col_plot_sel))]
    
    if(length(col_plot_se_sel) != 0){
      col_sel = c("length", as.character(unlist(col_plot_se_sel)))
      melt_plot_se <- melt.data.table(plot_dt[, ..col_sel],
                                      id.vars = c("length"),
                                      variable.name = "sample",
                                      value.name = "se"
      )[, sample := factor(sample,
                           levels = as.character(unlist(col_plot_se_sel)),
                           labels = names(col_plot_se_sel))]
      
      melt_plot_dt <- merge.data.table(melt_plot_dt, melt_plot_se,
                                       by = c("length", "sample"),
                                       all.x = TRUE, sort = FALSE)
    }
    
    plot <- ggplot(melt_plot_dt, aes(as.numeric(length), mean, fill = sample))
    
    if(plot_style != "mirrored") {
      plot <- plot + geom_bar(stat = "identity", position = position_dodge(0.9))
      if(length(col_plot_se_sel) != 0){
        plot <- plot + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = sample),
                                     width = 0.40, size = 1.1, na.rm = T, position = position_dodge(0.9))
      }
    } else {
      plot <- plot + geom_bar(stat = "identity")
      if(length(col_plot_se_sel) != 0){
        plot <- plot + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = sample),
                                     width = 0.35, size = 1.1, na.rm = T)
      }
    }
    
    plot <- plot + labs(x = "Read length", y = "Count (%)") +
      theme_bw(base_size = 18) +
      scale_x_continuous(limits = c(xmin - 0.5, xmax + 0.5),
                         breaks = seq(xmin + ((xmin) %% 2), xmax,
                                      by = max(c(1, floor((xmax - xmin)/7))))) + 
      theme(panel.grid.minor.x = element_blank())
    
    if(!is.null(colour)){
      plot <- plot + scale_fill_manual(labels = names(sample_l), values = as.character(colour)) +
        scale_color_manual(labels = names(sample_l), values = as.character(colour))
    }
    
    if(length(col_plot_sel) > 1 & plot_style != "split"){
      plot <- plot + theme(legend.position = c(0.98,1), legend.justification = c(1, 1),
                           legend.title = element_blank(), legend.background = element_blank())
    } else {
      plot <- plot + theme(legend.position = "none")
    }
    
    if(plot_style == "mirrored"){
      plot <- plot + scale_y_continuous(labels = abs)
    }
    
    if(plot_style == "split"){
      plot <- plot + facet_wrap(sample ~ ., ncol = 2) +
        theme(strip.background = element_blank())
    }
    
    if(multisamples == "average" & is.list(sample_l) &
       length(sample_l) == 1 && names(sample_l) != "sample_name"){
      plot <- plot + ggtitle(names(sample_l)) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    
    output[["plot"]] <- plot
    
  } else {
    if(!is.null(colour)){
      colour = colour[1]
    } else {
      colour = "gray40"
    }
    
    for(col_plot in paste0(as.character(unlist(col_plot_sel)), "_percentage")){
      plot <- ggplot() +
        geom_bar(data = final_dt, aes_string("length", col_plot),
                 stat = "identity", fill = colour) +
        labs(title = gsub("_percentage", "", col_plot), x = "Read length", y = "Count (%)") +
        theme_bw(base_size = 18) +
        scale_x_continuous(limits = c(xmin - 0.5, xmax + 0.5),
                           breaks = seq(xmin + ((xmin) %% 2), xmax,
                                        by = max(c(1, floor((xmax - xmin)/7))))) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        theme(panel.grid.minor.x = element_blank())
      
      output[[gsub("_percentage", "", paste0("plot_", col_plot))]] <- plot
    }
  }
  return(output)
}
