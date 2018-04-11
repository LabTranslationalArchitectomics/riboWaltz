#' Compute the percentage of P-sites per frame.
#'
#' For one or several samples this function computes the percentage of P-sites
#' falling on the three reading frames of the transcripts and generates a
#' barplot of the resulting values. This analysis is performed for the annotated
#' 5' UTR, coding sequence and 3' UTR, separately. It is possible to compute the
#' percentage of P-sites per frame using all the read lengths or to restrict the
#' analysis to a sub-range of read lengths.
#'
#' @param data A list of data tables from \code{\link{psite_info}}.
#' @param sample A character string vector specifying the name of the sample(s)
#'   of interest. By default this argument is NULL, meaning that all the samples
#'   in \code{data} are included in the analysis.
#' @param region Either "all" or a character string among "5utr", "cds", "3utr"
#'   specifying the regions of the transcript (5' UTR, CDS or 3' UTR,
#'   respectively) that must be included in the analysis. Default is "all",
#'   meaning that the all the regions are considered.
#' @param length_range Either "all", an integer or an integer vector. Default is
#'   "all", meaning that all the read lengths are included in the analysis.
#'   Otherwise, only the read lengths matching the specified value(s) are kept.
#' @param plot_title Any character string specifying the title of the plot. If
#'   "auto", the title of the plot reports the region specified by \code{region}
#'   (if any) and the length(s) of the reads used for generating the barplot.
#'   Default is NULL, meaning that no title will be added to the plot.
#' @return A list containing a ggplot2 object and a data table with the
#'   associated data.
#' @examples
#' data(reads_psite_list)
#'
#' ## Generate the barplot for all the read lengths
#' frame_whole <- frame_psite(reads_psite_list, sample = "Samp1")
#' frame_whole[["plot"]]
#'
#' ## Generate the barplot restricting the analysis to the coding sequence and
#' to the reads of 28 nucleotides
#' frame_sub <- frame_psite(reads_psite_list, sample = "Samp1", region = "cds",
#' length_range = 28)
#' frame_sub[["plot"]]
#' @import data.table
#' @import ggplot2
#' @export
frame_psite <- function(data, sample = NULL, region = "all", length_range = "all",
                        plot_title = NULL){
  if(length(sample) == 0) {
    sample <- names(data)
  }
  
  if(!identical(length_range, "all") & !inherits(length_range, "numeric") & !inherits(length_range, "integer")){
    cat("\n")
    warning("class of length_range is neither numeric nor integer. Set to default \"all\"\n")
    length_range = "all"
  }
  
  if(!identical(length_range, "all")){
    for(samp in sample){
      len_check <- unique(data[[samp]]$length)
      if(sum(length_range %in% len_check) == 0) {
        cat("\n")
        warning(sprintf("\"%s\" doesn't contain any reads of the specified lengths: sample removed\n", samp))
        sample <- sample[sample != samp]
      }
    }
  }
  
  if(length(sample) == 0){
    cat("\n")
    stop("none of the data tables in sample contains any reads of the specified lengths\n\n")
  }
  
  if(!region %in% c("all", "cds", "5utr", "3utr")){
    cat("\n")
    warning("region is invalid. Set to default \"all\"\n")
    region = "all"
  }
  
  length_temp <- vector()
  
  for (samp in sample) {
    
    if(identical(length_range, "all")){
      dt <- data[[samp]]
    } else {
      dt <- data[[samp]][length %in% length_range]
    }
    
    if (region == "all") {
      frame_dt <- dt[start_pos != 0 & stop_pos !=0
               ][, frame := psite_from_start %% 3
                 ][, list(count = .N), by = list(region = psite_region, frame)
                   ][, percentage := (count / sum(count)) * 100, by = region
                     ][is.na(percentage), percentage := 0
                       ][, region := factor(region, levels = c("5utr", "cds", "3utr"), labels = c("5' UTR", "CDS", "3' UTR"))]
    } else {
      frame_dt <- dt[psite_region == region
                     ][start_pos != 0 & stop_pos !=0
                       ][, frame := psite_from_start %% 3
                         ][, list(count = .N), by = frame
                           ][, percentage := (count / sum(count)) * 100
                             ][is.na(percentage), percentage := 0]
    }
    
    frame_dt[, sample := samp]
    
    if (exists("final_frame_dt")) {
      final_frame_dt <- rbind(final_frame_dt, frame_dt)
    } else {
      final_frame_dt <- frame_dt
    }
    
    length_temp <- unique(c(length_temp, data[[samp]]$length))
  }
  
  if(!identical(length_range, "all")){
    length_range <- sort(intersect(length_range, length_temp))
  } else {
    length_range <- sort(length_temp)
  }
  
  plot <- ggplot(final_frame_dt, aes(x = frame, y = percentage)) +
    geom_bar(stat = "identity") +
    theme_bw(base_size = 20) +
    labs(x = "Frame", y = "P-site signal (%)")
    
  if(region == "all") {
    plot <- plot + facet_grid(sample ~ region)
    final_frame_dt <- final_frame_dt[order(region, frame, sample)]
  } else {
    plot <- plot + facet_wrap( ~ sample, ncol = 3)
    final_frame_dt <- final_frame_dt[order(frame, sample)]
  }
  
  if(identical(plot_title, "auto")) {
    
    if(region == "all") {
      plottitle_region <- NULL
    } else {
      if(region == "5utr") { plottitle_region <- "Region: 5' UTR. " }
      if(region == "cds") { plottitle_region <- "Region: CDS. " }
      if(region == "3utr") { plottitle_region <- "Region: 3' UTR. " }
    }
    
    minlr <- min(length_range)
    maxlr <- max(length_range)
    
    if(minlr == maxlr) {
      plottitle_range <- paste0("Read length: ", minlr, " nts")
    } else {
      if(identical(length_range, minlr:maxlr) | identical(length_range, seq(minlr, maxlr, 1))){
        plottitle_range <- paste0("Read lengths: ", minlr, "-", maxlr, " nts")
      } else {
        nextl <- sort(length_range[c(which(diff(length_range) != 1), which(diff(length_range) != 1) + 1)])
        sep <- ifelse(nextl %in% length_range[which(diff(length_range) != 1)], ", ", "-")[-length(nextl)]
        if(1 %in% which(diff(length_range) == 1)){
          nextl <- c( length_range[1], nextl)
          sep <- c("-", sep)
        }
        if((length(length_range) - 1) %in% which(diff(length_range) == 1)){
          nextl <- c(nextl, length_range[length(length_range)])
          sep <- c(sep, "-")
        }
        sep <- c(sep, "")
        plottitle_range <- paste0("Read lengths: ", paste0(nextl, sep, collapse = ""), " nts")
      }
    }
    plottitle <- paste0(plottitle_region, plottitle_range)
    plot <- plot +
      labs(title = plottitle) +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    if(length(plot_title) != 0){
      plot <- plot +
        labs(title = plot_title) + 
        theme(plot.title = element_text(hjust = 0.5))
    }
  }
  
  ret_list<-list()
  ret_list[["dt"]] <- final_frame_dt
  ret_list[["plot"]] <- plot
  return(ret_list)
}

#' Compute the number of P-sites per frame stratified by read length.
#'
#' Similar to \code{\link{frame_psite}} but the results are stratified by the
#' length of the reads.
#'
#' @param data A list of data tables from \code{\link{psite_info}}.
#' @param sample A character string vector specifying the name of the sample(s)
#'   of interest. By default this argument is NULL, meaning that all the samples
#'   in \code{data} are included in the analysis.
#' @param region Either "all" or a character string among "5utr", "cds", "3utr"
#'   specifying the regions of the transcript (5' UTR, CDS or 3' UTR,
#'   respectively) that must be included in the analysis. Default is "all",
#'   meaning that the all the regions are considered.
#' @param cl An integer value in \emph{[1,100]} specifying the confidence level
#'   for restricting the analysis to a sub-range of read lengths. Default is 95.
#'   This parameter has no effect if \code{length_range} is specified.
#' @param length_range Either "all", an integer or an integer vector. Default is
#'   "all", meaning that all the read lengths are included in the analysis.
#'   Otherwise, only the read lengths matching the specified value(s) are kept.
#'   If specified, this parameter prevails over \code{cl}.
#' @param plot_title Any character string specifying the title of the plot. When
#'   "auto", the title of the plot reports the region specified by \code{region}
#'   (if any). Default is NULL, meaning that no title will be added to the plot.
#' @return A list containing a ggplot2 object and a data table with the
#'   associated data.
#' @examples
#' data(reads_psite_list)
#'
#' ## Generate the heatmap for all the read lengths
#' frame_len_whole <- frame_psite_length(reads_psite_list, sample = "Samp1")
#' frame_len_whole[["plot"]]
#'
#' ## Generate the heatmap for a sub-range of read lengths (the middle 90%) and
#' restricting the analysis to the coding sequence
#' frame_len_sub <- frame_psite_length(reads_psite_list, sample = "Samp1",
#' region = "cds", cl = 90)
#' frame_len_sub[["plot"]]
#' @import data.table
#' @import ggplot2
#' @export
frame_psite_length <- function(data, sample = NULL, region = "all", cl = 95,
                                  length_range = "all", plot_title = NULL){
  
  if(length(sample) == 0) {
    sample <- names(data)
  }
  
  if(!identical(length_range, "all") & !inherits(length_range, "numeric") & !inherits(length_range, "integer")){
    cat("\n")
    warning("class of length_range is neither numeric nor integer. Set to default \"all\"\n")
    length_range = "all"
  }
  
  if(!identical(length_range, "all")){
    for(samp in sample){
      len_check <- unique(data[[samp]]$length)
      if(sum(length_range %in% len_check) == 0) {
        cat("\n")
        warning(sprintf("\"%s\" doesn't contain any reads of the specified lengths: sample removed\n", samp))
        sample <- sample[sample != samp]
      }
    }
  }
  
  if(length(sample) == 0){
    cat("\n")
    stop("none of the data tables in sample contains any reads of the specified lengths\n\n")
  }
  
  if(!identical(length_range, "all")){
    minl <- min(length_range)
    maxl <- max(length_range)
  }

  if(!region%in%c("all", "cds", "5utr", "3utr")){
    cat("\n")
    warning("region is invalid. Set to default \"all\"\n")
    region = "all"
  }
  
  for (samp in sample) {
    
    if (region == "all") {
      
      dt <- data[[samp]][start_pos != 0 & stop_pos != 0
                         ][, frame := psite_from_start %% 3]
      
      if(identical(length_range, "all")){
        minl <- quantile(dt$length, (1 - cl/100) / 2)
        maxl <- quantile(dt$length, 1 - (1 - cl/100) / 2)
      }
      
      dt[, psite_region := factor(psite_region, levels = c("5utr", "cds", "3utr"))
         ][, frame := factor(frame, levels = c(0, 1, 2))
           ][, length := factor(length, levels = unique(length))]
      
      setkey(dt, length, psite_region, frame)
      frame_dt <- dt[CJ(levels(length), levels(psite_region), levels(frame)), list(count = .N), by = .EACHI
                     ][as.numeric(as.character(length)) %in% minl:maxl
                       ][, percentage := (count / sum(count)) * 100, by = list(length, psite_region)
                         ][is.na(percentage), percentage := 0
                           ][, psite_region := factor(psite_region, levels = c("5utr", "cds", "3utr"), labels = c("5' UTR", "CDS", "3' UTR"))]
      setnames(frame_dt, "psite_region", "region")

    } else {
      
      dt <- data[[samp]][psite_region == region
                         ][, frame := psite_from_start %% 3
                           ]
      
      if(length(length_range) == 0){
        minl <- quantile(dt$length, (1 - cl/100) / 2)
        maxl <- quantile(dt$length, 1 - (1 - cl/100) / 2)
      }
      
      dt[, frame := factor(frame, levels = c(0, 1, 2))
         ][, length := factor(length, levels = unique(length))]
      
      setkey(dt, length, frame)
      frame_dt <- dt[CJ(levels(length), levels(frame)), list(count = .N), by = .EACHI
                     ][as.numeric(as.character(length)) %in% minl:maxl
                       ][, percentage := (count / sum(count)) * 100, by = length
                         ][is.na(percentage), percentage := 0]
    }
    
    frame_dt$sample <- samp
    
    if (exists("final_frame_dt")) {
      final_frame_dt <- rbind(final_frame_dt, frame_dt)
    } else {
      final_frame_dt <- frame_dt
    }
  }
  
  minfp <- min(final_frame_dt$percentage)
  maxfo <- max(final_frame_dt$percentage)
  
  plot <- ggplot(final_frame_dt, aes(frame, as.numeric(as.character(length)))) +
    geom_tile(aes(fill = percentage)) +
    scale_fill_gradient("P-site signal (%)  ", low = "white", high = "#104ec1",
                        breaks = c(minfp, (maxfo - minfp) / 2, maxfo),
                        labels = c(round(minfp), round((maxfo - minfp) / 2), round(maxfo))) +
    labs(x = "Frame", y = "Read length") +
    theme_bw(base_size = 20) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(legend.position = "top") +
    scale_y_continuous(limits = c(minl - 0.5, maxl + 0.5), breaks = seq(minl + ((minl) %% 2), maxl, by = max(2, floor((maxl - minl) / 7))))
  
  if (region == "all") {
    plot <- plot + facet_grid(sample ~ region)
    final_frame_dt <- final_frame_dt[order(length, region, frame, sample)]
  } else {
    plot <- plot + facet_wrap( ~ sample, ncol = 3)
    final_frame_dt <- final_frame_dt[order(length, frame, sample)]
  }
  
  if(identical(plot_title, "auto") & !identical(region, "all")){
    if(region == "5utr") { plottitle <- "Region: 5' UTR" }
    if(region == "cds") { plottitle <- "Region: CDS" }
    if(region == "3utr") { plottitle <- "Region: 3' UTR" }
    plot <- plot +
      labs(title = plottitle) +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    if(length(plot_title) != 0 & !identical(plot_title, "auto")){
    plot <- plot +
      labs(title = plot_title) + 
      theme(plot.title = element_text(hjust = 0.5))
    }
  }
  
  ret_list<-list()
  ret_list[["dt"]] <- final_frame_dt
  ret_list[["plot"]] <- plot
  return(ret_list)
}
