#' Percentage of P-sites per reading frame.
#'
#' This function computes the percentage of P-sites falling in the three
#' possible translation reading frames and generates a bar plot of the resulting
#' values. It only handles annotated 5' UTRs, coding sequences and 3' UTRs,
#' separately.
#'
#' @param data List of data tables from \code{\link{psite_info}}.
#' @param sample Character string vector specifying the name of the sample(s) of
#'   interest. Default is NULL i.e. all samples in \code{data} are processed.
#' @param transcripts Character string vector listing the name of transcripts to
#'   be included in the analysis. Default is NULL i.e. all transcripts are used.
#' @param region Character string specifying the region(s) of the transcripts to
#'   be analysed. It can be either "5utr", "cds", "3utr" for 5' UTRs, CDSs and
#'   3' UTRs, respectively. Default is "all" i.e. all regions are considered.
#'   According to this parameter the bar plots are differently arranged to
#'   optimise the organization and the visualization of the data.
#' @param length_range Integer or an integer vector specyfying the read
#'   length(s) to be included in the analysis. Default is "all" i.e. all read
#'   lengths are used.
#' @param plot_title Character string specifying the title of the plot. If
#'   "auto", the title of the plot reports the region specified by \code{region}
#'   (if any) and the considered read length(s). Default is NULL i.e. no title
#'   is plotted.
#'   
#' @return A list containing a ggplot2 object ("plot") and the data table with
#'   the associated data ("dt").
#' @examples
#' data(reads_psite_list)
#'
#' ## Generate the bar plot for all read lengths:
#' frame_whole <- frame_psite(reads_psite_list, sample = "Samp1")
#'
#' ## Generate the bar plot restricting the analysis to coding sequences and
#' ## reads of 28 nucleotides:
#' frame_sub <- frame_psite(reads_psite_list, sample = "Samp1", region = "cds",
#' length_range = 28)
#' @import data.table
#' @import ggplot2
#' @export
frame_psite <- function(data, sample = NULL, transcripts = NULL, region = "all",
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
      
      if(length(transcripts) == 0) {
        dt <- data[[samp]]
      } else {
        dt <- data[[samp]][transcript %in% transcripts]
      }
      
      len_check <- unique(dt$length)
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
    
    if (length(transcripts) != 0) {
      dt <- dt[transcript %in% transcripts]
    }
    
    if (region == "all") {
      frame_dt <- dt[cds_start != 0 & cds_stop !=0
                     ][, frame := psite_from_start %% 3]
      
      setkey(frame_dt, psite_region, frame)
      frame_dt <- frame_dt[CJ(psite_region = c("5utr", "cds", "3utr"), frame = c(0,1,2)),
                           list(count = .N), by = .EACHI
                           ][, percentage := (count / sum(count)) * 100, by = psite_region
                             ][is.na(percentage), percentage := 0
                               ][, psite_region := factor(psite_region,
                                                          levels = c("5utr", "cds", "3utr"),
                                                          labels = c("5' UTR", "CDS", "3' UTR"))]
      setnames(frame_dt, "psite_region", "region")
    } else {
      frame_dt <- dt[psite_region == region
                     ][cds_start != 0 & cds_stop !=0
                       ][, frame := psite_from_start %% 3]
      
      setkey(frame_dt, frame)
      frame_dt <- frame_dt[CJ(frame = c(0,1,2)), list(count = .N), by = .EACHI
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
  final_frame_dt[, sample := factor(sample, levels = unique(sample))]
  
  if(!identical(length_range, "all")){
    length_range <- sort(intersect(length_range, length_temp))
  } else {
    length_range <- sort(length_temp)
  }
  
  plot <- ggplot(final_frame_dt, aes(x = frame, y = percentage)) +
    geom_bar(stat = "identity") +
    theme_bw(base_size = 22.5) +
    labs(x = "Frame", y = "P-site signal (%)") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  
  if(region == "all") {
    plot <- plot + facet_grid(sample ~ region)
    final_frame_dt <- final_frame_dt[order(sample, region, frame)]
  } else {
    plot <- plot + facet_wrap( ~ sample, ncol = ceiling(sqrt(length(sample))))
    final_frame_dt <- final_frame_dt[order(sample, frame)]
  }
  
  plot <- plot + theme(strip.background = element_blank())
  
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

#' Percentage of P-sites per reading frame stratified by read length.
#'
#' Similar to \code{\link{frame_psite}}, but the results are stratified by read
#' lengths and plotted as heatmaps.
#'
#' @param data List of data tables from \code{\link{psite_info}}.
#' @param sample Character string vector specifying the name of the sample(s) of
#'   interest. Default is NULL i.e. all samples in \code{data} are processed.
#' @param transcripts Character string vector listing the name of transcripts to
#'   be included in the analysis. Default is NULL i.e. all transcripts are used.
#' @param region Character string specifying the region(s) of the transcripts to
#'   be analysed. It can be either "5utr", "cds", "3utr" for 5' UTRs, CDSs and
#'   3' UTRs, respectively. Default is "all" i.e. all regions are considered.
#'   According to this parameter the heatmaps are differently arranged to
#'   optimise the organization and the visualization of the data.
#' @param cl Integer value in [1,100] specifying a confidence level for
#'   restricting the analysis to a sub-range of read lengths i.e. to the cl% of
#'   read lengths associated to the highest signals. Default is 100.
#'   This parameter has no effect if \code{length_range} is specified.
#' @param length_range Integer or an integer vector specyfying the read
#'   length(s) to be included in the analysis. Default is "all" i.e. all read
#'   lengths are used. If specified, this parameter prevails over \code{cl}.
#' @param plot_title Character string specifying the title of the plot. If
#'   "auto", the title of the plot reports the region specified by \code{region}
#'   (if any) and the considered read length(s). Default is NULL i.e. no title
#'   is displayed.
#' @param colour Character string specifying the colour of the plot. The colour
#'   scheme is as follow: tiles corresponding to the lowest signal are always
#'   white, tiles corresponding to the highest signal are of the specified
#'   colour and the progression between these two colours follows a linear
#'   gradient. Default is dark blue.
#' @return A list containing a ggplot2 object ("plot") and the data table with
#'   the associated data ("dt").
#' @examples
#' data(reads_psite_list)
#'
#' ## Generate the heatmap for all read lengths:
#' frame_len_whole <- frame_psite_length(reads_psite_list, sample = "Samp1")
#'
#' ## Generate the heatmap restricting the analysis to coding sequences and a 
#' ## sub-range of read lengths:
#' frame_len_sub <- frame_psite_length(reads_psite_list, sample = "Samp1",
#' region = "cds", cl = 90)
#' @import data.table
#' @import ggplot2
#' @export
frame_psite_length <- function(data, sample = NULL, transcripts = NULL,
                               region = "all", cl = 100, length_range = "all",
                               plot_title = NULL, colour = "#061b63"){
  
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
      
      if(length(transcripts) == 0) {
        dt <- data[[samp]]
      } else {
        dt <- data[[samp]][transcript %in% transcripts]
      }
      
      len_check <- unique(dt$length)
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
    
    if(length(transcripts) == 0) {
      dt <- data[[samp]]
    } else {
      dt <- data[[samp]][transcript %in% transcripts]
    }
    
    if (region == "all") {
      
      dt <- dt[cds_start != 0 & cds_stop != 0
               ][, frame := psite_from_start %% 3]
      
      if(identical(length_range, "all")){
        minl <- quantile(dt$length, (1 - cl/100) / 2)
        maxl <- quantile(dt$length, 1 - (1 - cl/100) / 2)
        length_range <- minl:maxl
      }
      
      dt[, psite_region := factor(psite_region, levels = c("5utr", "cds", "3utr"))
         ][, frame := factor(frame, levels = c(0, 1, 2))
           ][, length := factor(length, levels = unique(length))]
      
      setkey(dt, length, psite_region, frame)
      frame_dt <- dt[CJ(levels(length), levels(psite_region), levels(frame)), list(count = .N), by = .EACHI
                     ][as.numeric(as.character(length)) %in% length_range
                       ][, percentage := (count / sum(count)) * 100, by = list(length, psite_region)
                         ][is.na(percentage), percentage := 0
                           ][, psite_region := factor(psite_region, levels = c("5utr", "cds", "3utr"), labels = c("5' UTR", "CDS", "3' UTR"))]
      setnames(frame_dt, "psite_region", "region")

    } else {
      dt <- dt[psite_region == region
               ][, frame := psite_from_start %% 3]
      
      if(identical(length_range, "all")){
        minl <- quantile(dt$length, (1 - cl/100) / 2)
        maxl <- quantile(dt$length, 1 - (1 - cl/100) / 2)
        length_range <- minl:maxl
      }
      
      dt[, frame := factor(frame, levels = c(0, 1, 2))
         ][, length := factor(length, levels = unique(length))]
      
      setkey(dt, length, frame)
      frame_dt <- dt[CJ(levels(length), levels(frame)), list(count = .N), by = .EACHI
                     ][as.numeric(as.character(length)) %in% length_range
                       ][, percentage := (count / sum(count)) * 100, by = length
                         ][is.na(percentage), percentage := 0]
    }
    
    frame_dt[, sample := samp]
    
    if (exists("final_frame_dt")) {
      final_frame_dt <- rbind(final_frame_dt, frame_dt)
    } else {
      final_frame_dt <- frame_dt
    }
  }
  
  final_frame_dt[, sample := factor(sample, levels = unique(sample))]
  
  mins <- min(final_frame_dt$percentage)
  maxs <- max(final_frame_dt$percentage)
  
  plot <- ggplot(final_frame_dt, aes(frame, as.numeric(as.character(length)))) +
    geom_tile(aes(fill = percentage)) +
    scale_fill_gradient("P-site signal (%)  ", low = "white", high = colour,
                        breaks = c(mins, mins/2 + maxs/2, maxs),
                        labels = c(round(mins), round(mins/2 + maxs/2), round(maxs))) +
    labs(x = "Frame", y = "Read length") +
    theme_bw(base_size = 22.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = "top", legend.margin=margin(0,0,0,0), legend.box.margin=margin(5,0,-5,0)) +
    scale_y_continuous(limits = c(minl - 0.5, maxl + 0.5), breaks = seq(minl + ((minl) %% 2), maxl, by = max(2, floor((maxl - minl) / 7))))
  
  if (region == "all") {
    plot <- plot + facet_grid(sample ~ region)
    final_frame_dt <- final_frame_dt[order(sample, region, length, frame)]
  } else {
    plot <- plot + facet_wrap( ~ sample, ncol = ceiling(sqrt(length(sample))))
    final_frame_dt <- final_frame_dt[order(sample, length, frame)]
  }
  
  plot <- plot + theme(strip.background = element_blank())
  
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
