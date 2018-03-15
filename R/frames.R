#' Compute the percentage of P-sites per frame.
#'
#' For one or several samples this function computes the percentage of P-sites
#' falling on the three reading frames of the transcripts and generates a
#' barplot of the resulting values. This analysis is performed for the annotated
#' 5' UTR, coding sequence and 3' UTR, separately. It is possible to compute the
#' percentage of P-sites per frame using all the read lengths or to restrict the
#' analysis to a sub-range of read lengths.
#'
#' @param data A list of data frames from \code{\link{psite_info}}.
#' @param sample A character string vector specifying the name of the sample(s)
#'   of interest. By default this argument is NULL, meaning that all the samples
#'   in \code{data} are included in the analysis.
#' @param region Either "all" or a character string among "5utr", "cds", "3utr"
#'   specifying the regions of the transcript (5' UTR, CDS or 3' UTR,
#'   respectively) that must be included in the analysis. Default is "all",
#'   meaning that the all the regions are considered.
#' @param length_range Either "all", an integer or an integer vector. Default is
#'   "all", meaning that all the read lengths will be included in the analysis.
#'   Otherwise, only the read lengths matching the specified value(s) are kept.
#' @return A list containing a ggplot2 object and a data frame with the
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
#' length_range=28)
#' frame_sub[["plot"]]
#' @import ggplot2
#' @export
frame_psite<-function(data, sample=NULL, region="all", length_range="all"){
  if(length(sample) == 0) {
    sample <- names(data)
  }

  if(!region%in%c("all", "cds", "5utr", "3utr")){
    warning("region is invalid. Set to default \"all\"\n")
    region="all"
  }

  if(!identical(length_range, "all") & !inherits(length_range, "numeric") & !inherits(length_range, "integer")){
    warning("length_range is invalid. Set to default \"all\"\n")
    length_range = "all"
  }

  for (samp in sample) {
    if (region == "all") {

      if(length_range[1] == "all"){
        df <- data[[samp]]
      }else{
        df <- subset(data[[samp]], length%in%length_range)
      }
      
      df <- subset(df, start_pos!=0 & stop_pos!=0)
      df$frame <- df$psite_from_start %% 3
      frame_df <- as.data.frame(table(df$frame, factor(df$psite_region, levels = c("5utr", "cds", "3utr"), labels=c("5' UTR", "CDS", "3' UTR"))))
      colnames(frame_df) <- c("frame", "region", "count")
      frame_df$perc <- (frame_df$count/rep(by(frame_df$count, frame_df$region, function(x) sum(x)), each = 3)) * 100
    } else {
      df <- subset(data[[samp]], psite_region == region)
      df$frame <- df$psite_from_start %% 3
      frame_df <- as.data.frame(table(df$frame))
      colnames(frame_df) <- c("frame", "count")
      frame_df$perc <- (frame_df$count/sum(frame_df$count)) * 100
    }
    frame_df$samp <- samp

    if (exists("final_frame_df")) {
      final_frame_df <- rbind(final_frame_df, frame_df)
    } else {
      final_frame_df <- frame_df
    }
    rm(frame_df)
  }

  if(region == "all") {
    plottitle_region <- NULL
  } else {
    if(region == "5utr") { plottitle_region <- "Region: 5' UTR" }
    if(region == "cds") { plottitle_region <- "Region: CDS" }
    if(region == "3utr") { plottitle_region <- "Region: 3' UTR" }
  }

  if(length_range[1] == "all") {
    plottitle_range <- NULL
  } else{
    if(min(length_range) == max(length_range)) {
      plottitle_range <- paste("Read length: ", min(length_range), " nts", sep = "")
    } else {
      if(identical(length_range, min(length_range) : max(length_range)) | identical(length_range, seq(min(length_range), max(length_range), 1))){
        plottitle_range <- paste("Read lengths: ", min(length_range), "-", max(length_range), " nts", sep = "")
      } else {
        plottitle_range <- paste("Read lengths: ", paste(length_range, collapse=","), " nts", sep = "")
      }
    }
  }

  plot<-ggplot(final_frame_df, aes(x = frame, y = perc)) +
    geom_bar(stat = "identity") +
    theme_bw(base_size = 20) +
    labs(x = "Frame", y = "P-site signal (%)", title = paste(plottitle_region, plottitle_range, sep="; "))

  if(region == "all") {
    plot <- plot + facet_grid(samp ~ region)
  } else {
    plot <- plot + facet_wrap( ~ samp, ncol = 3)
  }

  ret_list<-list()
  ret_list[["df"]]<-final_frame_df
  ret_list[["plot"]]<-plot
  return(ret_list)
}

#' Compute the number of P-sites per frame stratified by read length.
#'
#' Similar to \code{\link{frame_psite}} but the results are stratified by the
#' length of the reads.
#'
#' @param data A list of data frames from \code{\link{psite_info}}.
#' @param sample A character string vector specifying the name of the sample(s)
#'   of interest. By default this argument is NULL, meaning that all the samples
#'   in \code{data} are included in the analysis.
#' @param region Either "all" or a character string among "5utr", "cds", "3utr"
#'   specifying the regions of the transcript (5' UTR, CDS or 3' UTR,
#'   respectively) that must be included in the analysis. Default is "all",
#'   meaning that the all the regions are considered.
#' @param cl An integer value in \emph{[1,100]} specifying the confidence level
#'   for restricting the analysis to a sub-range of read lengths. By default it is
#'   set to 100. This parameter has no effect if \code{length_range} is
#'   specified.
#' @param length_range Either "all", an integer or an integer vector. Default is
#'   "all", meaning that all the read lengths will be included in the analysis.
#'   Otherwise, only the read lengths matching the specified value(s) are kept.
#'   If specified, this parameter prevails over \code{cl}.
#' @return A list containing a ggplot2 object and a data frame with the
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
#' @import ggplot2
#' @export
frame_psite_length<-function(data, sample=NULL, region="all", cl=100, length_range=NULL){
  if(length(sample) == 0) {
    sample <- names(data)
  }

  if(!region%in%c("all", "cds", "5utr", "3utr")){
    warning("region is invalid. Set to default \"all\"\n")
    region="all"
  }

  for (samp in sample) {
    if (region == "all") {
      df <- data[[samp]]
      df <- subset(df, start_pos!=0 & stop_pos!=0)
      df$frame <- df$psite_from_start %% 3

      minl <- quantile(df$length, (1 - cl/100)/2)
      maxl <- quantile(df$length, 1 - (1 - cl/100)/2)
      if(length(length_range) != 0){
        if(!inherits(length_range, "numeric") & !inherits(length_range, "integer")){
          warning("length_range is invalid. Confidence interval is used\n")
        } else {
          minl <- min(length_range)
          maxl <- max(length_range)
        }
      }
      lenmax <- length(minl:maxl)

      frame_df <- as.data.frame(table(factor(df$length, levels = minl:maxl), df$frame, factor(df$psite_region, levels = c("5utr", "cds", "3utr"), labels=c("5' UTR", "CDS", "3' UTR"))))
      colnames(frame_df) <- c("length", "frame", "region", "count")
      sum_count <- as.vector(by(frame_df$count, frame_df[,c("length", "region")], function(x) sum(x)))
      frame_df$perc<-(frame_df$count / c(rep(sum_count[1:lenmax],3),rep(sum_count[(lenmax+1):(2 * lenmax)],3),rep(sum_count[(2 * lenmax+1):(3 * lenmax)],3))) * 100
    } else {
      df <- subset(data[[samp]], psite_region == region)
      df$frame <- df$psite_from_start %% 3

      minl <- quantile(df$length, (1 - cl/100)/2)
      maxl <- quantile(df$length, 1 - (1 - cl/100)/2)
      if(length(length_range) != 0){
        if(!inherits(length_range, "numeric") & !inherits(length_range, "integer")){
          warning("length_range is invalid. Confidence interval is used\n")
        } else {
          minl <- min(length_range)
          maxl <- max(length_range)
        }
      }
      lenmax <- length(minl:maxl)

      frame_df <- as.data.frame(table(factor(df$length, levels = minl:maxl), df$frame))
      colnames(frame_df) <- c("length", "frame", "count")
      sum_count <- as.vector(by(frame_df$count, frame_df$length, function(x) sum(x)))
      frame_df$perc<-(frame_df$count / rep(sum_count, 3)) * 100
    }
    frame_df$samp <- samp

    if (exists("final_frame_df")) {
      final_frame_df <- rbind(final_frame_df, frame_df)
    } else {
      final_frame_df <- frame_df
    }
    rm(frame_df)
  }

  final_frame_df[is.na(final_frame_df$perc),"perc"]<-0

  plot<-ggplot(final_frame_df, aes(frame, as.numeric(as.character(length)))) +
    geom_tile(aes(fill = perc)) +
    scale_fill_gradient("P-site signal (%) ", low = "white", high = "#104ec1",
                        breaks = c(min(final_frame_df$perc),(max(final_frame_df$perc)-min(final_frame_df$perc))/2,max(final_frame_df$perc)),
                        labels = c(round(min(final_frame_df$perc)), round((max(final_frame_df$perc)-min(final_frame_df$perc))/2),round(max(final_frame_df$perc)))) +
    labs(x="Frame",y="Read length") +
    theme_bw(base_size = 20) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(legend.position="top") +
    scale_y_continuous(limits=c(minl-0.5, maxl+0.5), breaks = seq(minl + ((minl) %% 2), maxl, by=max(2, floor((maxl-minl)/7))))
  if (region == "all") {
    plot <- plot + facet_grid(samp ~ region)
  } else {
    plot <- plot + facet_wrap( ~ samp, ncol = 3)
  }

  ret_list<-list()
  ret_list[["df"]]<-final_frame_df
  ret_list[["plot"]]<-plot
  return(ret_list)
}








