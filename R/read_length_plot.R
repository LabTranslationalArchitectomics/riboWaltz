#' Read length distributions.
#'
#' This function generates read length distributions, displayed as bar plots.
#' Multiple samples and replicates can be handled..
#'
#' @param data Either list of data tables or GRangesList object from
#'   \code{\link{bamtolist}}, \code{\link{bedtolist}},
#'   \code{\link{length_filter}} or \code{\link{psite_info}}.
#' @param sample Either character string, character string vector or named list
#'   of character string(s)/character string vector(s) specifying the name of
#'   the sample(s) and replicate(s) of interest. If a list is provided, each
#'   element of the list is considered as an independent sample associated with
#'   one ore multiple replicates. Multiple samples and replicates are handled
#'   and visualised according to \code{multisamples} and \code{plot_style}.
#' @param multisamples Either "average" or "independent". It specifies how to
#'   handle multiple samples and replicates stored in \code{sample}:
#'   * if \code{sample} is a character string vector and \code{multisamples} is
#'   set to "average" the elements of the vector are considered as replicates
#'   of one sample and a single bar plot is returned.
#'   * if \code{sample} is a character string vector and \code{multisamples} is
#'   set to "independent", each element of the vector is analysed independently
#'   of the others. The number of plots returned and their organization is
#'   specified by \code{plot_style}.
#'   * if \code{sample} is a list, \code{multisamples} must be set to "average".
#'   Each element of the list is analysed independently of the others, its
#'   replicates averaged and its name reported in the plot. The number of plots
#'   returned and their organization is specified by \code{plot_style}.
#'   Note: when this parameter is set to "average" the bar plot associated with
#'   each sample displays the length-specific mean signal computed across the
#'   replicates and the corresponding standard error is also reported. Default
#'   is "average".
#' @param plot_style Either "split", "facet", "dodge" or "mirror". It specifies
#'   how to organize and display multiple bar plots:
#'   * "split": one bar plot for each sample is returned as an independent
#'   ggplot object;
#'   * "facet": the bar plots are placed one next to the other, in
#'   independent boxes;
#'   * "dodge": all bar plots are displayed in one box and, for each length,
#'   samples are placed side by side.
#'   * "mirror": \code{sample} must be either a character string vector or
#'   a list of exactly two elements and the resulting bar plots are mirrored
#'   along the x axis.
#'   Default is "split".
#' @param scale_factors Either "auto", a named numeric vector or "none". It
#'   specifies how read length distributions should be scaled before merging
#'   multiple samples (if any):
#'   * "auto": each distribution is scaled so that the sum of all bars is 100.
#'   * named numeric vector: \code{scale_factors} must be the same length of
#'   unlisted \code{sample} and each scale factor must be named after the
#'   corresponding string in unlisted \code{sample}. No specific order is
#'   required. Each distribution is multiplied by the matching scale factor.
#'   * "none": no scaling is applied. 
#'   Default is "auto".
#' @param transcripts Character string vector listing the name of transcripts to
#'   be included in the analysis. Default is NULL, i.e. all transcripts are
#'   used.
#' @param length_range Integer or integer vector for restricting the plot to a
#'   chosen range of read lengths. Default is NULL, i.e. all read lengths are
#'   used. If specified, this parameter prevails over \code{cl}.
#' @param cl Integer value in [1,100] specifying a confidence level for
#'   restricting the plot to an automatically-defined range of read lengths. The
#'   new range is computed according to the most frequent read lengths, which
#'   accounts for the cl% of the sample and is defined by discarding the
#'   (100-cl)% of read lengths falling in the tails of the read lengths
#'   distribution. If multiple samples are analysed, a single range of read
#'   lengths is computed such that at least the cl% of all samples is
#'   represented. Default is 100.
#' @param colour Character string or character string vector specifying the
#'   colour of the bar plot(s). If \code{plot_style} is set to either "dodge" or
#'   "mirror", a colour for each sample is required. Default is NULL, i.e. the
#'   default R colour palette is used.
#' @return List containing: one or more ggplot object(s) and the data table with
#'   the corresponding x- and y-axis values ("plot_dt"); an additional data
#'   table with raw and scaled number of reads per length in each sample
#'   ("count_dt").
#' @examples
#' data(reads_list)
#' 
#' ## Generate fake samples and replicates
#' for(i in 2:6){
#'   samp_name <- paste0("Samp", i)
#'   set.seed(i)
#'   reads_list[[samp_name]] <- reads_list[["Samp1"]][sample(.N, 5000)]
#' }
#' 
#' ## Define the list of samples and replicate to use as input
#' input_samples <- list("S1" = c("Samp1", "Samp2"),
#'                       "S2" = c("Samp3", "Samp4", "Samp5"),
#'                       "S3" = c("Samp6"))
#' 
#' ## Generate the length distribution for a sub-range of read lengths:
#' example_length_dist <- rlength_distr(reads_list,
#'                                      sample = input_samples,
#'                                      multisamples = "average",
#'                                      plot_style = "facet",
#'                                      cl = 99,
#'                                      colour = c("#333f50", "#39827c", "gray70"))
#' @import data.table
#' @import ggplot2
#' @export
rlength_distr <- function(data, sample, multisamples = "average",
                          plot_style = "split", scale_factors = "auto",
                          transcripts = NULL, length_range = NULL, cl = 100,
                          colour = NULL) {
  
  if(class(data[[1]])[1] == "GRanges"){
    data_tmp <- list()
    for(i in unlist(sample)){
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
  
  if(is.numeric(scale_factors)) {
    if(!all(unlist(sample) %in% names(scale_factors))){
      cat("\n")
      stop("scale factor for one or more sample is missing\n\n")
    }
  }
  
  if(!(multisamples %in% c("average", "independent"))){
    cat("\n")
    warning("parameter multisamples must be either \"average\" or \"independent\".
            Set to default \"average\"\n", call. = FALSE)
    multisamples <- "average"
  }
  
  if(multisamples == "independent" & is.list(sample)) {
    cat("\n")
    warning("parameter multisamples is set to \"independent\" but parameter sample is a list:
            parameter multisamples will be coerced to default \"average\"\n", call. = FALSE)
    multisamples <- "average"
  }
  
  if(length(sample) != 2 & plot_style == "mirror") {
    cat("\n")
    warning("parameter plot_style is set to \"mirror\" but parameter sample is a list of dimension > 2:
            parameter plot_style will be coerced to default \"split\"\n", call. = FALSE)
    plot_style <- "split"
  }
  
  if(is.character(sample) & length(sample) == 1) {
    multisamples <- "independent"
    plot_style <- "split"
  }
  
  if(is.character(sample) & length(sample) > 1 & multisamples == "average") {
    sample <- list("Average" = sample)
    plot_style <- "split"
    cat("\n")
    warning("Default name of averaged samples is \"Average\":
            consider to use a named list of one element to provide a meaningful plot title\n", call. = FALSE)
  }

  if(is.list(sample) & length(sample) == 1){
    plot_style <- "split"
  }
  
  if(!(plot_style %in% c("mirror", "dodge", "split", "facet"))){
    cat("\n")
    warning("parameter plot_style must be either \"split\", \"facet\", \"dodge\" or \"mirror\".
            Set to default \"split\"\n", call. = FALSE)
    plot_style <- "split"
  }
 
  #define color vector
  if((length(colour) < length(sample)) & 
     ((plot_style %in% c("dodge", "mirror")) | 
      (plot_style %in% c("split", "facet") & length(colour) != 1))){
    
    if(length(colour) !=0){
      warning("Not enough colors specified:\n
            default ggplot color palette will be used\n", call. = FALSE)
    }
    
    default_gg_col <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    colour <- default_gg_col(length(sample))
    
  } else {
    if(plot_style %in% c("split", "facet") & length(colour) == 1){
      colour <- rep(colour, length(sample))
    }
  }
  
  #define length range taking into account all (unlisted) samples
  if(length(length_range) == 0){
    for(samp in as.character(unlist(sample))){
      if(length(transcripts) == 0) {
        dt <- data[[samp]]
      } else {
        dt <- data[[samp]][transcript %in% transcripts]
      }
      
      if(length(length_range) == 0){
        length_range <- seq(quantile(dt$length, (1 - cl/100)/2),
                            quantile(dt$length, 1 - (1 - cl/100)/2))
      } else {
        xmin <- min(min(length_range), quantile(dt$length, (1 - cl/100)/2))
        xmax <- max(max(length_range), quantile(dt$length, 1 - (1 - cl/100)/2))
        length_range <- seq(xmin, xmax)
      }
    }
  }
  
  xmin = min(length_range)
  xmax = max(length_range)
  
  # check if all samples have reads of the specified lengths
  # especially required if only one read length is passed
  if(length(length_range) != 0){
    if(is.list(sample)){
      samp_dt <- data.table(stack(sample))
      setnames(samp_dt, c("sample", "sample_l"))
    } else {
      samp_dt <- data.table("sample" = sample, "sample_l" = sample)
    }
    
    for(samp in samp_dt$sample){
      
      dt <- data[[samp]][cds_start != 0 & cds_stop !=0]
      
      if(length(transcripts) != 0) {
        dt <- dt[transcript %in% transcripts]
      }
      
      len_check <- unique(dt$length)
      if(sum(length_range %in% len_check) == 0) {
        cat("\n")
        warning(sprintf("\"%s\" doesn't contain any reads of the selected lengths: sample removed\n", samp), call. = FALSE)
        #select element of sample which include the sample to be removed (useful if sample is a list)
        sel_l_samp <- samp_dt[sample == samp, sample_l]
        #remove the sample from the list/vector
        if(is.list(samp)){
          sample[[sel_l_samp]] <- sample[[sel_l_samp]][sample[[sel_l_samp]] != samp]
        } else {
          sample <- sample[sample != samp]
        }
      }
    }
  }
  
  if(is.null(unlist(sample))){
    cat("\n")
    stop("none of the data tables listed in sample contains any reads of the specified lengths\n\n")
  }
  
  # compute count of reads of defined lengths and scale them 
  final_dt <- data.table()
  for(samp in as.character(unlist(sample))){
    if(length(transcripts) == 0) {
      dt <- data[[samp]]
    } else {
      dt <- data[[samp]][transcript %in% transcripts]
    }
    
    setkey(dt, length)
    dist_dt <- dt[CJ(length_range), .(count = .N), by = .EACHI]
    
    #scaling/normalization
    if(is.character(scale_factors) &  scale_factors[1] == "auto"){
      dist_dt[, scaled_count := (count / sum(count)) * 100]
      y_title <- "% reads"
    } else {
      y_title <- "# reads"
      if(is.numeric(scale_factors)){
        dist_dt[, scaled_count := count * scale_factors[samp]]
      } else {
        dist_dt[, scaled_count := count]
      }
    }
    
    dist_dt[, tmp_samp := samp]

    final_dt <- rbind(final_dt, dist_dt)
  }
  
  output <- list()
  output[["count_dt"]] <- copy(final_dt[, c("tmp_samp", "length", "count", "scaled_count")])
  if(is.character(scale_factors) & scale_factors[1] == "auto"){
    output[["count_dt"]][, scaled_count := NULL]
  }
  setnames(output[["count_dt"]], "tmp_samp", "sample")

  # compute mean and se of samples if a list is provided
  if(is.list(sample)){
    
      samp_dt <- data.table(stack(sample))
      setnames(samp_dt, c("tmp_samp", "sample"))
      
      # set names of samples as specified in parameter sample  
      final_dt <- merge.data.table(final_dt, samp_dt, sort = FALSE)[, tmp_samp := NULL]
      
      # compute mean and se
      plot_dt <- final_dt[, .(mean_scaled_count = mean(scaled_count),
                              se_scaled_count = sd(scaled_count/sqrt(.N))), by = .(length, sample)]
      
    if(any(lengths(sample) != 1)){
      output[["plot_dt"]] <- copy(plot_dt[, c("sample", "length", "mean_scaled_count", "se_scaled_count")])
      setnames(output[["plot_dt"]], c("length", "mean_scaled_count", "se_scaled_count"), c("x", "y", "y_se"))
    } else {
      output[["plot_dt"]] <- copy(final_dt[, c("sample", "length", "scaled_count")])
      setnames(output[["plot_dt"]], c("length", "scaled_count"), c("x", "y"))
    }

  } else {
    plot_dt <- final_dt[, sample := tmp_samp
                        ][, se_scaled_count := NA]
    setnames(plot_dt, "scaled_count", "mean_scaled_count")
    
    output[["plot_dt"]] <- copy(plot_dt[, c("sample", "length", "mean_scaled_count")])
    setnames(output[["plot_dt"]], c("length", "mean_scaled_count"), c("x", "y"))
  }
  
  plot_dt[, sample := factor(sample, levels = unique(sample))]

  oldw <- getOption("warn")
  options(warn=-1)
  
  if(plot_style == "split"){ 
    i <- 0
    for(samp in unique(plot_dt$sample)){ # generate a plot for each sample and store it
      i <- i + 1
      sel_col = colour[i]
      
      sub_plot_dt <- plot_dt[sample == samp]
      plot <- ggplot(sub_plot_dt, aes(as.numeric(length), mean_scaled_count)) +
        geom_bar(stat = "identity", fill = sel_col) +
        geom_errorbar(aes(ymin = mean_scaled_count - se_scaled_count,
                          ymax = mean_scaled_count + se_scaled_count),
                      width = 0.35, linewidth = 1.1, na.rm = T, color = sel_col,
                      show.legend = F) +
        labs(title = samp, x = "Read length", y = y_title) +
        theme_bw(base_size = 23) +
        scale_x_continuous(limits = c(xmin - 0.5, xmax + 0.5),
                           breaks = seq(xmin + ((xmin) %% 2), xmax,
                                        by = max(c(1, floor((xmax - xmin)/7))))) + 
        theme(panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5)) 
      
      output[[paste0("plot_", samp)]] <- plot
    }
  } else {
    if(plot_style == "mirror") {
      plot_dt[sample == unique(plot_dt$sample)[2], mean_scaled_count := -mean_scaled_count]
    }
    
    plot <- ggplot(plot_dt, aes(as.numeric(length), mean_scaled_count, fill = sample))

    if(plot_style %in% c("facet", "mirror")) {
      plot <- plot + geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean_scaled_count - se_scaled_count,
                          ymax = mean_scaled_count + se_scaled_count, color = sample),
                      width = 0.35, linewidth = 1.1, na.rm = T)
      if(identical(plot_style, "mirror")){
        plot <- plot + geom_hline(yintercept = 0, linetype = 2, color = "gray20")
      }
    } else {
      plot <- plot + geom_bar(stat = "identity", position = position_dodge(0.9)) +
        geom_errorbar(aes(ymin = mean_scaled_count - se_scaled_count,
                          ymax = mean_scaled_count + se_scaled_count, color = sample),
                      width = 0.35, linewidth = 1.1, na.rm = T, position = position_dodge(0.9))
    }

    plot <- plot + labs(x = "Read length", y = y_title) +
      theme_bw(base_size = 23) +
      scale_x_continuous(limits = c(xmin - 0.5, xmax + 0.5),
                         breaks = seq(xmin + ((xmin) %% 2), xmax,
                                      by = max(c(1, floor((xmax - xmin)/7))))) + 
      theme(panel.grid.minor.x = element_blank()) +
      scale_fill_manual(values = colour) + 
      scale_color_manual(values = colour) +
      scale_y_continuous(labels = abs)

    if(uniqueN(colour) > 1 & plot_style != "facet"){
      plot <- plot + theme(legend.position = c(0.98,1), legend.justification = c(1, 1),
                           legend.title = element_blank(), legend.background = element_blank())
    } else {
      plot <- plot + theme(legend.position = "none")
    }

    if(plot_style == "facet"){
      plot <- plot + facet_wrap(sample ~ ., ncol = ceiling(sqrt(length(sample)))) +
        theme(strip.background = element_blank())
    }

    output[["plot"]] <- plot

  } 
  
  options(warn = oldw)
  return(output)
}
