#' Read length filtering.
#'
#' This function provides multiple options for filtering the reads according to
#' their length. Read lengths to keep are either specified by the user or
#' automatichally selected on the basis of the trinucleotide periodicity of reads
#' mapping on the CDS.
#'
#' @param data List of data tables from \code{\link{bamtolist}},
#'   \code{\link{bedtolist}} or \code{\link{duplicates_filter}}.
#' @param sample Character string or character string vector specifying the name
#'   of the sample(s) to process. Default is NULL i.e. all samples are
#'   processed.
#' @param length_filter_mode Either "periodicity" or "custom". It specifies how
#'   read length selection should be performed. "periodicity": only read lengths
#'   satisfying a periodicity threshold (see \code{periodicity_threshold}) are
#'   kept. It ensures the removal of all reads with low or no periodicity;
#'   "custom": only read lengths specified by the user are kept (see
#'   \code{length_range}). Default is "periodicity".
#' @param periodicity_threshold Integer in [10, 100]. Only read lengths
#'   satisfying this threshold (i.e. a higher percentage of read extremities
#'   falls in one of the three reading frames along the CDS) are kept. This
#'   parameter is considered only if \code{length_filter_mode} is set to
#'   "periodicity". Default is 50.
#' @param length_range Integer or integer vector specifying one read
#'   length or a range of read lengths to keep, respectively. This parameter is
#'   considered only if \code{length_filter_mode} is set to "custom".
#' @param granges Logical value whether to return a GRangesList object. Default
#'   is FALSE i.e. a list of data tables is returned instead (the required input
#'   for \code{\link{duplicates_filter}}, \code{\link{psite}},
#'   \code{\link{psite_info}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}}).
#' @return A list of data tables or a GRangesList object.
#' @examples
#' data(reads_list)
#' 
#' ## Keep reads of length between 27 and 30 nucleotides (included):
#' filtered_list <- length_filter(reads_list, length_filter_mode = "custom",
#'                                length_range = 27:30)
#' 
#' ## Keep reads of lengths satisfying a periodicity threshold (70%):
#' filtered_list <- length_filter(reads_list, length_filter_mode = "periodicity",
#'                                periodicity_threshold = 70)
#' @import data.table
#' @export
length_filter <- function(data, sample = NULL,
                          length_filter_mode = "periodicity",
                          periodicity_threshold = 50,
                          length_range = NULL, granges = FALSE){
  
  check_sample <- setdiff(unlist(sample), names(data))
  if(length(check_sample) != 0){
    cat("\n")
    stop(sprintf("incorrect sample name(s): \"%s\" not found\n\n",
                 paste(check_sample, collapse = ", ")))
  }
  
  if(is.null(sample)){
    sample <- names(data)
  }
  
  if(!(length_filter_mode %in% c("custom", "periodicity"))){
    cat("\n")
    warning("parameter length_filter_mode must be either \"periodicity\" or \"custom\"\nset to default \"periodicity\"\n", call. = FALSE)
    length_filter_mode = "periodicity"
  }
  
  if(length_filter_mode == "custom" & !is.numeric(length_range)
     & !is.integer(length_range)){
    stop("length_range must be of class integer\n\n")
  }
  
  if(length_filter_mode == "periodicity" & ((!is.numeric(periodicity_threshold)
     & !is.integer(periodicity_threshold)) | periodicity_threshold < 10
     | periodicity_threshold > 100)){
    stop("periodicity_threshold must be an integer between 10 and 100 \n\n")
  }
  
  for(samp in sample) {
    cat(sprintf("processing %s\n", samp))
    dt <- data[[samp]]
    
    nreads <- nrow(dt)
    cat(sprintf("reads: %s M\n", format(round((nreads / 1000000), 2), nsmall = 2)))

    if(identical(length_filter_mode, "custom")) {
      dt <- dt[length %in% length_range]
      cat(sprintf("%s M  (%s %%) reads removed\n", 
                  format(round((nreads - nrow(dt))/ 1000000, 2), nsmall = 2), 
                  format(round(((nreads - nrow(dt))/nreads) * 100, 2), nsmall = 2) ))
      cat(sprintf("reads (kept): %s M\n\n", format(round((nrow(dt) / 1000000), 2), nsmall = 2)))
    } else {
      if(identical(length_filter_mode, "periodicity")){
        nreads <- nrow(dt)
        
        subdt5 <- dt[cds_start != 0 &
                       (end5 - cds_start) >= 0 &
                       (cds_stop - end5) >= 0]
        subdt5[, end5_frame := as.factor((end5 - cds_start) %% 3)]
        t_end5 <- subdt5[, .N, by = list(length, end5_frame)
                         ][, end5_perc := (N / sum(N)) * 100, by = length]
        keep_length5 <- unique(t_end5[end5_perc >= periodicity_threshold, length])
        
        subdt3 <- dt[cds_start != 0 &
                       (end3 - cds_start) >= 0 &
                       (cds_stop - end3) >= 0]
        subdt3[, end3_frame := as.factor((end3 - cds_start) %% 3)]
        t_end3 <- subdt3[, .N, by = list(length, end3_frame)
                         ][, end3_perc := (N / sum(N)) * 100, by = length]
        keep_length3 <- unique(t_end3[end3_perc >= periodicity_threshold, length])
        
        keep_length <- intersect(keep_length5, keep_length3)
        dt <- dt[length %in% keep_length]
        
        cat(sprintf("%s M  (%s %%) reads removed\n", 
                    format(round((nreads - nrow(dt))/ 1000000, 2), nsmall = 2), 
                    format(round(((nreads - nrow(dt))/nreads) * 100, 2), nsmall = 2) ))
        cat(sprintf("reads (kept): %s M\n\n", format(round((nrow(dt) / 1000000), 2), nsmall = 2)))
      }
    }

    if (granges == T || granges == TRUE) {
      dt <- GenomicRanges::makeGRangesFromDataFrame(dt,
                                                    keep.extra.columns = TRUE,
                                                    ignore.strand = TRUE,
                                                    seqnames.field = c("transcript"),
                                                    start.field = "end5",
                                                    end.field = "end3",
                                                    strand.field = "strand",
                                                    starts.in.df.are.0based = FALSE)
      GenomicRanges::strand(dt) <- "+"
    }
    
    data[[samp]] <- dt
  }
  
  if (granges == T || granges == TRUE) {
    data <- GenomicRanges::GRangesList(data)
  }
  
  return(data)
}
