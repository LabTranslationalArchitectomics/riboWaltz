#' Duplicates filtering.
#'
#' This function provides multiple options for remove duplicated reads: when two
#' or more reads are marked as duplicates, all of them are discarded but one.
#'
#' @param data List of data tables from \code{\link{bamtolist}},
#'   \code{\link{bedtolist}} or \code{\link{length_filter}}.
#' @param sample Character string or character string vector specifying the name
#'   of the sample(s) to process. Default is NULL i.e. all samples are
#'   processed.
#' @param extremity Either "both", "5end", "3end". It specifies the criterion to
#'   define which reads should be considered duplicates. Reads are marked as
#'   duplicates if they map on the same transcript and share: both the 5'
#'   estremity and the 3' extremity ("both"), only the 5' extremity ("5end"),
#'   only the same 3' extremity ("3end "). For "5end" and "3end", reads of
#'   different lengths can be marked as duplicates. See \code{keep} to choose
#'   which one should be kept.
#' @param keep Either "shortest" or "longest". It specifies wheter to keep the
#'   shortest or the longest read when duplicates display different lengths.
#'   This parameter is considered only if \code{extremity} is set to "5end" or
#'   "3end". Default is "shortest".
#' @param granges Logical value whether to return a GRangesList object. Default
#'   is FALSE i.e. a list of data tables is returned instead (the required input
#'   for \code{\link{length_filter}}, \code{\link{psite}},
#'   \code{\link{psite_info}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}}).
#' @return A list of data tables or a GRangesList object.
#' @examples
#' #generate an \emph{ad hoc} dataset:
#' library(data.table)
#' dt <- data.table(transcript = rep("ENSMUST00000000001.4", 6),
#'                  end5 = c(92, 92, 92, 94, 94, 95),
#'                  end3 = c(119, 119, 122, 122, 123, 123)
#'                  )[, length := end3 - end5 + 1
#'                    ][, cds_start := 14
#'                     ][, cds_stop := 1206]
#' example_reads_list <- list()
#' example_reads_list[["Samp_example"]] <- dt
#' 
#' ## Reads are duplicates if they share both the 5' estremity and the
#' ## 3' extremity:
#' filtered_list <- duplicates_filter(example_reads_list,
#'                                    extremity = "both")
#' 
#' ## Reads are duplicates if they only share the 5' estremity. Among duplicated 
#' ## reads we keep the shortes one:
#' filtered_list <- duplicates_filter(example_reads_list,
#'                                    extremity = "5end",
#'                                    keep = "shortest")
#' @import data.table
#' @export
duplicates_filter <- function(data, sample = NULL, extremity = "both",
                              keep = "shortest", granges = FALSE){
  
  check_sample <- setdiff(unlist(sample), names(data))
  if(length(check_sample) != 0){
    cat("\n")
    stop(sprintf("incorrect sample name(s): \"%s\" not found\n\n",
                 paste(check_sample, collapse = ", ")))
  }
  
  if(is.null(sample)){
    sample <- names(data)
  }
  
  if(!(extremity %in% c("both", "5end", "3end"))){
    cat("\n")
    warning("parameter extremity must be either \"both\", \"5end\" or \"3end\"\nset to default \"both\"\n", call. = FALSE)
    extremity = "both"
  }
  if(extremity %in% c("5end", "3end")){
    if(!(keep %in% c("shortest", "longest"))){
      cat("\n")
      warning("parameter keep must be either \"shortest\" or \"longest\"\nset to default \"shortest\"\n", call. = FALSE)
      keep = "shortest"
    }
  }

  for(samp in sample){
    cat(sprintf("processing %s\n", samp))
    dt <- data[[samp]]
    
    nreads <- nrow(dt)
    cat(sprintf("reads: %s M\n", format(round((nreads / 1000000), 2), nsmall = 3)))
    
    if(!is.null(extremity)){
      if(extremity == "both") {
        dt <- unique(dt, by = c("transcript" ,"end5", "end3"))
      } else {
        if(keep == "longest"){
          dt <- dt[order(transcript, end5, -end3)]
        }
        if(extremity == "5end") {
          dt <- unique(dt, by = c("transcript", "end5"))
        } else {
          dt <- unique(dt, by = c("transcript", "end3"))
        }
      }
    } 
    
    dt <- dt[order(transcript, end5, end3)]
    
    cat(sprintf("%s M  (%s %%) reads removed\n", 
                format(round((nreads - nrow(dt))/ 1000000, 2), nsmall = 3), 
                format(round(((nreads - nrow(dt))/nreads) * 100, 2), nsmall = 3) ))
    cat(sprintf("reads kept: %s M\n\n", format(round((nrow(dt) / 1000000), 2), nsmall = 3)))
    
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






