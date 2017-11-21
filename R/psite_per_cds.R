#' Compute the number of P-sites per coding sequence.
#'
#' Computes the number of P-sites in frame locatd within the coding sequence of 
#' each transcript for all the samples of the input list. It is possible to 
#' restrict the analysis to a subsequence of the coding region by exluding a 
#' chosen number of nucleotides at the beginiing and at the end of the CDS. Only
#' the transcripts associated to an annotated CDS are considered. The resulting
#' data frame report the name of the transript along with the length of the
#' considered region (in nucleotide) and the number of p-sites.
#'
#' @param data A list of data frames from \code{\link{psite_info}}.
#' @param annotation A data frame from \code{\link{create_annotation}}.
#' @param start_nts A positive integer specifying the number of nucleotides at 
#'   the beginning of the coding sequences to be exluded from the count. Default
#'   is 0.
#' @param stop_nts A positive integer specifying the number of nucleotides at 
#'   the end of the coding sequences to be exluded from the count. Default is 0.
#' @return A list of data frames.
#' @examples
#' data(reads_psite_list)
#' data(mm81cdna)
#'
#' ## Compute the number of p-sites in frame along the whole coding sequence.
#' psite_cds_list <- psite_per_cds(reads_psite_list, mm81cdna)
#' 
#' ## Compute the number of p-sites in frame along the coding sequence exluding
#' the first 15 nucleotides and the last 10 nucleotides.
#' psite_cds_list <- psite_per_cds(reads_psite_list, mm81cdna, start_nts = 15, stop_nts = 10)
#' @export
psite_per_cds <- function(data, annotation, start_nts = 0, stop_nts = 0) {
  names <- names(data)
  temp <- subset(annotation, l_cds > start_nts + stop_nts)[, c("transcript", "l_cds")]
  temp$l_cds <- temp$l_cds - start_nts - stop_nts
  colnames(temp) <- c("transcript", "l_region")
  psite_cds <- list()
  for (n in names) {
    cat(sprintf("processing %s\n", n))
    psite_cds[[n]] <- temp
    df <- subset(data[[n]], as.character(transcript) %in% as.character(temp$transcript) &
                   psite_from_start >= start_nts &
                   psite_from_stop <= -stop_nts &
                   psite_from_start %% 3 == 0)
    psite_cds[[n]]$psite_count <- table(factor(df$transcript,
                                               levels = as.character(psite_cds[[n]]$transcript)))
  }
  return(psite_cds)
}


