#' Annotation data table.
#'
#' This function generates transcript basic annotation data tables starting from
#' GTF files or TxDb objects. Annotation data tables include a column named
#' \emph{transcript} reporting the name of the reference transcripts and four
#' columns named \emph{l_tr}, \emph{l_utr5}, \emph{l_cds} and \emph{l_utr3}
#' reporting the length of the transcripts and of their annotated 5' UTRs, CDSs
#' and 3' UTRs, respectively. Please note: if a transcript region is not
#' annotated its length is set to 0.
#' 
#' @param gtfpath A character string specifying the path to a GTF file,
#'   including its name and extension. Please make sure the GTF file derives
#'   from the same release of the sequences used in the alignment step. Note
#'   that either \code{gtfpath} or \code{txdb} must be specified. Default is
#'   NULL.
#' @param txdb Character string specifying the TxDb annotation package to be
#'   loaded. If not already present in the system, it is automatically installed
#'   through the biocLite.R script (check
#'   \href{http://bioconductor.org/packages/release/BiocViews.html#___TxDb}{here}
#'    the list of available TxDb annotation packages). Please make sure the TxDb
#'   annotation package derives from the same release of the sequences used in
#'   the alignment step. Note that either \code{gtfpath} or \code{txdb} must be
#'   specified. Default is NULL.
#' @param dataSource Optional character string describing the origin of the GTF
#'   data file. This parameter is considered only if \code{gtfpath} is
#'   specified. For more information about this parameter please refer to the
#'   description of \emph{dataSource} of the
#'   \code{\link[GenomicFeatures]{makeTxDbFromGFF}} function included in the
#'   \code{GenomicFeatures} package.
#' @param organism Optional character string reporting the genus and species of
#'   the organism of the GTF data file. This parameter is considered only if
#'   \code{gtfpath} is specified. For more information about this parameter
#'   please refer to the description of \emph{organism} of the
#'   \code{\link[GenomicFeatures]{makeTxDbFromGFF}} function included in the
#'   \code{GenomicFeatures} package.
#' @return A data table.
#' @examples
#' ## gtf_file <- "path/to/GTF/file.GTF"
#' ## create_annotation(gtfpath = gtf_file, dataSource = "gencode6", organism = "Mus musculus")
#' @import data.table
#' @export
create_annotation <-  function(gtfpath = NULL, txdb = NULL, dataSource = NA, organism = NA) {
  
  if(length(gtfpath) != 0 & length(txdb) != 0){
    cat("\n")
    warning("gtfpath and txdb are both specified. Only gtfpath will be considered\n")
    txdb = NULL
  }
  
  if(length(gtfpath) == 0 & length(txdb) == 0){
    cat("\n")
    stop("neither gtfpath nor txdb is specified\n\n")
  }
  
  if(length(gtfpath) != 0){
    path_to_gtf <- gtfpath
    txdbanno <- GenomicFeatures::makeTxDbFromGFF(file = path_to_gtf, format = "gtf", dataSource = dataSource, organism = organism)
  } else {
    if(txdb %in% rownames(installed.packages())){
      library(txdb, character.only = TRUE)
    } else {
      source("https://bioconductor.org/biocLite.R")
      biocLite(txdb, suppressUpdates = TRUE)
      library(txdb, character.only = TRUE)
    }
    txdbanno <- get(txdb)
  }
  
  exon <- suppressWarnings(GenomicFeatures::exonsBy(txdbanno, by = "tx",use.names=T))
  utr5<- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(txdbanno,use.names=T))
  cds <- suppressWarnings(GenomicFeatures::cdsBy(txdbanno, by = "tx", use.names=T))
  utr3<- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(txdbanno,use.names=T))
  exon <- as.data.table(exon[unique(names(exon))])
  utr5 <- as.data.table(utr5[unique(names(utr5))])
  cds <- as.data.table(cds[unique(names(cds))])
  utr3 <-as.data.table(utr3[unique(names(utr3))])
  
  anno_df <- exon[, list(l_tr = sum(width)), by = list(transcript = group_name)]
  l_utr5 <- utr5[, list(l_utr5 = sum(width)), by = list(transcript = group_name)]
  l_cds <- cds[, list(l_cds = sum(width)), by = list(transcript = group_name)]
  l_utr3 <- utr3[, list(l_utr3 = sum(width)), by = list(transcript = group_name)]
  
  merge_allx <- function(x, y) merge(x, y, all.x=TRUE)
  anno_df  <-  Reduce(merge_allx, list(anno_df, l_utr5, l_cds, l_utr3))
  anno_df[is.na(anno_df)] <- 0
  
  return(anno_df)
}
