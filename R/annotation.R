#' Create an annotation data frame.
#'
#' Starting from a GTF file or a TxDb object this function generates a dada
#' frame containing a basic annotation of the transcripts. The data frame
#' includes a column named \emph{transcript} reporting the name of the reference
#' sequences and four columns named \emph{l_tr}, \emph{l_utr5}, \emph{l_cds} and
#' \emph{l_utr3} reporting the length of the transcripts and of their annotated
#' \emph{5' UTR}, \emph{CDS} and \emph{3' UTR}, respectively.
#' @param gtfpath A character string specifying the path to te GTF file,
#'   including its name and extension. Please make sure the GTF derives from the
#'   same release of the sequences used in the alignment step. Note that either
#'   \code{gtfpath} or \code{txdb} must be specified.
#' @param txdb A character string specifying the name of the annotation package
#'   for TxDb object(s) to be loaded. If it is not already present in your
#'   system, it will be installed through the biocLite.R script (check the list
#'   of TxDb annotation packages available in the Bioconductor repositories at
#'   http://bioconductor.org/packages/release/BiocViews.html#___TxDb )). Please
#'   make sure the annotation package derives from the same release of the
#'   sequences used in the alignment step. Note that either \code{gtfpath} or
#'   \code{txdb} must be specified.
#' @param dataSource An optional character string describing the origin of the
#'   GTF data file. For more information about this parameter please refer to
#'   the description of \emph{dataSource} of the
#'   \code{\link[GenomicFeatures]{makeTxDbFromGFF}} function included in the
#'   \code{GenomicFeatures} package.
#' @param organism A optional character string reporting the genus and species
#'   of the organism when \code{gtfpath} is specified.  For more information
#'   about this parameter please refer to the description of \emph{dataSource}
#'   of the \code{\link[GenomicFeatures]{makeTxDbFromGFF}} function included in
#'   the \code{GenomicFeatures} package.
#' @return A data frame.
#' @examples
#' gtf_file <- location_of_GTF_file
#' path_bed <- location_of_output_directory
#' bamtobed(gtfpath = gtf_file, dataSource = "gencode6", organism = "Mus musculus")
#' @export
create_annotation  <-  function(gtfpath = NULL, txdb = NULL, dataSource = NA, organism = NA) {
  
  if(length(gtfpath) != 0 & length(txdb) != 0){
    warning("gtfpath and txdb are both specified. Only gtfpath will be considered\n")
    txdb = NULL
  }
  
  if(length(gtfpath) == 0 & length(txdb) == 0){
    cat("\n")
    stop("\nERROR: neither gtfpath nor txdb is specified \n\n")
  }
  
  if(length(gtfpath) != 0){
    path_to_gtf <- gtfpath
    txdbanno <- GenomicFeatures::makeTxDbFromGFF(file=path_to_gtf, format="gtf", dataSource = dataSource, organism = organism)
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
  
  exon <- GenomicFeatures::exonsBy(txdbanno, by = "tx",use.names=T)
  utr5<- GenomicFeatures::fiveUTRsByTranscript(txdbanno,use.names=T)
  cds <- GenomicFeatures::cdsBy(txdbanno, by = "tx", use.names=T)
  utr3<- GenomicFeatures::threeUTRsByTranscript(txdbanno,use.names=T)
  
  anno_df <- data.frame("transcript"=names(exon), "l_tr" = sapply(exon,function(x) sum(width(x))))
  l_utr5<-data.frame("transcript"=names(utr5),"l_utr5"= sapply(utr5,function(x) sum(width(x))))
  l_cds<-data.frame("transcript"=names(cds),"l_cds"=sapply(cds,function(x) sum(width(x))))
  l_utr3<-data.frame("transcript"=names(utr3),"l_utr3"=sapply(utr3,function(x) sum(width(x))))
  
  merge_allx <- function(x, y) merge(x, y, all.x=TRUE)
  anno_df  <-  Reduce(merge_allx, list(anno_df, l_utr5, l_cds, l_utr3))
  anno_df[is.na(anno_df)] <- 0
  
  return(anno_df)
}




