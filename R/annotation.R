#' Create an annotation data frame.
#' 
#' Creates a data frame with a basic reference annotation of the transripts 
#' starting from a GTF file. The data frame will contain a column named 
#' \emph{transcript} reporting the name of the transcripts (as reported in the 
#' GTF, usually the ENST ID and version, dot separated) and four columns named 
#' \emph{l_tr}, \emph{l_utr5}, \emph{l_cds} and \emph{l_utr3} reporting the 
#' length of the whole transcript and of the annotated \emph{5' UTR}, the 
#' \emph{CDS} and the \emph{3' UTR}, respectively. Please note that the 
#' transcript names should coincide with those in the alignment file.
#' 
#' @param gtfpath A character string specifying the complete path to te GTF
#'   file, including its name and extension. Either \code{gtfpath} or
#'   \code{txdb} must be specified.
#' @param txdb A TxDb object storing transcript annotations. Either
#'   \code{gtfpath} or \code{txdb} must be specified.
#' @param dataSource A character string describing the origin of the data file.
#'   Please refer to the description of \emph{dataSource} of the
#'   \code{makeTxDbFromGFF} function included in the \code{GenomicFeatures}
#'   package.
#' @param organism A character string specifying the Genus and species of this
#'   organism. Please refer to the description of  \emph{organism} of the
#'   \code{makeTxDbFromGFF} function included in the \code{GenomicFeatures}
#'   package.
#' @return A data frame.
#' @examples
#' gtf_file <- location_of_GTF_file
#' path_bed <- location_of_output_directory
#' bamtobed(gtfpath = gtf_file, dataSource = "gencode6", organism = "Mus musculus")
#' @import GenomicFeatures
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
    txdbanno <- txdb
  }
  
  exon <- exonsBy(txdbanno, by = "tx",use.names=T)
  utr5<- fiveUTRsByTranscript(txdbanno,use.names=T)
  cds <- cdsBy(txdbanno, by = "tx", use.names=T)
  utr3<-threeUTRsByTranscript(txdbanno,use.names=T)
  
  anno_df <- data.frame("transcript"=names(exon), "l_tr" = sapply(exon,function(x) sum(width(x))))
  l_utr5<-data.frame("transcript"=names(utr5),"l_utr5"= sapply(utr5,function(x) sum(width(x))))
  l_cds<-data.frame("transcript"=names(cds),"l_cds"=sapply(cds,function(x) sum(width(x))))
  l_utr3<-data.frame("transcript"=names(utr3),"l_utr3"=sapply(utr3,function(x) sum(width(x))))
  
  merge_allx <- function(x, y) merge(x, y, all.x=TRUE)
  anno_df  <-  Reduce(merge_allx, list(anno_df, l_utr5, l_cds, l_utr3))
  anno_df[is.na(anno_df)] <- 0
  
  return(anno_df)
}




