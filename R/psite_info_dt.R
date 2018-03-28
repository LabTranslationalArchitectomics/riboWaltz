#' Update reads information according to the inferred P-sites.
#' 
#' Starting ftom the P-site position identfied by \code{\link{psite}}, this
#' function updates the data frames that contains information about the reads.
#' It attaches to the data frames 4 columns reporting the P-site position with
#' respect to the 1st nucleotide of the transcript, the start and the stop codon
#' of the annotated coding sequence (if any) and the region of the transcript
#' (5' UTR, CDS, 3' UTR) that includes the P-site. Please note: if a transcript
#' is not associated to any annotated CDS then the positions of the P-site from
#' both the start and the stop codon is set to NA. If either a FASTA file or a
#' BSgenome data package with the nucleotide  sequences is provided, an
#' additional column reporting the three nucleotides covered by the P-site is
#' attached.
#'
#' @param data A list of data tables from either \code{\link{bamtolist}} or
#'   \code{\link{bedtolist}}.
#' @param offset A data table from \code{\link{psite}}.
#' @param fastapath An optional character string specifying the path to the
#'   FASTA file used in the alignment step, including its name and extension.
#'   This file can contain reference nucleotide sequences either of a genome
#'   assembly or of all the transcripts (see \code{fasta_genome}). Please make
#'   sure the sequences derive from the same release of the annotation file used
#'   in the \code{\link{create_annotation}} function. Note: either
#'   \code{fastapath} or \code{bsgenome} is required to generate an
#'   additional column reporting the three nucletotides covered by the P-sites.
#'   Default is NULL.
#' @param fasta_genome A locigal value whether or not the FASTA file specified
#'   by \code{fastapath} contains nucleotide sequences of a genome assembly.
#'   FALSE means that the nucleotide sequences of all the transcripts are
#'   provided instead. When this parameter is TRUE (the default), an annotation
#'   object is required (see \code{gtfpath} and \code{txdb}).
#' @param bsgenome An optional character string specifying the name of the
#'   BSgenome data package containing the genome sequences to be loaded. If it
#'   is not already present in your system, it will be installed through the
#'   biocLite.R script (check the list of data packages available in the
#'   Bioconductor repositories for your version of R/Bioconductor by the
#'   \code{\link[BSgenome]{available.genomes}} function of the BSgenome
#'   package)). This parameter also requires an annotation object (see
#'   \code{gtfpath} and \code{txdb}). Please make sure the sequences included in
#'   the specified BSgenome data pakage are in agreement with the sequences used
#'   in the alignment step. Note: either \code{fastapath} or \code{bsgenome}
#'   is required to generate an additional column reporting the three
#'   nucletotides covered by the P-sites. Default is NULL.
#' @param gtfpath A character string specifying the path to te GTF file,
#'   including its name and extension. Please make sure the GTF derives from the
#'   same release of what is specified by \code{fastapath} or by
#'   \code{bsgenome}. Note that either \code{gtfpath} or \code{txdb} must be
#'   specified when the nucleotide sequences of a genome assembly are provided
#'   (see \code{fastapath} or \code{bsgenome}). Default is NULL.
#' @param txdb A character string specifying the name of the annotation package
#'   for TxDb object(s) to be loaded. If it is not already present in your
#'   system, it will be installed through the biocLite.R script (check the list
#'   of TxDb annotation packages available in the Bioconductor repositories at
#'   http://bioconductor.org/packages/release/BiocViews.html#___TxDb )). Please
#'   make sure the annotation package derives from the same release of what is
#'   specified by \code{fastapath} or by \code{bsgenome}. Note that either
#'   \code{gtfpath} or \code{txdb} must be specified when the nucleotide
#'   sequences of a genome assembly are provided (see \code{fastapath} or
#'   \code{bsgenome}). Default is NULL.
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
#' @param granges A logical value whether or not to return a GRangesList object.
#'   Default is FALSE, meaning that a list of data tables (the required input
#'   for the downstream analyses and graphical outputs provided by riboWaltz) is
#'   returned instead.
#' @return A list of data tables or a GRangesList object.
#' @examples
#' data(reads_list)
#' data(psite_offset)
#' data(mm81cdna)
#' 
#' reads_psite_list <- psite_info(reads_list, psite_offset)
#' @import dplyr
#' @import data.table
#' @export
psite_info_dt <- function(data, offset, fastapath = NULL, fasta_genome = TRUE,
                          bsgenome = NULL, gtfpath = NULL, txdb = NULL, 
                          dataSource = NA, organism = NA, granges = FALSE) {
                       
  if(((length(fastapath) != 0 & (fasta_genome == TRUE | fasta_genome == T)) |
      length(bsgenome) != 0) &
     length(gtfpath) == 0 & length(txdb) == 0){
    cat("\n")
    stop("\nERROR: annotation file not specified \n\n")
  }
  
  if(length(fastapath) != 0 & length(bsgenome) != 0){
    warning("fastapath and bsgenome are both specified. Only fastapath will be considered\n")
    bsgenome = NULL
  }
  
  if(length(gtfpath) != 0 & length(txdb) != 0){
    warning("gtfpath and txdb are both specified. Only gtfpath will be considered\n")
    txdb = NULL
  }
  
  if((length(gtfpath) != 0 | length(txdb) != 0) &
     ((length(fastapath) == 0 & length(bsgenome) == 0) |
      (length(fastapath) != 0 & (fasta_genome == FALSE | fasta_genome == F)))){
    warning("annotation file specified but no sequences from genome assembly provided\n")
  }
  
  if(length(gtfpath) != 0 | length(txdb) != 0){
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
  }
  
  if(length(fastapath) != 0 | length(bsgenome) != 0){
    if(length(fastapath) != 0) {
      if(fasta_genome == TRUE | fasta_genome == T){
        temp_sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
        exon <- GenomicFeatures::exonsBy(txdbanno, by="tx", use.names=TRUE)
        names(exon) <- paste(names(exon), c(1:length(names(exon))), sep="_._" )
        exon <-  as.data.frame(exon)
        sub_exon <- subset(exon, seqnames %in% names(temp_sequences))
        seq_df <- sub_exon %>% 
          group_by(group_name) %>%
          summarise(seq = paste(Biostrings::subseq(temp_sequences[as.character(seqnames)],
                                                   start = start,
                                                   end = end), collapse=""))
        sequences <- Biostrings::DNAStringSet(x = seq_df$seq)
        names(sequences) = unlist(lapply(strsplit(seq_df$group_name, "_._"), `[[`, 1))
      } else {
        sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
      }
    } else {
      if(bsgenome %in% installed.genomes()){
        library(bsgenome, character.only = TRUE)
      } else {
        source("http://www.bioconductor.org/biocLite.R")
        biocLite(bsgenome, suppressUpdates = TRUE)
        library(bsgenome, character.only = TRUE)
      }
    }
    sequences <- GenomicFeatures::extractTranscriptSeqs(get(bsgenome), txdbanno, use.names=T)
  }
  
  names <- names(data)
  for (n in names) {
    cat(sprintf("processing %s\n", n))
    df <- data[[n]]
    suboff <- offset[sample == n, .(length,adj_offset_from_3)]
    cat("1. adding p-site position\n")
    df[suboff,  on = 'length', psite := i.adj_offset_from_3]
    df[, psite := end3 - psite]
    setcolorder(df,c("transcript", "end5", "psite", "end3", "length", "start_pos", "stop_pos"))
    df[, psite_from_start := psite - start_pos
       ][stop_pos == 0, psite_from_start := 0]
    df[, psite_from_stop := psite - stop_pos
       ][stop_pos == 0, psite_from_stop := 0]
    cat("2. adding transcript region\n")
    df[, psite_region := "5utr"
       ][psite_from_start >= 0 & psite_from_stop <= 0, psite_region := "cds"
         ][psite_from_stop > 0, psite_region := "3utr"
           ][stop_pos == 0, psite_region := NA
             ]
    if(length(fastapath) != 0 | length(bsgenome) != 0){
      cat("3. adding nucleotide sequence\n\n")
      df[, psite_codon := as.character(subseq(sequences[as.character(df$transcript)],
                                              start = df$psite,
                                              end = df$psite + 2))]
    }
    
    if (granges == T || granges == TRUE) {
      df <- GenomicRanges::makeGRangesFromDataFrame(df,
                                                    keep.extra.columns=TRUE,
                                                    ignore.strand=TRUE,
                                                    seqnames.field=c("transcript"),
                                                    start.field="end5",
                                                    end.field="end3",
                                                    strand.field="strand",
                                                    starts.in.df.are.0based=FALSE)
      strand(df) <- "+"
    }
    
    data[[n]] <- df
  }
  
  if (granges == T || granges == TRUE) {
    data <- GenomicRanges::GRangesList(data)
  }
  
  return(data)
}
