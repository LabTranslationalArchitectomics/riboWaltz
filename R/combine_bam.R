#' Convert BAM files into BED files.
#'
#' Converts a list of BAM files into a list of BED files containing, for each
#' read, the name of the transcript on which it alignes, the position of its
#' fist and last nucleotide, its length and the associated strand.
#'
#' @param bamfolder A character string specifying the path to the BAM files.
#'   The function looks for BAM format file recursively starting from the
#'   specified folder.
#' @param bedfolder A character string specifying the (existing or not) location
#'   of the output directory i.e. where the BED files shuold be saved. By
#'   default this argument is NULL, which implies the folder is set as a
#'   subfolder of \code{bamfolder}, called \emph{bed}.
#' @examples
#' path_bam <- location_of_BAM_files
#' path_bed <- location_of_output_directory
#' bamtobed(bamfolder = path_bam, bedfolder = path_bed)
#' @export
bamtobed <- function(bamfolder, bedfolder = NULL) {
  if (length(bedfolder) == 0) {
    bedfolder <- paste(bamfolder, "/bed", sep = "")
  }
  if (!dir.exists(bedfolder)) {
    dir.create(bedfolder)
  }
  bamtobed <- paste("bam_folder=\"", bamfolder, "\" && out_folder=\"", bedfolder, "\" && bam_files=$(find $bam_folder -type \"f\" -name \"*.bam\" | sort) && for name in $bam_files; do outname=$(echo $name | rev | cut -d'/' -f 1 | cut -d '.' -f 2- | rev); printf \"from \\t\\t $name \\t\\t to \\t\\t $out_folder/$outname.bed\\n\"; bedtools bamtobed -i $name | cut -f1,2,3,6 | awk -F'\\t' '{OFS = \"\\t\"; $5=$4; $2=$2+1; $4=$3-$2+1; print}' > $out_folder/$outname.bed; done",
                    sep = "")
  system(bamtobed)
}

#' Convert a list of BED files into a list of data frames.
#'
#' Reads a list of BED files, converts each file into a data frame and combines
#' them into a list, so that each data frame in the list correspond to a sample.
#' Two additional columns are attached to the data frames, containing the
#' position of the start and the stop codons with respect to the beginning of
#' the transcript. Since the original BAM files come from an alingment on
#' transcripts, the reads associated to the negative strand should be present in
#' a low percentage, and they will be removed
#'
#' @param bedfolder A character string indicating the folder where the BED files
#'   are located.
#' @param annotation A data frame with a reference annotation of the transripts.
#'   It must contain at least five columns named \emph{transcript}, 
#'   \emph{transcript_type}, \emph{l_utr5}, \emph{l_cds} and \emph{l_utr3} 
#'   containing the name of the transcripts (the same as in the reference 
#'   transcriptome), its transcript type, the position of the first nucleotide 
#'   of the \emph{5' UTR}, the \emph{CDS} and the  \emph{3' UTR}, respectively. 
#'   No specific order is required.
#' @param list_name A character vector specifying the desired names to be 
#'   assigned to the data frames of the output list. Its length must be equal to
#'   the number of BED files within \code{bedfolder}. Pay attention to the order
#'   in which their are provided: the first string will be assigned to the first
#'   file, the second string to the second one and so on. By default this 
#'   argument is NULL, which implies the data frames will be named as the BED 
#'   files, leaving their path and extension out.
#' @return A list of data frames.
#' @examples
#' path_bed <- location_of_BED_files
#' annotation_df <- dataframe_with_transcript_annotation
#' bedtolist(bamfolder = path_bed, annotation = annotation_df)
#' @export
bedtolist <- function(bedfolder, annotation, list_name = NULL) {
  names <- list.files(path = bedfolder, pattern = ".bed")
  if (length(list_name) == 0) {
    list_name <- unlist((strsplit(names, ".bed")))
  } else {
    if (length(list_name) > length(names)) {
      cat("\n")
      stop("\nERROR: length of list_name greater than number of files\n\n")
    }
    if (length(list_name) < length(names)) {
      cat("\n")
      stop("\nERROR: length of list_name smaller than number of files\n\n")
    }
  }
  rownames(annotation) <- as.character(annotation$transcript)
  sample.reads.list <- list()
  i <- 0
  for (n in names) {
    i <- i + 1
    cat(sprintf("reading %s\n", n))
    sampname <- list_name[i]
    filename <- paste(bedfolder, n, sep = "/")
    df <- read.table(filename, header = FALSE, sep = "\t")
    colnames(df) <- c("transcript", "end5", "end3", "length", "strand")
    nreads <- nrow(df)
    cat(sprintf("reads: %f M\n", (nreads / 1e+06)))
    df <- subset(df, strand == "+")
    cat(sprintf("positive strand: %s %%\n", format(round((nrow(df)/nreads)*100, 2), nsmall = 2) ))
    cat(sprintf("negative strand: %s %%\n\n", format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
    df$start_pos <- annotation[as.character(df$transcript), "l_utr5"] + 1
    df$stop_pos <- annotation[as.character(df$transcript), "l_utr5"] + annotation[as.character(df$transcript), "l_cds"]
    df$start_pos <- ifelse(df$start_pos == 1 & df$stop_pos == 0, 0, df$start_pos)
    sample.reads.list[[sampname]] <- df[ , !(names(df) %in% "strand")]
  }
  return(sample.reads.list)
}
