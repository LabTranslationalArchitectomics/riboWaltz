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
#' position of the start and the stop codon of the CDS with respect to the 1st
#' nuclotide of the transcript. Please note: if a transcript doesn't present any
#' annotated CDS then the positions of both the start and the stop codon will be
#' set to 0. Since the original BAM files come from a transcriptome alingment,
#' the reads aligned to the negative strand should be present in a low
#' percentage, and they will be removed.
#'
#' @param bedfolder A character string indicating the path to the folder where
#'   the BED files are located.
#' @param annotation A data frame from \code{\link{create_annotation}}. Please
#'   note that the transcript names in the annotation data frame should coincide
#'   with those in the alignment file.
#' @param list_name A character vector specifying the desired names to be 
#'   assigned to the data frames in the output list. Its length must be equal to
#'   the number of BED files within \code{bedfolder}. Pay attention to the order
#'   their are provided: the first string will be assigned to the first
#'   file, the second string to the second one and so on. By default this 
#'   argument is NULL, which implies the data frames will be named after the BED 
#'   paths, leaving their path and extension out.
#' @return A list of data frames.
#' @examples
#' path_bed <- location_of_BED_files
#' annotation_df <- dataframe_with_transcript_annotation
#' bedtolist(bedfolder = path_bed, annotation = annotation_df)
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
  sample_reads_list <- list()
  i <- 0
  for (n in names) {
    i <- i + 1
    cat(sprintf("reading %s\n", n))
    sampname <- list_name[i]
    filename <- paste(bedfolder, n, sep = "/")
    df <- read.table(filename, header = FALSE, sep = "\t")
    colnames(df) <- c("transcript", "end5", "end3", "length", "strand")
    nreads <- nrow(df)
    cat(sprintf("reads (total): %f M\n", (nreads / 1e+06)))
    df <- subset(df, as.character(transcript) %in% rownames(annotation))
    cat(sprintf("reads (kept): %f M\n", (nreads / 1e+06)))
    df <- subset(df, strand == "+")
    cat(sprintf("positive strand: %s %%\n", format(round((nrow(df)/nreads)*100, 2), nsmall = 2) ))
    cat(sprintf("negative strand: %s %%\n\n", format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
    df$start_pos <- annotation[as.character(df$transcript), "l_utr5"] + 1
    df$stop_pos <- annotation[as.character(df$transcript), "l_utr5"] + annotation[as.character(df$transcript), "l_cds"]
    df$start_pos <- ifelse(df$start_pos == 1 & df$stop_pos == 0, 0, df$start_pos)
    sample_reads_list[[sampname]] <- df[ , !(names(df) %in% "strand")]
  }
  return(sample_reads_list)
}
