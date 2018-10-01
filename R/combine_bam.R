#' Convert BAM files into a list of data tables or into a GRangesList object.
#'
#' Reads one or several BAM files, converts each file into a data table and
#' combines them into a list. Alternatively, it returns a GRangesList i.e. a
#' list of GRanges objects. In both cases the data structure contains for each
#' read the name of the reference sequence (i.e. of the transcript) on which it
#' aligns, the leftmost and rightmost position of the read and its length. Two
#' additional columns are attached, reporting the leftmost and rightmost
#' position of the CDS of the reference sequence with respect to its 1st
#' nuclotide. Please note: if a transcript is not associated to any annotated
#' CDS then its start and the stop codon are set to 0.
#'
#' @param bamfolder A character string indicating the path to the folder
#'   containing the BAM files.
#' @param annotation A data table from \code{\link{create_annotation}}. Please
#'   make sure that the name of the reference sequences in the annotation data
#'   table coincides with those in the BAM files.
#' @param transcript_align A logical value whether or not the BAM files within
#'   \code{bamfolder} refers to a transcriptome alignment (intended as an
#'   alignment based on a reference FASTA of all the transcript sequences). When
#'   this parameter is TRUE (the default) no reads mapping on the negative
#'   strand should be present and they are therefore removed.
#' @param list_name A character string vector specifying the desired names for
#'   the data tables of the output list. Its length must coincides with the
#'   number of BAM files within \code{bamfolder}. Please pay attention to the
#'   order in which they are provided: the first string is assigned to the first
#'   file, the second string to the second one and so on. By default this
#'   argument is NULL, implying that the data tables are named after the name of
#'   the BAM file, leaving their path and extension out.
#' @param rm_version A logical value whether ot not to remove the version of the
#'   transcripts from the end of their ID, usually separated by a dot. This
#'   option might be useful to make the transcripts IDs in the BAM files match
#'   with those in the annotation table. Default is FALSE.
#' @param granges A logical value whether or not to return a GRangesList object.
#'   Default is FALSE, meaning that a list of data tables (the required input
#'   for \code{\link{length_filter}}, \code{\link{psite}} and
#'   \code{\link{psite_info}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}}) is returned instead.
#' @return A list of data tables or a GRangesList object.
#' @examples
#' ## path_bam <- "path/to/BAM/files"
#' ## annotation_dt <- datatable_with_transcript_annotation
#' ## bamtolist(bamfolder = path_bam, annotation = annotation_dt)
#' @import data.table
#' @export
bamtolist <- function(bamfolder, annotation, transcript_align = TRUE, 
                      list_name = NULL, rm_version = FALSE,
                      granges = FALSE) {
  names <- list.files(path = bamfolder, pattern = ".bam$")
  if (length(list_name) == 0) {
    list_name <- unlist((strsplit(names, ".bam")))
  } else {
    if (length(list_name) > length(names)) {
      cat("\n")
      stop("length of list_name greater than number of files\n\n")
    }
    if (length(list_name) < length(names)) {
      cat("\n")
      stop("length of list_name smaller than number of files\n\n")
    }
  }

  sample_reads_list <- list()
  i <- 0
  for (n in names) {
    i <- i + 1
    cat(sprintf("reading %s\n", n))
    sampname <- list_name[i]
    filename <- paste(bamfolder, n, sep = "/")
    dt <- as.data.table(GenomicAlignments::readGAlignments(filename))
    dt <- dt[,.(seqnames, start, end, width, strand)]
    names(dt) <- c("transcript", "end5", "end3", "length", "strand")
    
    if(rm_version == TRUE | rm_version == T){
      dt <- dt[, transcript := tstrsplit(transcript, ".", fixed = TRUE, keep = 1)]
    }

    nreads <- nrow(dt)
    cat(sprintf("reads: %s M\n", format(round((nreads / 1000000), 2), nsmall = 2)))
    dt <- dt[as.character(transcript) %in% as.character(annotation$transcript)]
    if(nreads != nrow(dt)){
      if(nrow(dt) == 0){
        stop("%s M  (%s %%) reads removed: referece transcript ID not found in the annotation table\n\n")
      } else{
        cat(sprintf("%s M  (%s %%) reads removed: referece transcript ID not found in the annotation table\n", 
                    format(round((nreads - nrow(dt)) / 1000000, 2), nsmall = 2), 
                    format(round(((nreads - nrow(dt)) / nreads) * 100, 2), nsmall = 2) )) 
      }
    }
    
    if(transcript_align == TRUE | transcript_align == T){
      dt <- dt[strand == "+"]
      cat(sprintf("%s M  (%s %%) reads removed: mapping on the negative strand\n", 
                  format(round((nreads - nrow(dt)) / 1000000, 2), nsmall = 2), 
                  format(round(((nreads - nrow(dt)) / nreads) * 100, 2), nsmall = 2) ))
    }
    
    dt[annotation, on = 'transcript', c("start_pos", "stop_pos") := list(i.l_utr5 + 1, i.l_utr5 + i.l_cds)]
    dt[start_pos == 1 & stop_pos == 0, start_pos := 0]
    dt[, strand := NULL]

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
    
    sample_reads_list[[sampname]] <- dt
  }
  
  if (granges == T || granges == TRUE) {
    sample_reads_list <- GenomicRanges::GRangesList(sample_reads_list)
  }
  
  return(sample_reads_list)
}

#' Convert BAM files into BED files.
#'
#' Converts one or several BAM files into a list of BED files containing for
#' each read the name of the reference sequence (i.e. of the transcript) on
#' which it aligns, the leftmost and rightmost position of the read, its length
#' and the associated strand. Please note: this function calls the
#' \code{\link{bamtobed}} utility of the BEDTools suite.
#'
#' @param bamfolder A character string specifying the path to the directory
#'   containing the BAM files. The function recursively looks for BAM format
#'   file starting from the specified folder.
#' @param bedfolder A character string specifying the (existing or not) location
#'   of the directory where the BED files should be stored. By default this
#'   argument is NULL, which implies the folder is set as a subdirectory of
#'   \code{bamfolder}, called \emph{bed}.
#' @examples
#' ## path_bam <- "location_of_BAM_files"
#' ## path_bed <- "location_of_output_directory"
#' ## bamtobed(bamfolder = path_bam, bedfolder = path_bed)
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

#' Convert BED files into a list of data tables or a GRangesList.
#'
#' Reads one or several BED files, converts each file into a data table and
#' combines them into a list. Alternatively, it returns a GRangesList i.e. a
#' list of GRanges objects. In both cases two additional columns are attached to
#' the data structures, reporting the leftmost and rightmost position of the CDS
#' of the reference sequence with respect to its 1st nuclotide. Please note: if
#' a transcript is not associated to any annotated CDS then its start and the
#' stop codon are set to 0.
#'
#' @param bedfolder A character string indicating the path to the folder
#'   containing the BED files from \code{\link{bamtobed}}.
#' @param annotation A data table from \code{\link{create_annotation}}. Please
#'   make sure that the name of the reference sequences in the annotation data
#'   table coincides with those in the BED files.
#' @param transcript_align A logical value whether or not the BED files within
#'   \code{bedfolder} refers to a transcriptome alignment (intended as an
#'   alignment based on a reference FASTA of all the transcript sequences). When
#'   this parameter is TRUE (the default) no reads mapping on the negative
#'   strand should be present and they are therefore removed.
#' @param list_name A character string vector specifying the desired names for
#'   the data tables of the output list. Its length must coincides with the
#'   number of BED files within \code{bedfolder}. Please pay attention to the
#'   order in which they are provided: the first string is assigned to the first
#'   file, the second string to the second one and so on. By default this
#'   argument is NULL, implying that the data tables are named after the name of
#'   the BED file, leaving their path and extension out.
#' @param rm_version A logical value whether ot not to remove the version of the
#'   transcripts from the end of their ID, usually separated by a dot. This
#'   option might be useful to make the transcripts IDs in the BED files match
#'   with those in the annotation table. Default is FALSE.
#' @param granges A logical value whether or not to return a GRangesList object.
#'   Default is FALSE, meaning that a list of data tables (the required input
#'   for \code{\link{length_filter}}, \code{\link{psite}} and
#'   \code{\link{psite_info}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}}) is returned instead.
#' @return A list of data tables or a GRangesList object.
#' @examples
#' ## path_bed <- "path/to/BED/files"
#' ## annotation_dt <- datatable_with_transcript_annotation
#' ## bedtolist(bedfolder = path_bed, annotation = annotation_dt)
#' @import data.table
#' @export
bedtolist <- function(bedfolder, annotation, transcript_align = TRUE,
                      list_name = NULL, rm_version = FALSE, 
                      granges = FALSE) {
  names <- list.files(path = bedfolder, pattern = ".bed$")
  if (length(list_name) == 0) {
    list_name <- unlist((strsplit(names, ".bed")))
  } else {
    if (length(list_name) > length(names)) {
      cat("\n")
      stop("length of list_name greater than number of files\n\n")
    }
    if (length(list_name) < length(names)) {
      cat("\n")
      stop("length of list_name smaller than number of files\n\n")
    }
  }

  sample_reads_list <- list()
  i <- 0
  for (n in names) {
    i <- i + 1
    cat(sprintf("reading %s\n", n))
    sampname <- list_name[i]
    filename <- paste(bedfolder, n, sep = "/")
    dt <- fread(filename, sep="\t", header = FALSE)
    names(dt) <- c("transcript", "end5", "end3", "length", "strand")
    
    if(rm_version == TRUE | rm_version == T){
      dt <- dt[, transcript := tstrsplit(transcript, ".", fixed = TRUE, keep = 1)]
    }
    
    nreads <- nrow(dt)
    cat(sprintf("reads: %s M\n", format(round((nreads / 1000000), 2), nsmall = 2)))
    dt <- dt[as.character(transcript) %in% as.character(annotation$transcript)]
    if(nreads != nrow(dt)){
      if(nrow(dt) == 0){
        stop("%s M  (%s %%) reads removed: referece transcript ID not found in the annotation table\n\n")
      } else{
        cat(sprintf("%s M  (%s %%) reads removed: referece transcript ID not found in the annotation table\n", 
                    format(round((nreads - nrow(dt)) / 1000000, 2), nsmall = 2), 
                    format(round(((nreads - nrow(dt)) / nreads) * 100, 2), nsmall = 2) )) 
      }
    }
    
    if(transcript_align == TRUE | transcript_align == T){
      dt <- dt[strand == "+"]
      cat(sprintf("%s M  (%s %%) reads removed: mapping on the negative strand\n", 
                  format(round((nreads - nrow(dt)) / 1000000, 2), nsmall = 2), 
                  format(round(((nreads - nrow(dt)) / nreads) * 100, 2), nsmall = 2) ))
    }
    
    dt[annotation, on = 'transcript', c("start_pos", "stop_pos") := list(i.l_utr5 + 1, i.l_utr5 + i.l_cds)]
    dt[start_pos == 1 & stop_pos == 0, start_pos := 0]
    dt[, strand := NULL]
    
    if(granges == T || granges == TRUE) {
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
    
    sample_reads_list[[sampname]] <- dt
  }
  
  if(granges == T || granges == TRUE) {
    sample_reads_list <- GenomicRanges::GRangesList(sample_reads_list)
  }
  
  return(sample_reads_list)
}
