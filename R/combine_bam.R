#' Convert BAM files into a list of data tables or into a GRangesList object.
#'
#' Reads one or several BAM files, converts each file into a data table and
#' combines them into a list. Alternatively, it returns a GRangesList i.e. a
#' list of GRanges objects. In both cases the data structure contains for each
#' read the name of the reference sequence (i.e. of the transcript) on which it
#' aligns, the leftmost and rightmost position of the read and its length. Two
#' additional columns are also attached, reporting the leftmost and rightmost
#' position of the CDS of the reference sequence with respect to its 1st
#' nuclotide. Please note: if a transcript is not associated to any annotated
#' CDS then its start and the stop codon are set to 0. Multiple options for
#' treating the read lengths are available.
#'
#' @param bamfolder A character string indicating the path to the folder
#'   containing the BAM files.
#' @param annotation A data table from \code{\link{create_annotation}}. Please
#'   make sure that the name of the reference sequences in the annotation data
#'   table coincides with those in the BAM files.
#' @param transcript_align A locigal value whether or not the BAM files within
#'   \code{bamfolder} refers to a transcriptome alignment (intended as an
#'   alignment based on a reference FASTA of all the transcript sequences). When
#'   this parameter is TRUE (the default) no reads mapping on the negative
#'   strand should be present and they are therefore removed.
#' @param length_filter_mode Either "none" (the default), "custom" or
#'   "periodicity". It specifies how to handle the selection of the read length.
#'   "none": all read lengths are included in the analysis; "custom": only read
#'   lengths specified by the user are included (see
#'   \code{length_filter_vector}); "periodicity": only read lengths satisfying a
#'   periodicity threshold (see \code{periodicity_threshold}) are included in
#'   the analysis. This mode enables the removal of all the reads without
#'   periodicity.
#' @param length_filter_vector An integer or an integer vector specifying either
#'   the read length or multiple read lengths to be kept, respectively. This
#'   parameter is considered only when \code{length_filter_mode} is set to
#'   "custom".
#' @param periodicity_threshold An integer in \emph{[10, 100]}. Only the read
#'   lengths satisfying this threshold (i.e. with a higher percentage of read
#'   extremities falling in the same frame) are kept. This parameter is
#'   considered only when \code{length_filter_mode} is set to "periodicity".
#'   Default is 50.
#' @param list_name A character string vector specifying the desired names for
#'   the data tables of the output list. Its length must coincides with the
#'   number of BAM files within \code{bamfolder}. Please pay attention to the
#'   order in which they are provided: the first string is assigned to the first
#'   file, the second string to the second one and so on. By default this
#'   argument is NULL, implying that the data tables are named after the name of
#'   the BAM file, leaving their path and extension out.
#' @param granges A logical value whether or not to return a GRangesList object.
#'   Default is FALSE, meaning that a list of data tables (the required input
#'   for \code{\link{psite}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}}) is returned instead.
#' @return A list of data tables or a GRangesList object.
#' @examples
#' ## path_bam <- "location_of_BAM_files"
#' ## annotation_dt <- datatable_with_transcript_annotation
#' ## bamtolist(bamfolder = path_bam, annotation = annotation_dt)
#' @import data.table
#' @export
bamtolist <- function(bamfolder, annotation, transcript_align = TRUE,
                      length_filter_mode = "none", length_filter_vector = NULL,
                      periodicity_threshold = 50, list_name = NULL,
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
  
  if(length_filter_mode == "custom" & !inherits(length_filter_vector, "numeric") & !inherits(length_filter_vector, "integer")){
    stop("length_filter_vector must be integer\n\n")
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

    if(transcript_align == TRUE | transcript_align == T){
      nreads <- nrow(dt)
      cat(sprintf("reads (total): %s M\n", format(round((nreads / 1e+06), 2), nsmall = 2)))
      dt <- dt[as.character(transcript) %in% as.character(annotation$transcript) & strand == "+"]
      cat(sprintf("%s M  (%s %%) reads removed: mapping on the negative strand\n", 
                  format(round((nreads - nrow(dt))/ 1e+06, 2), nsmall = 2), 
                  format(round(((nreads - nrow(dt))/nreads) * 100, 2), nsmall = 2) ))
    }
    
    dt[annotation, on = 'transcript', c("start_pos", "stop_pos") := list(i.l_utr5 + 1, i.l_utr5 + i.l_cds)]
    dt[start_pos == 1 & stop_pos == 0, start_pos := 0]
    
    if (identical(length_filter_mode, "custom")) {
      nreads <- nrow(dt)
      dt <- dt[length %in% length_filter_vector]
      cat(sprintf("%s M  (%s %%) reads removed: length_filter_mode applied\n", 
                  format(round((nreads - nrow(dt))/ 1e+06, 2), nsmall = 2), 
                  format(round(((nreads - nrow(dt))/nreads) * 100, 2), nsmall = 2) ))
      cat(sprintf("reads (kept): %s M\n\n", format(round((nrow(dt) / 1e+06), 2), nsmall = 2)))
    } else {
      if(identical(length_filter_mode, "periodicity")){
        nreads <- nrow(dt)
        
        subdt5 <- dt[start_pos != 0 &
                       (end5 - start_pos) >= 0 &
                       (stop_pos - end5) >= 0]
        subdt5[, end5_frame := as.factor((end5 - start_pos) %% 3)]
        t_end5 <- subdt5[, .N, by = list(length, end5_frame)
                         ][, end5_perc := (N / sum(N)) * 100, , by = length]
        keep_length5 <- unique(t_end5[end5_perc >= periodicity_threshold, length])
        
        subdt3 <- dt[start_pos != 0 &
                       (end3 - start_pos) >= 0 &
                       (stop_pos - end3) >= 0]
        subdt3[, end3_frame := as.factor((end3 - start_pos) %% 3)]
        t_end3 <- subdt3[, .N, by = list(length, end3_frame)
                         ][, end3_perc := (N / sum(N)) * 100, , by = length]
        keep_length3 <- unique(t_end3[end3_perc >= periodicity_threshold, length])
        
        keep_length <- intersect(keep_length5, keep_length3)
        dt <- dt[length %in% keep_length]
        
        cat(sprintf("%s M  (%s %%) reads removed: length_filter_mode applied\n", 
                    format(round((nreads - nrow(dt))/ 1e+06, 2), nsmall = 2), 
                    format(round(((nreads - nrow(dt))/nreads) * 100, 2), nsmall = 2) ))
        cat(sprintf("reads (kept): %s M\n\n", format(round((nrow(dt) / 1e+06), 2), nsmall = 2)))
      }
    }
    
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
#' and the associated strand. Please note that thus function calls the bamtobed
#' utility of the BEDTools suite.
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
#' stop codon are set to 0. Multiple options for treating the read lengths are
#' available.
#'
#' @param bedfolder A character string indicating the path to the folder
#'   containing the BED files.
#' @param annotation A data table from \code{\link{create_annotation}}. Please
#'   make sure that the name of the reference sequences in the annotation data
#'   table coincides with those in the BED files.
#' @param transcript_align A locigal value whether or not the BED files within
#'   \code{bedfolder} refers to a transcriptome alignment (intended as an
#'   alignment based on a reference FASTA of all the transcript sequences). When
#'   this parameter is TRUE (the default) no reads mapping on the negative
#'   strand should be present and they are therefore removed.
#' @param length_filter_mode Either "none" (the default), "custom" or
#'   "periodicity". It specifies how to handle the selection of the read length.
#'   "none": all read lengths are included in the analysis; "custom": only read
#'   lengths specified by the user are included (see
#'   \code{length_filter_vector}); "periodicity": only read lengths satisfying a
#'   periodicity threshold (see \code{periodicity_threshold}) are included in
#'   the analysis. This mode enables the removal of all the reads without
#'   periodicity.
#' @param length_filter_vector An integer or an integer vector specifying either
#'   the read length or multiple read lengths to be kept, respectively. This
#'   parameter is considered only when \code{length_filter_mode} is set to
#'   "custom".
#' @param periodicity_threshold An integer in \emph{[10, 100]}. Only the read
#'   lengths satisfying this threshold (i.e. with a higher percentage of read
#'   extremities falling in the same frame) are kept. This parameter is
#'   considered only when \code{length_filter_mode} is set to "periodicity".
#'   Default is 50.
#' @param list_name A character string vector specifying the desired names for
#'   the data tables of the output list. Its length must coincides with the
#'   number of BED files within \code{bedfolder}. Please pay attention to the
#'   order in which they are provided: the first string is assigned to the first
#'   file, the second string to the second one and so on. By default this
#'   argument is NULL, implying that the data tables are named after the name of
#'   the BED file, leaving their path and extension out.
#' @param granges A logical value whether or not to return a GRangesList object.
#'   Default is FALSE, meaning that a list of data tables (the required input
#'   for \code{\link{psite}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}}) is returned instead.
#' @return A list of data tables or a GRangesList object.
#' @examples
#' ## path_bed <- "location_of_BED_files"
#' ## annotation_dt <- datatable_with_transcript_annotation
#' ## bedtolist(bedfolder = path_bed, annotation = annotation_dt)
#' @import data.table
#' @export
bedtolist <- function(bedfolder, annotation, transcript_align = TRUE,
                      length_filter_mode = "none", length_filter_vector = NULL,
                      periodicity_threshold = 50, list_name = NULL,
                      granges = FALSE) {
  names <- list.files(path = bedfolder, pattern = ".bed")
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
  
  if(length_filter_mode == "custom" & !inherits(length_filter_vector, "numeric") & !inherits(length_filter_vector, "integer")){
    stop("'length_filter_vector' must be an integer.\n\n")
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

    if(transcript_align == TRUE | transcript_align == T){
      nreads <- nrow(dt)
      cat(sprintf("reads (total): %s M\n", format(round((nreads / 1e+06), 2), nsmall = 2)))
      dt <- dt[as.character(transcript) %in% as.character(annotation$transcript) & strand == "+"]
      cat(sprintf("%s M  (%s %%) reads removed: mapping on the negative strand\n", 
                  format(round((nreads - nrow(dt))/ 1e+06, 2), nsmall = 2), 
                  format(round(((nreads - nrow(dt))/nreads) * 100, 2), nsmall = 2) ))
    }
    
    dt[annotation, on = 'transcript', c("start_pos", "stop_pos") := list(i.l_utr5 + 1, i.l_utr5 + i.l_cds)]
    dt[start_pos == 1 & stop_pos == 0, start_pos := 0]
    
    if (identical(length_filter_mode, "custom")) {
      nreads <- nrow(dt)
      dt <- dt[length %in% length_filter_vector]
      cat(sprintf("%s M  (%s %%) reads removed: length_filter_mode applied\n", 
                  format(round((nreads - nrow(dt))/ 1e+06, 2), nsmall = 2), 
                  format(round(((nreads - nrow(dt))/nreads) * 100, 2), nsmall = 2) ))
      cat(sprintf("reads (kept): %s M\n\n", format(round((nrow(dt) / 1e+06), 2), nsmall = 2)))
    } else {
      if(identical(length_filter_mode, "periodicity")){
        nreads <- nrow(dt)
        
        subdt5 <- dt[start_pos != 0 &
                       (end5 - start_pos) >= 0 &
                       (stop_pos - end5) >= 0]
        subdt5[, end5_frame := as.factor((end5 - start_pos) %% 3)]
        t_end5 <- subdt5[, .N, by = list(length, end5_frame)
                         ][, end5_perc := (N / sum(N)) * 100, , by = length]
        keep_length5 <- unique(t_end5[end5_perc >= periodicity_threshold, length])
        
        subdt3 <- dt[start_pos != 0 &
                       (end3 - start_pos) >=0 &
                       (stop_pos - end3) >=0]
        subdt3[, end3_frame := as.factor((end3 - start_pos) %% 3)]
        t_end3 <- subdt3[, .N, by = list(length, end3_frame)
                         ][, end3_perc := (N / sum(N)) * 100, , by = length]
        keep_length3 <- unique(t_end3[end3_perc >= periodicity_threshold, length])
        
        keep_length <- intersect(keep_length5, keep_length3)
        dt <- dt[length %in% keep_length]
        
        cat(sprintf("%s M  (%s %%) reads removed: length_filter_mode applied\n", 
                    format(round((nreads - nrow(dt))/ 1e+06, 2), nsmall = 2), 
                    format(round(((nreads - nrow(dt))/nreads) * 100, 2), nsmall = 2) ))
        cat(sprintf("reads (kept): %s M\n\n", format(round((nrow(dt) / 1e+06), 2), nsmall = 2)))
      }
    }
    
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
