#' Convert BAM files into a list of data frames or into a GRangesList object.
#'
#' Reads one or several BAM files, converts each file into a data frame and
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
#' @param annotation A data frame from \code{\link{create_annotation}}. Please
#'   make sure that the name of the reference sequences in the annotation data
#'   frame coincides with those in the BAM files.
#' @param transcript_align A locigal value whether or not the BAM files within
#'   \code{bamfolder} refers to a transcriptome alignment (intended as an
#'   alignment based on a reference FASTA of all the transcript sequences). When
#'   this parameter is TRUE (the default) no reads mapping on the negative
#'   strand should be present and they are therefore removed.
#' @param filter Either "none" (the default), "custom" or "periodicity". It
#'   specifies how to handle the selection of the read length. "none": all read
#'   lengths are included in the analysis; "custom": only read lengths specified
#'   by the user are included (see \code{custom_range}); "periodicity": only
#'   read lengths satisfying a periodicity threshold (see \code{periodicity_th})
#'   are included in the analysis. This mode enables the removal of all the
#'   reads without periodicity.
#' @param custom_range An integer or an integer vector specifying either the
#'   read length or multiple read lengths to be kept, respectively. This
#'   parameter is considered only when \code{filter} is set to "custom".
#' @param periodicity_th An integer in \emph{[10, 100]}. Only the read lengths
#'   satisfying this threshold (i.e. with a higher percentage of read
#'   extremities falling in the same frame) are kept. This parameter is
#'   considered only when \code{filter} is set to "periodicity". Default is 50.
#' @param list_name A character string vector specifying the desired names for
#'   the data frames of the output list. Its length must coincides with the
#'   number of BAM files within \code{bamfolder}. Please pay attention to the
#'   order in which they are provided: the first string is assigned to the first
#'   file, the second string to the second one and so on. By default this
#'   argument is NULL, implying that the data frames are named after the name of
#'   the BAM file, leaving their path and extension out.
#' @param granges A logical value whether or not to return a GRangesList object.
#'   Default is FALSE, meaning that a list of data frames (the required input
#'   for \code{\link{psite}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}}) is returned instead.
#' @return A list of data frames or a GRangesList object.
#' @examples
#' path_bam <- "location_of_BAM_files"
#' annotation_df <- dataframe_with_transcript_annotation
#' bamtolist(bamfolder = path_bam, annotation = annotation_df)
#' @import dplyr
#' @export
bamtolist <- function(bamfolder, annotation, transcript_align = TRUE,
                      filter = "none", custom_range = NULL, periodicity_th = 50,
                      list_name = NULL, granges = FALSE) {
  names <- list.files(path = bamfolder, pattern = ".bam")
  if (length(list_name) == 0) {
    list_name <- unlist((strsplit(names, ".bam")))
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
  
  if(filter == "custom" & !inherits(custom_range, "numeric") & !inherits(custom_range, "integer")){
    stop("'custom_range' must be integer.\n\n")
  }
  
  rownames(annotation) <- as.character(annotation$transcript)
  sample_reads_list <- list()
  i <- 0
  for (n in names) {
    i <- i + 1
    cat(sprintf("reading %s\n", n))
    sampname <- list_name[i]
    filename <- paste(bamfolder, n, sep = "/")
    df <- as.data.frame(GenomicAlignments::readGAlignments(filename))
    df <- df[, c("seqnames", "start", "end", "width", "strand")]
    colnames(df) <- c("transcript", "end5", "end3", "length", "strand")
    
    if(transcript_align == TRUE | transcript_align == T){
      nreads <- nrow(df)
      cat(sprintf("reads (total): %f M\n", (nreads / 1e+06)))
      df <- subset(df, as.character(transcript) %in% rownames(annotation))
      df <- subset(df, strand == "+")
      cat(sprintf("positive strand: %s %%\n", 
                  format(round((nrow(df)/nreads)*100, 2), nsmall = 2) ))
      cat(sprintf("negative strand: %s %%\n", 
                  format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
      cat(sprintf("reads (kept): %f M\n\n", (nrow(df) / 1e+06)))
    }
    
    df$start_pos <- annotation[as.character(df$transcript), "l_utr5"] + 1
    df$stop_pos <- annotation[as.character(df$transcript), "l_utr5"] + 
      annotation[as.character(df$transcript), "l_cds"]
    df$start_pos <- ifelse(df$start_pos == 1 & df$stop_pos == 0, 0, df$start_pos)
    
    if (identical(filter, "custom")) {
      df <- subset(df, as.numeric(as.character(length)) %in% custom_range)
      cat(sprintf("%s (%s %%) reads have been removed\n\n", 
                  format(nreads - nrow(df), nsmall = 2), 
                  format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
    } else {
      if(identical(filter, "periodicity")){
        subdf5 <- subset(df, start_pos!=0 & end5 - start_pos >=0 & stop_pos - end5 >=0)
        t_temp5 <- subdf5 %>% 
          mutate(end5_frame = (end5 - start_pos) %% 3) %>%
          mutate(end5_frame = factor(end5_frame, levels = c("0", "1", "2")))
        t_end5 <- t_temp5 %>% 
          group_by(length, end5_frame) %>%
          summarise(end5_count = n()) %>%
          mutate(end5_perc = (end5_count / sum (end5_count)) * 100) %>%
          data.frame
        keep_length5 <- unique(subset(t_end5, end5_perc >= periodicity_th)$length)
        subdf3 <- subset(df, start_pos!=0 & end3 - start_pos >=0 & stop_pos - end3 >=0)
        t_temp3 <- subdf3 %>% 
          mutate(end3_frame = (end3 - start_pos) %% 3) %>%
          mutate(end3_frame = factor(end3_frame, levels = c("0", "1", "2")))
        t_end3 <- t_temp3 %>% 
          group_by(length, end3_frame) %>%
          summarise(end3_count = n()) %>%
          mutate(end3_perc = (end3_count / sum (end3_count)) * 100) %>%
          data.frame
        keep_length3 <- unique(subset(t_end3, end3_perc >= periodicity_th)$length)
        keep_length <- intersect(keep_length5, keep_length3)
        df <- subset(df, as.numeric(as.character(length)) %in% keep_length)
        cat(sprintf("%s (%s %%) reads have been removed\n\n", 
                    format(nreads - nrow(df), nsmall = 2), 
                    format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
      }
    }
    
    df <- df[, !(names(df) %in% "strand")]
    
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
    
    sample_reads_list[[sampname]] <- df
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
#' path_bam <- "location_of_BAM_files"
#' path_bed <- "location_of_output_directory"
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

#' Convert BED files into a list of data frames or a GRangesList.
#'
#' Reads one or several BED files, converts each file into a data frame and
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
#' @param annotation A data frame from \code{\link{create_annotation}}. Please
#'   make sure that the name of the reference sequences in the annotation data
#'   frame coincides with those in the BED files.
#' @param transcript_align A locigal value whether or not the BED files within
#'   \code{bedfolder} refers to a transcriptome alignment (intended as an
#'   alignment based on a reference FASTA of all the transcript sequences). When
#'   this parameter is TRUE (the default) no reads mapping on the negative
#'   strand should be present and they are therefore removed.
#' @param filter Either "none" (the default), "custom" or "periodicity". It
#'   specifies how to handle the selection of the read length. "none": all read
#'   lengths are included in the analysis; "custom": only read lengths specified
#'   by the user are included (see \code{custom_range}); "periodicity": only
#'   read lengths satisfying a periodicity threshold (see \code{periodicity_th})
#'   are included in the analysis. This mode enables the removal of all the
#'   reads without periodicity.
#' @param custom_range An integer or an integer vector specifying either the
#'   read length or multiple read lengths to be kept, respectively. This
#'   parameter is considered only when \code{filter} is set to "custom".
#' @param periodicity_th An integer in \emph{[10, 100]}. Only the read lengths
#'   satisfying this threshold (i.e. with a higher percentage of read
#'   extremities falling in the same frame) are kept. This parameter is
#'   considered only when \code{filter} is set to "periodicity". Default is 50.
#' @param list_name A character string vector specifying the desired names for
#'   the data frames of the output list. Its length must coincides with the
#'   number of BED files within \code{bedfolder}. Please pay attention to the
#'   order in which they are provided: the first string is assigned to the first
#'   file, the second string to the second one and so on. By default this
#'   argument is NULL, implying that the data frames are named after the name of
#'   the BED file, leaving their path and extension out.
#' @param granges A logical value whether or not to return a GRangesList object.
#'   Default is FALSE, meaning that a list of data frames (the required input
#'   for \code{\link{psite}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}}) is returned instead.
#' @return A list of data frames or a GRangesList object.
#' @examples
#' path_bed <- "location_of_BED_files"
#' annotation_df <- dataframe_with_transcript_annotation
#' bedtolist(bedfolder = path_bed, annotation = annotation_df)
#' @import dplyr
#' @export
bedtolist <- function(bedfolder, annotation, transcript_align = TRUE,
                      filter = "none", custom_range = NULL, periodicity_th = 50,
                      list_name = NULL, granges = FALSE) {
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
  
  if(filter == "custom" & !inherits(custom_range, "numeric") & !inherits(custom_range, "integer")){
    stop("'custom_range' must be an integer.\n\n")
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

    if(transcript_align == TRUE | transcript_align == T){
      nreads <- nrow(df)
      cat(sprintf("reads (total): %f M\n", (nreads / 1e+06)))
      df <- subset(df, as.character(transcript) %in% rownames(annotation))
      df <- subset(df, strand == "+")
      cat(sprintf("positive strand: %s %%\n", 
                  format(round((nrow(df)/nreads)*100, 2), nsmall = 2) ))
      cat(sprintf("negative strand: %s %%\n", 
                  format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
      cat(sprintf("reads (kept): %f M\n\n", (nrow(df) / 1e+06)))
    }
    
    df$start_pos <- annotation[as.character(df$transcript), "l_utr5"] + 1
    df$stop_pos <- annotation[as.character(df$transcript), "l_utr5"] + 
      annotation[as.character(df$transcript), "l_cds"]
    df$start_pos <- ifelse(df$start_pos == 1 & df$stop_pos == 0, 0, df$start_pos)
    
    if (identical(filter, "custom")) {
      df <- subset(df, as.numeric(as.character(length)) %in% custom_range)
      cat(sprintf("%s (%s %%) reads have been removed\n\n", 
                  format(nreads - nrow(df), nsmall = 2), 
                  format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
    } else {
      if(identical(filter, "periodicity")){
        subdf5 <- subset(df, start_pos!=0 & end5 - start_pos >=0 & stop_pos - end5 >=0)
        t_temp5 <- subdf5 %>% 
          mutate(end5_frame = (end5 - start_pos) %% 3) %>%
          mutate(end5_frame = factor(end5_frame, levels = c("0", "1", "2")))
        t_end5 <- t_temp5 %>% 
          group_by(length, end5_frame) %>%
          summarise(end5_count = n()) %>%
          mutate(end5_perc = (end5_count / sum (end5_count)) * 100) %>%
          data.frame
        keep_length5 <- unique(subset(t_end5, end5_perc >= periodicity_th)$length)
        subdf3 <- subset(df, start_pos!=0 & end3 - start_pos >=0 & stop_pos - end3 >=0)
        t_temp3 <- subdf3 %>% 
          mutate(end3_frame = (end3 - start_pos) %% 3) %>%
          mutate(end3_frame = factor(end3_frame, levels = c("0", "1", "2")))
        t_end3 <- t_temp3 %>% 
          group_by(length, end3_frame) %>%
          summarise(end3_count = n()) %>%
          mutate(end3_perc = (end3_count / sum (end3_count)) * 100) %>%
          data.frame
        keep_length3 <- unique(subset(t_end3, end3_perc >= periodicity_th)$length)
        keep_length <- intersect(keep_length5, keep_length3)
        df <- subset(df, as.numeric(as.character(length)) %in% keep_length)
        cat(sprintf("%s (%s %%) reads have been removed\n\n", 
                    format(nreads - nrow(df), nsmall = 2), 
                    format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
      }
    }
    
    df <- df[, !(names(df) %in% "strand")]
    
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
    
    sample_reads_list[[sampname]] <- df
  }
  
  if (granges == T || granges == TRUE) {
    sample_reads_list <- GenomicRanges::GRangesList(sample_reads_list)
  }
  
  return(sample_reads_list)
}
