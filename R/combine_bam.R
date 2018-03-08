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
#' percentage, and they will be removed. Additional filter based on the length
#' of the reads can be chosen by the user.
#'
#' @param bedfolder A character string indicating the path to the folder where
#'   the BED files are located.
#' @param annotation A data frame from \code{\link{create_annotation}}. Please
#'   note that the transcript names in the annotation data frame should coincide
#'   with those in the alignment file.
#' @param filter Either "none" (the default), "custom" or "periodicity". It
#'   indicates if and what filter must be applied to the data sets, i.e. which
#'   reads must be discharged. "none": all the reads are kept; "custom": the
#'   user specifies the read length or a range of read lengths to be kept
#'   (see \code{custom_range}); "periodicity": all the read lengths associated
#'   to a trinucleotide periodicity of either the 3' end or 5' end of the reads
#'   in any frame of the sequence are kept. It is possible to chose the
#'   threshold that "defines" the periodicity (see \code{periodicity_th}).
#' @param custom_range An integer or an integer vector specifying the read 
#'   length or a range of read lengths to be kept, respectively. This parameter
#'   is considered only if \code{filter} == "custom".
#' @param periodicity_th An integer in \emph{[10, 90]} used to detremine the
#'   existence of periodicity for any read length. Only the read lengths at
#'   which the percentage of read extremities falling in one frame is above the
#'   threshold are kept. This parameter is considered only if \code{filter} == 
#'   "periodicity". Default is 50.
#' @param list_name A character vector specifying the desired names to be 
#'   assigned to the data frames in the output list. Its length must be equal to
#'   the number of BED files within \code{bedfolder}. Pay attention to the order
#'   their are provided: the first string will be assigned to the first
#'   file, the second string to the second one and so on. By default this 
#'   argument is NULL, which implies the data frames will be named after the BED 
#'   paths, leaving their path and extension out.
#' @param granges A logical value whether or not to return a GRangesList 
#'   object. Default is FALSE, meaning that a list of data frames (the required
#'   input for \code{\link{psite}}, \code{\link{rends_heat}} and
#'   \code{\link{rlength_distr}) will be returned.
#' @return A list of data frames or a GRangesList object if \code{granges}
#'   == TRUE.
#' @examples
#' path_bed <- location_of_BED_files
#' annotation_df <- dataframe_with_transcript_annotation
#' bedtolist(bedfolder = path_bed, annotation = annotation_df)
#' @import dplyr
#' @export
bedtolist <- function(bedfolder, annotation, filter = "none", custom_range = NULL, 
                      periodicity_th = 50, list_name = NULL, granges = F) {
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
    nreads <- nrow(df)
    cat(sprintf("reads (total): %f M\n", (nreads / 1e+06)))
    df <- subset(df, as.character(transcript) %in% rownames(annotation))
    cat(sprintf("reads (kept): %f M\n", (nreads / 1e+06)))
    df <- subset(df, strand == "+")
    nreads <- nrow(df)
    cat(sprintf("positive strand: %s %%\n", 
                format(round((nrow(df)/nreads)*100, 2), nsmall = 2) ))
    cat(sprintf("negative strand: %s %%\n\n", 
                format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
    
    if (identical(filter, "custom")) {
      df <- subset(df, as.numeric(as.character(length)) %in% custom_range)
      cat(sprintf("%s (%s %%) reads have been removed\n\n", 
                  format(nreads - nrow(df), nsmall = 2), 
                  format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
    } else {
      if(identical(filter, "periodicity")){
        t_temp <- df %>% 
          mutate(end5_frame = end5 %% 3, end3_frame = end3 %% 3) %>%
          mutate(end5_frame = factor(end5_frame, levels = c("0", "1", "2"))) %>% 
          mutate(end3_frame = factor(end3_frame, levels = c("0", "1", "2"))) 
        t_end5 <- t_temp %>% 
          group_by(length, end5_frame) %>%
          summarise(end5_count = n()) %>%
          mutate(end5_perc = (end5_count / sum (end5_count)) * 100) %>%
          data.frame
        t_end3 <- t_temp %>% 
          group_by(length, end3_frame) %>%
          summarise(end3_count = n()) %>%
          mutate(end3_perc = (end3_count / sum (end3_count)) * 100) %>%
          data.frame
        t_final <- cbind(t_end5[,c(1,ncol(t_end5))], "end3_perc" = t_end3[,ncol(t_end3)])
        keep_length <- subset(t_final, end5_perc >= periodicity_th | end3_perc >= periodicity_th)$length
        df <- subset(df, as.numeric(as.character(length)) %in% keep_length)
        cat(sprintf("%s (%s %%) reads have been removed\n\n", 
                    format(nreads - nrow(df), nsmall = 2), 
                    format(round(((nreads - nrow(df))/nreads)*100, 2), nsmall = 2) ))
      }
    }
    
    df$start_pos <- annotation[as.character(df$transcript), "l_utr5"] + 1
    df$stop_pos <- annotation[as.character(df$transcript), "l_utr5"] + annotation[as.character(df$transcript), "l_cds"]
    df$start_pos <- ifelse(df$start_pos == 1 & df$stop_pos == 0, 0, df$start_pos)
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
  


  
  
 
  
  
  
