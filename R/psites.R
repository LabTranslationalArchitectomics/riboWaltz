#' Ribosome P-sites position within reads.
#'
#' This function identifies the exact position of the ribosome P-site within
#' each read, determined by the localisation of its first nucleotide (see
#' \code{Details}). It returns a data table containing, for all samples and read
#' lengths: i) the percentage of reads in the whole dataset, ii) the percentage
#' of reads aligning on the start codon (if any); iii) the distance of the
#' P-site from the two extremities of the reads before and after the correction
#' step; iv) the name of the sample. Optionally, this function plots a
#' collection of read length-specific occupancy metaprofiles displaying the
#' P-site offsets computed through the process.
#'
#' @param data List of data tables from \code{\link{bamtolist}},
#'   \code{\link{bedtolist}} or \code{\link{length_filter}}.
#' @param flanking Integer value specifying for the selected reads the minimum
#'   number of nucleotides that must flank the reference codon in both
#'   directions. Default is 6.
#' @param start Logical value whether to use the translation initiation site as
#'   reference codon. Default is TRUE. If FALSE, the second to last codon is
#'   used instead.
#' @param extremity Either "5end", "3end" or "auto". It specifies if the
#'   correction step should be based on 5' extremities ("5end") or 3'
#'   extremities ("3end"). Default is "auto" i.e. the optimal extremity is
#'   automatically selected.
#' @param plot Logical value whether to plot the occupancy metaprofiles
#'   displaying the P-site offsets computed in both steps of the algorithm.
#'   Default is FALSE.
#' @param plot_dir Character string specifying the directory where read
#'   length-specific occupancy metaprofiles shuold be stored. If the specified
#'   folder doesn't exist, it is automatically created. If NULL (the default),
#'   the metaprofiles are stored in a new subfolder of the working directory,
#'   called \emph{offset_plot}. This parameter is considered only if \code{plot}
#'   is TRUE.
#' @param plot_format Either "png" (the default) or "pdf". This parameter
#'   specifies the file format storing the length-specific occupancy
#'   metaprofiles. It is considered only if \code{plot} is TRUE.
#' @param cl Integer value in [1,100] specifying a confidence level for
#'   generating occupancy metaprofiles for to a sub-range of read lengths i.e.
#'   for the cl% of read lengths associated to the highest signals. Default is
#'   99. This parameter is considered only if \code{plot} is TRUE.
#' @param log_file Logical value whether to generate a plain text file, called
#'   \emph{best_offset.txt}, that reports the extremity used for the correction
#'   step and the best offset for each sample. Default is FALSE.
#' @param log_file_dir Character string specifying the directory where the log
#'   file shuold be saved. If the specified folder doesn't exist, it is
#'   automatically created. If NULL (the default), the file is stored in the
#'   working directory. This parameter is considered only if \code{log_file} is
#'   TRUE.
#' @details The P-site offset (PO) is defined as the distance between the
#'   extremities of a read and the first nucleotide of the P-site itself. The
#'   function processes all samples separately starting from reads mapping on
#'   the reference codon (either the start codon or the second to last codon,
#'   see \code{start}) of any annotated coding sequences. Read lengths-specific
#'   POs are inferred in two steps. First, reads mapping on the reference codon
#'   are grouped according to their length, each group corresponding to a bin.
#'   Reads whose extremities are too close to the reference codon are discarded
#'   (see \code{flanking}). For each bin temporary 5' and 3' POs are defined as
#'   the distances between the first nucleotide of the reference codon and the
#'   nucleotide corresponding to the global maximum found in the profiles of the
#'   5' and the 3' end at the left and at the right of the reference codon,
#'   respectively. After the identification of the P-site for all reads aligning
#'   on the reference codon, the POs corresponding to each length are assigned
#'   to each read of the dataset. Second, the most frequent temporary POs
#'   associated to the optimal extremity (see \code{extremity}) and the
#'   predominant bins are exploited as reference values for correcting the
#'   temporary POs of smaller bins. Briefly, the correction step defines for
#'   each length bin a new PO based on the local maximum, whose distance from
#'   the reference codon is the closest to the most frequent temporary POs. For
#'   further details please refer to the \strong{riboWaltz} article (available
#'   \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006169}{here}).
#' @return A data table.
#' @examples
#' data(reads_list)
#'
#' ## Compute the P-site offset automatically selecting the optimal read
#' ## extremity for the correction step and not plotting any metaprofile:
#' psite(reads_list, flanking = 6, extremity="auto")
#'
#' ## Compute the P-site offset specifying the extremity used in the correction
#' ## step and plotting the length-specific occupancy metaprofiles for a 
#' ## sub-range of read lengths (the middle 95%). The plots will be placed in 
#' ## the current working directory:
#' psite_offset <- psite(reads_list, flanking = 6, extremity = "3end", plot = TRUE, cl = 95)
#' @import data.table
#' @import ggplot2
#' @export
psite <- function(data, flanking = 6, start = TRUE, extremity = "auto",
                  plot = FALSE, plot_dir = NULL, plot_format="png", cl = 99,
                  log_file = FALSE, log_file_dir = NULL) {
  
  if(log_file == T | log_file == TRUE){
    if(length(log_file_dir) == 0){
      log_file_dir <- getwd()
    }
    if (!dir.exists(log_file_dir)) {
      dir.create(log_file_dir)
    }
    logpath <- paste0(log_file_dir, "/best_offset.txt")
    cat("sample\texremity\toffset(nts)\n", file = logpath)
  }
  
  names <- names(data)
  offset <- NULL
  for (n in names) { 
    cat(sprintf("processing %s\n", n))
    dt <- data[[n]]
    lev <- sort(unique(dt$length))
    if(start == T | start == TRUE){
      base <- 0
      dt[, site_dist_end5 := end5 - cds_start]
      dt[, site_dist_end3 := end3 - cds_start]
    } else {
      base <- -5
      dt[, site_dist_end5 := end5 - cds_stop - base]
      dt[, site_dist_end3 := end3 - cds_stop - base]
    }
    site_sub <- dt[site_dist_end5 <= -flanking & site_dist_end3 >= flanking - 1]
    minlen <- min(site_sub$length)
    maxlen <- max(site_sub$length)
    t <- table(factor(site_sub$length, levels = lev))
    
    # offset
    offset_temp <- data.table(length = as.numeric(as.character(names(t))), percentage = (as.vector(t)/sum(as.vector(t))) * 100)
    offset_temp[, around_site := "T"
                ][percentage == 0, around_site := "F"]
    tempoff <- function(v_dist){
      ttable <- sort(table(v_dist), decreasing = T)
      ttable_sr <- ttable[as.character(as.numeric(names(ttable))+1)]
      ttable_sl <- ttable[as.character(as.numeric(names(ttable))-1)]
      tsel <- rowSums(cbind(ttable > ttable_sr, ttable > ttable_sl), na.rm = T)
      return(as.numeric(names(tsel[tsel == 2][1])))
    }
    
    offset_temp5 <- site_sub[, list(offset_from_5 = tempoff(.SD$site_dist_end5)), by = length]
    offset_temp3 <- site_sub[, list(offset_from_3 = tempoff(.SD$site_dist_end3)), by = length]
    merge_allx <- function(x, y) merge(x, y, all.x = TRUE, by = "length")
    offset_temp  <-  Reduce(merge_allx, list(offset_temp, offset_temp5, offset_temp3))
    
    # adjusted offset
    adj_off <- function(dt_site, dist_site, add, bestoff){
      temp_v <- dt_site[[dist_site]]
      t <- table(factor(temp_v, levels = seq(min(temp_v) - 2, max(temp_v) + add)))
      t[1:2] <- t[3] + 1
      locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) == -2)]))) + 1
      adjoff <- locmax[which.min(abs(locmax - bestoff))]
      ifelse(length(adjoff) != 0, adjoff, bestoff)
    }
    
    best_from5_tab <- offset_temp[!is.na(offset_from_5), list(perc = sum(percentage)), by = offset_from_5
                                  ][perc == max(perc)]
    best_from3_tab <- offset_temp[!is.na(offset_from_5), list(perc = sum(percentage)), by = offset_from_3
                                  ][perc == max(perc)]
    
    if(extremity == "auto" &
       ((best_from3_tab[, perc] > best_from5_tab[, perc] &
         as.numeric(best_from3_tab[, offset_from_3]) <= minlen - 2) |
        (best_from3_tab[, perc] <= best_from5_tab[, perc] &
         as.numeric(best_from5_tab[, offset_from_5]) > minlen - 1)) |
       extremity == "3end"){
      best_offset <- as.numeric(best_from3_tab[, offset_from_3])
      line_plot <- "3end"
      adj_tab <- site_sub[, list(corrected_offset_from_3 = adj_off(.SD, "site_dist_end3", 0, best_offset)), by = length]
      offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "length")
      offset_temp[is.na(corrected_offset_from_3), corrected_offset_from_3 := best_offset
                  ][, corrected_offset_from_5 := -corrected_offset_from_3 + length - 1]
    } else {
      if(extremity == "auto" &
         ((best_from3_tab[, perc] <= best_from5_tab[, perc] &
           as.numeric(best_from5_tab[, offset_from_5]) <= minlen - 1) |
          (best_from3_tab[, perc] > best_from5_tab[, perc] &
           as.numeric(best_from3_tab[, offset_from_3]) > minlen - 2)) |
         extremity == "5end"){
        best_offset <- as.numeric(best_from5_tab[, offset_from_5])
        line_plot <- "5end"
        adj_tab <- site_sub[, list(corrected_offset_from_5 = adj_off(.SD, "site_dist_end5", 1, best_offset)), by = length]
        offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "length")
        offset_temp[is.na(corrected_offset_from_5), corrected_offset_from_5 := best_offset
                    ][, corrected_offset_from_5 := abs(corrected_offset_from_5)
                      ][, corrected_offset_from_3 := abs(corrected_offset_from_5 - length + 1)]
      }
    }
    
    cat(sprintf("best offset: %i nts from the %s\n", abs(best_offset), gsub("end", "' end", line_plot)))
    
    if(log_file == T | log_file == TRUE){
      cat(sprintf("%s\t%s\t%i\n", n, gsub("end", "'end", line_plot), abs(best_offset)), file = logpath, append = TRUE)
    }
    
    t <- table(factor(dt$length, levels = lev))
    offset_temp[!is.na(offset_from_5), offset_from_5 := abs(offset_from_5)
                ][, total_percentage := as.numeric(format(round((as.vector(t)/sum(as.vector(t))) * 100, 3), nsmall=4))
                  ][, percentage := as.numeric(format(round(percentage, 3), nsmall=4))
                    ][, sample := n]
    
    setcolorder(offset_temp, c("length", "total_percentage", "percentage", "around_site", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
    if(start == TRUE | start == T){
      setnames(offset_temp, c("length", "total_percentage", "start_percentage", "around_start", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
    } else {
      setnames(offset_temp, c("length", "total_percentage", "stop_percentage", "around_stop", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
    }
    
    # plot
    if (plot == T | plot == TRUE) {
      options(warn=-1)
      if (length(plot_dir) == 0) {
        dir <- getwd()
        plot_dir <- paste(dir, "/offset_plot", sep = "")
      }
      if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
      }
      minlen <- ceiling(quantile(site_sub$length, (1 - cl/100)/2))
      maxlen <- ceiling(quantile(site_sub$length, 1 - (1 - cl/100)/2))
      for (len in minlen:maxlen) {
        progress <- ceiling(((len + 1 - minlen)/(maxlen - minlen + 1)) * 25)
        cat(sprintf("\rplotting   %s\r", paste(paste(rep(c(" ", "<<", "-"),
                                                         c(25 - progress, 1, progress)), collapse = ""), " ", as.character(progress*4),
                                               "% ", paste(rep(c("-", ">>", " "), c(progress, 1, 25 - progress)), collapse = ""), sep = "")))
        site_temp <- dt[site_dist_end5 %in% seq(-len + 1, 0) & length == len]
        site_tab5 <- data.table(table(factor(site_temp$site_dist_end5, levels = (-len + 1) : (len))))
        site_temp <- dt[site_dist_end3 %in% seq(0, len - 2) & length == len]
        site_tab3 <- data.table(table(factor(site_temp$site_dist_end3, levels = (-len) : (len - 2))))
        setnames(site_tab5, c("distance", "reads"))
        setnames(site_tab3, c("distance", "reads"))
        site_tab5[, distance := as.numeric(as.character(site_tab5$distance))
                  ][, extremity := "5' end"]
        site_tab3[, distance := as.numeric(as.character(site_tab3$distance))
                  ][, extremity := "3' end"]
        final_tab <- rbind(site_tab5[distance <= 0], site_tab3[distance >= 0])
        final_tab[, extremity := factor(extremity, levels = c("5' end", "3' end"))]
        
        p <- ggplot(final_tab, aes(distance, reads, color = extremity)) +
          geom_line() +
          geom_vline(xintercept = seq(floor(min(final_tab$distance)/3) * 3, floor(max(final_tab$distance)/3) * 3, 3), linetype = 2, color = "gray90") +
          geom_vline(xintercept = 0, color = "gray50") +
          geom_vline(xintercept = - offset_temp[length == len, offset_from_5], color = "#D55E00", linetype = 2, size = 1.1) +
          geom_vline(xintercept = offset_temp[length == len, offset_from_3], color = "#56B4E9", linetype = 2, size = 1.1) +
          geom_vline(xintercept = - offset_temp[length == len, corrected_offset_from_5], color = "#D55E00", size = 1.1) +
          geom_vline(xintercept = offset_temp[length == len, corrected_offset_from_3], color = "#56B4E9", size = 1.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - len, xmax = -flanking , fill = "#D55E00", alpha = 0.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - 1 , xmax = len - flanking - 1, fill = "#56B4E9", alpha = 0.1) +
          labs(x = "Distance from start (nt)", y = "Number of read extremities", title = paste(n, " - length=", len, " nts", sep = ""), color= "Extremity") +
          theme_bw(base_size = 20) +
          scale_fill_discrete("") +
          theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5))
        
        if(line_plot == "3end"){
          p <- p + geom_vline(xintercept = best_offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best_offset - len + 1, color = "black", linetype = 3, size = 1.1)
        } else {
          p <- p + geom_vline(xintercept = best_offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best_offset + len - 1, color = "black", linetype = 3, size = 1.1)
        }
        
        p <- p + 
          scale_x_continuous(limits = c(min(final_tab$distance), max(final_tab$distance)),
                             breaks = seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5), 
                             labels = as.character(seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5) + base))
        
        subplot_dir <- paste(plot_dir, n, sep = "/")
        dir.create(subplot_dir)
        ggsave(paste(subplot_dir, "/", len, ".", plot_format, sep = ""), plot = p, width = 15, height = 5, units = "in")
      }
      cat(sprintf("\rplotting   %s\n",
                  paste(paste(rep(c(" ", "<<", "-"), c(25 - progress, 1, progress)), collapse = ""), " ", 
                        as.character(progress*4), "% ", 
                        paste(rep(c("-", ">>", " "), c(progress, 1, 25 - progress)), collapse = ""), sep = "")))
      options(warn=0)
    }
    
    dt[, c("site_dist_end5", "site_dist_end3") := NULL]
    offset <- rbind(offset, offset_temp)
  }
  return(offset)
}
                  

#' Update reads information according to the inferred P-sites.
#' 
#' This function provides additional reads information according to the position
#' of the P-site identfied by \code{\link{psite}}. It attaches to each data
#' table in a list four columns reporting i) the P-site position with respect to
#' the 1st nucleotide of the transcript, ii) the P-site position with respect to
#' the start and the stop codon of the annotated coding sequence (if any) and
#' iii) the region of the transcript (5' UTR, CDS, 3' UTR) that includes the
#' P-site. Please note: for transcripts not associated to any annotated CDS the
#' position of the P-site with respect to the start and the stop codon is set to
#' NA. Optionally, additional columns reporting the three nucleotides covered by
#' the P-site, the A-site and the E-site are attached, based on FASTA files or
#' BSgenome data packages containing the transcript nucleotide sequences.
#'
#' @param data List of data tables from \code{\link{bamtolist}},
#'   \code{\link{bedtolist}} or \code{\link{length_filter}}.
#' @param offset Data table from \code{\link{psite}}.
#' @param site Either "psite, "asite", "esite" or a combination of these
#'   strings. It specifies if additional column(s) reporting the three
#'   nucleotides covered by the ribosome P-site ("psite"), A-site ("asite") and
#'   E-site ("esite") should be added. Note: either \code{fastapath} or
#'   \code{bsgenome} is required for this purpose. Default is NULL.
#' @param fastapath Character string specifying the FASTA file used in the
#'   alignment step, including its path, name and extension. This file can
#'   contain reference nucleotide sequences either of a genome assembly or of
#'   all the transcripts (see \code{Details} and \code{fasta_genome}). Please
#'   make sure the sequences derive from the same release of the annotation file
#'   used in the \code{\link{create_annotation}} function. Note: either
#'   \code{fastapath} or \code{bsgenome} is required to generate additional
#'   column(s) specified by \code{site}. Default is NULL.
#' @param fasta_genome Logical value whether the FASTA file specified by
#'   \code{fastapath} contains nucleotide sequences of a genome assembly. If
#'   TRUE (the default), an annotation object is required (see \code{gtfpath}
#'   and \code{txdb}). FALSE implies the nucleotide sequences of all the
#'   transcripts is provided instead.
#' @param refseq_sep Character specifying the separator between reference
#'   sequences' name and additional information to discard, stored in the
#'   headers of the FASTA file specified by \code{fastapath} (if any). It might
#'   be required for matching the reference sequences' identifiers reported in
#'   the input list of data tables. All characters before the first occurrence
#'   of the specified separator are kept. Default is NULL i.e. no string
#'   splitting is performed.
#' @param bsgenome Character string specifying the BSgenome data package with
#'   the genome sequences to be loaded. If not already present in the system, it
#'   is automatically installed through the biocLite.R script (check the list of
#'   available BSgenome data packages by running the
#'   \code{\link[BSgenome]{available.genomes}} function of the BSgenome
#'   package). This parameter must be coupled with an annotation object (see
#'   \code{gtfpath} and \code{txdb}). Please make sure the sequences included in
#'   the specified BSgenome data pakage are in agreement with the sequences used
#'   in the alignment step. Note: either \code{fastapath} or \code{bsgenome} is
#'   required to generate additional column(s) specified by \code{site}. Default
#'   is NULL.
#' @param gtfpath Character string specifying the location of a GTF file,
#'   including its path, name and extension. Please make sure the GTF file and
#'   the sequences specified by \code{fastapath} or \code{bsgenome} derive from
#'   the same release. Note that either \code{gtfpath} or \code{txdb} is
#'   required if and only if nucleotide sequences of a genome assembly are
#'   provided (see \code{fastapath} or \code{bsgenome}). Default is NULL.
#' @param txdb Character string specifying the TxDb annotation package to be
#'   loaded. If not already present in the system, it is automatically installed
#'   through the biocLite.R script (check
#'   \href{http://bioconductor.org/packages/release/BiocViews.html#___TxDb}{here}
#'   the list of available TxDb annotation packages). Please make sure the TxDb
#'   annotation package and the sequences specified by \code{fastapath} or
#'   \code{bsgenome} derive from the same release. Note that either
#'   \code{gtfpath} or \code{txdb} is required if and only if nucleotide
#'   sequences of a genome assembly are provided (see \code{fastapath} or
#'   \code{bsgenome}). Default is NULL.
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
#' @param granges Logical value whether to return a GRangesList object. Default
#'   is FALSE i.e. a list of data tables (the required input for downstream
#'   analyses and graphical outputs provided by riboWaltz) is returned instead.
#' @details \strong{riboWaltz} only works for read alignments based on
#'   transcript coordinates. This choice is due to the main purpose of RiboSeq
#'   assays to study translational events through the isolation and sequencing
#'   of ribosome protected fragments. Most reads from RiboSeq are supposed to
#'   map on mRNAs and not on introns and intergenic regions. Nevertheless, BAM
#'   based on transcript coordinates can be generated in two ways: i) aligning
#'   directly against transcript sequences; ii) aligning against standard
#'   chromosome sequences, requiring the outputs to be translated in transcript
#'   coordinates. The first option can be easily handled by many aligners (e.g.
#'   Bowtie), given a reference FASTA file where each sequence represents a
#'   transcript, from the beginning of the 5' UTR to the end of the 3' UTR. The
#'   second procedure is based on reference FASTA files where each sequence
#'   represents a chromosome, usually coupled with comprehensive gene annotation
#'   files (GTF or GFF). The STAR aligner, with its option --quantMode
#'   TranscriptomeSAM (see Chapter 6 of its
#'   \href{http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf}{manual}),
#'    is an example of tool providing such a feature.
#' @return A list of data tables or a GRangesList object.
#' @examples
#' data(reads_list)
#' data(psite_offset)
#' data(mm81cdna)
#' 
#' reads_psite_list <- psite_info(reads_list, psite_offset)
#' @import data.table
#' @export
psite_info <- function(data, offset, site = NULL, fastapath = NULL, 
                       fasta_genome = TRUE, refseq_sep = NULL, bsgenome = NULL,
                       gtfpath = NULL, txdb = NULL, dataSource = NA,
                       organism = NA, granges = FALSE) {
  
  if(!(all(site %in% c("psite", "asite", "esite"))) & length(site) != 0){
    cat("\n")
    stop("parameter site must be either NULL, \"psite\", \"asite\", \"esite\" or a combination of the three strings \n\n")
  } else {
    if(length(site) != 0 & length(fastapath) == 0 & length(bsgenome) == 0){
      cat("\n")
      stop("parameter site is specified but both fastapath and bsgenome are missing \n\n")
    }
  }
  
  if(length(site) != 0){
    if(((length(fastapath) != 0 & (fasta_genome == TRUE | fasta_genome == T)) |
        length(bsgenome) != 0) &
       length(gtfpath) == 0 & length(txdb) == 0){
      cat("\n")
      stop("genome annotation file not specified (both GTF path and TxDb object are missing)\n\n")
    }
    
    if(length(fastapath) != 0 & length(bsgenome) != 0){
      cat("\n")
      warning("both fastapath and bsgenome are specified. Only fastapath will be considered\n")
      bsgenome = NULL
    }
    
    if(length(gtfpath) != 0 & length(txdb) != 0){
      cat("\n")
      warning("both gtfpath and txdb are specified. Only gtfpath will be considered\n")
      txdb = NULL
    }
    
    if((length(gtfpath) != 0 | length(txdb) != 0) &
       ((length(fastapath) == 0 & length(bsgenome) == 0) |
        (length(fastapath) != 0 & (fasta_genome == FALSE | fasta_genome == F)))){
      cat("\n")
      warning("a genome annotation file is specified but no sequences from genome assembly are provided\n")
    }
    
    if(length(gtfpath) != 0 | length(txdb) != 0){
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
    }
    
    if(length(fastapath) != 0 | length(bsgenome) != 0){
      if(length(fastapath) != 0) {
        if(fasta_genome == TRUE | fasta_genome == T){
          temp_sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
          if(length(refseq_sep) != 0){
            names(temp_sequences) <- tstrsplit(names(temp_sequences), refseq_sep, fixed = TRUE, keep = 1)[[1]]
          }
          exon <- suppressWarnings(GenomicFeatures::exonsBy(txdbanno, by = "tx", use.names = TRUE))
          exon <- as.data.table(exon[unique(names(exon))])
          sub_exon_plus <- exon[as.character(seqnames) %in% names(temp_sequences) & strand == "+"]
          sub_exon_minus <- exon[as.character(seqnames) %in% names(temp_sequences) & strand == "-"
                                 ][, new_end := Biostrings::width(temp_sequences[as.character(seqnames)]) - start + 1
                                   ][, new_start := Biostrings::width(temp_sequences[as.character(seqnames)]) - end + 1]
          
          seq_dt_plus <- sub_exon_plus[, nt_seq := "emp"
                                       ][, nt_seq := as.character(Biostrings::subseq(temp_sequences[as.character(seqnames)],
                                                                                     start = start,
                                                                                     end = end))
                                         ][, list(seq = paste(nt_seq, collapse = "")), by = group_name]
          
          revcompl_temp_sequences <- reverseComplement(temp_sequences)
          seq_dt_minus <- sub_exon_minus[, nt_seq := "emp"
                                         ][, nt_seq := as.character(Biostrings::subseq(revcompl_temp_sequences[as.character(seqnames)],
                                                                                       start = new_start,
                                                                                       end = new_end))
                                           ][, list(seq = paste(nt_seq, collapse = "")), by = group_name]
          
          sequences <- Biostrings::DNAStringSet(c(seq_dt_plus$seq, seq_dt_minus$seq))
          names(sequences) <- c(unique(sub_exon_plus$group_name), unique(sub_exon_minus$group_name))
        } else {
          sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
          if(length(refseq_sep) != 0){
            names(sequences) <- tstrsplit(names(sequences), refseq_sep, fixed = TRUE, keep = 1)[[1]]
          }
        }
      } else {
        if(bsgenome %in% installed.genomes()){
          library(bsgenome, character.only = TRUE)
        } else {
          source("http://www.bioconductor.org/biocLite.R")
          biocLite(bsgenome, suppressUpdates = TRUE)
          library(bsgenome, character.only = TRUE)
        }
        sequences <- GenomicFeatures::extractTranscriptSeqs(get(bsgenome), txdbanno, use.names=T)
      }
    }
  }
  
  names <- names(data)
  for (n in names) {
    cat(sprintf("processing %s\n", n))
    dt <- data[[n]]
    suboff <- offset[sample == n, .(length,corrected_offset_from_3)]
    cat("1. adding p-site position\n")
    dt[suboff,  on = 'length', psite := i.corrected_offset_from_3]
    dt[, psite := end3 - psite]
    setcolorder(dt,c("transcript", "end5", "psite", "end3", "length", "cds_start", "cds_stop"))
    dt[, psite_from_start := psite - cds_start
       ][cds_stop == 0, psite_from_start := 0]
    dt[, psite_from_stop := psite - cds_stop
       ][cds_stop == 0, psite_from_stop := 0]
    cat("2. adding transcript region\n")
    dt[, psite_region := "5utr"
       ][psite_from_start >= 0 & psite_from_stop <= 0, psite_region := "cds"
         ][psite_from_stop > 0, psite_region := "3utr"
           ][cds_stop == 0, psite_region := NA]
    if(length(site) != 0){
      cat("3. adding nucleotide sequence(s)\n")
      if("psite" %in% site){
        dt[, p_site_codon := as.character(Biostrings::subseq(sequences[as.character(dt$transcript)],
                                                             start = dt$psite,
                                                             end = dt$psite + 2))]
      }
      if("asite" %in% site){
        dt[, a_site_codon := as.character(Biostrings::subseq(sequences[as.character(dt$transcript)],
                                                             start = dt$psite + 3,
                                                             end = dt$psite + 5))]
      }
      if("esite" %in% site){
        dt[, e_site_codon := as.character(Biostrings::subseq(sequences[as.character(dt$transcript)],
                                                             start = dt$psite - 3,
                                                             end = dt$psite - 1))]
      }
    }
    
    setorder(dt, transcript, end5, end3)
    
    if (granges == T | granges == TRUE) {
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
    
    data[[n]] <- dt
  }
  
  if (granges == T | granges == TRUE) {
    data <- GenomicRanges::GRangesList(data)
  }
  
  return(data)
}
