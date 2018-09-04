#' Identify the ribosome P-site position within the reads.
#'
#' This function identifies within each read the position of the ribosome
#' P-site, determined by the localisation of its first nucleotide. The function
#' processes the samples separately starting from the reads aligning on the
#' reference codon (selected by the user between the start codon and the second
#' to last codon) of any annotated coding sequence. It then returns the position
#' of the P-site specifically inferred for all the read lengths. It also allows
#' to plot a collection of read length-specific occupancy metaprofiles
#' showing the P-sites offsets computed throughout the two steps of the
#' algorithm.
#'
#' @param data A list of data tables from \code{\link{bamtolist}},
#'   \code{\link{bedtolist}} or \code{\link{length_filter}}.
#' @param flanking An integer that specifies how many nucleotides, at least, of
#'   the reads mapping on the reference codon must flank the reference codon in
#'   both directions. Default is 6.
#' @param start A logical value whether ot not to compute the P-site offsets
#'   starting from the reads aligning on the translation initiation site. FALSE
#'   implies that the reads mapping on the last triplet before the stop codon
#'   are used instead. Default is TRUE.
#' @param extremity A character string specifing which extremity of the reads
#'   should be used in the correction step of the algorithm. It can be either
#'   "5end" or "3end" for the 5' and the 3' extremity, respectively. Default is
#'   "auto", meaning that the best extremity is automatically selected.
#' @param plot A logical value whether or not to plot the occupancy metaprofiles
#'   showing the P-sites offsets computed throughout the two steps of the
#'   algorithm. Default is FALSE.
#' @param plotdir A character string specifying the (existing or not) location
#'   of the directory where the occupancy metaprofiles shuold be stored. This
#'   parameter is considered only if \code{plot} is TRUE. By default this
#'   argument is NULL, which implies it is set as a subfolder of the working
#'   directory, called \emph{offset_plot}.
#' @param plotformat Either "png" (the default) or "pdf", this parameter
#'   specifies the file format of the generated metaprofiles. It is considered
#'   only if \code{plot} is TRUE.
#' @param cl An integer value in \emph{[1,100]} specifying the confidence level
#'   for restricting the generation of the occupancy metaprofiles to a sub-range
#'   of read lengths. By default it is set to 99. This parameter is considered
#'   only if \code{plot} is TRUE.
#' @return A data table.
#' @examples
#' data(reads_list)
#'
#' ## Compute the P-site offset automatically selecting the otimal read
#' ## extremity for the correction step and not plotting any metaprofile
#' psite(reads_list, flanking = 6, extremity="auto")
#'
#' ## Compute the P-site offset specifying the extremity used in the correction
#' ## step and plotting the metaprofiles only for a sub-range of read lengths (the
#' ## middle 95%). The plots will be placed in the current working directory.
#' psite_offset <- psite(reads_list, flanking = 6, extremity="3end", plot = TRUE, cl = 95)
#' @import data.table
#' @import ggplot2
#' @export
psite <- function(data, flanking = 6, start = TRUE, extremity="auto", plot = FALSE,
                  plotdir = NULL, plotformat="png", cl = 99) {
  names <- names(data)
  offset <- NULL
  for (n in names) { 
    cat(sprintf("processing %s\n", n))
    dt <- data[[n]]
    lev <- sort(unique(dt$length))
    
    if(start == T | start == TRUE){
      base <- 0
      dt[, site_dist_end5 := end5 - start_pos]
      dt[, site_dist_end3 := end3 - start_pos]
    } else {
      base <- -5
      dt[, site_dist_end5 := end5 - stop_pos - base]
      dt[, site_dist_end3 := end3 - stop_pos - base]
    }
    site_sub <- dt[site_dist_end5 <= -flanking & site_dist_end3 >= flanking - 1]
    minlen <- min(site_sub$length)
    maxlen <- max(site_sub$length)
    t <- table(factor(site_sub$length, levels = lev))
    
    # offset
    offset_temp <- data.table(length = as.numeric(as.character(names(t))), percentage = (as.vector(t)/sum(as.vector(t))) * 100)
    offset_temp[, around_site := "T"
                ][percentage == 0, around_site := "F"]
    offset_temp5 <- site_sub[, list(offset_from_5 = as.numeric(names(which.max(table(site_dist_end5))))), by = length]
    offset_temp3 <- site_sub[, list(offset_from_3 = as.numeric(names(which.max(table(site_dist_end3))))), by = length]
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
    
    best_from5_tab <- offset_temp[, list(perc = sum(percentage)), offset_from_5
                                  ][perc == max(perc)]
    best_from3_tab <- offset_temp[, list(perc = sum(percentage)), offset_from_3
                                  ][perc == max(perc)]

    if(extremity == "auto" &
       ((best_from3_tab[, perc] > best_from5_tab[, perc] &
         as.numeric(best_from3_tab[, offset_from_3]) <= minlen - 2) |
        (best_from3_tab[, perc] <= best_from5_tab[, perc] &
         as.numeric(best_from5_tab[, offset_from_5]) <= minlen - 1)) |
       extremity == "3end"){
      best_offset <- as.numeric(best_from3_tab[, offset_from_3])
      line_plot <- "from3"
      cat(sprintf("best offset: %i nts from the 3' end\n", best_offset))
      adj_tab <- site_sub[, list(adj_offset_from_3 = adj_off(.SD, "site_dist_end3", 0, best_offset)), by = length]
      offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "length")
      offset_temp[is.na(adj_offset_from_3), adj_offset_from_3 := best_offset
                  ][, adj_offset_from_5 := -adj_offset_from_3 + length - 1]
    } else {
      if(extremity == "auto" &
         ((best_from3_tab[, perc] <= best_from5_tab[, perc] &
           as.numeric(best_from5_tab[, offset_from_5]) <= minlen - 1) |
          (best_from3_tab[, perc] > best_from5_tab[, perc] &
           as.numeric(best_from3_tab[, offset_from_3]) > minlen - 2)) |
         extremity == "5end"){
        best_offset <- as.numeric(best_from5_tab[, offset_from_5])
        line_plot <- "from5"
        cat(sprintf("best offset: %i nts from the 5' end\n", -best_offset))
        adj_tab <- site_sub[, list(adj_offset_from_5 = adj_off(.SD, "site_dist_end5", 1, best_offset)), by = length]
        offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "length")
        offset_temp[is.na(adj_offset_from_5), adj_offset_from_5 := best_offset
                    ][, adj_offset_from_5 := abs(best_offset)
                      ][, adj_offset_from_3 := abs(adj_offset_from_5 - length + 1)]
      }
    }
    
    t <- table(factor(dt$length, levels = lev))
    offset_temp[!is.na(offset_from_5), offset_from_5 := abs(offset_from_5)
                ][, total_percentage := as.numeric(format(round((as.vector(t)/sum(as.vector(t))) * 100, 3), nsmall=4))
                  ][, percentage := as.numeric(format(round(percentage, 3), nsmall=4))
                    ][, sample := n]
                
    setcolorder(offset_temp, c("length", "total_percentage", "percentage", "around_site", "offset_from_5", "offset_from_3", "adj_offset_from_5", "adj_offset_from_3", "sample"))
    if(start == TRUE | start == T){
      setnames(offset_temp, c("length", "total_percentage", "start_percentage", "around_start", "offset_from_5", "offset_from_3", "adj_offset_from_5", "adj_offset_from_3", "sample"))
    } else {
      setnames(offset_temp, c("length", "total_percentage", "stop_percentage", "around_stop", "offset_from_5", "offset_from_3", "adj_offset_from_5", "adj_offset_from_3", "sample"))
    }

    # plot
    if (plot == T || plot == TRUE) {
      options(warn=-1)
      if (length(plotdir) == 0) {
        dir <- getwd()
        plotdir <- paste(dir, "/offset_plot", sep = "")
      }
      if (!dir.exists(plotdir)) {
        dir.create(plotdir)
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
          geom_vline(xintercept = - offset_temp[length == len, adj_offset_from_5], color = "#D55E00", size = 1.1) +
          geom_vline(xintercept = offset_temp[length == len, adj_offset_from_3], color = "#56B4E9", size = 1.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - len, xmax = -flanking , fill = "#D55E00", alpha = 0.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - 1 , xmax = len - flanking - 1, fill = "#56B4E9", alpha = 0.1) +
          labs(x = "Distance from start (nt)", y = "Number of read extremities", title = paste(n, " - length=", len, " nts", sep = ""), color= "Extremity") +
          theme_bw(base_size = 20) +
          scale_fill_discrete("") +
          theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5))
        
        if(line_plot == "from3"){
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
        
        subplotdir <- paste(plotdir, n, sep = "/")
        dir.create(subplotdir)
        ggsave(paste(subplotdir, "/", len, ".", plotformat, sep = ""), plot = p, width = 15, height = 5, units = "in")
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
#' Starting ftom the P-site position identfied by \code{\link{psite}}, this
#' function updates the data tables that contains information about the reads.
#' It attaches to the data tables 4 columns reporting the P-site position with
#' respect to the 1st nucleotide of the transcript, the start and the stop codon
#' of the annotated coding sequence (if any) and the region of the transcript
#' (5' UTR, CDS, 3' UTR) that includes the P-site. Please note: if a transcript
#' is not associated to any annotated CDS then the positions of the P-site from
#' both the start and the stop codon is set to NA. One or more additional
#' columns reporting the three nucleotides covered by the P-site, the A-site or
#' the E-site can be attached by providing either a FASTA file or a BSgenome
#' data package with the nucleotide sequences.
#'
#' @param data A list of data tables from \code{\link{bamtolist}},
#'   \code{\link{bedtolist}} or \code{\link{length_filter}}.
#' @param offset A data table from \code{\link{psite}}.
#' @param site Either NULL, "psite, "asite", "esite" or a vector with a
#'   combination of the three character strings. When this parameter is not NULL
#'   (the default), it specifies which of the column(s) reporting the three
#'   nucleotides covered by the P-site ("psite"), A-site ("asite") or E-site
#'   ("esite") must be added. Note: either \code{fastapath} or \code{bsgenome}
#'   is required to generate the additional column(s).
#' @param fastapath An optional character string specifying the path to the
#'   FASTA file used in the alignment step, including its name and extension.
#'   This file can contain reference nucleotide sequences either of a genome
#'   assembly or of all the transcripts (see \code{fasta_genome}). Please make
#'   sure the sequences derive from the same release of the annotation file used
#'   in the \code{\link{create_annotation}} function. Note: either
#'   \code{fastapath} or \code{bsgenome} is required to generate the additional
#'   column(s) specified by \code{site}. Default is NULL.
#' @param fasta_genome A logical value whether or not the FASTA file specified
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
#'   package). This parameter also requires an annotation object (see
#'   \code{gtfpath} and \code{txdb}). Please make sure the sequences included in
#'   the specified BSgenome data pakage are in agreement with the sequences used
#'   in the alignment step. Note: either \code{fastapath} or \code{bsgenome} is
#'   required to generate the additional column(s) specified by \code{site}.
#'   Default is NULL.
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
#' @import data.table
#' @export
psite_info <- function(data, offset, site = NULL, fastapath = NULL, 
                       fasta_genome = TRUE, bsgenome = NULL, gtfpath = NULL,
                       txdb = NULL, dataSource = NA, organism = NA,
                       granges = FALSE) {
  
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
    
    if((length(fastapath) != 0 | length(bsgenome) != 0)){
      if(length(fastapath) != 0) {
        if(fasta_genome == TRUE | fasta_genome == T){
          temp_sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
          exon <- suppressWarnings(GenomicFeatures::exonsBy(txdbanno, by = "tx", use.names=TRUE))
          exon <- as.data.table(exon[unique(names(exon))])
          sub_exon <- exon[seqnames %in% names(temp_sequences)]
          seq_dt <- sub_exon[, list(seq = paste(Biostrings::subseq(temp_sequences[as.character(seqnames)],
                                                                   start = start,
                                                                   end = end), collapse="")),
                             by = group_name]
          sequences <- Biostrings::DNAStringSet(seq_dt$seq)
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
        sequences <- GenomicFeatures::extractTranscriptSeqs(get(bsgenome), txdbanno, use.names=T)
      }
    }
  }
  
  names <- names(data)
  for (n in names) {
    cat(sprintf("processing %s\n", n))
    dt <- data[[n]]
    suboff <- offset[sample == n, .(length,adj_offset_from_3)]
    cat("1. adding p-site position\n")
    dt[suboff,  on = 'length', psite := i.adj_offset_from_3]
    dt[, psite := end3 - psite]
    setcolorder(dt,c("transcript", "end5", "psite", "end3", "length", "start_pos", "stop_pos"))
    dt[, psite_from_start := psite - start_pos
       ][stop_pos == 0, psite_from_start := 0]
    dt[, psite_from_stop := psite - stop_pos
       ][stop_pos == 0, psite_from_stop := 0]
    cat("2. adding transcript region\n")
    dt[, psite_region := "5utr"
       ][psite_from_start >= 0 & psite_from_stop <= 0, psite_region := "cds"
         ][psite_from_stop > 0, psite_region := "3utr"
           ][stop_pos == 0, psite_region := NA]
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
    
    dt <- dt[order(transcript, end5)]
    
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
    
    data[[n]] <- dt
  }
  
  if (granges == T || granges == TRUE) {
    data <- GenomicRanges::GRangesList(data)
  }
  
  return(data)
}
