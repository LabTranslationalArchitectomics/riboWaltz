#' Plot percentage of P-site per transcript region.
#' 
#' For each sample computes the percentage of P-sites falling in the three
#' regions of the transcripts (5' UTR, CDS and 3'UTR) and generates the
#' corresponding barplot. The function also calculates the percentage of region
#' length for the selected transcripts (reported in column "RNAs").
#' 
#' @param data A list of data frames from \code{\link{psite_info}}.
#' @param annotation A data frame with a reference annotation of the transripts.
#'   It must contain at least five columns named \emph{transcript}, 
#'   \emph{transcript_type}, \emph{l_utr5}, \emph{l_cds} and \emph{l_utr3} 
#'   containing the name of the transcripts (the same as in the reference 
#'   transcriptome), its transcript type, the position of the first nucleotide 
#'   of the \emph{5' UTR}, the \emph{CDS} and the  \emph{3' UTR}, respectively. 
#'   No specific order is required.
#' @param sample A character string vector specifying the name of the sample(s)
#'   of interest. By default this argument is NULL, meaning that all the data
#'   frames in \code{data} will be analysed.
#' @param transcripts A character string vector specifying the name of the 
#'   transcripts to be considered in the analysis. By default this argument is 
#'   NULL, which implies all the transcripts in \code{data} will be used. Either
#'   transcripts with no annotated \emph{5' UTR}, \emph{CDS} and
#'   \emph{3'UTR} or mRNAs not labeled as protein coding will be removed.
#' @param label A character string vector of the same length of \code{sample}
#'   specifying the name of the samples to be displaied in the plot. By default
#'   this argument is NULL i.e. the labels correspond to the sample names.
#' @param colour A character string vector of three elements specifying the
#'   colours of the bars (for \emph{5' UTR}, \emph{CDS} and \emph{3'UTR}
#'   respectively). By default is a grayscale.
#' @return A list containing a ggplot2 object, and a data frame with the
#'   associated data.
#' @examples
#' data(reads_list)
#' data(mm81cdna)
#'
#' reg_psite <- region_psite(reads_list, mm81cdna, sample = "Samp1")
#' reg_psite[["plot"]]
#' @import ggplot2
#' @export
region_psite <- function(data, annotation, sample = NULL, transcripts = NULL,
                       label = NULL, colour = c("gray70", "gray40", "gray10")) {
  rownames(annotation) <- as.character(annotation$transcript)
  
  l.transcripts <- rownames(annotation)[which(annotation$l_utr5 > 0 &
                                                annotation$l_cds >0 &
                                                annotation$l_cds %% 3 == 0 &
                                                annotation$l_utr3 > 0 &
                                                annotation$transcript_type == "protein_coding")]
  if (length(transcripts) == 0) {
    c.transcript <- l.transcripts
  } else {
    c.transcript <- intersect(l.transcripts, transcripts)
  }
  
  if(length(sample) == 0) {
    sample <- names(data)
  }
  
  if(length(label) == 0) {
    label <- sample
  }

  barplot.table <- as.data.frame(matrix("0", ncol = length(sample), nrow = 3))
  colnames(barplot.table) <- sample
  rownames(barplot.table) <- factor(c("5utr","cds","3utr"), levels = c("5utr","cds","3utr"))
  
  for(sam in sample) {
    df <- subset(data[[sam]], as.character(transcript) %in% c.transcript)
    barplot.table[, sam] <- as.data.frame(table(factor(df$psite_region, levels = c("5utr","cds","3utr"))))$Freq
  }
  
  norm.table <- t(t(barplot.table) / colSums(barplot.table)) * 100
  melt.table <- reshape::melt(norm.table)
  colnames(melt.table) <- c("region","sample","percentage")
  melt.table$class<-"mapped"
  melt.table$sample <- factor(melt.table$sample, levels = sample, labels = label)
  
  sub_anno <- subset(annotation, transcript %in% c.transcript)
  RNA_reg <- colSums(sub_anno[,c("l_utr5", "l_cds", "l_utr3")])
  RNA_reg_perc <- (RNA_reg / sum(RNA_reg)) *100
  
  RNA.table<-data.frame(region=factor(c("5utr","cds","3utr"), levels = c("5utr","cds","3utr")),
                        sample = rep("RNAs", 3),
                        percentage = RNA_reg_perc,
                        class = "rna")
  
  final_melt.table <- rbind(melt.table, RNA.table)
  final_melt.table$region <- factor(final_melt.table$region,
                                    levels = c("5utr", "cds", "3utr"),
                                    labels = c("5' UTR  ","CDS  ","3' UTR"))
  final_melt.table <- final_melt.table[order(final_melt.table$region),]
  
  bs <- 25
  bp <- ggplot(final_melt.table,aes(x = sample, y = percentage, fill = region)) +
    geom_bar(stat = "identity",color = "white",width = 0.65, alpha = 0.9, size = 0.025 * bs) +
    scale_fill_manual(name = "", values = colour) +
    theme_bw(base_size = bs) +
    theme(legend.position = "top", legend.key = element_blank()) +
    theme(axis.title.x = element_blank()) +
    scale_x_discrete(breaks = unique(final_melt.table$sample)) +
    facet_grid( . ~ class, scales = "free", space="free_x") +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    scale_y_continuous("P-sites (%)", sec.axis = sec_axis(~ . * 1 , name = "Length (%)"))
  
  output <- list()
  output[["plot"]] <- bp
  output[["df"]] <- final_melt.table
  return(output)
}









