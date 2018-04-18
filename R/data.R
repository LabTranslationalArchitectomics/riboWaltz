#' Annotation
#'
#' A dataset containing basic information about 109,712 mouse mRNA (using the
#' Ensembl v81 transcript annotation).
#'
#' @format A data table with 109,712 rows and 5 variables (the lengths are
#'   expressed in nucleotides):
#' \describe{
#'   \item{transcript}{Name of the transcript (ENST ID and version, dot
#'   separated)}
#'   \item{l_tr}{Length of the transcript}
#'   \item{l_utr5}{Length of the annotated 5' UTR (if any)}
#'   \item{l_cds}{Length of the annotated CDS (if any)}
#'   \item{l_utr3}{Length of the annotated 3' UTR (if any)}
#' }
"mm81cdna"


#' P-site offsets
#'
#' This dataset contains information on the offset computed by \code{\link{psite}}
#' starting from \code{\link{reads_list}}.
#'
#' @format A data table with 31 rows and 9 variables (the lengths and the
#'   distances are expressed in nucleotides):
#' \describe{
#'   \item{length}{Length of the read}
#'   \item{total_percentage}{Percentage of reads of the considered length in the
#'   whole dataset}
#'   \item{start_percentage}{Percentage of reads of the considered length
#'   aligning on the start codon (if any)}
#'   \item{around_start}{A logical value reporting whether at least one read of
#'   the specified length aligns on the start codon (T = yes, F = no)}
#'   \item{offset_from_5}{Temporary P-site offset from the 5' end of read
#'   (before the correction step)}
#'   \item{offset_from_3}{Temporary P-site offset from the 3' end of read
#'   (before the correction step)}
#'   \item{adj_offset_from_5}{P-site offset from the 5' end of read after the
#'   correction step}
#'   \item{adj_offset_from_3}{P-site offset from the 3' end of read after the
#'   correction step}
#'   \item{sample}{Name of the sample}
#' }
"psite_offset"


#' Reads information 
#'
#' This dataset contains details on mapping reads from BAM or BED files.
#'
#' @format A list of data tables with 1 object (named \emph{Samp1}) of 393,338
#'   rows and 6 variables (the lengths and the distances are expressed in
#'   nucleotides):
#' \describe{
#'   \item{transcript}{Name of the transcript (ENST ID and version, dot
#'   separated)}
#'   \item{end5}{Position of the 5' end of the read with respect to the first
#'   nuclotide of the transcript}
#'   \item{end3}{Position of the 3' end of the read with respect to the first
#'   nuclotide of the transcript}
#'   \item{length}{Length of the read}
#'   \item{start_pos}{Leftmost position of the annotated CDS with respect to the
#'   first nuclotide of the transcript}
#'   \item{stop_pos}{Rightmost position of the annotated CDS with respect to the
#'   first nuclotide of the transcript}
#' }
"reads_list"


#' P-sites and reads information 
#'
#' This dataset contains details on mapping reads after the identification of
#' the P-site and the update of \code{\link{reads_list}}.
#'
#' @format A list of data tables with 1 object (named \emph{Samp1}) of 393,338
#'   rows and 10 variables (the lengths and the distances are expressed in
#'   nucleotides):
#' \describe{
#'   \item{transcript}{Name of the transcript (ENST ID and version, dot
#'   separated)}
#'   \item{end5}{Position of the 5' end of the read with respect to the first
#'   nuclotide of the transcript}
#'   \item{psite}{Position of the P-site with respect to the first nuclotide of
#'   the transcript}
#'   \item{end3}{Position of the 3' end of the read with respect to the first
#'   nuclotide of the transcript}
#'   \item{length}{Length of the read}
#'   \item{start_pos}{Leftmost position of the CDS with respect to the first
#'   nuclotide of the transcript}
#'   \item{stop_pos}{Rightmost position of the CDS with respect to the first
#'   nuclotide of the transcript}
#'   \item{psite_from_start}{Position of the P-site with respect to the first
#'   nuclotide of the annotated CDS (if any)}
#'   \item{psite_from_stop}{Position of the P-site with respect to the last
#'   nuclotide of the annotated CDS (if any)}
#'   \item{psite_region}{Region of the transcript that includes the P-site
#'   (5utr, cds, 3utr)}
#' }
"reads_psite_list"