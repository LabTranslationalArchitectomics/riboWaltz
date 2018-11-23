#' Annotation
#'
#' A dataset containing basic information about 109,712 mouse mRNA (Ensembl v81
#' transcript annotation).
#'
#' @format A data table with 109,712 rows and 5 variables:
#' \describe{
#'   \item{transcript}{Name of the transcript (ENST ID and version, dot
#'   separated)}
#'   \item{l_tr}{Length of the transcript, in nucleotides}
#'   \item{l_utr5}{Length of the annotated 5' UTR (if any), in nucleotides}
#'   \item{l_cds}{Length of the annotated CDS (if any), in nucleotides}
#'   \item{l_utr3}{Length of the annotated 3' UTR (if any), in nucleotides}
#' }
"mm81cdna"


#' Reads information 
#'
#' An example dataset containing details on reads mapping on the mouse
#' transcriptome, generated from BAM or BED files. A subset of the original
#' dataset is provided, including only reads aligning on the translation
#' initiation site. Please contact the authors for more information.
#'
#' @format A list of data tables with 1 object (named \emph{Samp1}) of 393,338
#'   rows and 6 variables:
#' \describe{
#'   \item{transcript}{Name of the transcript (ENST ID and version, dot
#'   separated)}
#'   \item{end5}{Position of the 5' end of the read with respect to the first
#'   nuclotide of the transcript, in nucleotides}
#'   \item{end3}{Position of the 3' end of the read with respect to the first
#'   nuclotide of the transcript, in nucleotides}
#'   \item{length}{Length of the read, in nucleotides}
#'   \item{cds_start}{Leftmost position of the annotated CDS with respect to the
#'   first nuclotide of the transcript, in nucleotides}
#'   \item{cds_stop}{Rightmost position of the annotated CDS with respect to the
#'   first nuclotide of the transcript, in nucleotides}
#' }
"reads_list"


#' P-site offsets
#'
#' An example dataset containing length-specific ribosome P-site offsets as
#' returned by \code{\link{psite}} applied to \code{\link{reads_list}}.
#'
#' @format A data table with 31 rows and 9 variables:
#' \describe{
#'   \item{length}{Length of the read, in nucleotides}
#'   \item{total_percentage}{Percentage of reads of the considered length in the
#'   whole dataset}
#'   \item{start_percentage}{Percentage of reads of the considered length
#'   aligning on the start codon (if any)}
#'   \item{around_start}{A logical value whether at least one read of
#'   the considered length aligns on the start codon (T = yes, F = no)}
#'   \item{offset_from_5}{Temporary P-site offset from the 5' end of the read,
#'   in nucleotides (before the correction step)}
#'   \item{offset_from_3}{Temporary P-site offset from the 3' end of the read,
#'   in nucleotides (before the correction step)}
#'   \item{corrected_offset_from_5}{P-site offset from the 5' end of the read, in
#'   nucleotides (after the correction step)}
#'   \item{corrected_offset_from_3}{P-site offset from the 3' end of the read, in
#'   nucleotides (after the correction step)}
#'   \item{sample}{Name of the sample}
#' }
"psite_offset"


#' Reads details updated with P-site information
#'
#' An example dataset that combines details on reads mapping on the mouse
#' transcriptome (see \code{\link{reads_list}}) and length-specific ribosome
#' P-site offsets (see \code{\link{psite_offset}}), as returned by
#' \code{\link{psite_info}}.
#'
#' @format A list of data tables with 1 object (named \emph{Samp1}) of 393,338
#'   rows and 10 variables:
#' \describe{
#'   \item{transcript}{Name of the transcript (ENST ID and version, dot
#'   separated)}
#'   \item{end5}{Position of the 5' end of the read with respect to the first
#'   nuclotide of the transcript, in nucleotides}
#'   \item{psite}{Position of the P-site with respect to the first nuclotide of
#'   the transcript, in nucleotides}
#'   \item{end3}{Position of the 3' end of the read with respect to the first
#'   nuclotide of the transcript, in nucleotides}
#'   \item{length}{Length of the read, in nucleotides}
#'   \item{cds_start}{Leftmost position of the CDS with respect to the first
#'   nuclotide of the transcript, in nucleotides}
#'   \item{cds_stop}{Rightmost position of the CDS with respect to the first
#'   nuclotide of the transcript, in nucleotides}
#'   \item{psite_from_start}{Position of the P-site with respect to the first
#'   nuclotide of the annotated CDS (if any), in nucleotides}
#'   \item{psite_from_stop}{Position of the P-site with respect to the last
#'   nuclotide of the annotated CDS (if any), in nucleotides}
#'   \item{psite_region}{Region of the transcript that includes the P-site
#'   (5utr, cds, 3utr)}
#' }
"reads_psite_list"
