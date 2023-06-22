#' A function to turn non-corrected into corrected labelled peptide/protein fractions (Corr LPFs)
#'
#' This function allows users to input their non-corrected LPF matrix from "AnnotateProteins.R" and correct the fractional enrichment in individual peptides using the enrichment percentages in amino acid soluble pools from paired treatments. The function returns the corrected matrix, which can be directly used to calculate fractional synthesis rates.
#' @param InputPath Directory where the depurated qPCR matrix is stored.
#' @param InputPathTreatmentFile Directory where the position matrix is stored.
#' @keywords qPCR deltadeltaCt
#' @export
#' @examples
#'
#' ...

DeltaDeltaCtCc <- function(InputPath,
                           InputPathTreatmentFile) {

  qPCR_table <- read.table(InputPath, fill = TRUE)

  qPCR_pos_table <- read.table(InputPathTreatmentFile, fill = TRUE, header = T)





}


