#' Encode DNA sequences into a matrix
#'
#' This function converts DNA strings into a matrix representation,
#' where each nucleotide is split into columns for downstream analysis.
#'
#' @param dna_strings A character vector of DNA sequences.
#' @return A data frame with encoded DNA sequences.
#' @examples
#' dna_encoding(c("ATCG", "GGTA"))
#' @export
#Define function for DNA sequence encoding
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}
#' Predict m6A sites for multiple samples
#'
#' This function takes a feature data frame and predicts m6A probability
#' and classification status for each sample using a trained random forest model.
#'
#' @param ml_fit A trained random forest model.
#' @param feature_df A data frame containing features such as GC content,
#' RNA type, RNA region, exon length, distance to junction,
#' evolutionary conservation, and DNA 5-mer.
#' @param positive_threshold Numeric threshold (0–1) for classification, default 0.5.
#' @return The input data frame with additional columns:
#' \code{predicted_m6A_prob} and \code{predicted_m6A_status}.
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package="m6APrediction"))
#' prediction_multiple(ml_fit, example_df, positive_threshold = 0.6)
#' @export
#' @import randomForest
#Multiple
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  stopifnot(all(c("gc_content","RNA_type","RNA_region","exon_length",
                  "distance_to_junction","evolutionary_conservation","DNA_5mer") %in% colnames(feature_df)))

  feature_df$RNA_type <- factor(feature_df$RNA_type)
  feature_df$RNA_region <- factor(feature_df$RNA_region)

  infered_prob <- predict(ml_fit, feature_df, type="prob")[,2]
  feature_df$predicted_m6A_prob <- infered_prob
  feature_df$predicted_m6A_status <- ifelse(infered_prob > positive_threshold, "Positive", "Negative")
  return(feature_df)
}
#' Predict m6A sites for a single sample
#'
#' This function predicts m6A probability and classification status
#' for a single RNA sample by calling \code{prediction_multiple}.
#'
#' @param ml_fit A trained random forest model.
#' @param gc_content Numeric GC content.
#' @param RNA_type Character string specifying RNA type (e.g., "mRNA").
#' @param RNA_region Character string specifying RNA region (e.g., "CDS").
#' @param exon_length Numeric exon length.
#' @param distance_to_junction Numeric distance to junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer Character string representing DNA 5-mer.
#' @param positive_threshold Numeric threshold (0–1), default 0.5.
#' @return A named vector with predicted probability and status.
#' @examples
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#'                   exon_length = 12, distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8, DNA_5mer = "ATCGAT")
#' @export
#' @import randomForest
#Single
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  feature_df <- data.frame(gc_content, RNA_type, RNA_region, exon_length,
                           distance_to_junction, evolutionary_conservation, DNA_5mer,
                           stringsAsFactors = FALSE)
  return(prediction_multiple(ml_fit, feature_df, positive_threshold))
}
#
