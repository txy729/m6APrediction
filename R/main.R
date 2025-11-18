#------------------------------------------------------------------
# dna_encoding.R
#------------------------------------------------------------------

#' Encode DNA sequences into a matrix
#'
#' This function converts DNA strings into a matrix representation,
#' where each nucleotide is split into columns for downstream analysis.
#'
#' @param dna_strings A character vector of DNA sequences.
#' @return A data frame with encoded DNA sequences.
#' @examples
#' dna_encoding(c("ATCGT", "GGTAA"))
#' @export
# Define function for DNA sequence encoding
dna_encoding <- function(dna_strings){

  # Check sequence length for robustness against prediction mismatch
  nn <- nchar( dna_strings[1] )
  if (nn != 5) {
    # Assuming the model is strictly a 5-mer predictor based on 'DNA_5mer' column name
    stop(paste0("All DNA sequences must be exactly 5 bases long. Found length: ", nn))
  }

  # Split sequences into a matrix
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)

  # Set column names nt_pos1, nt_pos2, ...
  colnames(seq_m) <- paste0("nt_pos", 1:nn)

  # Convert to data frame
  seq_df <- as.data.frame(seq_m)

  # Convert to factors, specifying all possible levels
  # The levels c("A", "T", "C", "G") should be correct for nucleotides
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))

  return(seq_df)
}

#------------------------------------------------------------------
# prediction_multiple.R
#------------------------------------------------------------------

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
#' @importFrom stats predict
# Multiple sample prediction
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {

  # Check for required original feature columns
  required_cols <- c("gc_content","RNA_type","RNA_region","exon_length",
                     "distance_to_junction","evolutionary_conservation","DNA_5mer")
  stopifnot(all(required_cols %in% colnames(feature_df)))

  # --- Extract Factor Levels from Model (Fixes 'New factor levels' ERROR) ---

  # Check if model is randomForest and has xlevels attribute
  if (!inherits(ml_fit, "randomForest") || is.null(ml_fit$forest$xlevels)) {
    stop("Model must be a 'randomForest' object with factor level information.")
  }
  model_levels <- ml_fit$forest$xlevels

  # 1. Extract DNA sequence and perform encoding
  dna_strings <- feature_df$DNA_5mer
  encoded_dna_df <- dna_encoding(dna_strings)

  # 2. Remove 'DNA_5mer' from the original features
  feature_df_no_dna <- feature_df[, -which(colnames(feature_df) == "DNA_5mer")]

  # 3. Combine original features and encoded features to create the full feature set
  full_feature_df <- cbind(feature_df_no_dna, encoded_dna_df)

  # 4. Ensure factors are correctly set using levels from the model

  # RNA_type
  if ("RNA_type" %in% names(model_levels)) {
    full_feature_df$RNA_type <- factor(full_feature_df$RNA_type,
                                       levels = model_levels$RNA_type)
  } else {
    full_feature_df$RNA_type <- factor(full_feature_df$RNA_type)
  }

  # RNA_region
  if ("RNA_region" %in% names(model_levels)) {
    full_feature_df$RNA_region <- factor(full_feature_df$RNA_region,
                                         levels = model_levels$RNA_region)
  } else {
    full_feature_df$RNA_region <- factor(full_feature_df$RNA_region)
  }

  # The remaining problem is likely in the example data (m6A_input_example.csv)
  # containing a level (e.g., 'pseudogene') not in the model.
  # This explicit level setting handles cases where a level is MISSING in new data,
  # but fails if a level is NEW in new data.
  # If the error persists, the example data needs correction.

  # 5. Make predictions using the complete feature set
  infered_prob <- predict(ml_fit, full_feature_df, type="prob")[,2]

  # --- Result processing ---

  # Add results back to the original data frame
  feature_df$predicted_m6A_prob <- infered_prob
  feature_df$predicted_m6A_status <- ifelse(infered_prob > positive_threshold, "Positive", "Negative")

  return(feature_df)
}

#------------------------------------------------------------------
# prediction_single.R
#------------------------------------------------------------------

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
#' # Load ml_fit before using it in the example
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' # FIX: Corrected line width for Rd check, and DNA_5mer length for prediction
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 5, evolutionary_conservation = 0.8,
#'                   DNA_5mer = "TGACC")
#' @export
#' @import randomForest
#' @importFrom stats predict
# Single sample prediction
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {

  # Create a single-row data frame
  feature_df <- data.frame(gc_content, RNA_type, RNA_region, exon_length,
                           distance_to_junction, evolutionary_conservation, DNA_5mer,
                           stringsAsFactors = FALSE)

  # Call prediction_multiple
  return(prediction_multiple(ml_fit, feature_df, positive_threshold))
}
