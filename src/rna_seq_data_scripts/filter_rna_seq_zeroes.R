########################################################################
#
# Requires that RNA-Seq data have already been extracted (see
# script `extract_rna_seq_data.R`).
#
# Here, data are filtered by removing any covariates (columns) that have too
# many zeroes (as determined by the threshold). IMPORTANT NOTE: this operation
# should be performed only on the training set; whichever covariates are 
# removed can then be removed/ignored from any testing data. If this filtering
# is done on the entire dataset, this can introduce bias in the estimation of
# prediction error. 
#
########################################################################

filter_rna_seq_zeroes <- function(data_dir = "data/processed/", files,
                                  threshold = 0.8, training = NULL,
                                  patient_data_file = "data/patient_data") {
  cat(format(Sys.time()), "Filtering RNA-Seq Data...\n")
  cat("Running `filter_rna_seq_zeroes()` with threshold =", threshold, "and",
      "training indicies", ifelse(is.null(training), "not", ""), "provided.\n")
  
  cat("Filtering RNA-Seq data: removing columns if proportion of zeros is ", 
      "at least ", threshold * 100, "%.\n", sep = "")
  
  # Check that training indicies are provided. If not, use the default partition.
  if(is.null(training)) {
    cat("No training index provided.
        Default partition loaded from:", patient_data_file, "\n")
    data <- load_RData_file(patient_data_file)
    training <- which(data$training == "Training")
    rm(data)
  } else {
    cat("Training indicies provided:", training, "\n")
  }
  
  filter_data <- function(rna_seq_file_name, filtered_file_name) {
    cat("Loading RNA-seq data from `", paste(data_dir, rna_seq_file_name, sep = ""),
        "`.\n", sep = "")
    rna_seq_data <- load_RData_file(paste(data_dir, rna_seq_file_name, sep = ""))
    
    prop_of_zeroes <- apply(rna_seq_data[training, ], 2, function(x) mean(x == 0))
    index_high_prop_zeroes <- which(prop_of_zeroes >= threshold)
    
    cat("\t", length(index_high_prop_zeroes), " of ", ncol(rna_seq_data), 
        " columns removed with >= ", 100 * threshold, "% zeroes.\n", sep = "")
    if(length(index_high_prop_zeroes) > 0) {
      rna_seq_data <- rna_seq_data[, -index_high_prop_zeroes]
    }
    
    outfile <- paste(data_dir, filtered_file_name, sep = "")
    cat("Saving filtered saved to `", outfile, "`.\n\n", sep = "")
    save(rna_seq_data, file = outfile, compress = TRUE)
    
    rm(rna_seq_data)
  }
  
  for(file in files) {
    filter_data(file$rna_seq_file_name, file$filtered_file_name)
    gc()
  }
}
