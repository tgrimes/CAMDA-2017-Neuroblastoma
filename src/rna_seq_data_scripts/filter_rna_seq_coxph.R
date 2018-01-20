########################################################################
#
# Requires that RNA-Seq data have already been extracted (see
# script `extract_rna_seq_data.R`).
#
# Here, data are filtered by keeping only covariates (columns) that are found
# to be most highly associated with the response (as determined by Cox ph model). 
# IMPORTANT NOTE: this operation
# should be performed only on the training set; whichever covariates are 
# removed can then be removed/ignored from any testing data. If this filtering
# is done on the entire dataset, this can introduce bias in the estimation of
# prediction error. 
#
########################################################################


#TODO: input file name to load clinical values.
filter_rna_seq_coxph <- function(data_dir, files, clinical_data_file = NULL, 
                                 training = NULL, outcome = "OS",
                                 y = NULL, observed = NULL, use_hr_only = FALSE,
                                 patient_data_file = "data/patient_data",
                                 n_keep = 2000) {
  library(survival)
  
  cat(format(Sys.time()), "Filtering RNA-Seq Data...\n")
  cat("Running `filter_rna_seq_coxph()`\n")
  
  # Check that training indicies are provided. If not, use the default partition.
  if(is.null(training)) {
    cat("No training index provided.",
        "Default partition loaded from:", patient_data_file, "\n")
    data <- load_RData_file(patient_data_file)
    (training <- which(data$training == "Training"))
    rm(data)
  }else {
    cat("Training indicies provided:", training, "\n")
  }
  
  # Check if response `y`, censoring variable `observed`, and 
  # high risk indicator `hr` are provided.
  data <- load_RData_file(patient_data_file)
  if(!is.null(y)) {
    cat("Reponse `y` is provided:", y, "\n")
  } else {
    cat("No response variable provided.",
        "Loading survival times and censoring from", patient_data_file, "\n")
    if(toupper(outcome) == "OS") {
      cat('\t y <- os_days\n')
      y <- data$os_days[training]
      cat("\t", y, "\n")
    } else if(toupper(outcome) == "EFS") {
      cat('\t y <- efs_days\n')
      y <- data$efs_days[training]
      cat("\t", y, "\n")
    }
  }
  if(!is.null(observed)) {
    cat("Censoring indicator `observed` is provided:", observed, "\n")
  } else {
    cat("No censoring provided.",
        "Loading censoring from", patient_data_file, "\n")
    cat('\t observed <- death_from_disease == "Y"\n')
    observed <- (data$death_from_disease == "Y")[training]
    cat("\t", observed, "\n")
  }
  if(use_hr_only) {
    cat("Using high-risk patients only.\n")
    cat("Loading high-risk labels from", patient_data_file, "\n")
    cat('\t hr <- high_risk == "Y"\n')
    hr <- (data$high_risk == "Y")[training] # Use only high-risk observations.
    cat("\t", hr, "\n")
  } else {
    hr <- 1:length(y) # Use all observations.
  }
  rm(data)
  
  # Load clinical covariates, if file is provided.
  covariates <- NULL
  if(is.null(clinical_data_file)) {
    cat("No clinical data file provided.\n")
  } else {
    cat("Loading clinical data from:", clinical_data_file, "\n")
    covariates <- get_predictors(clinical_data_file)[training, ]
    if(ncol(covariates) > 1) {
      covariates <- covariates[hr, ]
    } else {
      covariates <- covariates[hr]
    }
  }
  
  cat("Reponse 'y' and censoring 'observed` are subset by high-risk index `hr`.\n")
  cat("y <- y[hr]\n")
  y <- y[hr]
  cat("\t", y, "\n")
  cat("observed <- observed[hr]\n")
  observed <- observed[hr]
  cat("\t", observed, "\n")
  
  filter_data <- function(intput_rna_seq_file_name, 
                          output_filtered_file_name,
                          output_pvalues_file_name) {
    cat("Loading RNA-seq data from `", 
        paste(data_dir, intput_rna_seq_file_name, sep = ""), "`.\n", sep = "")
    rna_seq_data <- load_RData_file(paste(data_dir, intput_rna_seq_file_name, sep = ""))
    
    # Run Cox PH on training set of RNA-seq data.
    coxph_p_value <- function(x) {
      fit <- coxph(Surv(y, observed) ~ cbind(covariates, x))
      fit_summary <- summary(fit)$coefficients
      pval <- fit_summary[nrow(fit_summary), ncol(fit_summary)]
      return(pval)
    }
    p_values <- apply((rna_seq_data[training, ])[hr, ], 2, coxph_p_value)
    cat("\t", sum(is.na(p_values)), "of", length(p_values), 
        "columns have NA p-values.\n")
    p_values[is.na(p_values)] <- 1
    
    n_keep <- min(n_keep, ncol(rna_seq_data))
    q <- sort(c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 
              round(n_keep / ncol(rna_seq_data), 4)))
    cat("\t", n_keep, "of", ncol(rna_seq_data), "columns are retained.\n")
    cat("\t p_values (quantile):", 
        paste(format(quantile(p_values, q), digits = 4), 
              paste("(", q, ")", sep = "")), "\n")
    
    top_genes <- order(p_values)[1:n_keep]
    rna_seq_data <- rna_seq_data[, top_genes]
    
    outfile <- paste(data_dir, output_filtered_file_name, "_", outcome, sep = "")
    cat("Saving filtered data to `", outfile, "`.\n", sep = "")   
    save(rna_seq_data, file = outfile, compress = TRUE)
    
    outfile <- paste(data_dir, output_pvalues_file_name, "_", outcome, sep = "")
    cat("Saving Cox PH p-values to `", outfile, "`.\n\n", sep = "")  
    save(p_values, file = outfile, compress = TRUE)
  }
  
  for(file in files) {
    filter_data(file[[1]], file[[2]], file[[3]])
    gc()
  }
}
