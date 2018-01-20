filter_clinical_variables <- function(
    keep_variables = c("age_group", "sex", "mycn_status"),
    file_name = "clinical_predictors",
    patient_data_file = "data/patient_data") {
  
  data_dir <- "data/processed/"
  
  cat("Running `script_save_clinical_predictors.R`\n")
  
  cat("Patient clinical variables loaded from:", patient_data_file, "\n")
  load(patient_data_file)
  
  # Clinical data from patient data. --------------------------------------------
  clinical_data <- cbind(patient_data$age_group == "> 18 months",
                         patient_data$sex == "M",
                         patient_data$mycn_status == "Amplified")
  clinical_data[which(is.na(clinical_data[, 3])), 3] <- FALSE # NA mycn status.
  clinical_data <- clinical_data * 1
  colnames(clinical_data) <- c("age_group", "sex", "mycn_status")
  cat("Including clinical_data: age group, sex, and MYCN status.\n")
  
  outfile <- paste(data_dir, "clinical_predictors", sep = "")
  cat("Clinical data saved to", outfile, "\n")
  save(clinical_data, file = outfile)
}