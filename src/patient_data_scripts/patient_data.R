create_patient_data <- function() {
  output_data_dir <- "data/"
  raw_data_dir <- "data/raw/"
  
  # Initialize output to text file. This records each use of the script.
  log_file <- paste("logs/patient_data_processing/", timestamp(), ".txt", sep = "")
  sink(log_file, append = TRUE)
  
  data1 <- read_meta_GSE49711(raw_data_dir)$data
  data2 <- read_meta_GSE62564(raw_data_dir)$data
  
  #Check that variables agree between data sets.
  all(data1$sample_id == data2$sample_id)
  all(data1$training == data2$training)
  all(data1$high_risk == data2$high_risk)
  all(data1$favorable == data2$favorable)
  all(data1$death_from_disease != data2$overall_survival)
  all(data1$age_at_diagnosis == data2$age_at_diagnosis)
  all(data1$progression == data2$event_free_survival)
  
  #Note, sex does not match up in both data sets.
  #Values from data1 are used; these align with the Su et al. paper.
  na_sex <- which(data2$sex == 'NA')
  which((data1$sex[-na_sex] == 'M') != (data2$sex[-na_sex] == 'M'))
  
  #Create binary variable for age.
  age_group <- factor(data1$age_at_diagnosis > 18*30)
  levels(age_group) <- c("< 18 months", "> 18 months")
  #Use data1 variables and id_ref, efs_days, and os_days from data2.
  patient_data <- cbind(data1, data2[, c(2, 6, 8)], age_group)
  patient_data <- patient_data[, c(1:2, 12, 3, 15, 4:11, 13:14)]
  rm(data1)
  rm(data2)
  
  outfile <- paste(output_data_dir, "patient_data", sep = "")
  cat("Patient data saved to", outfile, "\n")
  save(patient_data, file = outfile)
  
  
  # Close output to text file.
  cat("\n")
  sink()
}