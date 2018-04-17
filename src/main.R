# CAMDA: Neuroblastoma
cat("Initializing required packages and files.\n")
source("src/init_required_packages_and_files.R")

log_file <- paste("logs/main/", timestamp(), ".txt", sep = "")
open_log(log_file)
cat(format(Sys.time()), "Running `main.R`...\n")

# # Optionally: download raw data, store in data/ folder, and rerun 
# #             preprocessing steps:
# cat("Processing the patient data using `patient_data.R`.\n")
# create_patient_data(output_data_dir = "data/",
#                     raw_data_dir = "data/raw/")
# 
# cat("Processing the raw data using `run_processing.R`.\n")
# run_processing(processed_data_dir = "data/processed/",
#                raw_data_dir = "data/raw/",
#                threshold = 80,
#                gene_cards_file = "data/gene_cards")

seed <- 100


cat("Begin fitting models.\n")

load("data/patient_data")

# # Plot the patient data.
# plot_patient_data(patient_data)
# plot_survival_days(patient_data)

model_list = list(pls_model = pls_model,
                  spls_model = spls_model, 
                  lasso_model = lasso_model, 
                  elnet_model = elnet_model)

# Index for training samples.
training_index <- which(patient_data$training == "Training")

for(survival_type in c("OS", "EFS")) {
  # Reponse variable (time and censoring).
  if(tolower(survival_type) == "os") {
    y <- patient_data$os_days
    observed <- patient_data$death_from_disease == "Y"
  } else if(tolower(survival_type) == "efs") {
    y <- patient_data$efs_days
    observed <- patient_data$progression == "Y"
  } else {
    stop("Invalid survival type. Must be either 'OS' or 'EFS'")
  }
  
  # Initialize results list.
  all_results <- NULL
  
  x_file_list <- "data/processed/clinical_predictors"
  x <- get_predictors(x_file_list) # rna_data will be NULL.
  title <- paste(survival_type, ": Clinical Only", sep = "")
  result <- run_using(x, y, observed, training_index, 
                      list(simple_model = simple_model), seed)
  all_results <- c(all_results,
                   list(clinical_only = c(result, 
                                          list(x_file_list = x_file_list))))
  rm(x)
   
  x_file_list <- c(paste("data/processed/genes_filtered_coxph_", 
                         survival_type, sep = ""))
  x <- get_predictors(x_file_list)
  result <- run_using(x, y, observed, training_index, model_list, seed)
  all_results <- c(all_results,
                   list(genes = c(result, 
                                  list(x_file_list = x_file_list))))
  rm(x)
  
  x_file_list <- c(paste("data/processed/transcripts_filtered_coxph_", 
                         survival_type, sep = ""))
  x <- get_predictors(x_file_list)
  result <- run_using(x, y, observed, training_index, model_list, seed)
  all_results <- c(all_results,
                   list(transcripts = c(result, 
                                        list(x_file_list = x_file_list))))
  rm(x)
  
  x_file_list <- c(paste("data/processed/introns_filtered_coxph_", 
                         survival_type, sep = ""))
  x <- get_predictors(x_file_list)
  result <- run_using(x, y, observed, training_index, model_list, seed)
  all_results <- c(all_results,
                   list(introns = c(result, 
                                    list(x_file_list = x_file_list))))
  rm(x)
  
  x_file_list <- c(paste("data/processed/transcripts_filtered_coxph_", 
                         survival_type, sep = ""),
                   paste("data/processed/introns_filtered_coxph_", 
                         survival_type, sep = ""))
  x <- get_predictors(x_file_list)
  result <- run_using(x, y, observed, training_index, model_list, seed)
  all_results <- c(all_results,
                   list(transcripts_introns = c(result, 
                                                list(x_file_list = x_file_list))))
  rm(x)
  
  attr(all_results, "survival_type") <- survival_type
  attr(all_results, "seed") <- seed
  attr(all_results, "training_index") <- training_index
  save(all_results, file = paste("output/results/all_results", survival_type,
                                 timestamp()))
}

cat("Create ensembles and append to results.\n")
results_list <- list.files("output/results")
ensemble(paste("output/results/", results_list[[1]], sep = ""))
ensemble(paste("output/results/", results_list[[2]], sep = ""))

generate_figures_and_tables("output/results/",
                            "output/figures/",
                            "output/tables/")

close_log()
