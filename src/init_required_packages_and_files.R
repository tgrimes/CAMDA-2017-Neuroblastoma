init <- function() {
  
  # Get working directory. Required source files should be in the `src/` folder.
  dir <- getwd()
  src <- "src" # Folder name containing all source files.
  
  
  # Create any missing directories.
  dir.create("bin")
  dir.create("data/processed", recursive = TRUE)
  dir.create("docs")
  dir.create("logs/create_output", recursive = TRUE)
  dir.create("logs/ensemble", recursive = TRUE)
  dir.create("logs/rna_seq_processing", recursive = TRUE)
  dir.create("logs/patient_data_processing", recursive = TRUE)
  dir.create("logs/main", recursive = TRUE)
  dir.create("output/ensemble", recursive = TRUE)
  dir.create("output/figures/km", recursive = TRUE)
  dir.create("output/figures/measures", recursive = TRUE)
  dir.create("output/results", recursive = TRUE)
  dir.create("output/tables", recursive = TRUE)
  
  #Required files:
  packages <- c("dplyr", "survival", "maxLik", "truncnorm", "spls", "glmnet",
                "foreach", "doParallel", "RankAggreg", "survminer", "boot", 
                "km.ci", "Hmisc", "gdata")
  file_names <- c("ensemble.R", 
                  "fit.R", 
                  "imputation.R", 
                  "loss_functions.R",
                  "miscellaneous.R", 
                  "plots.R")
  file_directories <- c("models", 
                        "output_scripts",
                        "rna_seq_data_scripts", 
                        "patient_data_scripts")
  
  
  # Install any missing packages.
  if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
  }
  
  
  # Load any packages that are not already loaded.
  if(!all(packages %in% loadedNamespaces())) {
    missing_packages <- setdiff(packages, loadedNamespaces())
    sapply(missing_packages, require, character.only = TRUE)
  }
  
  
  # Load scripts from source files.
  for(file in file_names) {
    source(paste(dir, src, file, sep = "/"))
  }
  
  
  # Load scripts from other directories.
  source_dir <- function(dir) {
    for (name in list.files(dir)) {
      source(file.path(dir, name))
    }
  }
  for(directory in file_directories) {
    source_dir(paste(dir, src, directory, sep = "/"))
  }
}

init()
warnings()
rm(init)