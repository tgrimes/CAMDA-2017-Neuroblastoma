read_meta_GSE49711 <- function(raw_data_dir = "data/raw/") {
  
  get_file_data <- function(file_name) {
    cat("Processing patient data from:", file_name, "\n")
    
    f <- file(file_name)
    #Preprocess the file.
    open(f)
    n_row_head <- 0
    while(readLines(f, n = 1) != "") {
      n_row_head <- n_row_head + 1
    }
    n_row_meta_data <- 0 
    while(readLines(f, n = 1) != "!series_matrix_table_begin") {
      n_row_meta_data <- n_row_meta_data + 1
    } 
    n_row_sequencing_data <- 0
    while(readLines(f, n = 1) != "!series_matrix_table_end") {
      n_row_sequencing_data <- n_row_sequencing_data + 1
    }
    close(f)
    
    f <- file(file_name)
    #Extract the head, meta data, and sequcing data from the file.
    open(f)
    head <- as.matrix(readLines(f, n = n_row_head), ncol = 1)
    temp <- readLines(f, n = 1)
    if(temp != "") {
      stop(paste("Unexpected line after head:", temp))
    }
    
    meta_data <- t(read.table(f, sep = "\t", skip = 0, 
                              nrows = n_row_meta_data))
    colnames(meta_data) <- meta_data[1, ]
    meta_data <- as.data.frame(meta_data[-1, ])
    temp <- readLines(f, n = 1)
    if(temp != "!series_matrix_table_begin") {
      stop(paste("Unexpected line after meta data:", temp))
    }
    
    if(n_row_sequencing_data == 1) {
      sequencing_data <- read.table(f, sep = "\t", skip = 0, header = FALSE, 
                                    colClasses = "character", nrows = n_row_sequencing_data)
      colnames(sequencing_data) <- (sequencing_data[1, ])
      sequencing_data <- as.data.frame(sequencing_data[-1, ])
    } else {
      sequencing_data <- read.table(f, sep = "\t", skip = 0, header = TRUE, 
                                    nrows = n_row_sequencing_data - 1)
    }
    
    temp <- readLines(f, n = 1)
    if(temp != "!series_matrix_table_end") {
      stop(paste("Unexpected line after sequencing data:", temp))
    }
    close(f)
    
    file_data <- list(folder = folder, file = file,
                      head = head, meta_data = meta_data, 
                      sequencing_data = sequencing_data)
  }
  
  clean_meta_data <- function(meta_data) {
    char_columns <- which(colnames(meta_data) 
                          %in% 
                            c("!Sample_characteristics_ch1"))
    given_variables <- meta_data[, char_columns]
    given_variables[] <- lapply(given_variables, as.character)
    vars_of_interest <- c("dataset: ", "Sex: ", "age at diagnosis: ", "mycn status: ",
                          "high risk: ", "inss stage: ", "class label: ", "progression: ", 
                          "death from disease: ")
    var_names <- sapply(vars_of_interest, function(s) gsub(": ", "", s))
    variables <- as.data.frame(matrix(NA, nrow = nrow(given_variables), 
                                      ncol = length(vars_of_interest)))
    colnames(variables) <- c("training", "sex", "age_at_diagnosis", "mycn_status",
                             "high_risk", "inss_stage", "favorable", "progression",
                             "death_from_disease")
    for(i in 1:nrow(variables)) {
      for(j in 1:length(vars_of_interest)) {
        temp <- which(grepl(vars_of_interest[j], given_variables[i, ]))
        if(length(temp) > 0) {
          #jth variable was found for ith sample.
          value <-  gsub(vars_of_interest[j], "", given_variables[i, temp])
          variables[i, j] <- suppressWarnings(ifelse(!is.na(as.numeric(value)), 
                                                     as.numeric(value), value))
        }
      }
    }
    
    #Training/validation:
    variables[, 1] <- as.factor(variables[, 1])
    levels(variables[, 1]) <- c("Training", "Validation")
    #Sex:
    variables[, 2] <- as.factor(variables[, 2])
    levels(variables[, 2]) <- c("F", "M")
    #Mycn:
    variables[variables[, 4] == "N/A", 4] <- NA
    variables[, 4] <- as.factor(variables[, 4])
    levels(variables[, 4]) <- c("Normal", "Amplified")
    #High risk:
    variables[, 5] <- as.factor(variables[, 5])
    levels(variables[, 5]) <- c("N", "Y")
    #INSS:
    variables[, 6] <- as.factor(variables[, 6])
    levels(variables[, 6]) <- c("1", "2", "3", "4", "4S")
    #Favoriable (class)
    variables[, 7] <- as.factor(variables[, 7])
    levels(variables[, 7]) <- c("Favorable", "Unfavorable", "Neither")
    #Progression:
    variables[, 8] <- as.factor(variables[, 8])
    levels(variables[, 8]) <- c("N", "Y")
    #Death:
    variables[, 9] <- as.factor(variables[, 9])
    levels(variables[, 9]) <- c("N", "Y")
    
    
    variable_names <- c("Training set", "Sex", "Age at diagnosis(days)", 
                        "Mycn status", "High risk", "INSS stage", "Favoriable", 
                        "Progression", "Death from disease")
    variable_description <- c("dataset: Expression data set used for classification training (1) or validation (2)",
                              "Sex: M = male; F= female",
                              "age at diagnosis: the age in days at diagnosis",
                              "mycn status: Amplification status of the MYCN proto-oncogene (amplified = 1, no amplification = 0; no information = N/A)",
                              "high risk: Clinically considered as high-risk neuroblastoma (yes=1, no= 0)",
                              "INSS stage: disease stage according to International Neuroblastoma Staging System (INSS) (1, 2, 3, 4 and 4S)",
                              "class label: Maximally divergent disease courses - unfavorable (= 1): patient died despite intensive chemotherapy, favorable (=0): patient survived without chemotharapy for at least 1000 days post diagnosis; not applicable (N/A)",
                              "progression: Occurrence of a tumor progression event (yes=1; no=0)",
                              "death from disease: Occurrence of death from the disease (yes=1; no=0)")
    return(list(variables = variables, 
                names = variable_names,
                descriptions = variable_description))
  }
  
  all_file_data <- list()
  FLAG <- rep(FALSE, 2)
  
  #RNA-Seq Magic 
  folder <- paste(raw_data_dir, "RNA-SEQ_Magic_AceView_GSE49711/", sep = "")
  file_name <- paste(folder, "GSE49711_series_matrix.txt", sep = "")
  all_file_data[[1]] <- get_file_data(file_name)
  
  nrow(all_file_data[[1]]$meta_data)
  
  #The sample id's in this file do not have " [2]" appended.
  sample_id <- as.character(all_file_data[[1]]$meta_data[, 1])
  samples <- as.numeric(gsub("SEQC_NB", "", sample_id))
  id_ref <- colnames(all_file_data[[1]]$sequencing_data)[-1]
  
  #Check for any missing samples and verify samples are in order.
  FLAG[1] <- !(length((1:498)[!(1:498 %in% samples)]) == 0)
  FLAG[2] <- !all((1:498 - samples) == 0)
  
  data <- (clean_meta_data(all_file_data[[1]]$meta_data))
  variables <- data$variables
  var_names <- data$names
  
  #Add sample id and references to variables.
  variables <- cbind(sample_id, id_ref, variables)
  
  if(any(FLAG == TRUE)) {
    warning(paste("some flags were raised:", which(FLAG == TRUE)))
  }
  
  return(list(data = variables, 
              names = var_names,
              descriptions = data$descriptions,
              FLAG = FLAG))
}