read_meta_GSE62564 <- function(raw_data_dir = "data/raw/") {
  
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
                            c("!Sample_characteristics_ch1",
                              "!Sample_characteristics_ch2"))
    given_variables <- meta_data[, char_columns]
    given_variables[] <- lapply(given_variables, as.character)
    vars_of_interest <- c("training/validation: ", "Sex: ", "age: ", "efs day: ",
                          "efs bin: ", "os day: ", "os bin: ", "high risk: ", 
                          "a_efs_all: ", "b_os_all: ", "c_sex_all: ",
                          "d_fav_all: ", "e_efs_hr: ", "f_os_hr: ", "tissue: ")
    var_names <- sapply(vars_of_interest, function(s) gsub(": ", "", s))
    variables <- as.data.frame(matrix(NA, nrow = nrow(given_variables), 
                                      ncol = length(vars_of_interest)))
    colnames(variables) <- tolower(var_names)
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
    levels(variables[, 2]) <- c("F", "M", "NA")
    #Event free survival binary:
    variables[, 5] <- as.factor(variables[, 5])
    levels(variables[, 5]) <- c("N", "Y")
    #Overall survival binary: (death N/Y -> surival Y/N)
    variables[, 7] <- as.factor(variables[, 7])
    levels(variables[, 7]) <- c("Y", "N")
    variables[, 7] <- factor(variables[, 7], levels(variables[, 7])[c(2, 1)])
    #High risk: replace na values with 0 (not high risk).
    variables[is.na(variables[, 8]), 8] <- 0
    variables[, 8] <- as.factor(variables[, 8])
    levels(variables[, 8]) <- c("N", "Y")
    #Favoriable: 
    variables[is.na(variables[, 12]), 12] <- 2
    variables[, 12] <- as.factor(variables[, 12])
    levels(variables[, 12]) <- c("Favorable", "Unfavorable", "Neither")
    #Tissue:
    variables[, 15] <- as.factor(variables[, 15])
    levels(variables[, 15]) <- c("Neuroblastoma")
    
    variable_names <- c("Training set", "Sex", "Age (days)", "Event free survival (days)",
                        "Event free survival", "Overall survival (days)", 
                        "Overall survival", "High risk", 
                        "Event free survival (all)", "Overall survival (all)", "Sex (all)",
                        "Favorable", "Event free survival (HR)", "Overall survival (HR)", 
                        "Tissue")
    variable_description <- c("",
                              "",
                              "",
                              "EFS day: number of days for event free survival",
                              "EFS bin: binary classification of event free survival",
                              "OS day: number of days for overall survival",
                              "OS bin: binary classification of overall survival",
                              "High Risk: Indicating whether a sample belongs to high risk group or not",
                              "A_EFS_All: binary class label for event free survival for all samples",
                              "B_OS_All: binary class label for overall survival for all samples",
                              "C_SEX_All: binary class label for sex",
                              "D_FAV_All: binary class label for favorable and unfavorable samples",
                              "E_EFS_HR: binary class label for event free survival of High Risk group",
                              "F_OS_HR: binary class label for overall survival of High Risk group",
                              "")
    return(list(variables = variables, 
                names = variable_names,
                descriptions = variable_description))
  }
  
  all_file_data <- list()
  FLAG <- rep(FALSE, 12)
  
  #RNA-Seq Magic 
  folder <- paste(raw_data_dir, "RNA-SEQ_Magic_AceView_GSE49711/", sep = "")
  file_name <- paste(folder, "GSE62564_series_matrix.txt", sep = "")
  all_file_data[[1]] <- get_file_data(file_name)
  
  #Check that there are 498 rows in the data.
  FLAG[1] <- !(nrow(all_file_data[[1]]$meta_data) == 498)
  
  samples <- as.character(all_file_data[[1]]$meta_data[, 1])
  sample_id <- gsub(" \\[2\\]", "", samples)
  samples <- as.numeric(gsub("SEQC_NB", "", sample_id))
  
  #Check for any missing samples and verify samples are in order.
  FLAG[2] <- !(length((1:498)[!(1:498 %in% samples)]) == 0)
  FLAG[3] <- !all((1:498 - samples) == 0)
  
  data <- clean_meta_data(all_file_data[[1]]$meta_data)
  variables <- data$variables
  var_names <- data$names
  numeric_vars <- which(sapply(variables, is.numeric))
  summary(variables)
  
  #Check if high risk(8), e_efs_hr(13), and f_os_hr(14) are missing together (each has 65.66% NA).
  na_variables <- apply(variables, 2, is.na)
  FLAG[4] <- !all((variables[, 8] == 'N') == na_variables[, 13])
  FLAG[5] <- !all(na_variables[, 13] == na_variables[, 14])
  
  #Check that e_efs_hr(13) and f_os_hr(14) are just subsetting the HR samples.
  FLAG[6] <- !all((variables$e_efs_hr[!na_variables[, 13]] + 1) == 
                    as.numeric(variables$`efs bin`[!na_variables[, 13]]))
  FLAG[7] <- !all((variables$f_os_hr[!na_variables[, 14]] + 1) != 
                    as.numeric(variables$`os bin`[!na_variables[, 14]]))
  
  #Check if favoriable(12) and high risk(8) are mutually exclusive.
  all((!na_variables[, 8] & !na_variables[, 12]) == FALSE)
  #They are not; consider sample 492.
  #variables[492, ] #Is both favoriable and high risk.
  
  #Check if c_sex_all and sex are the same.
  FLAG[8] <- !all((as.numeric(variables$sex) == 3) == na_variables[, 11])
  FLAG[9] <- !all((as.numeric(variables$sex[!na_variables[, 11]]) == 2) == 
                    (variables[!na_variables[, 11], 11] == 0))
  FLAG[10] <- !all((as.numeric(variables$sex[!na_variables[, 11]]) == 1) == 
                    (variables[!na_variables[, 11], 11] == 1))
  
  #Check for other equivalencies.
  FLAG[11] <- !all(as.numeric(variables$`efs bin`) == (variables$a_efs_all + 1))
  FLAG[12] <- !all(as.numeric(variables$`os bin`) != (variables$b_os_all + 1))
  
  #Remove redundant variables.
  redundant_vars <- c(9, 10, 11, 13, 14, 15)
  variables <- variables[, -redundant_vars]
  var_names <- var_names[-redundant_vars]
  colnames(variables) <- c("training", "sex", "age_at_diagnosis", "efs_days",
                           "event_free_survival", "os_days", "overall_survival",
                           "high_risk", "favorable")
  
  
  id_ref <- colnames(all_file_data[[1]]$sequencing_data)[-1]
  variables <- cbind(sample_id, id_ref, variables)
  
  if(any(FLAG == TRUE)) {
    warning(paste("some flags were raised:", which(FLAG == TRUE)))
  }
  
  return(list(data = variables, 
              names = var_names,
              descriptions = data$descriptions[-redundant_vars],
              FLAG = FLAG))
}

# #Read in data from first folder.
# folder <- "aCGH_GSE25771/"
# 
# file_name <- paste(folder, "GSE25771-GPL2873_series_matrix.txt", sep = "")
# all_file_data[[1]] <- get_file_data(file_name)
# 
# file_name <- paste(folder, "GSE25771-GPL2879_series_matrix.txt", sep = "")
# all_file_data[[2]] <- get_file_data(file_name)
# 
# file_name <- paste(folder, "GSE25771-GPL4093_series_matrix.txt", sep = "")
# all_file_data[[3]] <- get_file_data(file_name)
# 
# file_name <- paste(folder, "GSE25771-GPL5477_series_matrix.txt", sep = "")
# all_file_data[[4]] <- get_file_data(file_name)
# 
# 
# sum(c(nrow(all_file_data[[1]]$meta_data), 
#       nrow(all_file_data[[2]]$meta_data),
#       nrow(all_file_data[[3]]$meta_data),
#       nrow(all_file_data[[4]]$meta_data)))
# 
# samples <- unique(c(levels(all_file_data[[1]]$meta_data[, 2]), 
#              levels(all_file_data[[2]]$meta_data[, 2]),
#              levels(all_file_data[[3]]$meta_data[, 2]),
#              levels(all_file_data[[4]]$meta_data[, 2])))
# samples <- sort(sapply(samples, function(s) as.numeric(gsub("GSM634", "", s))))
# missing_samples <- (58:190)[!(58:190 %in% samples)]
# all(samples %in% 58:190 == TRUE)
# 
# 
# 
# #Read in data from first folder.
# folder <- "aCGH_GSE35953/"
# 
# file_name <- paste(folder, "GSE35953-GPL2873_series_matrix.txt", sep = "")
# all_file_data[[5]] <- get_file_data(file_name)
# 
# sum(c(nrow(all_file_data[[5]]$meta_data)))
# 
# samples2 <- unique(c(levels(all_file_data[[5]]$meta_data[, 2])))
# samples2 <- sort(sapply(samples2, function(s) as.numeric(gsub("GSM878", "", s))))
# missing_samples2 <- (71:91)[!(71:91 %in% samples2)]
# all(samples2 %in% 71:91 == TRUE)
# 
# 
# #Read in data from first folder.
# folder <- "aCGH_GSE45480/"
# sub_folder <- "aCGH_44k_GSE45476/"
# 
# file_name <- paste(folder, sub_folder,
#                    "GSE45476-GPL2873_series_matrix.txt", sep = "")
# all_file_data[[6]] <- get_file_data(file_name)
# 
# file_name <- paste(folder, sub_folder, 
#                    "GSE45476-GPL8355_series_matrix.txt", sep = "")
# all_file_data[[7]] <- get_file_data(file_name)
# 
# file_name <- paste(folder, sub_folder, 
#                    "GSE45476-GPL16865_series_matrix.txt", sep = "")
# all_file_data[[8]] <- get_file_data(file_name)
# 
# sum(c(nrow(all_file_data[[6]]$meta_data), 
#       nrow(all_file_data[[7]]$meta_data),
#       nrow(all_file_data[[8]]$meta_data)))
# 
# samples3 <- unique(c(levels(all_file_data[[6]]$meta_data[, 2]), 
#                     levels(all_file_data[[7]]$meta_data[, 2]),
#                     levels(all_file_data[[8]]$meta_data[, 2])))
# samples3 <- sort(sapply(samples3, function(s) as.numeric(gsub("GSM1105", "", s))))
# missing_samples3 <- (159:234)[!(159:234 %in% samples3)]
# all(samples3 %in% 159:234 == TRUE)
# 
# #105k data.
# sub_folder <- "aCGH_105k_GSE45478/"
# 
# file_name <- paste(folder, sub_folder,
#                    "GSE45478_series_matrix.txt", sep = "")
# all_file_data[[9]] <- get_file_data(file_name)
# 
# sum(c(nrow(all_file_data[[9]]$meta_data)))
# samples4 <- unique(c(levels(all_file_data[[9]]$meta_data[, 2])))
# samples4 <- sort(sapply(samples4, function(s) as.numeric(gsub("GSM1105", "", s))))
# missing_samples4 <- (244:366)[!(244:366 %in% samples4)]
# all(samples4 %in% 244:366 == TRUE)
# 
# 
# 
# #Read in data from fourth folder.
# folder <- "aCGH_GSE56109/"
# 
# file_name <- paste(folder, "GSE56109-GPL2873_series_matrix.txt", sep = "")
# all_file_data[[9]] <- get_file_data(file_name)
# 
# file_name <- paste(folder, "GSE56109-GPL2879_series_matrix.txt", sep = "")
# all_file_data[[10]] <- get_file_data(file_name)
# 
# file_name <- paste(folder, "GSE56109-GPL4093_series_matrix.txt", sep = "")
# all_file_data[[11]] <- get_file_data(file_name)
# 
# sum(c(nrow(all_file_data[[9]]$meta_data), 
#       nrow(all_file_data[[10]]$meta_data),
#       nrow(all_file_data[[11]]$meta_data)))
# 
# samples4 <- unique(c(levels(all_file_data[[9]]$meta_data[, 2]), 
#                      levels(all_file_data[[10]]$meta_data[, 2]),
#                      levels(all_file_data[[11]]$meta_data[, 2])))
# samples4 <- sort(sapply(samples4, function(s) as.numeric(gsub("GSM1356", "", s))))
# missing_samples4 <- (9:48)[!(9:48 %in% samples4)]
# all(samples4 %in% 9:48 == TRUE)
# 
# 
# 
# 
# #3 samples. Consider filenames[1], filenames[4], and filenames[133]
# #Each have different ncol and nrow. 
# # Want logRatios (these are the values given in samples (*_series_matrix.csv)).
# 
# # path <- paste(dir, "/aCGH_GSE25771", sep = "")
# # filenames <- file.path(path, 
# #                        list.files(path, pattern = "*.gz"))
# # file <- gzfile(filenames[133])
# # open.connection(file)
# # keep_columns <- rep("NULL", 88)
# # keep_columns[c(9, 14)] <- NA
# # data <- read.csv(file = file, colClasses = NA,
# #                  skip = 9, header = TRUE, sep = "\t")
# # close.connection(file)
# 
# #High-resolution, genome-wide analysis of human metastatic neuroblastoma 
# # samples by array-Comparative Genomic Hybridization (aCGH)
# # - 6 files.
# 
# # path <- paste(dir, "/aCGH_GSE35953", sep = "")
# # filenames <- file.path(path, 
# #                        list.files(path, pattern = "*.gz"))
# # file <- gzfile(filenames[6])
# # open.connection(file)
# # keep_columns <- rep("NULL", 90)
# # keep_columns[c(11, 16)] <- NA
# # data <- read.csv(gzfile(filenames[1]), colClasses = keep_columns,
# #                  skip = 9, header = TRUE, sep = "\t")
# # close.connection(file)