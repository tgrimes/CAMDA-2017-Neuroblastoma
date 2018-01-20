load_RData_file <- function(file_name){
  #loads an RData file, and returns it.
  load(file_name)
  get(ls()[ls() != "file_name"])
}

load_RData_from_list <- function(file_name_list) {
  file_name_list <- as.list(file_name_list)
  n <- length(file_name_list)
  output_list <- vector("list", n)
  for(i in 1:n) {
    output_list[[i]] <- load_RData_file(file_name_list[[i]])
  }
  return(output_list)
}

#' @param data_file_list a string containing the file names for clincal
#'   predictor and/or RNA-Seq data to load. Can be a list or vector. 
get_predictors <- function(data_file_list) {
  if(is.null(data_file_list)) return(NULL)
  do.call("cbind", load_RData_from_list(data_file_list))
}

timestamp <- function() {
  format(Sys.time(), "%Y_%m_%d_%Hh_%Mm_%Ss")
}

open_log <- function(log_file, append = TRUE) {
  sink(log_file, append = append)
  cat("Log opened", format(Sys.time()), "\n\n")
}

close_log <- function() {
  cat("\n\nLog closed", format(Sys.time()))
  cat("\n\n***\nNotes:\n")
  sink()
}