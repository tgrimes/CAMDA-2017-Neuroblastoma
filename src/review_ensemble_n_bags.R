# See how enseble performs with B = 5, 15, ... 50.

source("src/init_required_packages_and_files.R")

for(n_bags in c(5, 10, 20, 30, 50)) {
  for(seed in c(1, 2, 3)) {
    results_list <- list.files("output/results")
    input_dir <- "output/results/"
    output_dir <- "output/results_bags/"
    ensemble(paste(input_dir, results_list[[2]], sep = ""),
             output_dir,
             n_bags,
             seed)
  }
}


results_list <- list.files(output_dir)
for(results_file in results_list[9]) {
  load(paste0(output_dir, results_file))
  outcome <- attr(all_results, "survival_type") 
  
  # File name for a given measure.
  get_file_name <- function(method_name) {
    paste(figures_dir, "km/", 
          outcome, "_LPS_", cutoff_years, "_", method_name, ".png", sep = "")
  }
  
  # Create figure for each model.
  for(i in 1:length(all_results)) {
    temp <- all_results
    temp[-i] <- NULL # Remove all but the ith result.
    plot.km(temp, 2)
  }
}
