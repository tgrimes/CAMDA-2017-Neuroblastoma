generate_figures_and_tables <- function(results_dir = "output/results/",
                                        figures_dir = "output/figures/",
                                        tables_dir = "output/tables/") {
  results_list <- list.files(results_dir)
  with_ensemble <- grep("ensemble", results_list)
  results_list <- as.list(paste(results_dir, 
                                results_list[with_ensemble], sep = ""))
  
  # Make KM figures.
  make_km_plots(results_list, figures_dir, cutoff_years = 2)
  
  # Make KM table.
  make_km_table(results_list, tables_dir, cutoff_years = 2)
  
  # Make performance measure plots.
  make_performance_plots(results_list, figures_dir, B = 1000)
}