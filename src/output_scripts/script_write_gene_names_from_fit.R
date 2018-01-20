get_gene_names <- function() {
  # all_results_file <- "output/results/all_results_OS_2017_09_20_14_09_10"
  # all_results_file <- "output/results/all_results_OS_2017_07_13_10_00_05"
  #load("output/results/all_results_EFS_2017_09_19_09_50_01")
  #load("output/results/all_results_OS_2017_06_09_00_38_00")
  #load("output/results/all_results_OS_2017_07_13_10_00_05")
  #load("output/results/all_results_OS_2017_09_19_14_50_37")
  load("output/results/all_results_OS_2017_09_20_14_09_10")
  
  
  for(model in c("genes", "transcripts", "transcripts_introns")) {
    for(method in c("spls_model", "elnet_model")) {
      fit <- all_results[[model]]$pred[[method]]$fit
      p <- sum(abs(fit$beta_hat) > 0) # Number of non-zero coefficients.
      p
      
      ordered_coef <- order(abs(fit$beta_hat), decreasing = TRUE)[1:p]
      beta_hat <- fit$beta_hat[ordered_coef]
      coef_sign <- sign(beta_hat)
      genes <- names(fit$mean_x)[ordered_coef]
      genes <- gsub("\\..*", "", genes)
      genes <- gsub("and", "\n", genes)
      genes <- unique(genes)
      chr_index <- grep("chr", genes)
      chr <- genes[chr_index]
      chr <- gsub(":.*", "", chr)
      if(length(chr_index) > 0) {
        genes <- genes[-chr_index]
      }
      write.table(genes, 
                  file = paste("gene_names_", 
                               model, "_", 
                               method, ".txt", sep = ""), 
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
      write.table(paste(genes, coef_sign, beta_hat, sep = "\t\t"), 
                  file = paste("gene_summary_", 
                               model, "_", 
                               method, ".txt", sep = ""), 
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }
  
}