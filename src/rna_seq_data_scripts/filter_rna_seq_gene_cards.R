filter_using_gene_cards <- function(data_dir = "data/processed/",
                                    files,
                                    gene_cards_file = "data/gene_cards") {
  load(gene_cards_file)
  symbols <- gene_cards$gene_symbol
  
  cat(format(Sys.time()), "Filtering RNA-Seq Data...\n")
  cat("Running `filter_using_gene_cards()`. ",
      "Gene cards loaded from ", gene_cards_file, ".\n", sep = "")
  
  filter_data <- function(rna_seq_file_name, filtered_file_name) {
    cat("Loading RNA-seq data from `", paste(data_dir, rna_seq_file_name, sep = ""),
        "`.\n", sep = "")
    rna_seq_data <- load_RData_file(paste(data_dir, rna_seq_file_name, sep = ""))
    names <- colnames(rna_seq_data)
    
    # Find matches in gene cards.
    exact_match <- NULL
    for(symbol in symbols) {
      symbol_search <- paste("(^|(and))", symbol, "(\\.|(and)|$)", sep = "")
      temp <- grep(symbol_search, names, TRUE)
      if(length(temp) > 0) {
        exact_match <- c(exact_match, temp)
      }
    }
    exact_match <- unique(exact_match)
    cat("\t", length(exact_match), "of", ncol(rna_seq_data), 
        "column names found in gene cards.\n")
    
    if(length(exact_match) > 0) {
      rna_seq_data <- rna_seq_data[, exact_match]
    }
    
    cat("Saving filtered saved to `", 
        paste(data_dir, filtered_file_name, sep = ""), "`.\n", sep = "")
    save(rna_seq_data, file = paste(data_dir, filtered_file_name, sep = ""))
    
    rm(rna_seq_data)
  }
  
  for(file in files) {
    filter_data(file$rna_seq_file_name, file$filtered_file_name)
    gc()
  }
}