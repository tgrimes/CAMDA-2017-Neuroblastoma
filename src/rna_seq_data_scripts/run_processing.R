run_processing <- function(processed_data_dir = "data/processed/",
                           raw_data_dir = "data/raw/",
                           threshold = 80,
                           gene_cards_file = "data/gene_cards") {
  training <- NULL
  
  # Initialize output to text file. This records the use of each script.
  log_file <- function(header) {
    paste("logs/rna_seq_processing/", header, "_", timestamp(), ".txt", sep = "")
  }
  
  open_log(log_file("extract"), append = TRUE)
  extract_rna_seq_data(processed_data_dir, raw_data_dir)
  close_log()
  
  open_log(log_file("clinical"), append = TRUE)
  filter_clinical_variables()
  close_log()
  
  files <- list(genes = 
                  list(rna_seq_file_name = "genes_ncbi", 
                       filtered_file_name = 
                         paste("genes_filtered_zeroes_", threshold, sep = "")),
                transcripts = 
                  list(rna_seq_file_name = "transcripts_ncbi",
                       filtered_file_name = 
                         paste("transcripts_filtered_zeroes_", threshold, sep = "")),
                introns = 
                  list(rna_seq_file_name = "introns_raw", 
                       filtered_file_name = 
                         paste("introns_filtered_zeroes_", threshold, sep = "")))
  open_log(log_file("zeroes"))
  filter_rna_seq_zeroes(processed_data_dir, files, threshold / 100, training)
  close_log()
  
  # Optionally: This part can be rewritten to skip filtering with gene cards database.
  files <- list(genes =
                  list(rna_seq_file_name =
                         paste("genes_filtered_zeroes_", threshold, sep = ""),
                       filtered_file_name = "genes_filtered_gene_cards"),
                transcripts =
                  list(rna_seq_file_name =
                         paste("transcripts_filtered_zeroes_", threshold, sep = ""),
                       filtered_file_name = "transcripts_filtered_gene_cards"))
  filter_using_gene_cards(processed_data_dir, files, gene_cards_file)
  
  # Filter using COX PH with/without clinical variables as covariates. 
  files <- list(genes = 
                  list(input_rna_seq_file_name = "genes_filtered_gene_cards",
                         #paste("genes_filtered_zeroes_", threshold, sep = ""), 
                       output_filtered_file_name = "genes_filtered_coxph",
                       output_pvalues_file_name = "genes_pvalues_coxph"),
                transcripts = 
                  list(input_rna_seq_file_name = "transcripts_filtered_gene_cards",
                         #paste("transcripts_filtered_zeroes_", threshold, sep = ""),
                       output_filtered_file_name = "transcripts_filtered_coxph",
                       output_pvalues_file_name = "transcripts_pvalues_coxph"),
                introns =
                  list(input_rna_seq_file_name =
                         paste("introns_filtered_zeroes_", threshold, sep = ""),
                       output_filtered_file_name = "introns_filtered_coxph",
                       output_pvalues_file_name = "introns_pvalues_coxph"))
  
  open_log(log_file("coxph_OS"))
  filter_rna_seq_coxph(processed_data_dir, files, outcome = "OS")
  close_log()
  
  open_log(log_file("coxph_EFS"))
  filter_rna_seq_coxph(processed_data_dir, files, outcome = "EFS")
  close_log()
}