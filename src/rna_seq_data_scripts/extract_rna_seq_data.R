########################################################################
#
# Extracts the RNA-Seq data from the raw files. No filtering is done.
#
########################################################################

extract_rna_seq_data <- function(output_data_dir = "data/processed/", 
                                 raw_data_dir = "data/raw/") {
  cat(format(Sys.time()), "RNA-Seq Data:\n")
  cat("Running `script_extract_rna_seq_data.R`\n")
  
  # Load RNA-Seq data.
  # There are three layers: genes, transcripts, and exon junctions.
  file_name <- paste(raw_data_dir,
                     "RNA-SEQ_Magic_AceView_GSE49711/",
                     "GSE49711_SEQC_NB_MAV_G_log2.20121127.txt", sep = "")
  g_data <- read.delim(file_name, stringsAsFactors = FALSE)
  g_data_ncbi <- g_data[g_data[, 2] != "", ]
  genes <- t(g_data_ncbi[, -(1:9)])
  colnames(genes) <- g_data_ncbi[, 1]
  genes <- as.data.frame(do.call(cbind, by(t(genes), 
                                           INDICES = colnames(genes), 
                                           FUN = colMeans)))
  save(genes, file = paste(output_data_dir, "genes_ncbi", sep = ""))
  cat("NCBI Gene data saved as `genes_ncbi` in", output_data_dir, "folder.\n")
  rm("g_data", "g_data_ncbi", "genes")
  gc()
  
  file_name <- paste(raw_data_dir,
                     "RNA-SEQ_Magic_AceView_GSE49711/",
                     "GSE49711_SEQC_NB_MAV_T_log2.20121127.txt", sep = "")
  t_data <- read.delim(file_name, stringsAsFactors = FALSE)
  t_data_ncbi <- t_data[!is.na(t_data[, 2]), ]
  transcripts <- t(t_data_ncbi[, -(1:4)])
  colnames(transcripts) <- t_data_ncbi[, 1]
  transcripts <- as.data.frame(do.call(cbind, by(t(transcripts), 
                                                 INDICES = colnames(transcripts), 
                                                 FUN = colMeans)))
  save(transcripts, file = paste(output_data_dir, "transcripts_ncbi", sep = ""))
  cat("NCBI Transcript data saved as `transcripts_ncbi` in", 
      output_data_dir, "folder.\n")
  rm("t_data", "t_data_ncbi", "transcripts")
  gc()
  
  file_name <- paste(raw_data_dir,
                     "RNA-SEQ_Magic_AceView_GSE49711/",
                     "GSE49711_SEQC_NB_MAV_J_log2.20121127.txt", sep = "")
  j_data <- read.delim(file_name, stringsAsFactors = FALSE)
  gc()
  introns <- t(j_data[, -(1:6)])
  intron_names <- j_data[, 1]
  rm(j_data)
  gc()
  colnames(introns) <- intron_names
  introns <- as.data.frame(do.call(cbind, by(t(introns), 
                                             INDICES = colnames(introns), 
                                             FUN = colMeans)))
  save(introns, file = paste(output_data_dir, "introns_raw", sep = ""))
  cat("Intron data saved as `introns_raw` in", output_data_dir, "folder.\n")
  rm("introns")
  gc()
}