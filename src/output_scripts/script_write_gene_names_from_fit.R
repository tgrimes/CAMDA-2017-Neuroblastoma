get_gene_sets <- function(results_dir = "output/results/",
                          tables_dir = "output/tables/") {

  # # Optional: get entrez id's for genes from GeneCards database.
  # # These can be used to set the background for pathway enrichment.
  # # Results saved as "background_gene_cards.txt"
  # library(biomaRt)
  # symbol_to_entrez <- function(df) {
  #   gene_info <- getBM(attributes = c("hgnc_symbol", "entrezgene"), 
  #                      mart = useMart(biomart = "ensembl",
  #                                     dataset = "hsapiens_gene_ensembl"))
  #   gene_info %>%
  #     group_by(hgnc_symbol) %>%
  #     summarise(entrezgene = first(entrezgene)) ->
  #     gene_info
  #   
  #   return(left_join(df, gene_info))
  # }
  # load("data/gene_cards")
  # colnames(gene_cards)[1] <- "hgnc_symbol"
  # gene_cards <- symbol_to_entrez(gene_cards)
  # gene_cards <- gene_cards[!is.na(gene_cards$entrezgene), ]
  # write.table(unique(gene_cards$entrezgene),
  #             file = paste0(tables_dir, "gene_sets/background_gene_cards.txt"),
  #             quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  
  
  
  results_list <- list.files(results_dir)
  with_ensemble <- grep("ensemble", results_list)
  results_list <- as.list(paste(results_dir, 
                                results_list[-with_ensemble], sep = ""))
  df <- tibble(beta_hat = numeric(),
               outcome = character(),
               model = character(),
               method = character())
  
  for(results_file in results_list) {
    load(results_file)
    outcome <- attr(all_results, "survival_type")
    outcome_gene_set <- NULL
    
    # File name for a given measure.
    get_file_name <- function(outcome, method, model) {
      paste(tables_dir, "gene_sets/", outcome, "_",
            method, "_", model, "_gene_set.txt", sep = "")
    }
    
    par(mfrow = c(3, 4))
    for(model in c("genes", "transcripts", "transcripts_introns")) {
      for(method in c("elnet_model")) {
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
        
        
        write.table(paste(genes, coef_sign, beta_hat, sep = "\t\t"), 
                    file = get_file_name(outcome, method, model), 
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        if(length(beta_hat) > 0 & p < length(fit$beta_hat) * 0.8) {
          outcome_gene_set <- c(outcome_gene_set, genes)
          df <- rbind(df,
                      tibble(beta_hat = beta_hat,
                             outcome = outcome,
                             model = model,
                             method = method))
        }
      }
    }
    
    write.table(unique(outcome_gene_set), 
                file = paste0(tables_dir, "gene_sets/", outcome, "_gene_set.txt"), 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  return(df)
}

df <- get_gene_sets()
df %>%
  group_by(outcome, model) %>%
  mutate(m = length(beta_hat),
         max = max(beta_hat),
         min = min(beta_hat)) %>%
  ungroup() ->
  df 

df %>%
  ggplot(aes(x = beta_hat)) +
  facet_grid(model ~ outcome, scales = "free") +
  geom_histogram() + 
  geom_text(aes(x = 0, y = m, label = m)) +
  geom_point(aes(x = max, y = -5), size = 2, shape = 17, color = "red") +
  geom_point(aes(x = min, y = -5), size = 2, shape = 17, color = "red")


  
