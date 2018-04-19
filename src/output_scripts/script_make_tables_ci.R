#library(survival)
#library(survminer)

get_km_ci <- function(all_results, cutoff_years = 2, hr = NULL, year = 2) {
  outcome <- attr(all_results, "survival_type")
  training_index <- attr(all_results, "training_index") # Needed for subsetting.
  rate <- 1 / 365 # Multiplier for days. 1/365 gives survival in year.
  rate_label <- "(years)"
  
  if(is.null(hr)) {
    # Subset data by high risk, favorable, INSS satage, etc. ----------------------
    load("data/patient_data")
    hr <- patient_data$high_risk[-training_index] == "Y"
    rm(patient_data)
  }
  
  output <- NULL
  
  result_names <- names(all_results)
  n_results <- length(all_results)
  
  for(k in 1:n_results) {
    results <- all_results[[k]]
    
    model_names <- gsub("_.*", "", names(results$pred))
    n_models <- length(results$pred)
    
    for(i in 1:n_models) {
      x <- results$pred[[i]]$data$estimate * rate
      y_test <- results$pred[[i]]$data$true * rate
      observed_test <- results$pred[[i]]$data$observed
      
      score <- x > cutoff_years
      pred_survival = c("< 2 years", ">= 2 years")[1 + (score * 1)]
      
      title <- paste("LPS", cutoff_years, sep = "")
      
      df <- data.frame(y = y_test[hr], 
                       observed = observed_test[hr],
                       pred_survival = pred_survival[hr])
      
      sdiff <- tryCatch(survdiff(Surv(y, observed) ~ pred_survival, data = df),
                        error = function(e) NULL)
      if(!is.null(sdiff)) {
        p_value <- 1 - pchisq(sdiff$chisq, 1)
      } else {
        p_value <- NA
      }
      
      strata_labels <- paste(c("LPS", "non-LPS"), cutoff_years, sep = "")
      sfit <- survfit(Surv(y, observed) ~ pred_survival, data = df)
      if(is.null(summary(sfit)$strata)) {
        if(length(unique(pred_survival)) != 1) warning("Expected single KM curve.")
        strata <- as.factor(rep(df$pred_survival[1], length(summary(sfit)$time)))
        levels(strata) <- ifelse(df$pred_survival[1] == "< 2 years", 
                                 strata_labels[1], strata_labels[2])
      } else {
        strata <- summary(sfit)$strata
        levels(strata) <- strata_labels
      }
      
      data_name_list <- c("clinical_only", "ensemble", "genes", "transcripts", 
                           "introns", "transcripts_introns")
      data_symbol_list <- c("null model", "ensemble", "G", "T", "I", "TI")
      
      method_name_list <- c("simple", "ensemble", "pls", "spls", "lasso", "elnet")
      method_symbol_list <- c("", "", "1", "2", "3", "4")
      
      data <- data_symbol_list[which(data_name_list %in% result_names[k])] 
      model <- method_symbol_list[which(method_name_list %in% model_names[i])]
      
      for(strat in strata_labels) {
        if(!(strat %in% unique(strata))) {
          output <- rbind(output, 
                          data.frame(prob = NA,
                                     se = NA,
                                     upper = NA,
                                     lower = NA,
                                     n = 0,
                                     strata = strat,
                                     cutoff_years = cutoff_years,
                                     model = model,
                                     data = data,
                                     outcome = outcome,
                                     p_value = NA))
        } else {
          index <- (summary(sfit)$time <= year) & (strata == strat)
          first <- last(which(index))
          index[-first] <- FALSE
          
          output <- rbind(output, 
                          data.frame(prob = summary(sfit)$surv[index],
                                     se = summary(sfit)$std.err[index],
                                     upper = summary(sfit)$upper[index],
                                     lower = summary(sfit)$lower[index],
                                     n = summary(sfit)$n[strat == unique(strata)],
                                     strata = strat,
                                     cutoff_years = cutoff_years,
                                     model = model,
                                     data = data,
                                     outcome = outcome,
                                     p_value = round(p_value, 6)))
        }
      }
      
      
    }
  }
  return(output)
}

make_km_table <- function(results_list, tables_dir, cutoff_years = 2) {
  n_results <- length(results_list)
  tab_list <- vector("list", length(results_list))
  
  for(i in 1:length(results_list)) {
    results_file <- results_list[[i]]
    load(results_file) # Loads `all_results`.
    
    outcome <- attr(all_results, "survival_type")
    training_index <- attr(all_results, "training")
    load("data/patient_data")
    hr <- patient_data$high_risk[-training_index] == "Y"
    #hr <- patient_data$inss_stage[-training_index] == "4"
    #hr <- patient_data$mycn_status[-training_index] == "Amplified"
    rm(patient_data)
    
    tab <- get_km_ci(all_results, cutoff_years, hr)
    
    # Format p-values, prob, and se.
    tab$p_value <- ifelse(tab$p_value < 0.001,
                          rep("< 0.001", nrow(tab)),
                          as.character(round(tab$p_value, 3)))
    tab$p_value <- paste("$", tab$p_value, "$", sep = "")
    tab$prob <- round(tab$prob, 2)
    tab$se <- round(tab$se, 2)
    
    # Create adjacent columns for LPS and non-LPS results.
    labels <- paste(c("LPS", "non-LPS"), cutoff_years, sep = "")
    tab <- cbind(tab[tab$strata == labels[1], c(10, 1, 2, 5)],
                 tab[tab$strata == labels[2], c(1, 2, 5, 9, 8, 11)])
    
    # Replace underscores with spaces.
    colnames(tab) <- gsub("_", " ", colnames(tab))
    tab$data <- gsub("_", " ", tab$data)
    tab$model <- gsub("_", " ", tab$model)
    
    tab_list[[i]] <- tab
    # Save the table. Formatting is for Latex.
    file_name <- paste(tables_dir, outcome, "_LPS_", cutoff_years,
                       "_summary.txt", sep = "")
    write.table(tab, file = file_name, quote = FALSE,
                sep = " \t& ", eol = " \\\\\n", row.names = FALSE)
  }
  
  # Combine the two tables.
  tab <- rbind(tab_list[[2]], tab_list[[1]])
  index <- rep(1:nrow(tab_list[[2]]), each = 2)
  index[seq(2, length(index), 2)] <- index[seq(2, length(index), 2)] + nrow(tab_list[[2]])
  tab <- tab[index, ]
  file_name <- paste(tables_dir, "all_LPS_", cutoff_years,
                     "_summary.txt", sep = "")
  write.table(tab, file = file_name, quote = FALSE,
              sep = " \t& ", eol = " \\\\\n", row.names = FALSE)
}
