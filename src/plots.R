# library(Hmisc) #For rcorr().

plot_patient_data <- function(data, save_plot = TRUE) {
  original_mfrow <- par()$mfrow
  
  if(save_plot) {
    png(paste("output/figures/patient_data_histograms.png", sep = ""), 
        2000, 2000*(9/10), res = 350)
  }
  par(mfrow = c(3, 4))
  for(i in 5:ncol(data)) {
    var_name <- colnames(data)[i]
    if(is.numeric(data[, i])) {
      hist(data[, i], 
           main = paste(100 * round(sum(sapply(data[, i], is.na)) / 
                                      length(data[, i]), 4),
                        "% na", sep = ""),
           xlab = paste("[", i, "] ", var_name, sep = ""))
    } else {
      plot(data[, i], 
           main = paste(100 * round(sum(sapply(data[, i], is.na)) / 
                                      length(data[, i]), 4),
                        "% na", sep = ""),
           xlab = paste("[", i, "] ", var_name, sep = ""))
    }
  }
  if(save_plot) {
    dev.off()
  }
  
  # # Spearman's correlation between each pair of patient variables.
  # panel.dist <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  #   usr <- par("usr"); on.exit(par(usr))
  #   par(usr = c(0, 1, 0, 1))
  #   r <- round(rcorr(x, y, type = "spearman")[[1]][1, 2], 3)
  #   txt <- format(c(r, 0.123456789), digits=digits)[1]
  #   txt <- paste(prefix, txt, sep="")
  #   if(missing(cex.cor)) cex.cor <- 5/strwidth(txt)
  #   text(0.5, 0.5, txt, cex = max(sqrt(1 * abs(r)*4), 0.75))
  # }
  # pairs(~ ., data = data[, -(1:4)],
  #       upper.panel = panel.dist,
  #       pch = 16, main = "Spearman correlations")
  
  if(save_plot) {
    png(paste("output/figures/ef_survival_boxplots.png", sep = ""), 
        2000, 2000*(9/10), res = 350)
  }
  par(mfrow = c(2, 3))
  plot(data$progression, data$efs_days,
       xlab = "progression", ylab = "event free survival")
  plot(data$high_risk, data$efs_days,
       xlab = "high risk", ylab = "event free survival")
  plot(data$favorable, data$efs_days,
       xlab = "favorable", ylab = "event free survival")
  plot(data$age_group, data$efs_days,
       xlab = "age", ylab = "event free survival")
  plot(data$inss_stage, data$efs_days,
       xlab = "INSS", ylab = "event free survival")
  plot(data$mycn_status, data$efs_days,
       xlab = "mycn", ylab = "event free survival")
  if(save_plot) {
    dev.off()
  }
  
  if(save_plot) {
    png(paste("output/figures/os_survival_boxplots.png", sep = ""), 
        2000, 2000*(9/10), res = 350)
  }
  par(mfrow = c(2, 3))
  plot(data$death_from_disease, data$os_days, 
       xlab = "death from disease", ylab = "overall survival")
  plot(data$high_risk, data$os_days, 
       xlab = "high risk", ylab = "overall survival")
  plot(data$favorable, data$os_days, 
       xlab = "favorable", ylab = "overall survival")
  plot(data$age_group, data$os_days, 
       xlab = "age", ylab = "overall survival")
  plot(data$inss_stage, data$os_days, 
       xlab = "INSS", ylab = "overall survival")
  plot(data$mycn_status, data$os_days, 
       xlab = "mycn", ylab = "overall survival")
  if(save_plot) {
    dev.off()
  }
  
  if(save_plot) {
    png(paste("output/figures/survival_censoring_histograms.png", sep = ""), 
        2000, 2000*(9/10), res = 350)
  }
  par(mfrow = c(2, 2))
  hist(data$os_days[data$death_from_disease == 'Y'], 
       main = "No censoring", xlab = "overall survival")
  hist(data$os_days[data$death_from_disease == 'N'], 
       main = "Right-censored", xlab = "overall survival")
  hist(data$efs_days[data$progression == 'Y'], 
       main = "No censoring", xlab = "event-free survival")
  hist(data$efs_days[data$progression == 'N'], 
       main = "Right-censored", xlab = "event-free survival")
  if(save_plot) {
    dev.off()
  }
  
  par(mfrow = original_mfrow)
}

plot_survival_days <- function(data, save_plot = TRUE) {
  original_mfrow <- par()$mfrow
  
  if(save_plot) {
    png(paste("output/figures/survival_censoring_line_plot.png", sep = ""), 
        2000, 2000*(9/10), res = 350)
  }
  
  par(mfrow = c(1, 3))
  max_age <- max(data$age_at_diagnosis + data$os_days)
  n <- nrow(data)
  plot(x = seq(0, max_age, length.out = n), y = 1:n, 
       xlab = "patient age (in days)", ylab = "", type = 'n', yaxt = "n",
       main = "OS \n death from disease (red) \n right-censored (black)")
  death <- data$death_from_disease == 'Y'
  start <- data$age_at_diagnosis
  end <- start + data$os_days
  for(i in 1:n) {
    lines(x = c(start[i], end[i]), y = rep(i, 2), 
          type = "l", lwd = 1, 
          col = ifelse(death[i], "red", "black"))
  }
  
  os_less_5 <- data$os_days < 365 * 5 + 1
  plot(x = seq(0, max_age, length.out = n), y = 1:n,
       xlab = "patient age (in days)", ylab = "", type = 'n', yaxt = "n",
       main = paste("OS < 5 years \n",
                    "death (", sum(os_less_5 & death), "/", sum(os_less_5), ") (red) \n",
                    "censored (", sum(os_less_5 & !death), "/", sum(os_less_5), ") (black)", sep = ""))
  start <- data$age_at_diagnosis
  end <- start + data$os_days
  for(i in 1:n) {
    if(os_less_5[i]) {
      lines(x = c(start[i], end[i]), y = rep(i, 2), 
            type = "l", lwd = 1, 
            col = ifelse(death[i], "red", "black"))
    }
  }
  
  os_more_5 <- data$os_days >= 365 * 5 + 1
  plot(x = seq(0, max_age, length.out = n), y = 1:n,
       xlab = "patient age (in days)", ylab = "", type = 'n', yaxt = "n",
       main = paste("OS >= 5 years \n",
                    "death (", sum(os_more_5 & death), "/", sum(os_more_5), ") (red) \n",
                    "censored (", sum(os_more_5 & !death), "/", sum(os_more_5), ") (black)", sep = ""))
  start <- data$age_at_diagnosis
  end <- start + data$os_days
  for(i in 1:n) {
    if(os_more_5[i]) {
      lines(x = c(start[i], end[i]), y = rep(i, 2), 
            type = "l", lwd = 1, 
            col = ifelse(death[i], "red", "black"))
    }
  }
  
  if(save_plot) {
    dev.off()
  }
  
  par(mfrow = original_mfrow)
}
