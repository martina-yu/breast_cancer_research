library(writexl)

library(dplyr)
library(matrixStats)

library("survival")
library("survminer")
library(ggplot2)

### -------------- varibale definition -------------- ###

# getwd()
data_survival <- read.csv("./W2_DFA&SA_on_clinic_data/placebo_both_clinic_gene_data_for_training_after_merged.csv")
output_dir <- "/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W5_Feature_CoxPH"
p_val_threshold <- 0.05
n_total_cols <- ncol(data_survival)
n_clinical_cols <- 10
gene_columns_start <- 2
gene_columns_end <- n_total_cols - 10
gene_cols_index <- 2:(n_total_cols - n_clinical_cols)
# gene_cols_index <- 2:1000
gene_expression_data <- data_survival[, gene_cols_index] 
clinical_data <- data_survival[, -gene_cols_index]

TPM <- c(0.8, 0.9)
MAD <- c(0.2, 0.25, 0.3, 0.5, 0.6)
parameter_grid <- expand.grid(TPM = TPM, MAD = MAD)
num_patients <- nrow(gene_expression_data)
results_list <- list()
set.seed(42)

### ----------------- function definition ------------------ ###

# return_results <- function(tpm, mad){
#   num_genes <- sample(50:500, 1)
#   return(list(num_of_genes <- num_genes))
# }

uniCoxFunction <- function(gene_names, data_survival, p_val_threshold){
  results_list <- list()
  for (gene_name in gene_names) {
    
    formula_str <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", gene_name, "`")
    cox_formula <- as.formula(formula_str)
    
    cox_fit <- tryCatch({
      coxph(cox_formula, data = data_survival)
    }, error = function(e) {
      message(paste("At", gene_name, "running model wrongly：", e$message))
      return(NULL)
    })
    
    if (!is.null(cox_fit)) {
      cox_summary <- summary(cox_fit)
      
      unicox_result <- data.frame(
        gene = gene_name,
        p_value = cox_summary$coefficients[1, "Pr(>|z|)"],
        hazard_ratio = cox_summary$coefficients[1, "exp(coef)"],
        lower_95 = cox_summary$conf.int[1, "lower .95"],
        upper_95 = cox_summary$conf.int[1, "upper .95"]
      )
      results_list[[gene_name]] <- unicox_result
    }
  }
  
  final_results <- do.call(rbind, results_list)
  significant_genes <- final_results[final_results$p_value < p_val_threshold, ]
  number_of_significant_genes <- nrow(significant_genes)
  # print(paste("Totally", number_of_significant_genes, "genes are significant for P < ", p_val_threshold))
  return(list(final_results=final_results, significant_genes=significant_genes))
}

multiCoxFunction <- function(uni_sig_genes, data_survival, p_val_threshold){
  result_list <- list()
  uni_sig_genes <- uni_sig_genes$significant_genes
  significant_gene_names <- uni_sig_genes$gene
  all_predictors <- significant_gene_names
  # clinical_variable_names <- c("`Menospausal.Status`", "`Risk.Group`", "Age")
  # all_predictors <- c(significant_gene_names, clinical_variable_names)
  formula_predictors <- paste(all_predictors, collapse = " + ")
  formula_str <- paste("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~", formula_predictors)
  multivariate_formula <- as.formula(formula_str)
  multivariate_cox_fit <- coxph(multivariate_formula, data = data_survival)
  
  multicox_summary <- summary(multivariate_cox_fit)
  multicox_coef_table <- as.data.frame(multicox_summary$coefficients) ## coef exp(coef)  se(coef)   z  Pr(>|z|)
  multicox_conf_int_table <- as.data.frame(multicox_summary$conf.int) ## exp(coef) exp(-coef) lower .95 upper .95
  multicox_result <- cbind(multicox_coef_table, multicox_conf_int_table)
  
  multicox_result$p.adjusted.FDR <- p.adjust(multicox_result$`Pr(>|z|)`, method = "fdr")
  final_results <- multicox_result[multicox_result$p.adjusted.FDR < p_val_threshold, ]
  # print(paste('multicox regression < p threshold:', nrow(final_results)))
  significant_genes <- row.names(final_results)
  # print(paste('row names of significant_genes:', paste(significant_genes, collapse = ',')))
  return(list(multicox_final_results=final_results, multicox_significant_genes=significant_genes))
}

# formula_str <- paste("Surv(time, status) ~", formula_predictors)
# print(multivariate_formula)
# summary(multivariate_cox_fit)

# multivariate_significant_results <- cox_coef_table[cox_coef_table$`Pr(>|z|)` < 0.05,]
# print(nrow(multivariate_significant_results))
# 
# mul_sign_rslt <- multivariate_significant_results_fdr[order(multivariate_significant_results_fdr$`p.adjusted.FDR`), ]
# mul_sign_rslt_zero <- mul_sign_rslt[mul_sign_rslt$p.adjusted.FDR < 1e-303, ]
# all_significant_variable_names <- row.names(mul_sign_rslt_zero)


### ----------------- loop for different tpm and mad ------------------ ###

for (i in 1:nrow(parameter_grid)) {
  tpm <- parameter_grid$TPM[i]
  mad <- parameter_grid$MAD[i]
  print(paste('TPM:', tpm,'|','MAD:', mad))
  
  threshold <- tpm * num_patients
  genes_to_keep_low_expr <- colSums(gene_expression_data < 1) < threshold
  gene_expression_data_filtered1 <- gene_expression_data[, genes_to_keep_low_expr]
  print(paste('after TPM:',ncol(gene_expression_data_filtered1)))
  
  mad_values <- apply(gene_expression_data_filtered1, 2, mad, na.rm = TRUE)
  mad_sorted <- sort(mad_values, decreasing = TRUE)
  
  num_genes_to_keep <- floor(length(mad_sorted) * mad)
  
  genes_to_keep_mad <- names(mad_sorted[1:num_genes_to_keep]) ## keeped gene names
  gene_expression_data_filtered2 <- gene_expression_data_filtered1[, genes_to_keep_mad]
  
  print(paste('after TPM and MAD:',ncol(gene_expression_data_filtered2)))
  
  gene_expression_log2 <- log2(gene_expression_data_filtered2 + 1)
  data_survival <- cbind(clinical_data, gene_expression_log2)
  gene_names <- genes_to_keep_mad
  
  uni_results_genes <- uniCoxFunction(gene_names, data_survival, p_val_threshold)
  # uni_results_df <- uniCoxFunction(gene_names)
  print(paste('unicox regression kept:', nrow(uni_results_genes$significant_genes), 'genes.'))
  
  multi_results_genes <-multiCoxFunction(uni_results_genes, data_survival, p_val_threshold)
  print(paste('multicox regression kept:', length(multi_results_genes$multicox_significant_genes), 'genes.'))

  final_kept_genes <- multi_results_genes$multicox_significant_genes
  
  if (length(final_kept_genes) > 0) {
    current_result <- data.frame(
      TPM_param = tpm,
      MAD_param = mad,
      GeneName = final_kept_genes
    )
    results_list[[i]] <- current_result
    
  }
}

all_selected_genes_df <- bind_rows(results_list)
head(all_selected_genes_df)

gene_stability <- all_selected_genes_df %>%
  group_by(GeneName) %>%
  summarise(Frequency = n(), .groups = 'drop') %>%
  arrange(desc(Frequency))
head(gene_stability)

# high_freq_gene <- gene_stability[gene_stability$Frequency == 10,]
# high_freq_gene <- gene_stability[gene_stability$Frequency == 9,]
# high_freq_gene <- gene_stability[gene_stability$Frequency == 8,]
# high_freq_gene <- gene_stability[gene_stability$Frequency == 7,]
# high_freq_gene <- gene_stability[gene_stability$Frequency == 6,]
# high_freq_gene <- gene_stability[gene_stability$Frequency == 5,]
high_freq_gene <- gene_stability
high_freq_gene_name <- high_freq_gene$GeneName

sig_file_name <- "gene_stability.xlsx"
sig_file_name <- file.path(output_dir, sig_file_name)
write_xlsx(gene_stability, path = sig_file_name)


# 
# one_gene_name <- 'BHLHE40'
# formula_str_manual <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", one_gene_name, "`")
# cox_formula_manual <- as.formula(formula_str_manual)
# cox_fit_manual <- coxph(cox_formula_manual, data = data_survival)
# 
# summary_manual <- summary(cox_fit_manual)
# print(summary_manual)

### ----------------- Step 4: Plotting significant genes ----------------- 

gene_name_list <- c()
sum_p_val <- 0
p_val_data <- list()

for (gene_name in high_freq_gene_name) {
  quartile_breaks <- quantile(data_survival[[gene_name]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  data_survival$quartile_group <- cut(data_survival[[gene_name]],
                                      breaks = c(-Inf, quartile_breaks, Inf),
                                      labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
                                      include.lowest = TRUE)
  # print(table(data_survival$tercile_group))
  
  # log-rank testing
  survdiff_quartiles <- survdiff(Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ quartile_group, data = data_survival)
  p_value <- 1 - pchisq(survdiff_quartiles$chisq, df = length(survdiff_quartiles$n) - 1)
  p_val_data <- p_value
  if (p_value < p_val_threshold){
    sum_p_val <- sum_p_val + 1
    gene_name_list <- c(gene_name_list, gene_name)
  }
  fit_quartiles <- survfit(Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ quartile_group, data = data_survival)

  # p <- ggsurvplot(
  #   fit_quartiles,
  #   data = data_survival,
  #   pval = TRUE,
  #   risk.table = TRUE,
  #   conf.int = TRUE,
  #   risk.table.col = "strata",
  #   palette = c("#E7B800", "#2E9FDF", "steelblue", "darkred"),
  #   ggtheme = theme_bw(),
  #   title = paste("Survival Curve for", gene_name),
  #   legend.title = "Expression Quartile",
  #   legend.labs = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"))
  # 
  # print(p)

}

print(paste('sum of significant genes:', sum_p_val))
print(paste('p-value significant gene names:', paste(gene_name_list, collapse = ',')))

genes_df <- data.frame(Gene_Names = gene_name_list)
name_file_name <- "gene_name_list.xlsx"
name_file_name <- file.path(output_dir, name_file_name)
write_xlsx(genes_df, path = name_file_name)
