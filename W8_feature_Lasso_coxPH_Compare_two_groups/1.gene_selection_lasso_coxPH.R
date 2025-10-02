library(writexl)

library(dplyr)
library(matrixStats)

library("survival")
library(glmnet)
library("survminer")
library(ggplot2)

### -------------- varibale definition -------------- ###

getwd()
data_survival <- read.csv("/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W2_DFA&SA_on_clinic_data/placebo_both_clinic_gene_data_for_training_after_merged.csv")
output_dir <- "/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W8_feature_Lasso_coxPH_Compare_two_groups"
p_val_threshold <- 0.05
p_threshold <- 0.05
n_total_cols <- ncol(data_survival)
n_clinical_cols <- 10
gene_columns_start <- 2
gene_columns_end <- n_total_cols - 10
gene_cols_index <- 2:(n_total_cols - n_clinical_cols)
# gene_cols_index <- 2:1000
gene_expression_data <- data_survival[, gene_cols_index] 
clinical_data <- data_survival[, -gene_cols_index]

TPM <- 0.8
MAD <- 0.2
# parameter_grid <- expand.grid(TPM = TPM, MAD = MAD)
num_patients <- nrow(gene_expression_data)
results_list <- list()
set.seed(42)

result_list <- list()
hr_threshold_high <- 1.5
hr_threshold_low <- 1/1.5

n <- 55

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
  print(paste("Totally", number_of_significant_genes, "genes are significant for P < ", p_val_threshold))
  return(list(final_results=final_results, significant_genes=significant_genes))
}

### ----------------- loop for different tpm and mad ------------------ ###


tpm <- TPM
mad <- MAD
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

### ----------------- Lasso-CoxPH for feature selection ------------------ ###

survival_object <- Surv(data_survival$`Time.to.event.if.any..days.`, data_survival$`IDFS.Event`)
x_matrix <- as.matrix(gene_expression_log2)
cv_fit <- cv.glmnet(x = x_matrix, 
                    y = survival_object, 
                    family = "cox", 
                    alpha = 1)
best_lambda <- cv_fit$lambda.min

selected_genes_coefs <- coef(cv_fit, s = best_lambda)
selected_genes <- selected_genes_coefs[selected_genes_coefs[, 1] != 0, ]
selected_gene_names <- names(selected_genes)

print("Selected genes by Lasso-Coxph:")
print(selected_gene_names)
print("Number of selected genes:")
print(length(selected_gene_names))

final_formula_predictors <- paste(selected_gene_names, collapse = " + ")
final_formula_str <- paste("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~", final_formula_predictors)
final_multivariate_formula <- as.formula(final_formula_str)
final_multivariate_cox_fit <- coxph(final_multivariate_formula, data = data_survival)
final_multicox_summary <- summary(final_multivariate_cox_fit)

coef_table <- as.data.frame(final_multicox_summary$coefficients)
conf_int_table <- as.data.frame(final_multicox_summary$conf.int)

final_results_df <- cbind(coef_table, conf_int_table)
final_results_df$p.adjusted.FDR <- p.adjust(final_results_df$`Pr(>|z|)`, method = "fdr")
final_results_df$gene_name <- rownames(final_results_df)
final_results_df <- final_results_df[, c("gene_name", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "p.adjusted.FDR", "lower .95", "upper .95")]
print(final_results_df)

file_name <- paste0("Result_Lasso_with_Multicox.xlsx")
full_file_path <- file.path(output_dir, file_name)
write_xlsx(final_results_df, path = full_file_path)


### ----------------- ------------------ ###

data_survival <- cbind(clinical_data, gene_expression_log2)
gene_names <- genes_to_keep_mad

uni_results_genes <- uniCoxFunction(gene_names, data_survival, p_val_threshold)
print(paste('unicox regression kept:', nrow(uni_results_genes$significant_genes), 'genes.'))


# multiCoxFunction <- function(uni_sig_genes, data_survival, p_val_threshold){

uni_sig_genes <- uni_results_genes
uni_sig_genes <- uni_sig_genes$significant_genes
significant_gene_names <- uni_sig_genes$gene
all_predictors <- significant_gene_names

formula_predictors <- paste(all_predictors, collapse = " + ")
formula_str <- paste("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~", formula_predictors)
multivariate_formula <- as.formula(formula_str)
multivariate_cox_fit <- coxph(multivariate_formula, data = data_survival) 

### Warning message: In coxph.fit(X, Y, istrat, offset, init, control, weights = weights,
### :Ran out of iterations and did not converge

multicox_summary <- summary(multivariate_cox_fit)
multicox_coef_table <- as.data.frame(multicox_summary$coefficients) ## coef exp(coef)  se(coef)   z  Pr(>|z|)
multicox_conf_int_table <- as.data.frame(multicox_summary$conf.int) ## exp(coef) exp(-coef) lower .95 upper .95
multicox_result <- cbind(multicox_coef_table, multicox_conf_int_table)

multicox_result$p.adjusted.FDR <- p.adjust(multicox_result$`Pr(>|z|)`, method = "fdr")

#### --------------- by just P-value ------------------- #### 
# multicox_result_sorted <- multicox_result[order(multicox_result$p.adjusted.FDR), ]

multicox_result_filtered <- subset(multicox_result,
                                   `p.adjusted.FDR` < p_threshold &
                                     (`exp(coef)` > hr_threshold_high | `exp(coef)` < hr_threshold_low))
multicox_result_sorted <- multicox_result_filtered[order(multicox_result_filtered$`p.adjusted.FDR`), ]

multicox_result_top_n <- head(multicox_result_sorted, n)
multicox_result_top_n_with_names <- multicox_result_top_n
multicox_result_top_n_with_names$gene_name <- row.names(multicox_result_top_n_with_names)

file_name <- paste0("Multicox_results_top_", nrow(multicox_result_top_n_with_names), "_filtered_by_HR_P.xlsx")
full_file_path <- file.path(output_dir, file_name)
write_xlsx(multicox_result_top_n_with_names, path = full_file_path)

  # final_results <- multicox_result[multicox_result$p.adjusted.FDR < p_val_threshold, ]
  # print(paste('multicox regression < p threshold:', nrow(final_results)))
#   significant_genes <- row.names(multicox_result_top50)
#   print(paste('row names of significant_genes:', paste(significant_genes, collapse = ',')))
#   
#   file_name <- "Multicox_results_top_20.xlsx"
#   full_file_path <- file.path(output_dir, file_name)
#   write_xlsx(multicox_result_top50, path = full_file_path)
#   
#   return(list(multicox_final_results=multicox_result_top50, multicox_significant_genes=significant_genes))
# }

# formula_str <- paste("Surv(time, status) ~", formula_predictors)
# print(multivariate_formula)
# summary(multivariate_cox_fit)

# multivariate_significant_results <- cox_coef_table[cox_coef_table$`Pr(>|z|)` < 0.05,]
# print(nrow(multivariate_significant_results))
# 
# mul_sign_rslt <- multivariate_significant_results_fdr[order(multivariate_significant_results_fdr$`p.adjusted.FDR`), ]
# mul_sign_rslt_zero <- mul_sign_rslt[mul_sign_rslt$p.adjusted.FDR < 1e-303, ]
# all_significant_variable_names <- row.names(mul_sign_rslt_zero)




### ----------------- Step 4: Plotting significant genes ----------------- 

gene_name_list <- c()
sum_p_val <- 0
p_val_data <- list()

for (gene_name in row.names(multicox_result_top50_with_names)) {
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

}

print(paste('sum of significant genes:', sum_p_val))
print(paste('p-value significant gene names:', paste(gene_name_list, collapse = ',')))

genes_df <- data.frame(Gene_Names = gene_name_list)
name_file_name <- "gene_name_list.xlsx"
name_file_name <- file.path(output_dir, name_file_name)
write_xlsx(genes_df, path = name_file_name)
