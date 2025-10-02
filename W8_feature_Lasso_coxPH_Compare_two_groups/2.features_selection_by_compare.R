library(writexl)

library(dplyr)
library(matrixStats)

library("survival")
library(glmnet)
library("survminer")
library(ggplot2)

getwd()
data_survival <- read.csv("/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W2_DFA&SA_on_clinic_data/placebo_both_clinic_gene_data_for_training_after_merged.csv")
output_dir <- "/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W8_feature_Lasso_coxPH_Compare_two_groups/"

### -------------- varibale definition -------------- ###

p_val_threshold <- 0.05
n_total_cols <- ncol(data_survival)
n_clinical_cols <- 10

gene_columns_start <- 2
gene_columns_end <- n_total_cols - 10
gene_cols_index <- 2:(n_total_cols - n_clinical_cols)
gene_expression_data <- data_survival[, gene_cols_index] 
clinical_data <- data_survival[, -gene_cols_index]

clinical_data_filtered <- clinical_data[, c("Time.to.event.if.any..days.", "IDFS.Event")]

compare_at_time_point <- function(data_survival, gene_expression_data, time_in_days) {
  # Create a grouping variable
  data_survival <- clinical_data
  gene_expression_data <- gene_expression_log2
  time_in_days <- 365
  
  group <- rep("Censor", nrow(data_survival))
  group[data_survival$Time.to.event.if.any..days. <= time_in_days & data_survival$IDFS.Event == 1] <- "Event"
  
  # Remove patients who are censored before the time point
  valid_patients <- data_survival$Time.to.event.if.any..days. >= time_in_days | data_survival$IDFS.Event == 1
  gene_data_subset <- gene_expression_data[valid_patients, ]
  group_subset <- group[valid_patients]
  
  if(length(unique(group_subset)) < 2) {
    print(paste("Not enough groups to compare at", time_in_days, "days."))
    return(NULL)
  }
  
  # Perform Wilcoxon test for each gene
  p_values <- apply(gene_data_subset, 2, function(x) {
    wilcox.test(x ~ group_subset)$p.value
  })
  
  # Calculate Fold Change
  # Fold Change: A组平均值/B组平均值，通常使用log_2(Fold_Change)以便symmetrical
  # 设置一个阈值（例如，Fold Change > 1.5 或 < 1/1.5）。
  fold_changes <- apply(gene_data_subset, 2, function(x) {
    mean_event <- mean(x[group_subset == "Event"])
    mean_censor <- mean(x[group_subset == "Censor"])
    return(mean_event / mean_censor)
  })
  
  # Combine results into a dataframe
  results <- data.frame(
    Gene = names(p_values),
    p_value = p_values,
    FoldChange = fold_changes
  )
  
  # Sort by p-value and then by absolute log2(FoldChange)
  results$log2FoldChange <- log2(results$FoldChange)
  results <- results %>%
    arrange(p_value, desc(abs(log2FoldChange)))
  
  return(results)
}

TPM <- 0.8
MAD <- 0.2
num_patients <- nrow(gene_expression_data)
results_list <- list()
set.seed(42)

result_list <- list()
hr_threshold_high <- 1.5
hr_threshold_low <- 1/1.5

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

hr_threshold_high <- 1.2
hr_threshold_low <- 1/1.2

# Run the analysis for 1, 3, and 5 years (365, 1095, 1825 days)
results_1_year <- compare_at_time_point(clinical_data, gene_expression_log2, 365)
results_1_year <-
results_3_year <- compare_at_time_point(clinical_data, gene_expression_log2, 1095)
results_5_year <- compare_at_time_point(clinical_data, gene_expression_log2, 1825)

all_results <- list(
  "Results_1_Year" = results_1_year,
  "Results_3_Year" = results_3_year,
  "Results_5_Year" = results_5_year
)

file_name <- "Results_by_Groups.xlsx"
file_path <- paste0(output_dir, file_name)
write_xlsx(all_results,file_path)

results_1_year_filtered <- results_1_year %>%
  filter(p_value < p_val_threshold & (FoldChange > hr_threshold_high | FoldChange < hr_threshold_low))
results_3_year_filtered <- results_3_year %>%
  filter(p_value < p_val_threshold & (FoldChange > hr_threshold_high | FoldChange < hr_threshold_low))
results_5_year_filtered <- results_5_year %>%
  filter(p_value < p_val_threshold & (FoldChange > hr_threshold_high | FoldChange < hr_threshold_low))

all_results_filtered <- list(
  "results_1_year_filtered" = results_1_year_filtered,
  "results_3_year_filtered" = results_3_year_filtered,
  "results_5_year_filtered" = results_5_year_filtered
)

file_name <- "Results_by_Groups_filtered.xlsx"
file_path <- paste0(output_dir, file_name)
write_xlsx(all_results_filtered,file_path)

# Example: Get the top 20 genes for each time point
top_genes_1y <- head(results_1_year, 20)
top_genes_3y <- head(results_3_year, 20)
top_genes_5y <- head(results_5_year, 20)