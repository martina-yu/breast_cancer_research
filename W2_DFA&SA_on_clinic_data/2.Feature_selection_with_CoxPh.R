# install.packages("survival")
# install.packages("survminer")

# install.packages("dplyr")
# install.packages("matrixStats")
# install.packages("writexl")

library(writexl)

library(dplyr)
library(matrixStats)

library("survival")
library("survminer")
library(ggplot2)

## ----------------- Step 1: Preprocess -----------------
getwd()
data_survival <- read.csv("./W2_DFA&SA_on_clinic_data/placebo_both_clinic_gene_data_for_training_after_merged.csv")
p_val_threshold <- 0.05
output_dir <- "/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W2_DFA&SA_on_clinic_data/survival_results"

# tail(colnames(data_survival), -10)
## 低表达 TPM < 1, less than 20%, no expression
## 变异度 MAD IQR，排序，Top 50%，25%
## 标准化 log2+1，z-score

n_total_cols <- ncol(data_survival)
n_clinical_cols <- 10

gene_columns_start <- 2
# gene_columns_end <- 21 ## mini trial
gene_columns_end <- n_total_cols - 10

## extract gene expression matrix & clinic information
gene_cols_index <- 2:(n_total_cols - n_clinical_cols)
# gene_mini_index <- 2:21 ## mini trial

# gene_expression_data <- data_survival[, gene_mini_index] ## mini trial
gene_expression_data <- data_survival[, gene_cols_index] 
clinical_data <- data_survival[, -gene_cols_index]
# 
## --------- 1.1 screen genes with TPM >=1, less than 20%(30 person) ---------
num_patients <- nrow(gene_expression_data)
threshold <- 0.80 * num_patients
genes_to_keep_low_expr <- colSums(gene_expression_data < 1) < threshold
gene_expression_data_filtered1 <- gene_expression_data[, genes_to_keep_low_expr]
print(ncol(gene_expression_data_filtered1))

## --------- 1.2 Mean Absolution Deviation, keep highest 25% ---------
mad_values <- apply(gene_expression_data_filtered1, 2, mad, na.rm = TRUE)
mad_sorted <- sort(mad_values, decreasing = TRUE)

num_genes_to_keep <- floor(length(mad_sorted) * 0.6) 
# num_genes_to_keep <- floor(length(mad_sorted) * 0.3)
# num_genes_to_keep <- floor(length(mad_sorted) * 0.5) ## set up as 50%
# num_genes_to_keep <- floor(length(mad_sorted) * 0.25)

genes_to_keep_mad <- names(mad_sorted[1:num_genes_to_keep]) ## keeped gene names
gene_expression_data_filtered2 <- gene_expression_data_filtered1[, genes_to_keep_mad]

print(ncol(gene_expression_data_filtered2))

## --------- 1.3 Normalization ---------

gene_expression_log2 <- log2(gene_expression_data_filtered2 + 1) ### log

# gene_expression_zscore <- as.data.frame(scale(gene_expression_log2)) ### z-score
# data_survival <- cbind(clinical_data, gene_expression_zscore) ### z-score

data_survival <- cbind(clinical_data, gene_expression_log2)
gene_names <- genes_to_keep_mad

## ----------------- Step 2: Result Storage -----------------
cox_results <- data.frame(
  gene = character(),
  p_value = numeric(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  stringsAsFactors = FALSE
)

subgrp_results_pre <- data.frame(
  gene = character(),
  p_value = numeric(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  stringsAsFactors = FALSE
)

subgrp_results_post <- data.frame(
  gene = character(),
  p_value = numeric(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  stringsAsFactors = FALSE
)

subgrp_results_age_lt_50 <- data.frame(
  gene = character(),
  p_value = numeric(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  stringsAsFactors = FALSE
)

subgrp_results_age_ge_50 <- data.frame(
  gene = character(),
  p_value = numeric(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  stringsAsFactors = FALSE
)

subgrp_results_rg_1 <- data.frame(
  gene = character(),
  p_value = numeric(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  stringsAsFactors = FALSE
)

subgrp_results_rg_2 <- data.frame(
  gene = character(),
  p_value = numeric(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  stringsAsFactors = FALSE
)

subgrp_results_rg_3 <- data.frame(
  gene = character(),
  p_value = numeric(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  stringsAsFactors = FALSE
)

subgrp_results_rg_4 <- data.frame(
  gene = character(),
  p_value = numeric(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  stringsAsFactors = FALSE
)


clinical_cox_results <- data.frame(
  variable = character(),
  hazard_ratio = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# ----------------- Step 3: Feature Selection -----------------
# 单变量（可以加亚组分析）、多变量
# 其他方法 Lasso-Cox

# --------- 3.1  Univariate Cox Regression ---------
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
    
    p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
    hr <- cox_summary$coefficients[1, "exp(coef)"] #### lower,upper 95%
    
    lower_95 <- cox_summary$conf.int[1, "lower .95"]
    upper_95 <- cox_summary$conf.int[1, "upper .95"]

    cox_results <- rbind(cox_results, data.frame(
      gene = gene_name,
      p_value = p_val,
      hazard_ratio = hr,
      lower_95 = lower_95,
      upper_95 = upper_95
    ))
  }
}

# ----------------- 3.2: Multiple testing correction -----------------

# --------- final results with Cox model---------

cox_results$adj_p_value <- p.adjust(cox_results$p_value, method = "BH") ## Benjamini-Hochberg，FDR
final_results <- cox_results[order(cox_results$adj_p_value), ] ## Benjamini-Hochberg，FDR

head(final_results, 20)

# output_dir <- "/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W2_DFA&SA_on_clinic_data/survival_results"

# file_name <- "all_cox_results_withFDR.xlsx" ## Benjamini-Hochberg，FDR
# file_name <- "all_cox_results_withFDR_25_MAD.xlsx" ## MAD 25%
# full_file_path <- file.path(output_dir, file_name)
# write_xlsx(final_results, path = full_file_path)

# --------- significant results ---------

# significant_genes <- final_results[final_results$adj_p_value < p_val_threshold, ] ## Benjamini-Hochberg，FDR
significant_genes <- final_results[final_results$p_value < p_val_threshold, ]
number_of_significant_genes <- nrow(significant_genes)
print(paste("Totally", number_of_significant_genes, "genes are significant for P < ", p_val_threshold))

# sig_file_name <- "significant_genes_cox_results.xlsx"
# sig_file_name <- "significant_genes_cox_results_25_MAD.xlsx" ## MAD 25%
# sig_gene_file_path <- file.path(output_dir, sig_file_name)
# write_xlsx(significant_genes, path = sig_gene_file_path)

# --------- Plotting Survival curves ---------

top_20_genes <- head(significant_genes$gene, 20)
plot_outpit_dir <- file.path(getwd(), "W2_DFA&SA_on_clinic_data", "survival_plots")

# for (gene_name in top_20_genes) {
#   
#   ### --------- group by median gene expression ---------
#   # median_expr <- median(data_survival[[gene_name]], na.rm = TRUE)
#   # data_survival$gene_group <- ifelse(data_survival[[gene_name]] > median_expr, "High", "Low")
#   # 
#   # ggsurvplot(fit_terciles,
#   #            data = data_survival,
#   #            pval = TRUE,
#   #            risk.table = TRUE,
#   #            title = "Survival Curve")
#   
#   # km_fit <- survfit(Surv(Time.to.event.if.any..days., IDFS.Event) ~ gene_group, data = data_survival)
#   
#   ### --------- group by quatile gene expression ---------
#   # gene_name <- "GGA3"
#   quartile_breaks <- quantile(data_survival[[gene_name]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
#   
#   data_survival$quartile_group <- cut(data_survival[[gene_name]],
#                                       breaks = c(-Inf, quartile_breaks, Inf),
#                                       labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
#                                       include.lowest = TRUE)
#   
#   # table(data_survival$tercile_group)# 50, 50, 50
#   
#   #perform log-rank test
#   survdiff_quartiles <- survdiff(Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ quartile_group, data = data_survival)
#   print(survdiff_quartiles)
#   
#   # fit K-M curve
#   fit_quartiles <- survfit(Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ quartile_group, data = data_survival)
# 
#   # survive plot
#   p <- ggsurvplot(
#     fit_quartiles,
#     data = data_survival,
#     pval = TRUE,
#     risk.table = TRUE,
#     conf.int = TRUE,
#     risk.table.col = "strata",
#     palette = c("#E7B800", "#2E9FDF", "steelblue", "darkred"),
#     ggtheme = theme_bw(),
#     title = paste("Survival Curve for", gene_name),
#     legend.title = "Expression Quartile",
#     legend.labs = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"))
#   
#   print(p)
#   
#   # file_name <- paste0(gene_name, "_survival_curve.png")
#   # file_path <- file.path(output_dir, file_name)
#   # ggsave(file_path, plot = p$plot, width = 8, height = 6)
#   
# }


# --------- 3.3: Subgroup Analysis ---------

pre_menospausal_data <- subset(data_survival, `Menospausal.Status` == "Pre-Menopausal")
post_menospausal_data <- subset(data_survival, `Menospausal.Status` == "Post-Menopausal")

print(nrow(pre_menospausal_data))
print(nrow(post_menospausal_data))

# --------- 3.3.1: Subgroup Analysis pre menopausal status ---------

## pre menopausal status
for (gene_name in gene_names) {
  
  formula_str <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", gene_name, "`")
  cox_formula <- as.formula(formula_str)
  
  cox_fit <- tryCatch({
    coxph(cox_formula, data = pre_menospausal_data)
  }, error = function(e) {
    message(paste("At", gene_name, "running model wrongly：", e$message))
    return(NULL)
  })
  
  if (!is.null(cox_fit)) {
    cox_summary <- summary(cox_fit)
    
    p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
    hr <- cox_summary$coefficients[1, "exp(coef)"] #### lower,upper 95%
    
    lower_95 <- cox_summary$conf.int[1, "lower .95"]
    upper_95 <- cox_summary$conf.int[1, "upper .95"]
    
    subgrp_results_pre <- rbind(subgrp_results_pre, data.frame(
      gene = gene_name,
      p_value = p_val,
      hazard_ratio = hr,
      lower_95 = lower_95,
      upper_95 = upper_95
    ))
  }
}

subgrp_results_pre$adj_p_value <- p.adjust(subgrp_results_pre$p_value, method = "BH")
final_pre_results <- subgrp_results_pre[order(subgrp_results_pre$adj_p_value), ]
file_name <- "Unicox_results_preMeno.xlsx"
full_file_path <- file.path(output_dir, file_name)
write_xlsx(final_pre_results, path = full_file_path)

# --------- 3.3.2: Subgroup Analysis post menopausal status ---------

### post menopausal status
for (gene_name in gene_names) {
  
  formula_str <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", gene_name, "`")
  cox_formula <- as.formula(formula_str)
  
  cox_fit <- tryCatch({
    coxph(cox_formula, data = post_menospausal_data)
  }, error = function(e) {
    message(paste("At", gene_name, "running model wrongly：", e$message))
    return(NULL)
  })
  
  if (!is.null(cox_fit)) {
    cox_summary <- summary(cox_fit)
    
    p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
    hr <- cox_summary$coefficients[1, "exp(coef)"] #### lower,upper 95%
    
    lower_95 <- cox_summary$conf.int[1, "lower .95"]
    upper_95 <- cox_summary$conf.int[1, "upper .95"]
    
    subgrp_results_post <- rbind(subgrp_results_post, data.frame(
      gene = gene_name,
      p_value = p_val,
      hazard_ratio = hr,
      lower_95 = lower_95,
      upper_95 = upper_95
    ))
  }
}

subgrp_results_post$adj_p_value <- p.adjust(subgrp_results_post$p_value, method = "BH")
final_post_results <- subgrp_results_post[order(subgrp_results_post$adj_p_value), ]
# head(final_post_results, 20)
file_name <- "Unicox_results_postMeno.xlsx"
full_file_path <- file.path(output_dir, file_name)

write_xlsx(final_post_results, path = full_file_path)

# --------- 3.3.3: Subgroup Analysis post & pre significant genes---------

sig_pre <- subgrp_results_pre[subgrp_results_pre$p_value < p_val_threshold, ]
sig_post <- subgrp_results_post[subgrp_results_post$p_value < p_val_threshold, ]
num_sig_pre <- nrow(sig_pre)
num_sig_post <- nrow(sig_post)
print(paste("On Pre-Menopausal status: Totally", num_sig_pre, "genes are significant for P < ", p_val_threshold))
print(paste("On Post-Menopausal status: Totally", num_sig_post, "genes are significant for P < ", p_val_threshold))

# --------- 3.3.4: Subgroup Analysis: compare and find the intersection ---------

name_sig_pre <- subgrp_results_pre$gene[subgrp_results_pre$p_value < p_val_threshold]
name_sig_post <- subgrp_results_pre$gene[subgrp_results_post$p_value < p_val_threshold]

overlapping_genes <- intersect(name_sig_pre, name_sig_post)
print(paste("significant on both：", length(overlapping_genes)))

unique_pre_genes <- setdiff(name_sig_pre, name_sig_post)
print(paste("significant only on pre group：", length(unique_pre_genes)))

unique_post_genes <- setdiff(name_sig_post, name_sig_pre)
print(paste("significant only on post group：", length(unique_post_genes)))

library(VennDiagram)

venn.diagram(
  x = list(name_sig_pre, name_sig_post),
  category.names = c("Pre-Menopausal", "Post-Menopausal"),
  filename = 'Venn_Diagram_Subgroups.png',
  output = TRUE
)

# --------- 3.3.5: Subgroup Analysis Age on 50 groups ---------

age_lt_50_data <- subset(data_survival, Age < 50)
age_ge_50_data <- subset(data_survival, Age >= 50)

print(nrow(age_ge_50_data))
print(nrow(age_lt_50_data))

for (gene_name in gene_names) {
  
  formula_str <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", gene_name, "`")
  cox_formula <- as.formula(formula_str)
  
  cox_fit <- tryCatch({
    coxph(cox_formula, data = age_lt_50_data)
  }, error = function(e) {
    message(paste("At", gene_name, "running model wrongly：", e$message))
    return(NULL)
  })
  
  if (!is.null(cox_fit)) {
    cox_summary <- summary(cox_fit)
    
    p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
    hr <- cox_summary$coefficients[1, "exp(coef)"] #### lower,upper 95%
    
    lower_95 <- cox_summary$conf.int[1, "lower .95"]
    upper_95 <- cox_summary$conf.int[1, "upper .95"]
    
    subgrp_results_age_lt_50 <- rbind(subgrp_results_age_lt_50, data.frame(
      gene = gene_name,
      p_value = p_val,
      hazard_ratio = hr,
      lower_95 = lower_95,
      upper_95 = upper_95
    ))
  }
}

subgrp_results_age_lt_50$adj_p_value <- p.adjust(subgrp_results_age_lt_50$p_value, method = "BH")
final_post_results <- subgrp_results_age_lt_50[order(subgrp_results_age_lt_50$adj_p_value), ]
file_name <- "Unicox_results_age_lessthan_50.xlsx"
full_file_path <- file.path(output_dir, file_name)
write_xlsx(final_post_results, path = full_file_path)

for (gene_name in gene_names) {
  
  formula_str <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", gene_name, "`")
  cox_formula <- as.formula(formula_str)
  
  cox_fit <- tryCatch({
    coxph(cox_formula, data = age_ge_50_data)
  }, error = function(e) {
    message(paste("At", gene_name, "running model wrongly：", e$message))
    return(NULL)
  })
  
  if (!is.null(cox_fit)) {
    cox_summary <- summary(cox_fit)
    
    p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
    hr <- cox_summary$coefficients[1, "exp(coef)"] #### lower,upper 95%
    
    lower_95 <- cox_summary$conf.int[1, "lower .95"]
    upper_95 <- cox_summary$conf.int[1, "upper .95"]
    
    subgrp_results_age_ge_50 <- rbind(subgrp_results_age_ge_50, data.frame(
      gene = gene_name,
      p_value = p_val,
      hazard_ratio = hr,
      lower_95 = lower_95,
      upper_95 = upper_95
    ))
  }
}

subgrp_results_age_ge_50$adj_p_value <- p.adjust(subgrp_results_age_ge_50$p_value, method = "BH")
final_post_results <- subgrp_results_age_ge_50[order(subgrp_results_age_ge_50$adj_p_value), ]
file_name <- "Unicox_results_age_greater_equal_50.xlsx"
full_file_path <- file.path(output_dir, file_name)
write_xlsx(final_post_results, path = full_file_path)

sig_lt <- subgrp_results_age_lt_50[subgrp_results_age_lt_50$p_value < p_val_threshold, ]
sig_ge <- subgrp_results_age_ge_50[subgrp_results_age_ge_50$p_value < p_val_threshold, ]
num_sig_lt <- nrow(sig_lt)
num_sig_ge <- nrow(sig_ge)
print(paste("On age < 50: Totally", num_sig_lt, "genes are significant for P < ", p_val_threshold))
print(paste("On age >= 50: Totally", num_sig_ge, "genes are significant for P < ", p_val_threshold))

name_sig_lt <- subgrp_results_age_lt_50$gene[subgrp_results_age_lt_50$p_value < p_val_threshold]
name_sig_ge <- subgrp_results_age_ge_50$gene[subgrp_results_age_ge_50$p_value < p_val_threshold]

overlapping_genes <- intersect(name_sig_ge, name_sig_lt)
print(paste("significant on both：", length(overlapping_genes)))

unique_pre_genes <- setdiff(name_sig_lt, name_sig_ge)
print(paste("significant only on age < 50 group：", length(unique_pre_genes)))

unique_post_genes <- setdiff(name_sig_ge, name_sig_lt)
print(paste("significant only on age >= 50 group：", length(unique_post_genes)))

# library(VennDiagram)

venn.diagram(
  x = list(name_sig_lt, name_sig_ge),
  category.names = c("Age < 50", "Age >= 50"),
  filename = 'Venn_Diagram_age.png',
  output = TRUE
)

# --------- 3.3.6: Subgroup Analysis on risk groups ---------

risk_group_1_data <- subset(data_survival, `Risk.Group` == "1")
risk_group_2_data <- subset(data_survival, `Risk.Group` == "2")
risk_group_3_data <- subset(data_survival, `Risk.Group` == "3")
risk_group_4_data <- subset(data_survival, `Risk.Group` == "4")


nrow(risk_group_4_data)

for (gene_name in gene_names) {
  
  formula_str <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", gene_name, "`")
  cox_formula <- as.formula(formula_str)
  
  cox_fit <- tryCatch({
    coxph(cox_formula, data = risk_group_1_data)
  }, error = function(e) {
    message(paste("At", gene_name, "running model wrongly：", e$message))
    return(NULL)
  })
  
  if (!is.null(cox_fit)) {
    cox_summary <- summary(cox_fit)
    
    p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
    hr <- cox_summary$coefficients[1, "exp(coef)"] #### lower,upper 95%
    
    lower_95 <- cox_summary$conf.int[1, "lower .95"]
    upper_95 <- cox_summary$conf.int[1, "upper .95"]
    
    subgrp_results_rg_1 <- rbind(subgrp_results_rg_1, data.frame(
      gene = gene_name,
      p_value = p_val,
      hazard_ratio = hr,
      lower_95 = lower_95,
      upper_95 = upper_95
    ))
  }
}

subgrp_results_rg_1$adj_p_value <- p.adjust(subgrp_results_rg_1$p_value, method = "BH")
final_rg1 <- subgrp_results_rg_1[order(subgrp_results_rg_1$adj_p_value), ]

file_name <- "Unicox_results_rg1.xlsx"
full_file_path <- file.path(output_dir, file_name)
write_xlsx(final_rg1, path = full_file_path)

sig_rg1 <- subgrp_results_rg_1[subgrp_results_rg_1$p_value < p_val_threshold, ]
num_sig_rg1 <- nrow(sig_rg1)

print(paste("On risk group 1: Totally", num_sig_rg1, "genes are significant for P < ", p_val_threshold))


for (gene_name in gene_names) {
  
  formula_str <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", gene_name, "`")
  cox_formula <- as.formula(formula_str)
  
  cox_fit <- tryCatch({
    coxph(cox_formula, data = risk_group_2_data)
  }, error = function(e) {
    message(paste("At", gene_name, "running model wrongly：", e$message))
    return(NULL)
  })
  
  if (!is.null(cox_fit)) {
    cox_summary <- summary(cox_fit)
    
    p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
    hr <- cox_summary$coefficients[1, "exp(coef)"] #### lower,upper 95%
    
    lower_95 <- cox_summary$conf.int[1, "lower .95"]
    upper_95 <- cox_summary$conf.int[1, "upper .95"]
    
    subgrp_results_rg_2 <- rbind(subgrp_results_rg_2, data.frame(
      gene = gene_name,
      p_value = p_val,
      hazard_ratio = hr,
      lower_95 = lower_95,
      upper_95 = upper_95
    ))
  }
}

subgrp_results_rg_2$adj_p_value <- p.adjust(subgrp_results_rg_2$p_value, method = "BH")
final_rg2 <- subgrp_results_rg_2[order(subgrp_results_rg_2$adj_p_value), ]
# file_name <- "Unicox_results_rg2.xlsx"
# full_file_path <- file.path(output_dir, file_name)
# write_xlsx(final_rg2, path = full_file_path)

sig_rg2 <- subgrp_results_rg_2[subgrp_results_rg_2$p_value < p_val_threshold, ]
num_sig_rg2 <- nrow(sig_rg2)
print(paste("On risk group 2: Totally", num_sig_rg2, "genes are significant for P < ", p_val_threshold))

for (gene_name in gene_names) {
  
  formula_str <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", gene_name, "`")
  cox_formula <- as.formula(formula_str)
  
  cox_fit <- tryCatch({
    coxph(cox_formula, data = risk_group_3_data)
  }, error = function(e) {
    message(paste("At", gene_name, "running model wrongly：", e$message))
    return(NULL)
  })
  
  if (!is.null(cox_fit)) {
    cox_summary <- summary(cox_fit)
    
    p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
    hr <- cox_summary$coefficients[1, "exp(coef)"] #### lower,upper 95%
    
    lower_95 <- cox_summary$conf.int[1, "lower .95"]
    upper_95 <- cox_summary$conf.int[1, "upper .95"]
    
    subgrp_results_rg_3 <- rbind(subgrp_results_rg_3, data.frame(
      gene = gene_name,
      p_value = p_val,
      hazard_ratio = hr,
      lower_95 = lower_95,
      upper_95 = upper_95
    ))
  }
}

subgrp_results_rg_3$adj_p_value <- p.adjust(subgrp_results_rg_3$p_value, method = "BH")
final_rg3 <- subgrp_results_rg_3[order(subgrp_results_rg_3$adj_p_value), ]
# file_name <- "Unicox_results_rg3.xlsx"
# full_file_path <- file.path(output_dir, file_name)
# write_xlsx(final_rg3, path = full_file_path)

sig_rg3 <- subgrp_results_rg_3[subgrp_results_rg_3$p_value < p_val_threshold, ]
num_sig_rg3 <- nrow(sig_rg3)
print(paste("On risk group 3: Totally", num_sig_rg3, "genes are significant for P < ", p_val_threshold))


for (gene_name in gene_names) {
  
  formula_str <- paste0("Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `", gene_name, "`")
  cox_formula <- as.formula(formula_str)
  
  cox_fit <- tryCatch({
    coxph(cox_formula, data = risk_group_4_data)
  }, error = function(e) {
    message(paste("At", gene_name, "running model wrongly：", e$message))
    return(NULL)
  })
  
  if (!is.null(cox_fit)) {
    cox_summary <- summary(cox_fit)
    
    p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
    hr <- cox_summary$coefficients[1, "exp(coef)"] #### lower,upper 95%
    
    lower_95 <- cox_summary$conf.int[1, "lower .95"]
    upper_95 <- cox_summary$conf.int[1, "upper .95"]
    
    subgrp_results_rg_4 <- rbind(subgrp_results_rg_4, data.frame(
      gene = gene_name,
      p_value = p_val,
      hazard_ratio = hr,
      lower_95 = lower_95,
      upper_95 = upper_95
    ))
  }
}

subgrp_results_rg_4$adj_p_value <- p.adjust(subgrp_results_rg_4$p_value, method = "BH")
final_rg4 <- subgrp_results_rg_4[order(subgrp_results_rg_4$adj_p_value), ]
# file_name <- "Unicox_results_rg4.xlsx"
# full_file_path <- file.path(output_dir, file_name)
# write_xlsx(final_rg4, path = full_file_path)

sig_rg4 <- subgrp_results_rg_4[subgrp_results_rg_4$p_value < p_val_threshold, ]
num_sig_rg4 <- nrow(sig_rg4)
print(paste("On risk group 4: Totally", num_sig_rg4, "genes are significant for P < ", p_val_threshold))

# top_gene_pre <- final_pre_results[1, ]
# top_gene_post <- final_post_results[1, ]
# 
# print(paste("Pre-Menopausal，smallest adj_p_value gene：", top_gene_pre$gene))
# print(paste("Post-Menopausal，smallest adj_p_value gene：", top_gene_post$gene))
# 
# top_gene_pre$subgroup <- "Premenopausal"
# top_gene_post$subgroup <- "Postmenopausal"
# plot_data <- rbind(top_gene_pre, top_gene_post)
# print(plot_data)

# library(ggplot2)
# ggplot(plot_data, aes(y = subgroup, x = hazard_ratio)) +
#   geom_point(aes(x = hazard_ratio), size = 4, shape = 15) +
#   geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), height = 0.2) +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
#   labs(title = paste("IDFS HR (95% CI) for", top_gene_pre$gene, "by Menopausal Status"),
#        y = NULL,
#        x = "Hazard Ratio (HR)") +
#   scale_x_log10(limits = c(0.2, 5), breaks = c(0.2, 0.5, 1, 2, 5)) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_blank())

# --------- 3.4 Multivariate Cox Regression Analysis ---------

print(nrow(significant_genes))
significant_gene_names <- significant_genes$gene
clinical_variable_names <- c("`Menospausal.Status`", "`Risk.Group`", "Age")
all_predictors <- c(significant_gene_names, clinical_variable_names)
length(all_predictors)

formula_predictors <- paste(all_predictors, collapse = " + ")
formula_str <- paste("Surv(Time.to.event.if.any..days., IDFS.Event) ~", formula_predictors)
multivariate_formula <- as.formula(formula_str)
print(multivariate_formula)

multivariate_cox_fit <- coxph(multivariate_formula, data = data_survival)
summary(multivariate_cox_fit)

cox_summary <- summary(multivariate_cox_fit)
cox_coef_table <- as.data.frame(cox_summary$coefficients)
cox_conf_int_table <- as.data.frame(cox_summary$conf.int)

multivariate_significant_results <- cox_coef_table[cox_coef_table$`Pr(>|z|)` < 0.05,]
print(nrow(multivariate_significant_results))

cox_coef_table$p.adjusted.FDR <- p.adjust(cox_coef_table$`Pr(>|z|)`, method = "fdr")
multivariate_significant_results_fdr <- cox_coef_table[cox_coef_table$p.adjusted.FDR < 0.05, ]
print(nrow(multivariate_significant_results_fdr))

# print(row.names(multivariate_significant_results))

# write_xlsx(list(
#   all_results = cox_coef_table,
#   significant_results = multivariate_significant_results,
#   conf_int = cox_conf_int_table
# ), path = file.path(output_dir, "multivariate_significant_cox_summary_25_MAD.xlsx"))


### --------------- volcano data check --------------- ----------------

volcano_data <- multivariate_significant_results %>%
  mutate(
    log2_hr = log2(`exp(coef)`),
    neg_log10_p = -log10(`Pr(>|z|)`)
  )

hr_threshold <- 1.5
log2_hr_threshold <- log2(hr_threshold)

ggplot(volcano_data, aes(x = log2_hr, y = neg_log10_p)) +
  geom_point(aes(color = `Pr(>|z|)` < p_val_threshold), alpha = 0.6) +
  geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(log2_hr_threshold, -log2_hr_threshold), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot of Multivariate Cox Results",
       x = "log2(Hazard Ratio)",
       y = "-log10(P-value)",
       color = "Significant") +
  theme_minimal() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  theme(legend.position = "top")

highly_significant_genes <- volcano_data %>%
  filter(neg_log10_p > 300)
highly_significant_genes_names <- row.names(highly_significant_genes)
print(highly_significant_genes_names)

# check outliers

for (gene_name in highly_significant_genes_names) {
  
  median_expr <- median(data_survival[[gene_name]], na.rm = TRUE)
  data_survival$gene_group <- ifelse(data_survival[[gene_name]] > median_expr, "High", "Low")
  
  km_fit <- survfit(Surv(Time.to.event.if.any..days., IDFS.Event) ~ gene_group, data = data_survival)
  
  p <- ggsurvplot(
    km_fit,
    data = data_survival,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = TRUE,
    risk.table.col = "strata",
    palette = c("#E7B800", "#2E9FDF"),
    ggtheme = theme_bw(),
    title = paste("Survival Curve for", gene_name),
    legend.title = "Gene Expression"
  )
  print(p)
  
}

# # ----------------- Step 4: Plotting significant genes ----------------- 

# multivariate_significant_results_fdr <- multivariate_significant_results_fdr[order(multivariate_significant_results_fdr$`Pr(>|z|)`), ]
# all_significant_variable_names <- row.names(multivariate_significant_results_fdr)

mul_sign_rslt <- multivariate_significant_results_fdr[order(multivariate_significant_results_fdr$`p.adjusted.FDR`), ]
mul_sign_rslt_zero <- mul_sign_rslt[mul_sign_rslt$p.adjusted.FDR < e-303, ]
all_significant_variable_names <- row.names(multivariate_significant_results_fdr)

write_xlsx(multivariate_significant_results_fdr, path = file.path(output_dir, "multivariate_for_sign_genes_25_MAD.xlsx"))

write.csv(data.frame(gene = all_significant_variable_names),
          file = file.path(output_dir, "significant_gene_names.csv"),
          row.names = FALSE)

for (gene_name in all_significant_variable_names) {
  
  quartile_breaks <- quantile(data_survival[[gene_name]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  data_survival$quartile_group <- cut(data_survival[[gene_name]],
                                      breaks = c(-Inf, quartile_breaks, Inf),
                                      labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
                                      include.lowest = TRUE)
  # print(table(data_survival$tercile_group))
  
  # log-rank testing
  survdiff_quartiles <- survdiff(Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ quartile_group, data = data_survival)
  # print(survdiff_quartiles)
  
  fit_quartiles <- survfit(Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ quartile_group, data = data_survival)
  
  p <- ggsurvplot(
    fit_quartiles,
    data = data_survival,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = TRUE,
    risk.table.col = "strata",
    palette = c("#E7B800", "#2E9FDF", "steelblue", "darkred"),
    ggtheme = theme_bw(),
    title = paste("Survival Curve for", gene_name),
    legend.title = "Expression Quartile",
    legend.labs = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"))
  
  # print(p)
  
  file_name <- paste0(gene_name, "_survival_curve.png")
  file_path <- file.path(plot_outpit_dir, file_name)
  ggsave(file_path, plot = p$plot, width = 8, height = 6)

}

# # ----------------- Step 5: Univariate Cox regression on clinic events ----------------- 


# ------- Age --------

fit_age <- coxph(Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ Age, data = data_survival)
sum_age <- summary(fit_age)
clinical_cox_results <- rbind(clinical_cox_results, data.frame(
  variable = "Age",
  hazard_ratio = sum_age$conf.int[1, "exp(coef)"],
  lower_95 = sum_age$conf.int[1, "lower .95"],
  upper_95 = sum_age$conf.int[1, "upper .95"],
  p_value = sum_age$coefficients[1, "Pr(>|z|)"]
))

# ------- Risk Group --------

data_survival$`Risk.Group` <- as.factor(data_survival$`Risk.Group`)

fit_risk_group <- coxph(Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `Risk.Group`, data = data_survival)
sum_risk_group <- summary(fit_risk_group)

detailed_results <- as.data.frame(sum_risk_group$coefficients)
conf_int_results <- as.data.frame(sum_risk_group$conf.int)

detailed_results$hazard_ratio <- conf_int_results$`exp(coef)`
detailed_results$lower_95 <- conf_int_results$`lower .95`
detailed_results$upper_95 <- conf_int_results$`upper .95`

print(detailed_results)

# -------- Menopausal Status ----------

fit_meno <- coxph(Surv(`Time.to.event.if.any..days.`, `IDFS.Event`) ~ `Menospausal.Status`, data = data_survival)
sum_meno <- summary(fit_meno)
p_value_meno <- sum_meno$sctest["pvalue"]
for (i in 1:nrow(sum_meno$conf.int)) {
  row_name <- row.names(sum_meno$conf.int)[i]
  clinical_cox_results <- rbind(clinical_cox_results, data.frame(
    variable = paste("Menopausal Status:", row_name),
    hazard_ratio = sum_meno$conf.int[i, "exp(coef)"],
    lower_95 = sum_meno$conf.int[i, "lower .95"],
    upper_95 = sum_meno$conf.int[i, "upper .95"],
    p_value = sum_meno$coefficients[i, "Pr(>|z|)"]
  ))
}

print(clinical_cox_results)
file_name <- "clinical_unicox_results.xlsx"
full_file_path <- file.path(output_dir, file_name)
write_xlsx(clinical_cox_results, path = full_file_path)
