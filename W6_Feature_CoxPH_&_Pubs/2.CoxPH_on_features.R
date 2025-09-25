library(survival)
library(caret)

library(ggplot2)
library(survminer)
library(dplyr)

library(survcomp)
library(glmnet)
library(openxlsx)
library(caret)

getwd()
setwd("/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W6_Feature_CoxPH_&_Pubs")

n <- 40
k_folds <- 4

X_train <- read.csv(paste0("X_train_",n,"_HR_FDR.csv"))
X_test <- read.csv(paste0("X_test_",n,"_HR_FDR.csv"))
y_train <- read.csv(paste0("y_train_",n,"_HR_FDR.csv"))
y_test <- read.csv(paste0("y_test_",n,"_HR_FDR.csv"))

feature_names <- read.xlsx(paste0('./by_HR_FDR_P/Multicox_results_top_',n,"_filtered_by_HR_P.xlsx"), 1)
feature_names <- feature_names$gene_name

y_train$event <- as.logical(y_train$event)
y_train$event <- as.numeric(y_train$event)
y_test$event <- as.logical(y_test$event)
y_test$event <- as.numeric(y_test$event)

train_df <- cbind(X_train, y_train)
test_df <- cbind(X_test, y_test)

combined_data <- rbind(train_df, test_df)
folds <- createFolds(combined_data$event, k = k_folds, list = TRUE, returnTrain = FALSE)

formula <- as.formula(paste("Surv(time, event) ~ ", 
                            paste(feature_names, collapse = " + ")))
# formula <- Surv(time, event) ~ .
# coxph_model <- coxph(formula, data = train_df, control = coxph.control(iter.max = 50))
coxph_model <- coxph(formula, data = train_df)
plot(cox.zph(coxph_model))

# Create an empty vector to store C-index results from each fold
c_indices <- c()

for (i in 1:k_folds) {
  validation_indices <- folds[[i]]

  temp_train_df <- combined_data[-validation_indices, ]
  temp_validation_df <- combined_data[validation_indices, ]

  coxph_model_cv <- coxph(formula, data = temp_train_df)

  predicted_logits_cv <- predict(coxph_model_cv, newdata = temp_validation_df, type = 'lp')

  c_index_result_cv <- concordance(Surv(temp_validation_df$time, temp_validation_df$event) ~ predicted_logits_cv)

  c_indices <- c(c_indices, c_index_result_cv$concordance)
}

average_c_index <- mean(c_indices)
sd_c_index <- sd(c_indices)

print(paste("Average C-index:", round(average_c_index, 3)))
print(paste("Standard deviation of C-index:", round(sd_c_index, 3)))







risk_scores <- predict(coxph_model, newdata = test_df, type = "risk")
print(head(risk_scores))
quartile_breaks <- quantile(risk_scores, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
test_df$quartile_group <- cut(risk_scores,
                              breaks = c(-Inf, quartile_breaks, Inf),
                              labels = c("Q1 (Low Risk)", "Q2", "Q3", "Q4 (High Risk)"),
                              include.lowest = TRUE)
fit <- survfit(Surv(time, event) ~ quartile_group, data = test_df)

ggsurvplot(
  fit,
  data = test_df,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.title = "Risk Quartile",
  legend.labs = c("Q1 (Low Risk)", "Q2", "Q3", "Q4 (High Risk)"),
  palette = "jco",
  ggtheme = theme_bw(),
  title = "Kaplan-Meier Curves by Risk Quartile"
)

