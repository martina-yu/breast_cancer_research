library(survival)
library(caret)

library(ggplot2)
library(survminer)
library(dplyr)

library(survcomp)
library(glmnet)
library(openxlsx)

getwd()
setwd("/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W5_Feature_CoxPH")
X_train <- read.csv("X_train.csv")
X_test <- read.csv("X_test.csv")
y_train <- read.csv("y_train.csv")
y_test <- read.csv("y_test.csv")

feature_names <- read.xlsx('gene_name_list.xlsx', 1)
feature_names <- feature_names$Gene_Names

y_train$event <- as.logical(y_train$event)
y_train$event <- as.numeric(y_train$event)
y_test$event <- as.logical(y_test$event)
y_test$event <- as.numeric(y_test$event)

train_df <- cbind(X_train, y_train)
test_df <- cbind(X_test, y_test)

feature_names <- feature_names[1:20]
formula <- as.formula(paste("Surv(time, event) ~ ", 
                            paste(feature_names, collapse = " + ")))
# formula <- Surv(time, event) ~ .
# coxph_model <- coxph(formula, data = train_df, control = coxph.control(iter.max = 50))
coxph_model <- coxph(formula, data = train_df)
plot(cox.zph(coxph_model))

for (gene in feature_names){
  boxplot(train_df[[gene]] ~ train_df$event,
          xlab = "Survival Event (0=Censor, 1=Event)",
          ylab = paste("Gene Expression", gene))
}

# table(train_df$Menospausal.Status_Pre.Menopausal, train_df$event)

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
  # linetype = "strata",
  legend.title = "Risk Quartile",
  legend.labs = c("Q1 (Low Risk)", "Q2", "Q3", "Q4 (High Risk)"),
  palette = "jco",
  ggtheme = theme_bw(),
  title = "Kaplan-Meier Curves by Risk Quartile"
)

predicted_logits <- predict(coxph_model, newdata = test_df, type = 'lp')

c_index_result <- concordance(Surv(test_df$time, test_df$event) ~ predicted_logits)
print(c_index_result$concordance)

