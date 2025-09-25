library(survival)
library(caret)

library(ggplot2)
library(survminer)
library(dplyr)

library(survcomp)
library(glmnet)

getwd()

setwd("/Users/yuzimeng/Desktop/CBB/Yale/Lajos_Lab/W2_DFA&SA_on_clinic_data")
X_train <- read.csv("X_train.csv")
X_test <- read.csv("X_test.csv")
y_train <- read.csv("y_train.csv")
y_test <- read.csv("y_test.csv")

y_train$event <- as.logical(y_train$event)
y_train$event <- as.numeric(y_train$event)
y_test$event <- as.logical(y_test$event)
y_test$event <- as.numeric(y_test$event)

train_df <- cbind(X_train, y_train)
test_df <- cbind(X_test, y_test)

features <- colnames(X_train)[1:50]
formula <- as.formula(paste("Surv(time, event) ~ ", 
                            paste(features, collapse = " + ")))
#formula <- Surv(time, event) ~ .
coxph_model <- coxph(formula, data = train_df)

predicted_logits <- predict(coxph_model, newdata = test_df, type = "risk")
print(head(predicted_logits))
