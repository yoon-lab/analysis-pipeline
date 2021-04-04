---
title: "Dr. Choi_Generalized toxicity prediction model_SciRep 2018"
author: "Xiao Xiao"
date: "2021/4/4"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Load necessary packages and load data

```{r load data}
setwd("E:/coding materials/Dr Choi data")
library(clusterSim)
library(caret)
library(openxlsx)
cjsdata<-read.xlsx("Cjs.xlsx",sheet=2)
```

### Convert character values into factor & Normalization & Split data

### Different normalization methods are used here to fit different modelling methods 

```{r convert}
word<- c("assay","cell_line","cell_species","cell_origin","cell_type","toxicity")
for ( i in word) {
cjsdata[,i] <- as.factor(cjsdata[,i]) }
 
 preObj <- cjsdata[,c(1:8,14,15)] 
 
 # z-score normalization
preObj_zsocre <- data.Normalization(preObj, type = "n1", normalization = "column")
cjs.pre_zscore <- data.frame(preObj_zsocre, cjsdata[,c(9:13,16)]) 
 
 inTrain_zscore <- createDataPartition(y = cjs.pre_zscore$toxicity, p = 0.6, list = FALSE)
 train_zscore <- cjs.pre_zscore[inTrain_zscore,] 
 test_zscore <- cjs.pre_zscore[-inTrain_zscore,] 
 
 # log 10 preprocessing
 preObj_log <- log10(abs(preObj))
 cjs.pre_log <- data.frame(preObj_log, cjsdata[,c(9:13,16)])
 
  inTrain_log <- createDataPartition(y = cjs.pre_log$toxicity, p = 0.6, list = FALSE)
 train_log <- cjs.pre_zscore[inTrain_log,] 
 test_log <- cjs.pre_zscore[-inTrain_log,] 
 
 
 # min-max normalization
 normalize <- function(x)
{
    return((x- min(x)) /(max(x)-min(x)))
 }

```

## Random Forest (use z-score preprocessed data as an example)

```{r random forest}
mod_rf <- train(toxicity~., data=train_zscore, method='rf', metric='Accuracy')
pred_rf <- predict(mod_rf, test_zscore) 
confusionMatrix(pred_rf,as.factor(test_zscore$toxicity))
```

## Generalized Linear Model (Use log transformed data as an example)

```{r GLM}
options(warn=-1)   
mod_fit <- train(toxicity ~.,  data=train_log, method="glm", family="binomial")
pred_glm <- predict(mod_fit, newdata=test_log)
confusionMatrix(pred_glm,as.factor(test_log$toxicity))
```


## Support Vector Machine

```{R SVM}
mod_svm <- train(toxicity ~., data =train_zscore, method = "svmLinear",preProcess = c("center", "scale"),tuneLength = 10)
pred_svm <- predict(mod_svm, newdata = test_zscore)
confusionMatrix(pred_svm,as.factor(test_zscore$toxicity))
```

## Neural Network

```{R NN}
mod_nnet <- train(toxicity ~., data = train_zscore, method = "nnet")
pred_nnet <- predict(mod_nnet, newdata = test_zscore)
confusionMatrix(pred_nnet,as.factor(test_zscore$toxicity))
```
                                   

