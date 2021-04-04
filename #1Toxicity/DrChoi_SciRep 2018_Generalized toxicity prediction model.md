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

```{r random forest, echo=TRUE}
mod_rf <- train(toxicity~., data=train_zscore, method='rf', metric='Accuracy')
pred_rf <- predict(mod_rf, test_zscore) 
confusionMatrix(pred_rf,as.factor(test_zscore$toxicity))
```

Confusion Matrix and Statistics
          Reference
Prediction nonToxic Toxic
  nonToxic      193     3
  Toxic           3    30
               Accuracy : 0.9738          
                 95% CI : (0.9438, 0.9903)
    No Information Rate : 0.8559          
    P-Value [Acc > NIR] : 1.689e-09       
                  Kappa : 0.8938          
 Mcnemar's Test P-Value : 1               
            Sensitivity : 0.9847          
            Specificity : 0.9091          
         Pos Pred Value : 0.9847          
         Neg Pred Value : 0.9091          
             Prevalence : 0.8559          
         Detection Rate : 0.8428          
   Detection Prevalence : 0.8559          
      Balanced Accuracy : 0.9469          
       'Positive' Class : nonToxic 

## Generalized Linear Model (Use log transformed data as an example)

```{r GLM, echo=TRUE}
options(warn=-1)   
mod_fit <- train(toxicity ~.,  data=train_log, method="glm", family="binomial")
pred_glm <- predict(mod_fit, newdata=test_log)
confusionMatrix(pred_glm,as.factor(test_log$toxicity))
```

Confusion Matrix and Statistics
          Reference
Prediction nonToxic Toxic
  nonToxic      192    11
  Toxic           4    22
               Accuracy : 0.9345          
                 95% CI : (0.8943, 0.9629)
    No Information Rate : 0.8559          
    P-Value [Acc > NIR] : 0.0001641       
                  Kappa : 0.7088          
 Mcnemar's Test P-Value : 0.1213353       
            Sensitivity : 0.9796          
            Specificity : 0.6667          
         Pos Pred Value : 0.9458          
         Neg Pred Value : 0.8462          
             Prevalence : 0.8559          
         Detection Rate : 0.8384          
   Detection Prevalence : 0.8865          
      Balanced Accuracy : 0.8231          
       'Positive' Class : nonToxic 



## Support Vector Machine

```{R SVM, echo=TRUE}
mod_svm <- train(toxicity ~., data =train_zscore, method = "svmLinear",preProcess = c("center", "scale"),tuneLength = 10)
pred_svm <- predict(mod_svm, newdata = test_zscore)
confusionMatrix(pred_svm,as.factor(test_zscore$toxicity))
```

Confusion Matrix and Statistics
          Reference
Prediction nonToxic Toxic
  nonToxic      192    18
  Toxic           4    15
               Accuracy : 0.9039          
                 95% CI : (0.8582, 0.9388)
    No Information Rate : 0.8559          
    P-Value [Acc > NIR] : 0.020010        
                  Kappa : 0.5271          
 Mcnemar's Test P-Value : 0.005578        
            Sensitivity : 0.9796          
            Specificity : 0.4545          
         Pos Pred Value : 0.9143          
         Neg Pred Value : 0.7895          
             Prevalence : 0.8559          
         Detection Rate : 0.8384          
   Detection Prevalence : 0.9170          
      Balanced Accuracy : 0.7171          
       'Positive' Class : nonToxic  


## Neural Network

```{R NN}
mod_nnet <- train(toxicity ~., data = train_zscore, method = "nnet")
pred_nnet <- predict(mod_nnet, newdata = test_zscore)
confusionMatrix(pred_nnet,as.factor(test_zscore$toxicity))
```

Confusion Matrix and Statistics
          Reference
Prediction nonToxic Toxic
  nonToxic      192     6
  Toxic           4    27
               Accuracy : 0.9563          
                 95% CI : (0.9212, 0.9789)
    No Information Rate : 0.8559          
    P-Value [Acc > NIR] : 7.413e-07       
                  Kappa : 0.8184          
 Mcnemar's Test P-Value : 0.7518          
            Sensitivity : 0.9796          
            Specificity : 0.8182          
         Pos Pred Value : 0.9697          
         Neg Pred Value : 0.8710          
             Prevalence : 0.8559          
         Detection Rate : 0.8384          
   Detection Prevalence : 0.8646          
      Balanced Accuracy : 0.8989          
       'Positive' Class : nonToxic 
       
       

