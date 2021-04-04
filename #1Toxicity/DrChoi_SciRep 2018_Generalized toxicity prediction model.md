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

There are altogether 15 descriptors, 5 of which are categorical and 11 are numerical. First we change the categorical descriptors into factors. Then we normalize the remaining numerical descriptors using z-score and log transformation approach, as the paper indicates.

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
The paper uses four different algorithms for modeling. For each algorithm different normalization methods are used.

## Random Forest (use z-score preprocessed data as an example)
First is the random forest. Normalization method for this algorithm is z-score.

```{r random forest, echo=TRUE}
mod_rf <- train(toxicity~., data=train_zscore, method='rf', metric='Accuracy')
pred_rf <- predict(mod_rf, test_zscore) 
confusionMatrix(pred_rf,as.factor(test_zscore$toxicity))
```

These are the corresponding results.

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
Second is Generalized linear model. As the paper indicates, preprocessing method for this algorithm is log transformation, which minimizes the skewness of the data.

```{r GLM, echo=TRUE}
mod_fit <- train(toxicity ~.,  data=train_log, method="glm", family="binomial")
pred_glm <- predict(mod_fit, newdata=test_log)
confusionMatrix(pred_glm,as.factor(test_log$toxicity))
```
These are the results.

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
Third is the support vector machine, with z-score normalization method.

```{R SVM, echo=TRUE}
mod_svm <- train(toxicity ~., data =train_zscore, method = "svmLinear",preProcess = c("center", "scale"),tuneLength = 10)
pred_svm <- predict(mod_svm, newdata = test_zscore)
confusionMatrix(pred_svm,as.factor(test_zscore$toxicity))
```
Results:

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
Last is the neural Network, with z-score normalization method.

```{R NN}
mod_nnet <- train(toxicity ~., data = train_zscore, method = "nnet")
pred_nnet <- predict(mod_nnet, newdata = test_zscore)
confusionMatrix(pred_nnet,as.factor(test_zscore$toxicity))
```

Results:

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
       

In summary, this paper introduces a generalized toxicity modeling approach that could be used for all sorts of nanomaterials and all sorts of descriptors. It also introduces four modeling algorithms included in a single package (caret), which is convenient to use and easy to evaluate.

## Reference
Towards a generalized toxicity prediction model for oxide nanomaterials using integrated data from different sources. Choi, Jang-Sik et al. Sci Rep. 2018 Apr 17;8(1):6110. doi: 10.1038/s41598-018-24483-z.