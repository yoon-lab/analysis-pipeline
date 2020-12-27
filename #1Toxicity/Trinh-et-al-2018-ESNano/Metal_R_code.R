library(openxlsx) 
library(clusterSim) 
library(caret) 
library(randomForest) 
library(PresenceAbsence)

#Read Dataset C sheet as DataFrame 
DataFile <- file.choose(); DataFolder <- dirname(DataFile)
DataMetal <- read.xlsx(DataFile, sheet = 1, startRow = 1, colNames = TRUE,
                       rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = FALSE,
                       namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
score_org <- read.xlsx(DataFile, sheet = 1, colNames = TRUE) 

#Convert character values into factor 
for (i in c(3:5,18,19,21,23,26)) { 
  score_org[,i] <- as.factor(score_org[,i])
} 
rm(i) 


# Normalize numeric data 
# n1 means standardization, Standardize all numeric values 
for (i in c(6,9,12,15,22,24)) { 
  score_org[,i] <- as.numeric(score_org[,i])
} 
rm(i) 

preObj <- score_org[,c(6,9,12,15,22,24)] 
preObj <- data.Normalization(preObj, type = "n1", normalization = "column")

# Combine normalized data into dataset 
score.pre <- data.frame(preObj, score_org[,c(3:5,18,19,21,23,26)]) 

# Split data into training and test sets 
# createDataPartition function does startified sampling 
inTrain <- createDataPartition(y = score.pre$Toxicity, p = 0.6, list = FALSE) 
train <- score.pre[inTrain,] 
test <- score.pre[-inTrain,] 

# Fit a model on the training set to predict toxicity based on all other attributes using random forest algorithm 
mod <- randomForest(Toxicity ~ ., data = train, importance = T, proximity = T)

# Apply the model on test set to predict toxicity 
pred <- predict(mod, test) 
pred.prob <-predict(mod, test, type="prob") 

# Show the confusion matrix 
table(test$Toxicity, pred) 
(accuracy <- sum(pred==test$Toxicity) / nrow(test) * 100) 

# Build ROC plot & calculate AUC 
score_org$prob <- ifelse(score_org$Toxicity == "Nontoxic",0,1) 
ObsPred <- data.frame(Material=score_org[-inTrain,]$Chemical.composition, 
                      Probability=score_org[-inTrain,]$prob, 
                      Prediction=pred.prob[,2]) 

auc.roc.plot(ObsPred, which.model = 1, model.names = "All.attributes", 
             color = FALSE, line.type = "solid", threshold = 101, 
             xlab = "1 - Specificity (false positives)", 
             ylab = "Sensitivity (true positive)", main = "(f)", 
             find.auc = TRUE, opt.thresholds = FALSE, opt.methods = 9)

#Plot error rate as the number of tree increase 
plot(mod) 
legend("top", colnames(mod$err.rate), col = 1:3, cex = 0.8, fill = 1:3) 

#value important 
importance(mod) 
varImpPlot(mod) 

#margin graph 
plot(margin(mod, train$Toxicity)) 
