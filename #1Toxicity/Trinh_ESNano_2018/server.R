#############################################
# Install packages if they are not installed yet
# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("caret", "randomForest", "data.table", "shiny", "shinythemes", "shinydashboard")
ipak(packages)

# Load library
library(caret) 
library(randomForest)
library(data.table)
library(shiny)
library(shinythemes)
library(shinydashboard)
#############################################

#Load data and trained model
load("www/train.RData")
load("www/test.RData")
load("www/RFmodel.RData")
load("www/Performance.RData")


server <- function(input, output) {
  
  output$Performance <- renderTable({
    Performance
  }, digits = 2)
  
  output$SummaryInput <- renderTable({
    data.frame("Descriptor" = c("Metal NPs","Core size (nm)", "Hydrodynamic diameter (nm)", "Zeta potential (mV)","Cell line","Toxic Assay","Exposure time (h)","Concentration (ug/L)"),
               "Values" = c(input$Metal,input$CoreSize,input$HSize, input$Zeta,input$Cell_line,input$Assay,input$Time,input$Dose))
  }, digits = 2)
  
  output$Prediction <- renderTable({
    table1 <- data.frame("Toxicity" = "UNKNOWN",
                         "Dose" = input$Dose,
                         "Assay" = input$Assay,
                         "Time" = input$Time,
                         "Species" = "Human",
                         "Cancer" = 1,
                         "Cell_Tissue" = "Lung",
                         "Cell_line" = input$Cell_line,
                         "SSA" = 20,
                         "Zeta" = input$Zeta,
                         "HSize" = input$HSize,
                         "CoreSize" = input$CoreSize,
                         "Coating" = "None",
                         "Shape" = "Sphere",
                         "Metal" = input$Metal)
    
    table2 <- as.data.frame(predict(RFmodel, table1))
    
    data.frame("Observed.Toxicity" = table1[1,1], "Predicted.Toxicity" = table2[1,1])
    
    
  })
}

