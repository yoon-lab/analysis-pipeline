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
load("www/mod.RData")


server <- function(input, output) {
  
  output$SummaryInput <- renderTable({
    data.frame("Descriptor" = c("CoreSize (nm)","HSize (nm)", "Surface charge (mV)", "Surface area (m2/g)","Hsf (eV)","Ec (eV)","Ev (eV)","Meo (eV)","Exposure Time","Exposure dose (ug/mL)","Assay","Cell species","Cell origin","Cell type"),
               "Values" = c(input$CoreSize,input$HSize,input$`Surface charge`, input$`Surface area`,input$Hsf,input$Ec,input$Ev,input$MeO,input$`Exposure Time`, input$`Exposure Dose`,input$Assay,input$`Cell species`,input$`Cell origin`,input$`Cell type`))
  },digits=2)
  
  output$Prediction <- renderTable({  
    table1 <- data.frame("Coresize" = input$CoreSize,
                         "Hsize"  = input$HSize,
                         "Surface_charge" = input$`Surface charge`,
                         "Surface_area"  = input$`Surface area`,
                         "Hsf"= input$Hsf,
                         "Ec"  = input$Ec,
                         "Ev" = input$Ev,
                         "MeO"  = input$MeO,
                         "Exposure_time"  = input$`Exposure Time`,
                         "Exposure_dose" = input$`Exposure Dose`,
                         "Assay" = input$Assay,
                         "Cell_species" = input$`Cell species`,
                         "Cell_origin" = input$`Cell origin`,
                         "Cell_type" = input$`Cell type`,
                         "Toxicity" = "UNKNOWN" )
    table1 <- rbind(train[1,],table1)
   # colnames(table1) <- colnames(train)
    table2 <- as.data.frame(predict(mod, table1))
    data.frame("Observed.Toxicity" = table1[2,15], "Predicted.Toxicity" = table2[2,1])
    
  })
 
  
}
