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

ui <- dashboardPage(
  dashboardHeader(title = "Toxicity Classification of Oxide Nanomaterials", titleWidth = 350),
  dashboardSidebar(
    width = 350,
    sidebarMenu(
      menuItem("Introduction", tabName = "Introduction", icon = icon("dashboard")),
      menuItem("Model", tabName = "Model", icon = icon("th"))
    )
  ),
  dashboardBody(
    tabItems(
      # Model tab content
      tabItem(tabName = "Model",
              # Boxes need to be put in a row (or column)
              fluidRow(
                box(
                  height = 800,
                  title = "Input properties of NPs:",
                  
                  sliderInput("CoreSize", "Core diameter (nm):", 1, 500,25),
                  sliderInput("HSize", "Hydrodynamic diameter (nm):", 1, 1500, 50),
                  sliderInput("Surface charge", "Zeta potential (mV):", -60, 70, 0),
                  sliderInput("Surface area", "surface area (m2/g):", 0, 500, 50),
                  sliderInput("Hsf", "enthalpy of formation:", -30, 0, -5),
                  sliderInput("Ec", "conduction band:", -10, 0, -5),
                  sliderInput("Ev", "valence band energies:", -12, -3, -5),
                  sliderInput("MeO", "electronegativity:", 2, 10, 5),
                  sliderInput("Exposure Time", "Exposure time (h):", 2, 80, 5),
                  sliderInput("Exposure Dose", "Exposure dose (ug/mL):", 0, 1500, 50),
                  selectInput("Assay", "Assay:",
                              c("MTT" = "MTT",
                                "CellTiter-Glo" = "CellTiter-Glo",
                                "CCK-8" = "CCK-8",
                                "RTCA"="RTCA",
                                "LDH" = "LDH",
                                "Alamar blue" = "Alamar blue",
                                "Annexiv V/PI staining" = "Annexiv V/PI staining",
                                "ATP" = "ATP",
                                "MTS" = "MTS",
                                "CyQuant Assay" = "CyQuant Assay",
                                "Trypan blue" = "Trypan blue",
                                "Luminometric assay" = "Luminometric assay",
                                "PrestoBlue" = "PrestoBlue",
                                "Calcein AM" = "Calcein AM",
                                "WST" = "WST")),
                  selectInput("Cell species", "Cell species:",
                              c("Human" = "Human",
                                "Hamster" = "Hamster",
                                "Mouse" = "Mouse"
                                )),
                  selectInput("Cell origin", "Cell origin:",
                              c("Blood" = "Blood",
                                "Lung" = "Lung",
                                "Mesothelium" = "Mesothelium",
                                "Skin" = "Skin",
                                "Breast" = "Breast",
                                "nose" = "nose",
                                "Liver" = "Liver",
                                "Cervix" = "Cervix",
                                "Bone" = "Bone",
                                "Adrenal gland" = "Adrenal gland",
                                "Adipose tissue" = "Adipose tissue"
                                )),
                  selectInput("Cell type", "Cell type:",
                              c("Normal" = "Normal",
                                "Cancer" = "Cancer"
                              ))
                  
                  
                  
              ) ,
                
                
                box(height = 120, title = "Predicted toxicity:", tableOutput("Prediction")),
                
                box(height = 600, title = "Summary of input:", tableOutput("SummaryInput"))
              )
      ),
      
      # Introduction tab content
      tabItem(tabName = "Introduction",
              h2("Toxicity Classification of Oxide Nanomaterials"),
              
              column(
                br(),
                p("Dataset for model development was published in:", strong("My Kieu Ha (2018). Toxicity Classification of Oxide Nanomaterials: Effects of Data Gap Filling and PChem Score-based Screening Approaches. Sci Rep 8, 3141."), style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                br(),
               
                width=8),
      )
    )
  )
)