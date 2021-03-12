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


ui <- dashboardPage(
  dashboardHeader(title = "QSAR metal NPs", titleWidth = 350),
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
                  height = 500,
                  title = "Input properties of NPs:",
                  selectInput("Metal", "Metal NPs:",
                              c("Au" = "Au",
                                "Ag" = "Ag")),
                  selectInput("Shape", "Shape:",
                              c("Nanorod" = "Nanorod",
                                "Sphere" = "Sphere",
                                "Hollow" = "Hollow")),
                  sliderInput("CoreSize", "Core diameter (nm):", 1, 100, 10),
                  sliderInput("HSize", "Hydrodynamic diameter (nm):", 1, 300, 50),
                  sliderInput("Zeta", "Zeta potential (mV):", -20, 20, 0)
                  
                ),
                box(
                  height = 500,
                  title = "Input experimental conditions:",
                  selectInput("Cell_line", "Cell line:",
                              c("HeLa" = "HeLa",
                                "HepG2" = "HepG2",
                                "BEAS-2B" = "BEAS-2B",
                                "A549" = "A549")),
                  selectInput("Assay", "Toxicity assay method:",
                              c("MTS" = "MTS",
                                "MTT" = "MTT",
                                "AlamarBlue" = "Alamar Blue",
                                "NRU" = "NRU")),
                  sliderInput("Time", "Exposure time (h):", 1, 96, 10),
                  sliderInput("Dose", "Concentration (ug/L):", 1, 10^3, 50),
                ),
                
                box(
                  height = 500,
                  title = "Model performance:",
                  tableOutput("Performance")
                ),
                
                
                box(height = 120, title = "Predicted toxicity:", tableOutput("Prediction")),
                
                box(height = 360, title = "Summary of input:", tableOutput("SummaryInput"))
              )
      ),
      # Introduction tab content
      tabItem(tabName = "Introduction",
              h2("A web-based app for predicting cytotoxicity of metal nanoparticles (i.e., Ag and Au)"),
              
              column(
                br(),
                p("Dataset for model development was published in:", strong("Trinh et al.  Environmental Science: Nano 5.8 (2018): 1902-1910."),style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                br(),
                tags$img(src="article.png",width="716px",height="606px"),
                width=8),
      )
    )
  )
)