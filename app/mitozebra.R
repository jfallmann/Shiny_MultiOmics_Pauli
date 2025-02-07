suppressPackageStartupMessages({
  library(InteractiveComplexHeatmap)
  library(ComplexHeatmap)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyverse)
  library(bslib)
  library(EnhancedVolcano)
  library(shiny)
  library(DT)
  library(shinydashboard)
  library(shinyWidgets)
  library(plotly)
  
})

pdf(file = NULL)

source("ui.R")
source("server.R")

shinyApp(ui=ui, server=server)