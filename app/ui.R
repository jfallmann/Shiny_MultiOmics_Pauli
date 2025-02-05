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



ui <- dashboardPage(
  skin = "blue", 
  dashboardHeader(title = "MitoZebra"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Exploration (Mt and Total)", tabName = "raw_data", icon = icon("table")),
      menuItem("Mt-Proteomics DEP", tabName = "volcano_table", icon = icon("chart-bar")),
      menuItem("Interactive figure 4A (Heatmap)", tabName = "heatmap", icon = icon("th")),
      menuItem("Downloads", tabName = "downloads", icon = icon("download")),
      actionBttn("help", "Help", icon=icon("question-circle"), size="sm", style = "bordered"), 
      h2(),
      actionBttn("legend", "Embryo Stage Legend", icon=icon("fish-fins"), size="xs", style="minimal")
      
    )
  ),
  dashboardBody(
    tabItems(
      # Tab 1: Raw Data Exploration
      tabItem(tabName = "raw_data",
              fluidRow(
                box(title="Search Bar", width=5, solidHeader = TRUE, status="primary", style="overflow-x: auto",
                    textInput("genes", "Search your gene of interest:", value="", placeholder = "e.g. mdh2,hspe1"),
                    selectInput("scale", "Select scale:", choices = c("Continous", "Log-scale"), selected = "Continous"),
                    h3(),
                    actionBttn("help1", "Help", icon = icon("question-circle"), color = "default",style="jelly", size="md")),
                box(                  
                  title = "Plots", width = 12, solidHeader = TRUE, status = "info", style="overflow-x: auto;",
                  h3(),
                  fluidRow(
                    box(title="Mt Proteomics", width=12, solidHeader = T, status="primary", style="overflow-x: auto",
                        plotOutput("s1_line")
                    )
                  ), 
                  fluidRow(
                    box(title="Mt Phosphoproteomics", width=12, solidHeader = T, status="primary", style="overflow-x: auto",
                        plotlyOutput("s3_line")
                    )
                  ),
                  fluidRow(
                    box(title="Total Proteomics", width=12, solidHeader = T, status="primary", style="overflow-x: auto",
                        plotOutput("tot_exp_line")
                    )
                  ),
                  fluidRow(
                    box(title="Total Gene Expression", width=12, solidHeader = T, status="primary", style="overflow-x: auto",
                        plotOutput("rnaseq_line")
                    )
                  )
                )
              ),
              fluidRow(
                box(
                  title = "Tables", width = 12, solidHeader = TRUE, status = "warning",
                  h3("Mt Proteomics Table"),
                  DTOutput("s1"),
                  h3("Mt Phosphoproteomics Table"),
                  DTOutput("s3"), 
                  h3("Total Proteomics Table"), 
                  DTOutput("tot_prot"),
                  h3("Total Expression Table"),
                  DTOutput("tot_exp_tab")
                  
                )
              )
              
      ),
      # Tab 2: Volcano & Table
      tabItem(tabName = "volcano_table",
              fluidRow(
                box(
                  title = "Settings", width = 4, solidHeader = TRUE, status = "primary",
                  virtualSelectInput("comparison", "Select Comparison:", 
                              choices = c("0hpf - 2hpf" = "Oocyte.1_vs_cell64.128.2", 
                                          "2hpf - 4hpf" = "cell64.128.2_vs_Sphere.3",  
                                          "4hpf - 6hpf" = "Sphere.3_vs_Shield.4", 
                                          "6hpf - 8hpf" = "Shield.4_vs_per75Epiboly.5", 
                                          "8hpf - 24hpf" = "per75Epiboly.5_vs_dpf1.6"),
                              selected = "Oocyte.1_vs_cell64.128.2"),
                  textInput("volcano_search", "Search for Gene:", value = "", placeholder = "Enter gene name"),
                  chooseSliderSkin("Flat", color="#3d8dbc"),
                  sliderInput("pcutoff2", "FDR Threshold:", 
                              min = 0, max = 10, value = 1.3, step = 0.1),
                  sliderInput("fccutoff2", "FC Threshold:", 
                              min = 0, max = 10, value = 1, step = 0.1),

                  h5(),
                  actionBttn("help2", "Help", icon = icon("question-circle"), color = "default",style="jelly", size="md")
                  
                ),
                box(
                  title = "Volcano Plot", width = 6, solidHeader = TRUE, status = "info",
                  plotlyOutput("volcano2", height = "500px"),
                  actionBttn("color_pick", "Change the palette", icon = icon("palette"), style="material-circle", size="md")
                )
              ),
              fluidRow(
                box(
                  title = "Significant Genes Table", width = 12, solidHeader = TRUE, status = "warning",
                  DTOutput("significant", height = "400px")
                )
              )
      ),
      # Tab 3: Heatmap & Volcano
      tabItem(tabName = "heatmap",
              fluidRow(
                box(
                  title = "Heatmaps", width = 12, solidHeader = TRUE, status = "primary",
                  fluidRow(
                    column(
                      width = 6,
                      originalHeatmapOutput("ht", title = NULL, width = 550, height = 450)
                    ),
                    column(
                      width = 6,
                      subHeatmapOutput("ht", title = NULL, width = 550, height = 450)
                    )
                  ),
                  h5(),
                  actionBttn("help3", "Help", icon = icon("question-circle"), color = "default",style="jelly", size="md")
                  
                )
              ),
              fluidRow(
                box(
                  title = "Gene Table", width = 12, solidHeader = TRUE, status = "warning",
                  DTOutput("htBrush", height = "300px")
                )
              )
      ),
      # Tab 4: Download possibilities
      tabItem(tabName = "downloads", 
              fluidRow(
                box(
                  title = "Downloads",
                  downloadBttn("downloadS1", "Download Mt Proteomics (S1)", size="sm"),
                  h5(),
                  downloadBttn("downloadS2", "Download Mt Proteomics (S2)", size="sm"),
                  h5(),
                  downloadBttn("downloadS3", "Download Mt Phosphoproteomics (S3)", size="sm"),
                  h5(),
                  downloadBttn("downloadTotProt", "Download Total Proteomics", size = "sm"), 
                  h5(), 
                  downloadBttn("downloadRNAseq", "Download Total Expression", size = "sm")
                  
                )
              )
              
      )
    )
  )
)

