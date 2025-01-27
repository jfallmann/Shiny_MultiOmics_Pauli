
suppressPackageStartupMessages({
  library(InteractiveComplexHeatmap)
  library(ComplexHeatmap)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(bslib)
  library(EnhancedVolcano)
  library(shiny)
  library(DT)
  library(shinydashboard)
  library(shinyWidgets)
  library(plotly)
})

pdf(file = NULL)

#### DATA FILES READING ####

s1 <- fread("/srv/data/TableS1_fixed.tsv")
s2 <- fread("/srv/data/TableS2_fixed.tsv")
s3 <- fread("/srv/data/TableS3_fixed.tsv")
mat <- read.delim("/srv/data/Fig4_AC_QExHFX2.padj05.lfc1.3in1c_2in2c.limma.deg.ch_kmeans.k5.tsv", row.names = 1)
clusters <- mat[,19]
clusters <- factor(clusters, levels=c(1,2,3,4,5))
col_order <- c("Oocyte.1.1", "Oocyte.1.2", "Oocyte.1.3", "cell64.128.2.1", "cell64.128.2.2", "cell64.128.2.3", "Sphere.3.1", "Sphere.3.2", "Sphere.3.3", "Shield.4.1", "Shield.4.2", "Shield.4.3", "per75Epiboly.5.1", "per75Epiboly.5.2", "per75Epiboly.5.3", "dpf1.6.1", "dpf1.6.2", "dpf1.6.3")
volc_data <- read.delim("/srv/data/AC_QExHFX2.padj05.lfc0.any2.limma.deg.1.tsv")
uscs_data <- read.delim("/srv/data/mart_export.txt")

############################



ui <- dashboardPage(
  skin = "blue", 
  dashboardHeader(title = "Data Visualization App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Raw Data Exploration", tabName = "raw_data", icon = icon("table")),
      menuItem("Volcano & Table", tabName = "volcano_table", icon = icon("chart-bar")),
      menuItem("Heatmap & Volcano", tabName = "heatmap", icon = icon("th"))
    )
  ),
  dashboardBody(
    tabItems(
      # Tab 1: Raw Data Exploration
      tabItem(tabName = "raw_data",
              fluidRow(
                box(
                  title = "Gene Search", width = 4, solidHeader = TRUE, status = "primary",
                  textInput("genes", "Gene Expression Visualization (comma separated):", 
                            value = "", placeholder = "Gene1,Gene2,... (e.g. mdh2,hspe1)"),
                  textInput("gene", "Phosphorylation Visualization:", 
                            value = "", placeholder = "Search for a single gene (e.g. ctnnd1 )"),
                  selectInput("mean", "Choose Data to Display:", 
                              choices = c("Average", "Raw Data"), selected = "Raw Data"),
                  div(
                    style = "padding: 10px; background-color: #f9f9f9; border-radius: 5px; margin: 10px; color: #333;",
                    "For further analysis, you can click on the 'Name' column of every table to go to the corresponding genomic region in the UCSC browser.
                     You can also click on the 'Gene_stable_ID' to go to the corresponding Ensembl page."
                  )
                ),
                box(
                  title = "Supplementary Data", width = 8, solidHeader = TRUE, status = "info", style="overflow-x: auto;",
                  h4("Supplementary Data 1"),
                  plotOutput("s1_line"),
                  DTOutput("s1"),
                  h4("Supplementary Data 2"),
                  DTOutput("s2"),
                  h4("Supplementary Data 3"),
                  plotlyOutput("s3_line"),
                  DTOutput("s3")
                )
              )
      ),
      # Tab 2: Volcano & Table
      tabItem(tabName = "volcano_table",
              fluidRow(
                box(
                  title = "Settings", width = 4, solidHeader = TRUE, status = "primary",
                  selectInput("comparison", "Select Comparison:", 
                              choices = c("Oocyte - cell64.128" = "Oocyte.1_vs_cell64.128.2", 
                                          "cell64.128 - Sphere" = "cell64.128.2_vs_Sphere.3",  
                                          "Sphere - Shield" = "Sphere.3_vs_Shield.4", 
                                          "Shield - per75Epiboly" = "Shield.4_vs_per75Epiboly.5", 
                                          "per75Epiboly - dpf1" = "per75Epiboly.5_vs_dpf1.6"),
                              selected = "Oocyte.1_vs_cell64.128.2"),
                  sliderInput("pcutoff2", "FDR Threshold:", 
                              min = 0, max = 10, value = 1.3, step = 0.1),
                  sliderInput("fccutoff2", "FC Threshold:", 
                              min = 0, max = 10, value = 1, step = 0.1),
                  div(
                    style = "padding: 10px; background-color: #f9f9f9; border-radius: 5px; margin: 10px; color: #333;",
                    "For further analysis, you can click on the 'UCSC_Browser' column of every table to go to the corresponding genomic region in the UCSC genome browser.
                     You can also click on the 'Ensembl_Link' to go to the corresponding ENSEMBL page."
                  )
                ),
                box(
                  title = "Volcano Plot", width = 8, solidHeader = TRUE, status = "info",
                  plotlyOutput("volcano2", height = "400px")
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
                      originalHeatmapOutput("ht", title = NULL)
                    ),
                    column(
                      width = 6,
                      subHeatmapOutput("ht", title = NULL)
                    )
                  )
                )
              ),
              fluidRow(
                box(
                  title = "Gene Table", width = 12, solidHeader = TRUE, status = "primary",
                  DTOutput("htBrush", height = "300px")
                )
              )
      )
    )
  )
)

server <- function(input, output, session){
  
  highlight <- NULL
  
  brush_action = function(df, output) {
    row_index = unique(unlist(df$row_index))
    column_index = unique(unlist(df$column_index))
    row_label = unique(unlist(df$row_label))
    output$htBrush <- renderDT({
      if (!is.null(row_label)){
        brush<-data.frame(row_label)
        colnames(brush)<-c("Accession_ID")
        brush <- brush %>%
          separate(Accession_ID, into=c("Accession_ID", "Gene_ID"), sep="\\|")
        brush$Accession_ID <- paste0(
          '<a href="https://www.ensembl.org/id/', brush$Accession_ID, '" target="_blank">', 
          brush$Accession_ID, 
          '</a>'
        )
        
        datatable(brush, escape = F)
      }
    })
  }
  
  
  s1_by_gene <- reactive({
    gene <- strsplit(input$genes, ",")[[1]]
    gene <- trimws(gene)
    if (length(gene)==0){
      return(NULL)
    }
    subset(s1, Name %in% gene)
  })
  
  output$s1_line <- renderPlot({
    req(input$genes)
    data <- s1_by_gene()
    data_filter <- data[, c("Name", "Egg_1", "Egg_2" ,"Egg_3", "2_1"  , "2_2",   "2_3", "4_1" , "4_2" , "4_3" , "6_1" ,"6_2"  , "6_3" , "8_1" ,  "8_2"  ,"8_3", "24_1", "24_2",  "24_3" )]
    
    data_long <- data_filter %>%
      pivot_longer(cols = -Name, names_to = "Stage", values_to = "Value") %>%
      mutate(StageGroup = gsub("_.*", "", Stage)) %>% 
      group_by(StageGroup) %>%
      mutate(StageMean = mean(Value, na.rm = TRUE)) %>% 
      ungroup()
    
    data_long$Stage <- factor(data_long$Stage, levels=c("Egg_1", "Egg_2" ,"Egg_3", "2_1"  , "2_2",   "2_3", "4_1" , "4_2" , "4_3" , "6_1" ,"6_2"  , "6_3" , "8_1" ,  "8_2"  ,"8_3", "24_1", "24_2",  "24_3"))
    data_long$StageGroup <- factor(data_long$StageGroup, levels=c("Egg", "2", "4", "6", "8",  "24"))
    
    if (input$mean == "Average") {
      data_avg <- data_long %>%
        group_by(StageGroup, Name) %>%
        summarize(StageMean = mean(Value, na.rm = TRUE), .groups = "drop")
      
      plot <- ggplot(data_avg, 
                     aes(x = StageGroup, y = StageMean, color = as.factor(Name), group = Name)) +
        geom_line(size = 1) +
        geom_point() +
        labs(
          x = "Stage",
          y = "Expression",
          title = paste("Time Line showing the expression of genes", input$genes)
        ) +
        theme_minimal()
    } else if (input$mean == "Raw Data") {
      plot <- ggplot(data_long, 
                     aes(x = Stage, y = Value, color = as.factor(Name), group = Name)) +
        geom_line(size = 1) +
        geom_point() +
        labs(
          x = "Stage",
          y = "Expression",
          title = paste("Time Line showing the expression of genes", input$genes),
          color="Genes"
        ) +
        theme_minimal()
    } 
    return(plot)
  })
  
  output$s1 <- renderDT({
    
    req(input$genes)
    data<-s1_by_gene()
    

    colnames(uscs_data) <- c("Gene_stable_ID", "chr", "pos1", "pos2")
    data <- merge(data, uscs_data, by="Gene_stable_ID", all.x=T)
    
    data_cols <- data[, c("Description", "Name", "Gene_stable_ID", "Transcript_stable_ID", "Protein_stable_ID", "Transcript_name", "Gene_Synonym", "Gene_type", "Gene_description", "chr", "pos1", "pos2")]
    
    data_cols$Gene_stable_ID <- paste0(
      '<a href="https://www.ensembl.org/id/', data_cols$Gene_stable_ID, '" target="_blank">', 
      data_cols$Gene_stable_ID, 
      '</a>'
    )
    data_cols$Name <- paste0(
      '<a href="https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=danRer11&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr
',
      data_cols$chr,'%3A',data_cols$pos1,'-',data_cols$pos2  ,'" target="_blank">', data_cols$Name, '</a>'
    )
    data_cols <- data_cols[, c("Description", "Name", "Gene_stable_ID", "Transcript_stable_ID", "Protein_stable_ID", "Transcript_name", "Gene_Synonym", "Gene_type", "Gene_description")]
    
    datatable(data_cols, escape = F)
  })
  
#https://ucsc.vbc.ac.at/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=joerg.fallmann&hgS_otherUserSessionName=danRer11_Ribo_RNA_Nastya&nonVirtPosition=&position=chr'

  
  
  s2_by_gene <- reactive({
    gene <- strsplit(input$genes, ",")[[1]]
    gene <- trimws(gene)
    if (length(gene)==0){
      return(NULL)
    }
    subset(s2, Name %in% gene)
  })
  
  
  output$s2 <- renderDT({
    req(input$genes)
    data<-s2_by_gene()
    
    colnames(uscs_data) <- c("Gene_stable_ID", "chr", "pos1", "pos2")
    data <- merge(data, uscs_data, by="Gene_stable_ID", all.x=T)
    
    data_cols <- data[, c("GeneID", "Name", "ProteinID", "Gene_stable_ID", "Transcript_stable_ID", "Protein_stable_ID", "Transcript_name", "Gene_Synonym", "Gene_type", "Gene_description", "chr", "pos1", "pos2")]
    data_cols$Gene_stable_ID <- paste0(
      '<a href="https://www.ensembl.org/id/', data_cols$Gene_stable_ID, '" target="_blank">', 
      data_cols$Gene_stable_ID, 
      '</a>'
    )
    data_cols$Name <- paste0(
      '<a href="https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=danRer11&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr',
      data_cols$chr,'%3A',data_cols$pos1,'-',data_cols$pos2  ,'" target="_blank">', data_cols$Name, '</a>'
    )
    data_cols <- data_cols[, c("GeneID", "Name", "Gene_stable_ID", "ProteinID", "Transcript_stable_ID", "Protein_stable_ID", "Transcript_name", "Gene_Synonym", "Gene_type", "Gene_description")]
    
    datatable(data_cols, escape = F)
  })
  
  s3_by_gene <- reactive({
    if (input$gene == ""){
      return(NULL)
    }
    subset(s3, Gene_name == input$gene)
  })
  
  output$s3_line <- renderPlotly({
    req(input$gene)
    data <- s3_by_gene()
    
    temporal_data <- data[, c("Gene_name", "2_1", "2_2", "2_3", "4_1", "4_2", "4_3", "6_1", "6_2", "6_3", "24_1", "24_2", "24_3", "pepG")]
    temporal_data <- temporal_data %>%
      filter(grepl("p[1-9]", pepG, ignore.case = TRUE))
    
    row <- temporal_data[, 2:14]
    row_long <- row %>%
      pivot_longer(cols = -pepG, names_to = "Stage", values_to = "Value") %>%
      mutate(StageGroup = gsub("_.*", "", Stage)) # Extract stage group (e.g., 2, 4, 6, 24)
    
    row_long$Stage <- factor(row_long$Stage, levels = c("2_1", "2_2", "2_3", "4_1", "4_2", "4_3", "6_1", "6_2", "6_3", "24_1", "24_2", "24_3"))
    row_long$StageGroup <- factor(row_long$StageGroup, levels = c("2", "4", "6", "24"))
    
    if (input$mean == "Average") {
      row_avg <- row_long %>%
        group_by(StageGroup, pepG) %>%
        summarize(StageMean = mean(Value, na.rm = TRUE), .groups = "drop")
      
      plot <- ggplot(row_avg, aes(x = StageGroup, y = StageMean, color = as.factor(pepG), group = pepG)) +
        geom_line(size = 1) +
        geom_point() +
        labs(
          x = "Stage",
          y = "Expression",
          title = paste("Phosphosite time line for gene", input$gene),
          color="Phosphosites"
        ) +
        theme_minimal()
      plot <- plot + aes(text =pepG)
    } else if (input$mean == "Raw Data") {
      # Plot raw data
      plot <- ggplot(row_long, aes(x = Stage, y = Value, color = as.factor(pepG), group = pepG)) +
        geom_line(size = 1) +
        geom_point() +
        labs(
          x = "Stage",
          y = "Expression",
          title = paste("Phosphosite time line for gene", input$gene),
          color = "Phosphosites"
        ) +
        theme_minimal()
      plot <- plot + aes(text =pepG)
    }
    
    return(ggplotly(plot, tooltip="text")) 
  })
  
  output$s3 <- renderDT({
    req(input$gene)
    data<-s3_by_gene()
    
    colnames(uscs_data) <- c("Gene_stable_ID", "chr", "pos1", "pos2")
    data <- merge(data, uscs_data, by="Gene_stable_ID", all.x=T)
    data_cols <- data[, c("pepG", "ProtID", "Modifications", "Gene_stable_ID", "Gene_name", "foldChange_4h", "foldChange_6h", "foldChange_1d", "chr", "pos1", "pos2")]
    colnames(data_cols)[colnames(data_cols)=="Gene_name"] <- "Name"
    data_cols <- data_cols %>%
      filter(grepl("p1|p2", pepG, ignore.case = TRUE))
    data_cols$Gene_stable_ID <- paste0(
      '<a href="https://www.ensembl.org/id/', data_cols$Gene_stable_ID, '" target="_blank">', 
      data_cols$Gene_stable_ID, 
      '</a>'
    )
    data_cols$Name <- paste0(
      '<a href="https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=danRer11&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr',
      data_cols$chr,'%3A',data_cols$pos1,'-',data_cols$pos2  ,'" target="_blank">', data_cols$Name, '</a>'
    )
    data_cols <- data_cols[, c("pepG", "ProtID", "Modifications", "Gene_stable_ID", "Name", "foldChange_4h", "foldChange_6h", "foldChange_1d")]
    datatable(data_cols, escape = F)
  })
  
  
  filtered_v_data_reactive <- reactive({
    filtered_v_data$significant <- with(filtered_v_data, 
                                        ifelse(PValue < 10^(-input$pcutoff) & abs(logFC) > input$fccutoff, "Both",
                                               ifelse(PValue < 10^(-input$pcutoff), "p-value",
                                                      ifelse(abs(logFC) > input$fccutoff, "FC", "Not significant"))))
    return(filtered_v_data)
  })
  
  
  
  
  data_ht <- reactive({
    mat <- mat
  })
  
  observe({
    req(data_ht()) 
    mat <- data_ht() 
    
    clusters <- mat[,19]
    clusters <- factor(clusters, levels=c("1","2","3","4","5"))
    mat <- mat[,1:18]
    mat <- as.matrix(mat) # Convert to matrix
    ht <- draw(Heatmap(mat, show_row_dend = F, show_row_names = F, split = clusters, cluster_row_slices = F , column_order = col_order))
    makeInteractiveComplexHeatmap(
      input, output, session,
      ht,
      heatmap_id = "ht", 
      brush_action = brush_action
    )
  })
  
  
  
  significant_genes <- reactive({
    
    uscs_data <- read.delim("/srv/data/mart_export.txt")
    colnames(uscs_data) <- c("gene_id", "chr", "pos1", "pos2")
    volc_data <- merge(volc_data, uscs_data, by="gene_id", all.x=T)
    
    sel=input$comparison
    concp <- paste0(sel,"_adjp")
    columnp <- volc_data[[concp]]
    concf <- paste0(sel,"_LFC")
    columnf <- volc_data[[concf]]
    gene <- volc_data$oID
    ensembl_ids <- volc_data$gene_id
    chr <- volc_data$chr
    posi <- volc_data$pos1
    posf <- volc_data$pos2
    
    df <- data.frame(gene, ensembl_ids, chr, columnp, columnf, posi, posf)
    df <- na.omit(df)
    df <- subset(df, !(columnp==0 & columnf==0))
    colnames(df) <- c("gene","Ensembl_Link", "UCSC_Browser",  "FDR","FC", "Initial_Pos", "Final_Pos")
    
    pCutoff <- 10^(-input$pcutoff2)
    FCcutoff <- input$fccutoff2
    sigGenes <- df %>%
      filter(FDR < pCutoff & abs(FC) > FCcutoff)
    
    sigGenes$Ensembl_Link <- paste0(
      '<a href="https://www.ensembl.org/id/', sigGenes$Ensembl_Link, '" target="_blank">', 
      sigGenes$Ensembl_Link, 
      '</a>'
    )
    sigGenes$UCSC_Browser <- paste0(
      '<a href="https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=danRer11&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr',
      sigGenes$UCSC_Browser,'%3A',sigGenes$Initial_Pos,'-',sigGenes$Final_Pos  ,'" target="_blank">UCSC Browser</a>'
    )
    sigGenes <- sigGenes %>% select("gene", "Ensembl_Link", "UCSC_Browser", "FDR", "FC")
    return(sigGenes)
  })
  
  output$volcano2 <- renderPlotly({
    
    sel=input$comparison
    concp <- paste0(sel,"_adjp")
    columnp <- volc_data[[concp]]
    concf <- paste0(sel,"_LFC")
    columnf <- volc_data[[concf]]
    gene <- volc_data$oID
    
    df <- data.frame(gene,columnp, columnf)
    df <- na.omit(df)
    df <- subset(df, !(columnp==0 & columnf==0))
    colnames(df) <- c("gene", "FDR","FC")
    volc <- EnhancedVolcano(df, lab = df$gene, x = "FC", y = "FDR", pCutoff = 10^(-input$pcutoff2), FCcutoff = input$fccutoff2, selectLab = c("a"), title = "")
    
    volc <- volc + aes(text =gene)
    final <- ggplotly(volc + aes(x = FC, y = -log10(FDR)), tooltip="text")
    
    return(final)
  })
  
  output$significant <- renderDT({
    gene_data <- significant_genes()  
    datatable(gene_data, escape = F)  
  })
}

shinyApp(ui, server)
