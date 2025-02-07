
s1 <- fread("/srv/data/TableS1_fixed.tsv")
s2 <- fread("/srv/data/TableS2_fixed.tsv")
s3 <- fread("/srv/data/TableS3_fixed.tsv")
mat <- read.delim("/srv/data/Fig4_AC_QExHFX2.padj05.lfc1.3in1c_2in2c.limma.deg.ch_kmeans.k5.tsv", row.names = 1)
clusters <- mat[,19]
clusters <- factor(clusters, levels=c(1,2,3,4,5))
col_order <- c("Oocyte.1.1", "Oocyte.1.2", "Oocyte.1.3", "cell64.128.2.1", "cell64.128.2.2", "cell64.128.2.3", "Sphere.3.1", "Sphere.3.2", "Sphere.3.3", "Shield.4.1", "Shield.4.2", "Shield.4.3", "per75Epiboly.5.1", "per75Epiboly.5.2", "per75Epiboly.5.3", "dpf1.6.1", "dpf1.6.2", "dpf1.6.3")
volc_data <- read.delim("/srv/data/AC_QExHFX2.padj05.lfc0.any2.limma.deg.1.tsv")
uscs_data <- read.delim("/srv/data/mart_export.txt")
total_expression <- read.delim("/srv/data/TMT_total_protein.tsv.gz")
rnaseq <- read.delim("/srv/data/CountTable_TMM_normalized_CPM.tsv")


server <- function(input, output, session){
  
  highlight <- NULL
  
  suppressWarnings({ 
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
        brush <- brush %>%
          separate(Gene_ID, into = c("Gene_Name", "Protein_ID"), sep = "\\.")
        brush$Gene_Name <- str_to_title(brush$Gene_Name)
        
        datatable(brush, escape = F)
      }
    })
  }
  })
  
  
  s1_by_gene <- reactive({
    gene <- strsplit(input$genes, ",")[[1]]
    gene <- str_to_lower(trimws(gene))
    if (length(gene)==0){
      return(NULL)
    }
    subset(s1,  str_to_lower(Name) %in% gene | str_to_lower(Gene_stable_ID) %in% gene)
  })
  
  
  s2_by_gene <- reactive({
    gene <- strsplit(input$genes, ",")[[1]]
    gene <- str_to_lower(trimws(gene))
    if (length(gene)==0){
      return(NULL)
    }
    subset(s2, str_to_lower(Name) %in% gene | str_to_lower(Gene_stable_ID) %in% gene)
  })
  
  output$s1_line <- renderPlot({
    validate(need(input$genes, "Please select one or more protein(s) of interest"))
    
    
    # Get data
    data <- s1_by_gene()
    data_filter <- data[, c("Name", "Egg_1", "Egg_2" ,"Egg_3", "2_1"  , "2_2",   "2_3", "4_1" , "4_2" , "4_3" , 
                            "6_1" ,"6_2"  , "6_3" , "8_1" ,  "8_2"  ,"8_3", "24_1", "24_2",  "24_3" )]
    colnames(data_filter) <- c("Name", "0 hpf_1", "0 hpf_2" ,"0 hpf_3", "2 hpf_1"  , "2 hpf_2",   "2 hpf_3", "4 hpf_1" , "4 hpf_2" , "4 hpf_3" , 
                               "6 hpf_1" ,"6 hpf_2"  , "6 hpf_3" , "8 hpf_1" ,  "8 hpf_2"  ,"8 hpf_3", "24 hpf_1", "24 hpf_2",  "24 hpf_3" )
    # Convert data to long format
    data_long <- data_filter %>%
      pivot_longer(cols = -Name, names_to = "Stage", values_to = "Expression") %>%
      mutate(StageGroup = gsub("_.*", "", Stage)) 
    
    # Convert factors for correct ordering
    data_long$StageGroup <- factor(data_long$StageGroup, 
                                   levels=c("0 hpf", "2 hpf", "4 hpf", "6 hpf", "8 hpf", "24 hpf"))
    colnames(data_long)
    # Compute mean expression per stage group
    data_avg <- data_long %>%
      group_by(StageGroup, Name) %>%
      summarize(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop")
    
    # Generate the plot
    if (input$scale == "Log-scale") {
      plot <- ggplot() +
        # Plot individual replicates as dots (aligned vertically)
        geom_point(data = data_long, 
                   aes(x = StageGroup, y = Expression, color = as.factor(str_to_title(Name))), 
                   size = 3, alpha = 0.6) +  # Removed jitter to ensure vertical alignment
        
        # Plot the mean expression line
        geom_line(data = data_avg, 
                  aes(x = StageGroup, y = MeanExpression, color = as.factor(str_to_title(Name)), group = str_to_title(Name)), 
                  linewidth = 1.2) +
        
        # Add mean expression points (triangles)
        geom_point(data = data_avg, 
                   aes(x = StageGroup, y = MeanExpression, color = as.factor(str_to_title(Name))), 
                   size = 4, shape = 17) +
        # Labels and theme
        labs(
          x = "Stage",
          y = "log Expression",
          title = paste("Protein Expression Time Line", str_to_title(input$genes)),
          color = "Proteins"
        ) +
        theme_minimal() +
        theme(text = element_text(size = 14)) +
        scale_y_log10()
    } else {
    plot <- ggplot() +
      # Plot individual replicates as dots (aligned vertically)
      geom_point(data = data_long, 
                 aes(x = StageGroup, y = Expression, color = as.factor(str_to_title(Name))), 
                 size = 3, alpha = 0.6) +  # Removed jitter to ensure vertical alignment
      
      # Plot the mean expression line
      geom_line(data = data_avg, 
                aes(x = StageGroup, y = MeanExpression, color = as.factor(str_to_title(Name)), group = str_to_title(Name)), 
                linewidth = 1.2) +
      
      # Add mean expression points (triangles)
      geom_point(data = data_avg, 
                 aes(x = StageGroup, y = MeanExpression, color = as.factor(str_to_title(Name))), 
                 size = 4, shape = 17) +
      # Labels and theme
      labs(
        x = "Stage",
        y = "Expression",
        title = paste("Time Line of Protein Expression", str_to_title(input$genes)),
        color = "Proteins"
      ) +
      theme_minimal() +
      theme(text = element_text(size = 14))
    
    }
    
    return(plot)
  })
  output$s1 <- renderDT({
    
    req(input$genes)
    data<-s1_by_gene()
    data2 <- s2_by_gene()
    
    colnames(uscs_data) <- c("Gene_stable_ID", "Gene_name","chr", "pos1", "pos2")
    data <- merge(data, uscs_data, by="Gene_stable_ID", all.x=T)
    data <- cbind(data, data2)
    data_cols <- data[, c("Description", "Name", "Gene_stable_ID", "Transcript_stable_ID", "Protein_stable_ID", "Transcript_name", "Gene_Synonym", "Gene_type", "Gene_description", "chr", "pos1", "pos2",
                          "Egg_1", "Egg_2", "Egg_3", "2_1", "2_2", "2_3", "4_1", "4_2", "4_3", "6_1", "6_2", "6_3", "8_1", "8_2", "8_3", "24_1", "24_2", "24_3", 
                          "LFC_2hpf/Egg", "LFC_4hpf/Egg", "LFC_6hpf/Egg", "LFC_8hpf/Egg" ,"LFC_24hpf/Egg","Padj_2hpf/Egg", "Padj_4hpf/Egg", "Padj_6hpf/Egg", "Padj_8hpf/Egg", "Padj_24hpf/Egg")]
    data_cols$Gene_stable_ID <- paste0(
      '<a href="https://www.ensembl.org/id/', data_cols$Gene_stable_ID, '" target="_blank">', 
      data_cols$Gene_stable_ID, 
      '</a>'
    )
    data_cols$Name <- paste0(
      '<a href="https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=danRer11&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr',
      data_cols$chr,'%3A',data_cols$pos1,'-',data_cols$pos2  ,'" target="_blank">', str_to_title(data_cols$Name), '</a>'
    )

    data_cols <- data_cols[, c( "Name", "Gene_stable_ID", "Transcript_stable_ID", "Protein_stable_ID", 
                               "2_1", "2_2", "2_3", "4_1", "4_2", "4_3", "6_1", "6_2", "6_3", "8_1", "8_2", "8_3", "24_1", "24_2", "24_3", 
                               "LFC_2hpf/Egg", "LFC_4hpf/Egg", "LFC_6hpf/Egg", "LFC_8hpf/Egg" ,"LFC_24hpf/Egg","Padj_2hpf/Egg", "Padj_4hpf/Egg", "Padj_6hpf/Egg", "Padj_8hpf/Egg", "Padj_24hpf/Egg")]
    
    datatable(data_cols, escape = F, extensions = "Buttons", options = list(scrollX=T, paging=T, dom='tB', buttons = c("copy","csv", "excel")))
  })
  
  #https://ucsc.vbc.ac.at/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=joerg.fallmann&hgS_otherUserSessionName=danRer11_Ribo_RNA_Nastya&nonVirtPosition=&position=chr'
  
  

  
  s3_by_gene <- reactive({
    if (input$genes == ""){
      return(NULL)
    }
    gene <- strsplit(input$genes, ",")[[1]]
    gene <- str_to_lower(trimws(gene))
    subset(s3,  str_to_lower(Gene_name) %in% gene | str_to_lower(Gene_stable_ID) %in% gene)
  })
  
  output$s3_line <- renderPlotly({
    validate(need(input$genes, "Please select one or more protein(s) of interest"))
    data <- s3_by_gene()
    data$Gene_pepG <- paste(str_to_title(data$Gene_name), data$pepG, sep="_")
    temporal_data <- data[, c("Gene_name", "2_1", "2_2", "2_3", "4_1", "4_2", "4_3", "6_1", "6_2", "6_3", "24_1", "24_2", "24_3", "Gene_pepG")]
    colnames(temporal_data) <- c("Gene_name", "2 hpf_1", "2 hpf_2", "2 hpf_3", "4 hpf_1", "4 hpf_2", "4 hpf_3", "6 hpf_1", "6 hpf_2", "6 hpf_3", "24 hpf_1", "24 hpf_2", "24 hpf_3", "Gene_pepG")
    temporal_data <- temporal_data %>%
      filter(grepl("p[1-9]", Gene_pepG, ignore.case = TRUE))
    
    row <- temporal_data[, 2:14]
    row_long <- row %>%
      pivot_longer(cols = -Gene_pepG, names_to = "Stage", values_to = "Value") %>%
      mutate(StageGroup = gsub("_.*", "", Stage))
    
    row_long$Stage <- factor(row_long$Stage, levels = c("2_1", "2_2", "2_3", "4_1", "4_2", "4_3", "6_1", "6_2", "6_3", "24_1", "24_2", "24_3"))
    row_long$StageGroup <- factor(row_long$StageGroup, levels = c("2 hpf", "4 hpf", "6 hpf", "24 hpf"))
    
    row_avg <- row_long %>%
      group_by(StageGroup, Gene_pepG) %>%
        summarize(StageMean = mean(Value, na.rm = TRUE), .groups = "drop")
    
    if ( input$scale == "Log-scale") {
      
      plot <- ggplot() +
        geom_point(data=row_long,
                   aes(x=StageGroup, y=Value, color=as.factor(Gene_pepG)),
                   size=2, alpha=0.6) +
        geom_line(data=row_avg,
                  aes(x=StageGroup, y=StageMean, color=as.factor(Gene_pepG), group=Gene_pepG),
                  linewidth=1)+
        geom_point(data=row_avg,
                   aes(x=StageGroup, y=StageMean, color=as.factor(Gene_pepG)),
                   size=4, shape=17) +
        labs(
          x="Stage",
          y="Expression",
          title=paste("Time line of Phosphopeptide Expression", str_to_title(input$genes)),
          color="Phosphopeptides"
        ) +
        theme_minimal() +
        theme(text = element_text(size=14)) +
        scale_y_log10()
        
      plot <- plot + aes(text =Gene_pepG)
      
    } else {
      
      plot <- ggplot() +
        geom_point(data=row_long,
                   aes(x=StageGroup, y=Value, color=as.factor(Gene_pepG)),
                   size=2, alpha=0.6) +
        geom_line(data=row_avg,
                  aes(x=StageGroup, y=StageMean, color=as.factor(Gene_pepG), group=Gene_pepG),
                  linewidth=1)+
        geom_point(data=row_avg,
                   aes(x=StageGroup, y=StageMean, color=as.factor(Gene_pepG)),
                   size=4, shape=17) +
        labs(
          x="Stage",
          y="Expression",
          title=paste("Time line Phosphopeptide Expression", str_to_title(input$genes)),
          color="Phosphopeptides"
        ) +
        theme_minimal() +
        theme(text = element_text(size=14))
      
      plot <- plot + aes(text =Gene_pepG)
      
    }
    
    return(ggplotly(plot, tooltip="text")) 
  })
  
  output$s3 <- renderDT({
    req(input$genes)
    data<-s3_by_gene()
    data$Gene_pepG <- paste(str_to_title(data$Gene_name), data$pepG, sep="_")
    colnames(uscs_data) <- c("Gene_stable_ID","name", "chr", "pos1", "pos2")
    data <- merge(data, uscs_data, by="Gene_stable_ID", all.x=T)
    data_cols <- data[, c("Gene_pepG", "ProtID", "Modifications", "Gene_stable_ID", "Gene_name", "foldChange_4h", "foldChange_6h", "foldChange_1d", "chr", "pos1", "pos2", "2_1", "2_2", "2_3", "4_1", "4_2", "4_3", "6_1", "6_2", "6_3", "24_1", "24_2", "24_3")]
    colnames(data_cols)[colnames(data_cols)=="Gene_name"] <- "Name"
    data_cols <- data_cols %>%
      filter(grepl("p1|p2", Gene_pepG, ignore.case = TRUE))
    data_cols$Gene_stable_ID <- paste0(
      '<a href="https://www.ensembl.org/id/', data_cols$Gene_stable_ID, '" target="_blank">', 
      data_cols$Gene_stable_ID, 
      '</a>'
    )
    data_cols$Name <- paste0(
      '<a href="https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=danRer11&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr',
      data_cols$chr,'%3A',data_cols$pos1,'-',data_cols$pos2  ,'" target="_blank">', str_to_title(data_cols$Name), '</a>'
    )
    data_cols <- data_cols[, c("Gene_pepG", "ProtID", "Modifications", "Gene_stable_ID", "Name", "foldChange_4h", "foldChange_6h", "foldChange_1d", "2_1", "2_2", "2_3", "4_1", "4_2", "4_3", "6_1", "6_2", "6_3", "24_1", "24_2", "24_3")]
    datatable(data_cols, escape = F, extensions = "Buttons", options = list(scrollX=T, paging=T, dom='tB', buttons = c("copy","csv", "excel")))
  })
  
  
  filtered_v_data_reactive <- reactive({
    filtered_v_data$significant <- with(filtered_v_data, 
                                        ifelse(PValue < 10^(-input$pcutoff) & abs(logFC) > input$fccutoff, "Both",
                                               ifelse(PValue < 10^(-input$pcutoff), "p-value",
                                                      ifelse(abs(logFC) > input$fccutoff, "FC", "Not significant"))))
    return(filtered_v_data)
  })
  
  tot_exp_by_gene <- reactive({
    if (input$genes == ""){
      return(NULL)
    }
    gene <- strsplit(input$genes, ",")[[1]]
    gene <- str_to_lower(trimws(gene))
    subset(total_expression,  str_to_lower(Genes) %in% gene | str_to_lower(GeneID) %in% gene)
  })
  
  output$tot_exp_line <- renderPlot({
    validate(need(input$genes, "Please select one or more protein(s) of interest"))
    data <- tot_exp_by_gene()
    
    data_filter <- data[,c("Genes","Egg", "X2hpf", "X4hpf", "X6hpf", "X1dpf")]
    colnames(data_filter) <- c("Genes", "0 hpf", "2 hpf", "4 hpf", "6 hpf", "24 hpf")
    tot_prot_long <- data_filter %>% 
      pivot_longer(cols = -Genes, names_to = "Stage", values_to = "Expression") %>%
      mutate(StageGroup = gsub("_.*","",Stage))
    
    tot_prot_long$StageGroup <- factor(tot_prot_long$StageGroup, levels=c("0 hpf", "2 hpf", "4 hpf", "6 hpf", "24 hpf"))
    
    if (input$scale == "Log-scale") {
      plot <- ggplot(tot_prot_long, aes(x = StageGroup, y = Expression, group = str_to_title(Genes), color = str_to_title(Genes))) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 4, shape = 17) +
        labs(
          x = "Stage",
          y = "Expression",
          title = paste("Total proteomics normalized to ribosomal proteins of", paste(str_to_title(input$genes), collapse = ", ")),
          color = "Proteins"
        ) +
        theme_minimal() +
        theme(text = element_text(size = 14))+
        scale_y_log10()
    } else {
      plot <- ggplot(tot_prot_long, aes(x = StageGroup, y = Expression, group = str_to_title(Genes), color = str_to_title(Genes))) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 4, shape = 17) +
        labs(
          x = "Stage",
          y = "Expression",
          title = paste("Total proteomics normalized to ribosomal proteins of", paste(str_to_title(input$genes), collapse = ", ")),
          color = "Proteins"
        ) +
        theme_minimal() +
        theme(text = element_text(size = 14))
    
    }
    
    return(plot)
  })
  
  output$tot_prot <- renderDT({
    req(input$genes)
    data <- tot_exp_by_gene()
    
    colnames(uscs_data) <- c("Gene_ID", "Genes", "chr", "pos1", "pos2")
    data <- merge(data, uscs_data, by="Genes", all.x=T)
    data_cols <- data[, c("Gene_ID", "Genes", "Accession", "Protein_group_accessions", "Description", "chr", "pos1", "pos2", "Egg", "X1hpf", "X16cell", "X2hpf", "X3hpf", "Oblong", "X4hpf", "X6hpf", "X1dpf", "Bridge")]
    data_cols$Gene_ID <- paste0(
      '<a href="https://www.ensembl.org/id/', data_cols$Gene_ID, '" target="_blank">', 
      data_cols$Gene_ID, 
      '</a>'
    )
    data_cols$Genes <- paste0(
      '<a href="https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=danRer11&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr',
      data_cols$chr,'%3A',data_cols$pos1,'-',data_cols$pos2  ,'" target="_blank">', data_cols$Genes, '</a>'
    )
    data_cols <- data_cols[, c("Gene_ID", "Genes", "Accession", "Protein_group_accessions", "Description" , "Egg", "X1hpf", "X16cell", "X2hpf", "X3hpf", "Oblong", "X4hpf", "X6hpf", "X1dpf", "Bridge")]
    
    datatable(data_cols, escape=F, extensions = "Buttons", options = list(scrollX=T, paging=T, dom='tB', buttons = c("copy","csv", "excel")))
  })
  
  
  tot_by_gene <- reactive({
    if (input$genes == ""){
      return(NULL)
    }
    gene <- strsplit(input$genes, ",")[[1]]
    gene <- str_to_lower(trimws(gene))
    subset(rnaseq,  str_to_lower(Gene.name) %in% gene | str_to_lower(GeneID) %in% gene)
  })
  
  
  output$rnaseq_line <- renderPlot({
    validate(need(input$genes, "Please select genes of interest"))
    data <- tot_by_gene()
    
    data_filter <- data[,c("Gene.name","X1hpf", "X2hpf", "X4hpf", "X8hpf", "X28hpf")]
    colnames(data_filter) <- c("Gene.name", "1 hpf", "2 hpf", "4 hpf", "8 hpf", "28 hpf")
    tot_exp_long <- data_filter %>% 
      pivot_longer(cols = -Gene.name, names_to = "Stage", values_to = "Expression") %>%
      mutate(StageGroup = gsub("_.*","",Stage))
    
    tot_exp_long$StageGroup <- factor(tot_exp_long$StageGroup, levels=c("1 hpf", "2 hpf", "4 hpf", "8 hpf", "28 hpf"))
    
    if (input$scale == "Log-scale") {
      plot <- ggplot(tot_exp_long, aes(x = StageGroup, y = Expression, group = Gene.name, color = Gene.name)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 4, shape = 17) +
        labs(
          x = "Stage",
          y = "TMM normalized CPM",
          title = paste("Total expression (RNA) of", paste(input$genes, collapse = ", ")),
          color = "Genes"
        ) +
        theme_minimal() +
        theme(text = element_text(size = 14))+
        scale_y_log10()
    } else {
      plot <- ggplot(tot_exp_long, aes(x = StageGroup, y = Expression, group = Gene.name, color = Gene.name)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 4, shape = 17) +
        labs(
          x = "Stage",
          y = "TMM normalized CPM",
          title = paste("Total expression (RNA) of", paste(input$genes, collapse = ", ")),
          color = "Genes"
        ) +
        theme_minimal() +
        theme(text = element_text(size = 14))
      
    }
    
    return(plot)
  })
  
  output$tot_exp_tab <- renderDT({
    req(input$genes)
    data <- tot_by_gene()
    
    colnames(uscs_data) <- c("GeneID", "Gene_name","chr", "pos1", "pos2")
    data <- merge(data, uscs_data, by="GeneID", all.x=T)
    data_cols <- data[, c("GeneID", "Gene.name", "Gene.type", "Gene.Synonym", "chr", "pos1", "pos2", "X1hpf", "X2hpf", "X3hpf", "X4hpf", "X8hpf", "X10hpf", "X28hpf")]
    data_cols$GeneID <- paste0(
      '<a href="https://www.ensembl.org/id/', data_cols$GeneID, '" target="_blank">', 
      data_cols$GeneID, 
      '</a>'
    )
    data_cols$Gene.name <- paste0(
      '<a href="https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=danRer11&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr',
      data_cols$chr,'%3A',data_cols$pos1,'-',data_cols$pos2  ,'" target="_blank">', data_cols$Gene.name, '</a>'
    )
    data_cols <- data_cols[, c("GeneID", "Gene.name", "Gene.type", "Gene.Synonym", "X1hpf", "X2hpf", "X3hpf", "X4hpf", "X8hpf", "X10hpf", "X28hpf")]
    
    datatable(data_cols, escape=F, extensions = "Buttons", options = list(scrollX=T, paging=T, dom='tB', buttons = c("copy","csv", "excel")))
  })
  
  
  data_ht <- reactive({
    mat <- mat
  })
  suppressWarnings({
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
  })
  
  
  
  significant_genes <- reactive({
    
    uscs_data <- read.delim("/srv/data/mart_export.txt")
    colnames(uscs_data) <- c("gene_id", "name", "chr", "pos1", "pos2")
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
    
    df <- df %>%
      separate(gene, into=c("Protein_ID", "Gene"), sep="\\|")
    df$Gene <- str_to_title(df$Gene)
    
    highlight_gene <- str_split(input$volcano_search, ",", simplify = TRUE) %>%
      str_trim() %>%
      str_to_title()
    
    col_custom <- rep("grey", nrow(df))
    
    # Apply color rules
    
    col_custom[df$FDR < 10^(-input$pcutoff2)] <- if (is.null(input$color_P)) "dodgerblue" else input$color_P
    col_custom[abs(df$FC) > input$fccutoff2] <- if (is.null(input$color_FC)) "limegreen" else input$color_FC
    col_custom[df$FDR < 10^(-input$pcutoff2) & abs(df$FC) > input$fccutoff2] <- if (is.null(input$color_both)) "hotpink" else input$color_both
    col_custom[df$Gene %in% highlight_gene] <- "black"
    
    # Handle NA values (if needed)
    col_custom[is.na(col_custom)] <- "black"
    
    
    names(col_custom)[col_custom == "black"] <- input$volcano_search
    names(col_custom)[col_custom == if (is.null(input$color_both)) "hotpink" else input$color_both] <- "P value and FC"
    names(col_custom)[col_custom == if (is.null(input$color_P)) "dodgerblue" else input$color_P] <- "P value"
    names(col_custom)[col_custom == if (is.null(input$color_FC)) "limegreen" else input$color_FC] <- "FC"
    names(col_custom)[col_custom == "grey"] <- "Non Significant"
    
    volc <- EnhancedVolcano(df, 
                            lab = df$Gene,
                            selectLab = c("a"),
                            x = "FC", 
                            y = "FDR", 
                            pCutoff = 10^(-input$pcutoff2), 
                            FCcutoff = input$fccutoff2, 
                            title = "",
                            colAlpha = 0.75,
                            pointSize = ifelse(
                              df$Gene %in% highlight_gene, 3, 1
                            ),
                            colCustom = col_custom,
                            legendPosition = 'none'
    )
    volc <- volc + aes(text =Gene) + labs(color = "Legend")   
    
    final <- ggplotly(volc + aes(x = FC, y = -log10(FDR)), tooltip="text")
    
    return(final)
  })
  
  output$significant <- renderDT({
    gene_data <- significant_genes() 
    gene_data <- gene_data %>%
      separate(gene, into=c("Protein_ID", "Gene"), sep="\\|")
    gene_data$Gene <- str_to_title(gene_data$Gene)
    datatable(gene_data, escape = F)  
  })
  
  output$downloadS1 <- downloadHandler(
    filename = function() {
      paste("Suplementary_data_1.csv")
    }, 
    content = function(file) {
      write.csv(s1, file, row.names = T)
    }
    
  )
  
  output$downloadS2 <- downloadHandler(
    filename = function() {
      paste("Suplementary_data_2.csv")
    }, 
    content = function(file) {
      write.csv(s2, file, row.names = T)
    }
    
  )
  output$downloadS3 <- downloadHandler(
    filename = function() {
      paste("Suplementary_data_3.csv")
    }, 
    content = function(file) {
      write.csv(s3, file, row.names = T)
    }
    
  )

  output$downloadTotProt <- downloadHandler(
    filename = function() {
      paste("Total_Proteomics.csv")
    }, 
    content = function(file) {
      write.csv(total_expression, file, row.names = T)
    }
    
  )
  
  output$downloadRNAseq <- downloadHandler(
    filename = function() {
      paste("RNA_Expression_TMM_normalized_CPM.csv")
    }, 
    content = function(file) {
      write.csv(rnaseq, file, row.names = T)
    }
    
  )
  
  help_general <- paste(readLines("/srv/docs/general_help.txt"), collapse = "\n")
  
  observeEvent(input$help, {
    showModal(modalDialog(
      title = "General Information",
      HTML(help_general),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  help_raw <- paste(readLines("/srv/docs/help_tab1.txt"), collapse = "\n")
  
  observeEvent(input$help1, {
    showModal(modalDialog(
      title = "",
      HTML(help_raw),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  help_volc <- paste(readLines("/srv/docs/help_tab2.txt"), collapse = "\n")
  
  observeEvent(input$help2, {
    showModal(modalDialog(
      title = "",
      HTML(help_volc),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  help_ht <- paste(readLines("/srv/docs/help_tab3.txt"), collapse = "\n")
  
  observeEvent(input$help3, {
    showModal(modalDialog(
      title = "",
      HTML(help_ht),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  
  observeEvent(input$legend, {
    showModal(modalDialog(
      div(style="text-align: center;",
      renderImage({
        list(src = "/srv/data/Timeline_Scetch.png",
             contentType="image/png",
             width=750,
             height=400)
      }, deleteFile = F)),
      size="l",
      easyClose = T,
      footer = NULL
    ))

  })
  

  
  observeEvent(input$color_pick, {
    showModal(
      modalDialog(
        title = "Change the palette",
        colorPickr("color_both", "Select a color for `Significant in P-val and FC`", selected="#FF69B4", preview = T, hue = T),
        colorPickr("color_FC", "Select a color for `Significant in FC`", selected="#32CD32", preview = T, hue = T),
        colorPickr("color_P", "Select a color for `Significant in P-val`", selected="#1E90FF", preview = T, hue = T),
        footer = tagList(
          modalButton("Close")
        )
      )
    )
  })
  
  
  
  
}
