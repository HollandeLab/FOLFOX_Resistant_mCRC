##---- rownames annotData must be the same as colnames expr data.  
plotHeatmapExpr <-
  function(genes,
    exprData = logCPM,
    dataType = "logCPM",
    orgTum = "Org",
    treatmentOrg = "Control",
    responseOrg = c("Sensitive", "Resistant"),
    treatmentTum = NULL,
    responseTum = NULL,
    annotData = dgeFiltNorm$samples,
    annotColumns = c("PatientID", "Tissue", "Response"),
    geneAnnot = NULL, ### one column data frame with "Genes" as column name
    geneCol = NULL,
    clustRows = T,
    clustCols = T,
    showRowNames = T,
    showColNames = T,
    textSize = 10,
    plotTitle = "Heatmap") {
    
    ##-------------------------- Subset data for the genes and samples 
    ##                            then scale the expression matrix:
    
    gg <- genes[genes %in% rownames(exprData)]
    if (length(gg) == 0) {
      stop("Selected geens do not present in the expression data")
    }
    if(! all(rownames(annotData) == colnames(exprData))){
      stop("Make sure annotation data is in the same order as the expression data and has the same row names as the column names in the expression data")
    }
    
    if(orgTum == "Org"){
      subsetSample <- annotData$Treatment == treatmentOrg & annotData$Response  %in% responseOrg
    }
    if(orgTum == "Tum"){
      subsetSample <- annotData$Treatment == treatmentTum & annotData$Response_Tum  %in% responseTum
    }
    
    annotData <- annotData[subsetSample, , drop = F]
    ddgene <- exprData[gg, rownames(annotData), drop = F]
    ddgene <- t(scale(t(ddgene)))
    
    minExpr <- min(ddgene)
    maxExpr <- max(ddgene)
    
    ##-------------------------- Set up colours for th eheatmap:
    ##---- patients' colours
    colfunc <-
      colorRampPalette(c(
        brewer.pal(11, 'BrBG')[-c(1, 6, 11)],
        brewer.pal(11, 'PRGn')[-c(1, 6, 11)],
        brewer.pal(11, 'RdYlBu')
      ))
    
    patientCol  <- colfunc(length(unique(annotData$PatientID)))
    names(patientCol) <-
      unique(as.character(annotData$PatientID))[order(unique(as.character(annotData$PatientID)))]
    
    ##---- arrange samples in org data and assign treatment colour for org samples
    if (orgTum == "Org") {
      annotData <-
        annotData %>% arrange(Treatment_Response, PatientID, Tissue)
      ddgene <- ddgene[, annotData$UniqueIDs]
      
      treatmentCol <- brewer.pal(9, "Set1")[c(3, 4)]
      names(treatmentCol) <- c("Control", "FOLFOX")
    }
    
    ##---- arrange samples in tum data and assign treatment colour for tum samples
    if (orgTum == "Tum") {
      annotData <-
        annotData %>% arrange(Response_Tum, Treatment, PatientID, Tissue)
      ddgene <- ddgene[, annotData$UniqueIDs]
      
      treatmentCol <- brewer.pal(9, "Oranges")[c(2, 4, 6, 8, 9)]
      names(treatmentCol) <-
        c(
          "CAPOX",
          "FOLFOX",
          "FOLFOX + Bevacizumab",
          "FOLFIRI + Cetuximab",
          "FOLFOXIRI + Bevacizumab"
        )
    }
    
    respCol <- c("pink", "pink3", "darkred")
    names(respCol) <-
      c("Sensitive", "Semi-sensitive", "Resistant")
    
    respTumCol <- c("pink", "pink3", "darkred", "darkgoldenrod")
    names(respTumCol) <- c("CR", "PR", "PD", "Relapse")
    
    tissueCol <- c("gray20", "gray85")
    names(tissueCol) <- c("LT", "CoT")
    
    genderCol <- c("aquamarine", "aquamarine4")
    names(genderCol) <- c("Male", "Female")
    
    chemoCol  <- brewer.pal(8, "Blues")[c(3, 7, 5)]
    names(chemoCol) <- c("Naive", "Distant", "Treated")
    
    runCol <- brewer.pal(9, "Paired")[7:8]
    names(runCol) <- c("Run 1", "Run 2")
    
    siteCol <- brewer.pal(9, "Greens")[c(3, 7, 5)]
    names(siteCol) <- names(table(dge$samples$Site))
    
    synchMetachCol  <-
      c(brewer.pal(8, "Set2")[3], brewer.pal(8, "Set3")[3])
    names(synchMetachCol) <- names(table(dge$samples$SynchMetach))
    
    currentCol <-  list(
      PatientID = patientCol,
      Tissue = tissueCol,
      Gender = genderCol,
      Site = siteCol,
      SynchMetach = synchMetachCol,
      Treatment = treatmentCol,
      ChemoStatus = chemoCol,
      Run = runCol,
      Response_Tum = respTumCol,
      Response = respCol
    )
    
    ss <-
      sapply(names(currentCol), function(x)
        list(
          list(
            legend_direction = "horizontal",
            nrow = 2,
            title_gp = gpar(fontsize = textSize + 2, fontface = "italic"),
            labels_gp = gpar(fontsize = textSize), 
            grid_height = unit(0.8, "cm")
          )
        ))
    ss <- ss[annotColumns]
    
    sampleAnnotHm <-
      HeatmapAnnotation(df = annotData[, annotColumns, drop = F],
        col = currentCol,
        annotation_name_side = "left", 
        annotation_legend_param = ss,
        ## Added this:
        annotation_name_gp = gpar(fontsize = textSize))
    
    if(!is.null(geneAnnot)){
      geneAnnot <- geneAnnot[gg, , drop = F]
      geneHmAnn <- HeatmapAnnotation(df = geneAnnot,
        col =   list(Genes = geneCol),
        which = "row", name = "\nGenes")
      
      # geneHm <- Heatmap(geneAnnot,
      #   name = "Genes",
      #   col = geneCol,
      #   width = unit(5, "mm"), 
      #   show_row_names = showRowNames)
      
      expr_hm <- ComplexHeatmap::Heatmap(
        ddgene,
        col = circlize::colorRamp2(c(minExpr, maxExpr), c("yellow3" , "navy")),
        # name = "logExpr",
        column_title = plotTitle,
        cluster_rows = clustRows,
        cluster_columns = clustCols,
        show_row_names = showRowNames,
        show_column_names = showColNames,
        row_names_gp = gpar(fontsize = textSize),
        column_names_gp = gpar(fontsize = textSize),
        
        # column_title_gp = gpar(fontsize = textSize), 
        heatmap_legend_param = list(
          title_position = "topcenter",
          title = "Scaled log(CPM)",
          title_gp = gpar(fontsize =  textSize * 1.2, fontface = "italic"),
          labels_gp = gpar(fontsize = textSize),
          color_bar = "continuous", 
          legend_direction = "horizontal",
          legend_width = unit(5, "cm")
          # legend_height = unit(3, "cm")
        ),
        top_annotation = sampleAnnotHm,
        right_annotation = geneHmAnn,
        # left_annotation = geneHmAnn,
        na_col = "gray90"
      )
      
    } else{
      expr_hm <- ComplexHeatmap::Heatmap(
        ddgene,
        col = circlize::colorRamp2(c(minExpr, maxExpr), c("yellow3" , "navy")),
        # name = "logExpr",
        column_title = plotTitle,
        cluster_rows = clustRows,
        cluster_columns = clustCols,
        show_row_names = showRowNames,
        show_column_names = showColNames,
        row_names_gp = gpar(fontsize = textSize),
        column_names_gp = gpar(fontsize = textSize),
        # column_title_gp = gpar(fontsize = textSize), 
        heatmap_legend_param = list(
          title_position = "lefttop",
          title = "Scaled log(CPM)", 
          title_gp = gpar(fontsize = textSize * 1.2, fontface = "italic"),
          labels_gp = gpar(fontsize = textSize),
          color_bar = "continuous", 
          legend_direction = "horizontal", 
          legend_width = unit(5, "cm")),
        top_annotation = sampleAnnotHm, 
        na_col = "gray90"
      ) 
      
    }
    
    ComplexHeatmap::draw(
      expr_hm,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      legend_title_gp = gpar(fontsize = textSize, fontface = "italic")
    )
    
  }
