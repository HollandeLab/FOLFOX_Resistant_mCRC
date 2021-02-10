
library(matrixStats)
library(ggrepel)
library(data.table)
library(dplyr)
library(RColorBrewer)

RLE_ggplot <- function(exprData = exprData,
  annot = annot,
  largeData = T,
  textSize = 1.5,
  ylim = NULL,
  mainTitle = "RLE plot",
  cols =
    c(brewer.pal(8, "Dark2")[-5],
      brewer.pal(10, "Paired"),
      brewer.pal(12, "Set3"))
  
){

  if(!is.null(annot)){

    colnames(annot) <- "Annotation"
    # annot <- annot[order(annot$Annotation), , drop = F]
    # exprData <- exprData[, row.names(annot)]
    
    ## If the annotation data is not already factor, change it to be a factor
    if(!is.factor(annot$Annotation)){
      annot$Annotation <- factor(annot$Annotation,
        levels = unique(annot$Annotation), ordered = TRUE)
    }

    
    # cols <- c(brewer.pal(8, "Dark2")[-5], brewer.pal(10, "Paired"), brewer.pal(12, "Set3"))
  } else{
    annot = data.frame(Annotation = cbind(rep("", ncol(exprData))))
    row.names(annot) <- colnames(exprData)
    cols <- brewer.pal(8, "Accent")[2]
  }

  rle <- exprData - rowMedians(exprData)
  rleLong <- melt(rle)

  ## Here for merging data sets I use data table which is MUCH faster than data frame:
  annot$Var2 <- row.names(annot)
  rleLong <- merge(data.table::data.table(rleLong),
                   data.table::data.table(annot),
                   by = "Var2")


  ## Calculate the median of the RLE boxplots:
  rleLong <- rleLong %>%
    group_by(Var2) %>%
    mutate(MedRLE = median(value)) %>%
    ungroup() %>%
    data.frame()
  
  ## Make the order of the samples according to the annotation
  rleLong <- rleLong[order(rleLong$Annotation), ]
  rleLong$Var2 <- factor(rleLong$Var2, levels = unique(rleLong$Var2))
  
  ## Remove Whiskers of the boxplots if the sample size is very large
  if(largeData){
    p <- ggplot(rleLong, aes(x = Var2, y = value, fill = Annotation)) +
      geom_boxplot(outlier.shape = NA,
                   alpha = 0.4,
                   coef = 0
        # fatten = NULL
        ) +  ## remove median
      scale_x_discrete(name = paste0("Samples (n = ", ncol(exprData), ")")) +
      theme_classic() +
      theme(axis.text.x = element_blank())
  }else{
    p <- ggplot(rleLong, aes(
      x = Var2,
      y = value, 
      fill = Annotation)) +
      geom_boxplot(outlier.shape = NA,
                   alpha = 0.4
        # fatten = NULL
        ) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, size = rel(textSize)))
  }

  p <- p +
    geom_point(data = rleLong[! duplicated(rleLong$Var2), ],
               aes(x = Var2, y = MedRLE, fill = Annotation),
               size = 2,
               shape = 21,
               colour = "black",
               lwd = 2) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(name = "RLE") +
    scale_x_discrete(name = paste0("Samples (n = ", ncol(exprData), ")")) +
    geom_hline(yintercept = 0, col = "darkred", lwd = 1) +
    ggtitle(mainTitle) +
    theme(
      # axis.title.x = textSize,
      # axis.title = element_text(size = rel(textSize)),

      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),

      axis.text.y = element_text(angle = 0, size = rel(textSize)),

      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(size = rel(textSize)),

      axis.line = element_line(colour = "black", size = 1.5),
      # axis.ticks = element_line(),
      legend.position = "right",
      legend.direction = "vertical",

      legend.margin = margin(unit(0, "cm")),
      legend.title = element_text(size = rel(textSize * 0.8), face="italic"),
      legend.text = element_text(size = rel(textSize * 0.8)),
      plot.title = element_text(
        face = "bold",
        size = rel(textSize),
        hjust = 0.5
      )
    )
  if( ! largeData){
    p <- p + 
      theme(axis.text.x = element_text(angle = 90, size = rel(textSize)))
  }

  if (!is.null(ylim)){
    p <- p + scale_y_continuous(name = "RLE",
      limits = ylim)
  }

  print(p)

}





