
##############################################################################################################
################################### Author: Momeneh (Sepideh) Foroutan  ######################################
##############################################################################################################

## Before running this code:
## 1. Perform SVD: for svd genes should be in rows and samples in cols:
## s<- svd(apply(data , 1, function(x) scale(x, scale=FALSE , center = TRUE)))  

## 2. Read in the clin data having batch or subtypes info; the order of teh samples in clinical data should be the same in expr data
## match_indx<- match(colnames(expr), clin$SampleID)

##==================================================== INPUTS:

## expr: expression data
## clin: clinical data containing info regarding biology and batch-- NOTE: It has to have a column called "SampleID"
## group_to_test: column names from the clin data based on which we colour the samples in PCA plot
## data: a character vector defining the name of the output PCA plot 
## svdResult: result of svd
## labeled: If TRUE, then SampleIDs will be labeled on teh PC1 vs PC2
## group_to_shape: column name in the clin data to shape the points on teh PC1 vs PC2
## out_path: a path to where the output should be saved

##===================================================== OUTPUT:
## Several PCA plots 

#############################################################################################################

PCA_plot <- function(expr,
  clin,
  group_to_test,
  # data = "data",
  svdResult = s,
  plot_cex = 1,
  legend_cex = 1,
  labeled = FALSE,
  group_to_shape = NULL,
  cols =
    c(brewer.pal(8, "Dark2")[-5],
      brewer.pal(10, "Paired"),
      brewer.pal(12, "Set3"))
  # out_path
  ) {
  ## both arguments should be character.
  
  ## Calculate percentage of the variations for each PCs
  percent <- svdResult$d ^ 2 / sum(svdResult$d ^ 2) * 100
  labs <- sapply(seq_along(percent), function(i)
    paste("PC ", i, " (", round(percent[i], 2), "%)", sep = ""))
  
  if (length(grep("SampleID", colnames(clin))) == 0) {
    stop(
      "clin data HAS TO have a column named 'SampleID' in order to match sample IDs in expression and clin data"
    )
    
  }
  
  match_indx <- match(colnames(expr), clin$SampleID)
  
  if (!is.null(group_to_shape)) {
    group_shape <- clin[, group_to_shape][match_indx]
  }
  
  
  # pdf(paste0(out_path, "PCA_1_2_3_", data, ".pdf"),
  #   width = 12,
  #   height = 12)
  par(mar = c(6, 6, 2, 1.5))
  par(mfrow = c(2, 2))
  
  for (g in 1:length(group_to_test)) {
    group <- clin[, group_to_test[g]][match_indx]
    group <- as.factor(group)
    
    
    # cols <-
    #   c(brewer.pal(8, "Dark2")[-5],
    #     brewer.pal(10, "Paired"),
    #     brewer.pal(12, "Set3"))
    
    x <- 20
    y <- 20
    
    if (!is.null(group_to_shape)) {
      plot(
        svdResult$u[, 1:2],
        col = cols[as.numeric(group)],
        pch = c(11, 19, 17)[as.factor(group_shape)] ,
        las = 1 ,
        bty = "l" ,
        cex = plot_cex ,
        cex.lab = 1.6 ,
        cex.axis = 1.4,
        xlab = labs[1],
        ylab = labs[2],
        mgp = c(3.5, .7, 0)
      )
      
      if (labeled) {
        text(
          svdResult$u[, 1:2],
          labels = clin$SampleID,
          cex = 1,
          pos = 3
        )
      }
      
      plot(
        svdResult$u[, c(1, 3)],
        col = cols[as.numeric(group)],
        pch = c(11, 19, 17)[as.factor(group_shape)] ,
        las = 1 ,
        bty = "l" ,
        cex = plot_cex ,
        cex.lab = 1.6 ,
        cex.axis = 1.4,
        xlab = labs[1],
        ylab = labs[3],
        mgp = c(3.5, .7, 0)
      )
      
      plot(
        svdResult$u[, c(2, 3)],
        col = cols[as.numeric(group)] ,
        pch = c(11, 19, 17)[as.factor(group_shape)] ,
        las = 1 ,
        bty = "l" ,
        cex = plot_cex ,
        cex.lab = 1.6 ,
        cex.axis = 1.4,
        xlab = labs[2],
        ylab = labs[3],
        mgp = c(3.5, .7, 0)
      )
      
    } else{
      plot(
        svdResult$u[, 1:2],
        col = cols[as.numeric(group)],
        pch = 19,
        las = 1 ,
        bty = "l" ,
        cex = plot_cex ,
        cex.lab = 1.6 ,
        cex.axis = 1.4,
        xlab = labs[1],
        ylab = labs[2],
        mgp = c(3.5, .7, 0)
      )
      
      if (labeled) {
        text(
          svdResult$u[, 1:2],
          labels = clin$SampleID,
          cex = 1,
          pos = 3
        )
      }
      
      plot(
        svdResult$u[, c(1, 3)],
        col = cols[as.numeric(group)],
        pch = 19 ,
        las = 1 ,
        bty = "l" ,
        cex = plot_cex ,
        cex.lab = 1.6 ,
        cex.axis = 1.4,
        xlab = labs[1],
        ylab = labs[3],
        mgp = c(3.5, .7, 0)
      )
      
      plot(
        svdResult$u[, c(2, 3)],
        col = cols[as.numeric(group)] ,
        pch = 19 ,
        las = 1 ,
        bty = "l" ,
        cex = plot_cex ,
        cex.lab = 1.6 ,
        cex.axis = 1.4,
        xlab = labs[2],
        ylab = labs[3],
        mgp = c(3.5, .7, 0)
      )
    }
    
    
    x <- 20
    y <- 20
    plot(
      x,
      y,
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n",
      main = group_to_test[g],
      cex.main = 1.5,
      col = "white",
      axes = F
    )
    legend(
      "topleft" ,
      legend = levels(group),
      cex = legend_cex,
      pch = 19 ,
      col = cols[seq(along = levels(group))],
      bty = "n", 
      ncol= 2
    )
    
    if (!is.null(group_to_shape)) {
      legend(
        "bottomleft" ,
        legend = levels(group),
        cex = legend_cex,
        pch = c(11, 19, 17)[unique(as.factor(group_shape))],
        bty = "n"
      )
    }
    
  }
  
  # dev.off()
}
