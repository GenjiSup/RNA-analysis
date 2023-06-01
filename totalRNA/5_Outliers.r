library(tidyverse)
library(magrittr)
library(cowplot)
library(DESeq2)
require("DESeq2")
require("edgeR")
library(knitr)
library(DT)
library(plotly)
library(PCAtools)
library(dplyr)

DESIGN<- "Compound"
#needs loading of metadata
#input is the dds
check_outliers_vst = function(input, DESIGN, OUTPUT.DIR, VAR.Threshold) {  
  OutNames<-NULL
  vsd = vst(input, blind = T) %>% assay()
  
  p = PCAtools::pca(vsd, 
                    metadata = DE_Design, 
                    removeVar = 0.1)
  
  scree_plot = PCAtools::screeplot(p, axisLabSize = 10, titleLabSize = 22)
  PCA = PCAtools::biplot(p, x = 'PC1', y = 'PC2', lab = rownames(p$metadata),
                         drawConnectors = F, colby = DESIGN, 
                         legendPosition = "bottom", legendLabSize = 13, legendIconSize = 6.0) + coord_fixed(ratio = 1)
  PCA
  ggsave(paste0(OUTPUT.DIR, '/', '_pca.png')) #saves the last plot
  
  cowplot::plot_grid(scree_plot, PCA, ncol = 2, nrow = 1) #helps with plot formatting
  ggsave(paste0(OUTPUT.DIR, '/', '_scree-pca.png'))
  
  # get replicates group
  #  levels_compare = p$metadata$Dose %>% droplevels() %>% unique() #original code Martha
  levels_compare = unique(c(condition1, condition2))
  
  # find PC with var >= VAR.Threshold
  pc_check = p$variance[which(p$variance >= VAR.Threshold)]
  
  # if there are no PC >= VAR.Threshold%, print message and stop
  if (length(pc_check) < 1) {
    print(paste0("No PC contributes to ", VAR.Threshold, "% or more of total variance"))
    
  } else {
    
    # check outliers in all PC that contribute for at least VAR.Threshold% of total variance
    for (x in seq(1, length(pc_check))) {
      
      print(pc_check[x])
      # max distance along PC
      max_pc = abs(max(p$rotated[x]) - min(p$rotated[x]))
      
      for (i in seq(1, length(levels_compare))) {
        # name of all replicates
        all_reps = p$metadata %>% dplyr::filter(.[, 1] %in% levels_compare[i]) %>% rownames()
        
        # PC values for all replicates
        all_reps_value = p$rotated %>% 
          .[rownames(.) %in% all_reps, x, drop = F] %>% 
          as.data.frame() %>% 
          .[order(.[, 1], decreasing = T), , drop = F]
        
        # distance on PC for all replicates. Pairwise and between closest
        all_reps_dist = all_reps_value %>% 
          dplyr::select(1) %>% 
          dist %>%  # get distance
          as.matrix() %>% 
          .[row(.) == col(.) + 1]  # get diagonal
        
        # check if any distance is >= VAR.Threshold% of `max_pc`
        # max_pc : var(PC) = all_reps_dist : x%
        # TRUE = outlier 
        print(paste0("max(Variation_between_samples) = ", max(
          (as.numeric(p$variance[x]) * all_reps_dist)/ max_pc)
        ))
        outliers = ( (as.numeric(p$variance[x]) * all_reps_dist) / max_pc ) > VAR.Threshold
        
        # Obtaining names of outliers (Marcha's addition)
        outlierNames1 = NULL
        HalfPoint = as.integer(length(outliers)/2)
        check_first_half = sum(outliers[1:HalfPoint])
        if (check_first_half > 0) {
          TRUEpos = NULL
          for (pos in 1:HalfPoint) { if (outliers[pos]==TRUE){TRUEpos<-c(TRUEpos, pos)} }
          outlierNames1 = row.names(all_reps_value)[1:max(TRUEpos)]  
        }
        
        outlierNames2 = NULL
        check_second_half = sum(outliers[(HalfPoint+1):length(outliers)])
        if (check_second_half>0) {
          TRUEpos = NULL
          for (pos in (HalfPoint+1):length(outliers)) {	if (outliers[pos]==TRUE){TRUEpos<-c(TRUEpos, pos)} }
          outlierNames2 = row.names(all_reps_value)[(min(TRUEpos)+1):(length(outliers)+1)]  
        }
        
        outlierNames = c(outlierNames1, outlierNames2)
        
        print(paste0(levels_compare[i], " contained ", length(outlierNames), " outliers: "))
        print(outlierNames)
        OutNames = c(OutNames, outlierNames)
      }
    }
  }
  OutNames<<-OutNames
  All_OutNames<<-c(All_OutNames, OutNames)
}






## Use the function

outputdir <- ("C:/Users/carlo/OneDrive/Documenti/University/Intership/RNA")

# Normalize data 
dds <- DESeqDataSetFromMatrix(countData = round(samples), colData = as.data.frame(DE_Design), 
                              design = as.formula(paste0("~", DESIGN)))
dds <- DESeq(dds, quiet=TRUE) #"Get some coffee, next step will take a while"
norm_data <<- counts(dds,normalized=TRUE)

vsd = vst(dds, blind = T) %>% assay()
p = PCAtools::pca(vsd, 
                  metadata = DE_Design, 
                  removeVar = 0.1)
PCA<-PCAtools::biplot(p, x = 'PC1', y = 'PC2', lab = rownames(p$metadata),
                      drawConnectors = F, colby = DESIGN, 
                      legendPosition = "bottom", legendLabSize = 13, legendIconSize = 6.0) + coord_fixed(ratio = 1)



# first outlier check
All_OutNames<-NULL

VAR.Threshold<- 20
print(paste0("Remove outliers with >", VAR.Threshold, "% variance between samples"))

check_outliers_vst(dds, DESIGN, outputdir, VAR.Threshold)

if (length(OutNames>0)) {
  for (name in OutNames) {
    Remove_Pcode<-paste0(unlist(strsplit(name, "_"))[1], "_")
    DE_Design <- DE_Design [c(grep(Remove_Pcode,rownames(DE_Design), invert=TRUE)),]
    samples <- samples[, rownames(DE_Design) ]
  }
  print("These outliers (and matched samples) have been removed from the dataset. Check for more outliers")
} else {print("No outliers found")}
