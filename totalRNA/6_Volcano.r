# VOLCANO PLOT OF GENES
library(tidyverse)
library(magrittr)
library(cowplot)
#use the res text files from the RNA_DESeq2 script
setwd("C:/Users/carlo/OneDrive/Documenti/University/Intership/RNA/Fixed")
res_AZAvsCYC <- read.delim("res_AZAvsCYC.txt", header = TRUE, sep = "\t")

# Create a new column in 'res' indicating significant DEGs
   res_AZAvsCYC$isDEG <- ifelse(res_AZAvsCYC$padj < 0.01, "DEG", "Not DEG")
  
    # Define the desired range for the y-axis
     yAxisRange <- c(0, 80)  # Set the desired range here
    # Define the desired range for the x-axis
     xAxisRange <- c(-15, 15)  # Set the desired range here
     
       # Create the volcano plot with custom y-axis range
       volcano_plot <- ggplot(res_AZAvsCYC, aes(x = log2FoldChange, y = -log10(padj), color = isDEG)) +
           geom_point(size = 0.8) +
           scale_color_manual(values = c("DEG" = "red", "Not DEG" = "black")) +
           labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", title = "", color = "") +
           ylim(yAxisRange)  # Set the y-axis range
           xlim(xAxisRange)

       
         # Display the plot
         print(volcano_plot)
