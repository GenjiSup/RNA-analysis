#VENN diagram making
library(ggplot2)
library(ggvenn)
setwd("C:/Users/carlo/OneDrive/Documenti/University/Intership/miRNA/Results/Therapeutic")
data <- read.csv("miRNA.csv")
y <- list()
AZA <- data[, "AZA"]
CYC <- data[, "CYC"]
AZAvsCYC <- data[, "AZAvsCYC"]
y <- list(AZA = AZA, CYC = CYC, AZAvsCYC = AZAvsCYC)
ggvenn(y,
fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
stroke_size = 0.5, set_name_size = 4
)
