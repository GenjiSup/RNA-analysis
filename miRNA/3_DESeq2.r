---
title: "miRNA DESeq2"
author: "Carlo Alberto Zani"
date: "2023-04-26"
output: html_document
---

```{r, setup, include = FALSE}
################################
# Load R libraries					  #
################################ 
library(DESeq2)
require("DESeq2")
require("edgeR")
library(knitr)
library(DT)
library(plotly)
```

```{r echo=T, include=TRUE}
# Specify which groups need to be compared 
Samp4compare<- "Azathioprine" 
Cont4compare<- "Cyclosporine"
DESIGN<- "Compound"

# Set thresholds for Differential Expression (based on R-ODAF)
minCoverage <- 5000000
MinCount<- 1
pAdjValue<- 0.01 
```


``` {r echo=F, include=FALSE}
# Set file locations				
#analysisID <- paste(AO, GEO_GSE, sep="_")
#sampledir <- paste("C:/A_DESeq2/", GEO_GSE, "_Data/", sep="") 
#if(!dir.exists(sampledir)) {dir.create(sampledir)}
#outputdir <- paste(sampledir, GEO_GSE, "_Output/", sep="")
#if(!dir.exists(outputdir)) {dir.create(outputdir)}


#This comma delimited file contains at least 2 columns: NAME (sample names identical to the column names of sampleData) and Group (needs to identify to which group the sample belongs -> Disease/Control, ..., ExperimentalGroup & ControlGroup)

# Load input files 
setwd("C:/Users/carlo/OneDrive/Documenti/University/Intership/miRNA")

# Load gene table
miRNA_counts <- read.table("miRNA_counts.txt", header=TRUE, sep="\t")
# Find columns that contain the word "Tox"
tox_cols <- grep("Tox", colnames(miRNA_counts))
# Remove columns with the word "Tox"
miRNA_counts_filtered <- miRNA_counts[,-tox_cols]
# Save filtered gene table
write.table(miRNA_counts_filtered, "miRNA_counts_filtered.txt", sep="\t", row.names=FALSE)
# Load filtered gene table
sampleData <- read.delim("miRNA_counts_filtered.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE, quote="\"", row.names=1)

# Load metadata table
metadata <- read.table("metadata_miRNA.txt", header=TRUE, sep=" ")
# Filter samples with "toxic" status
metadata_filtered <- metadata[metadata$Dose != "Tox", ]
# Save filtered metadata table
write.table(metadata_filtered, "metadata_filtered.txt", sep=" ", row.names=FALSE)
DESeqDesign <- read.delim("metadata_filtered.txt", sep=" ", stringsAsFactors=FALSE, header=TRUE,  quote="\"", row.names="Name")

```

### **Overview of samples and groups**
```{r echo=F, results='asis'}
datatable(DESeqDesign)
```



### **Data processing:**

##### *Data clean-up: replace NA with 0*
``` {r echo=F, include=T}
ZeroDetected<-count(sampleData[ is.na(sampleData) ])
sampleData[ is.na(sampleData) ] <- 0 
print(paste0(ZeroDetected, "x replacement of NA values"))
```
```{r echo=FALSE, include=TRUE}
sampleData <- sampleData[,abs(log2(colSums(sampleData)/mean(colSums(sampleData)))) < 2]
DESeqDesign <- DESeqDesign[rownames(DESeqDesign) %in% colnames(sampleData),]
```

##### *Remove samples with total readcount < threshold (1M)*
``` #{r echo=F, include=T}
Keep<-ncol(sampleData[,(colSums(sampleData)> 1000000)])
Remove<-ncol(sampleData[,(colSums(sampleData)< 1000000)])
print(paste0("From the total of ", ncol(sampleData), " samples, ", Remove, " samples had to be removed due to low sequencing depth (<1M) -> ", Keep, " samples remaining"))
sampleData<- sampleData[,(colSums(sampleData)> 1000000)]
DESeqDesign <- DESeqDesign[rownames(DESeqDesign) %in% colnames(sampleData),]
```


``` {r echo=F, include=F}
# DEFINE FUNCTION for creating PCA
############################################
Make_PCA<-function(DATA,TITLE,FILENAME,LEGEND_LOCATION){
#DATA = norm_data, TITLE= title or PCA, FILENAME=x(without .png), LABLES = T/F, LEGEND_LOCATION= topright, topleft, ...
#X11(type='cairo')
	pc <- prcomp(t(DATA)) 
		PCvariance <- pc$sdev^2/sum(pc$sdev^2)
		PC1var<-signif((PCvariance[1]*100), 2)
		PC2var<-signif((PCvariance[2]*100), 2)
	Colors<-matrix(data=NA, ncol=ncol(DATA), nrow=1)
		colnames(Colors)<-colnames(DATA)
		pch<-NULL
	for (k in 1:ncol(Colors)) {
	  if (DESeqDesign[colnames(Colors)[k], DESIGN] == Cont4compare) {Colors[1,k] <- "dodgerblue" 
			pch<- c(pch, 21)
		} else if (DESeqDesign[colnames(Colors)[k], DESIGN] == Samp4compare) {Colors[1,k] <- "red" 
			pch<- c(pch, 21)
		} else {Colors[1,k] <- "grey" 
		  pch<- c(pch, 21)
		}
	}
	c<-t(Colors)
	
	#PCA plot of the data with colors
	plot( pc$x[ , 1:2 ], col=1, bg=c , cex=1, pch=pch, main=paste("PCA", TITLE, sep=" "), xlab=paste0("PC1: ", PC1var, "% variance"), ylab=paste0("PC2: ", PC2var, "% variance"))
	legend(LEGEND_LOCATION, legend=c(Cont4compare, Samp4compare, "other group(s)"), col=c("dodgerblue", "red", "grey"), pch=16)

	#to add labels to the plot
#text(pc$x[ , 1], pc$x[ ,2 ], rownames(c),pos=1, cex = 0.6)
#savePlot(filename = paste0(FILENAME, "_WithLables.png"), type = "png")
} # Make_PCA function defined

Make_PCA_SampleID<-function(DATA, TITLE, FILENAME, LEGEND_LOCATION){
#DATA = norm_data, TITLE= title or PCA, FILENAME=x(without .png), LABLES = T/F, LEGEND_LOCATION= topright, topleft, ...
#X11(type='cairo')
	pc <- prcomp(t(DATA)) 
		PCvariance <- pc$sdev^2/sum(pc$sdev^2)
		PC1var<-signif((PCvariance[1]*100), 2)
		PC2var<-signif((PCvariance[2]*100), 2)
	Colors<-matrix(data=NA, ncol=ncol(DATA), nrow=1)
		colnames(Colors)<-colnames(DATA)
		pch<-NULL
	for (k in 1:ncol(Colors)) {
	  if (DESeqDesign[colnames(Colors)[k], DESIGN] == Cont4compare) {Colors[1,k] <- "dodgerblue" 
			pch<- c(pch, 21)
		} else if (DESeqDesign[colnames(Colors)[k], DESIGN] == Samp4compare) {Colors[1,k] <- "red" 
			pch<- c(pch, 21)
		} else {Colors[1,k] <- "grey" 
		  pch<- c(pch, 21)
		}
	}
	c<-t(Colors)
	
	#PCA plot of the data with colors
	plot( pc$x[ , 1:2 ], col=1, bg=c , cex=1, pch=pch, main=paste("PCA", TITLE, sep=" "), xlab=paste0("PC1: ", PC1var, "% variance"), ylab=paste0("PC2: ", PC2var, "% variance"))
	legend(LEGEND_LOCATION, legend=c(Cont4compare, Samp4compare, "other group(s)"), col=c("dodgerblue", "red", "grey"), pch=16)
	#to add labels to the plot
text(pc$x[ , 1], pc$x[ ,2 ], rownames(c),pos=1, cex = 0.6)
} # Make_PCA function defined
```

##### **PCA plot to exclude the post-processing outliers replicates**
```{r echo=F, results=F,out.width="95%",fig.width=6,fig.height=6, fig.align='center',cache=F}
#Make_PCA(sampleData, paste0( " RawData ", Samp4compare, " vs ", Cont4compare), paste0("_RawData_", Samp4compare, "_vs_", Cont4compare), "topright")
```


```{r echo=F, results=F,out.width="95%",fig.width=6,fig.height=6, fig.align='center',cache=F}
#Make_PCA_SampleID(sampleData, paste0(analysisID, " RawData ", Samp4compare, " vs ", Cont4compare), paste0(outputdir, analysisID, "_RawData_", Samp4compare, "_vs_", Cont4compare), "topright")
```

``` {r echo=F, include=T}
#DESeqDesign <- DESeqDesign [c(grep("SRR9883048",rownames(DESeqDesign), invert=TRUE)),]
#DESeqDesign <- DESeqDesign [c(grep("SRR9883049",rownames(DESeqDesign), invert=TRUE)),]
#  sampleData <- sampleData[, rownames(DESeqDesign) ]
```  

##### **PCA plot after exclusion**
```{r echo=F, fig.align='center', fig.height=6, fig.width=6, cache=FALSE, out.width="95%", results=F}
Make_PCA(sampleData, paste0(" RawData ", Samp4compare, " vs ", Cont4compare), paste0("_RawData_", Samp4compare, "_vs_", Cont4compare), "bottomright")
```

```{r echo=F, results=F,out.width="95%",fig.width=6,fig.height=6, fig.align='center',cache=F}
#Make_PCA_SampleID(sampleData, paste0(analysisID, " RawData ", Samp4compare, " vs ", Cont4compare), paste0(outputdir, analysisID, "_RawData_", Samp4compare, "_vs_", Cont4compare), "right")
```

############################################
## Differential expression analysis: DESeq2 
############################################

``` {r echo=F, include=T}
x=1 #Rmarkdown should only run 1 comparison at a time (to create nice output)
#for (x in 1:length(Samp4compare)){	## for all comparisons to be done	

  condition1<- Cont4compare[x]	    		
  condition2<- Samp4compare[x]  
  
  DE_Design <- matrix(data=NA, ncol=2)
  DE_Design <- DESeqDesign [c(grep(condition1,DESeqDesign[,DESIGN[x]]), grep(condition2,DESeqDesign[,DESIGN[x]])),]
  samples <- sampleData[, rownames(DE_Design) ]
  
  ###########
  print(paste("Comparison: ", condition2, " vs ", condition1, ":"))		
```


``` {r echo=F, include=F}
  samples <- samples[rowSums(samples)!=0,]
  colnames(samples)<-NULL
  dds <- DESeqDataSetFromMatrix(countData = round(samples), colData = as.data.frame(DE_Design), design = as.formula(paste0("~", DESIGN[x])))
```


##### *Executing dds step*
``` {r echo=F, include=T} 
  dds <- DESeq(dds, quiet=TRUE)
 
# Filtering genes with low readcounts: 
# 75% of at least 1 group need to be above MinCount CPM

  SampPerGroup<-table(DE_Design[,DESIGN])
  kable(SampPerGroup, caption = "Amount of samples per group")

  Counts<-counts(dds, normalized=TRUE)
  CPMdds<-cpm(counts(dds, normalized=TRUE))
  
  Filter <- matrix(data=NA, ncol=3, nrow= nrow(Counts))
  rownames(Filter) <- rownames(Counts)
  colnames(Filter) <- c("Low readcounts","quantile","spike")
  
  for (gene in 1:nrow(dds)) {
    
    CountsPass<-NULL
    for (group in 1:length(SampPerGroup)) { 
      sampleCols<-grep(dimnames(SampPerGroup)[[1]][group],DE_Design[,DESIGN])
      Check<-sum(CPMdds[gene,sampleCols] >= MinCount)>= 0.75*SampPerGroup[group]
      CountsPass<-c(CountsPass, Check)
    }
    
    if ( sum(CountsPass) > 0 ) {Filter[gene,1] <- 1 }	else { Filter[gene,1] <- 0 }
    
  }	

    compte <- Counts[Filter[,1] == 1,]
  Filter <- Filter[rownames(Filter) %in% rownames(compte),]
print("Executing dds step ... Done")
```
##### *Filtering genes with low readcounts: 75% of at least 1 group need to be above "MinCount" CPM*
``` {r echo=F, include=T} 
  print(paste("low readcount filtering removed ",nrow(dds)- nrow(Filter)," genes from the ",nrow(dds)," assessed. ", nrow(Filter)," genes remaining",sep=""))
```
##### *Obtaining the DESeq2 results*
``` {r echo=F, include=T}  
  # compute the DEGs on the genes passing the Relevance condition
  
  res <- results(dds[rownames(compte),], contrast=c(DESIGN[x], condition2, condition1), pAdjustMethod= 'fdr')
  
  setwd("C:/Users/carlo/OneDrive/Documenti/University/Intership/miRNA/AZAvsCYC")
  FileName<-paste(condition2,"vs",condition1, "FDR", pAdjValue, sep="_")
  
  #Save output tables		
  norm_data <<- counts(dds[rownames(compte)],normalized=TRUE) 
  DEsamples <<- subset(res,res$padj < pAdjValue)	
  
  print(paste("A total of ", nrow(DEsamples), " DEGs were selected (before filtering)"))
  
  DECounts <- compte[rownames(compte) %in% rownames(DEsamples),]
  Filter <- Filter[rownames(Filter) %in% rownames(DECounts),]
```

##### *Filtering for relevance: Check median against third quantile*
#####     AND
##### *Filtering: the presence of a spurious spike*
``` {r echo=F, include=F}   
  for (gene in 1:nrow(DECounts)) {
    
    # Check the median against third quantile
    quantilePass <-NULL
    sampleColsg1 <- grep(dimnames(SampPerGroup)[[1]][1],DE_Design[,DESIGN])
    sampleColsg2 <- grep(dimnames(SampPerGroup)[[1]][2],DE_Design[,DESIGN])
    
    Check <- median(DECounts[gene,sampleColsg1]) > quantile(DECounts[gene,sampleColsg2], 0.75)[[1]]
    quantilePass <-c(quantilePass, Check)
    Check <- median(DECounts[gene,sampleColsg2]) > quantile(DECounts[gene,sampleColsg1], 0.75)[[1]]
    quantilePass <-c(quantilePass, Check)
    
    if ( sum(quantilePass) > 0 ) {Filter[gene,2] <- 1 }	else { Filter[gene,2] <- 0 }
    
    # Check for spurios spike (not biologically relevant)
    spikePass <- NULL
    for (group in 1:length(SampPerGroup)) { 
      sampleCols<-grep(dimnames(SampPerGroup)[[1]][group],DE_Design[,DESIGN])
      if (max(DECounts[gene,sampleCols]) ==0) {Check <- FALSE} else {
        Check <- (max(DECounts[gene,sampleCols])/sum(DECounts[gene,sampleCols])) >= 1.4*(SampPerGroup[group])^(-0.66)
        spikePass<-c(spikePass, Check)
      }
    }
    if ( sum(spikePass) > 1 ) {Filter[gene,3] <- 0 }	else { Filter[gene,3] <- 1 }
    
  }		
  
  # extract the final list of DEGs
  
  DECounts_real <- DEsamples[rowSums(Filter) == 3 ,]
  DECounts_no_quant <- DEsamples[Filter[,2] == 0 ,]
  DECounts_spike <- DEsamples[Filter[,3] == 0 ,]
```

``` {r echo=F, include=T}   
  print(paste("A total of ",nrow(DECounts_real), " DEGs were selected, after ",nrow(DECounts_no_quant)," genes(s) removed by the quantile rule and ", nrow(DECounts_spike)," gene(s) with a spike",sep=""))
```  

### **Save the normalized counts and the list of DEGs**
``` {r echo=F, include=T} 
write.table(norm_data,file=paste0(FileName, "_Norm_Data.csv"), sep="\t", quote=FALSE)
  print(paste0("Normalized count data has been stored in: ", FileName, "_Norm_Data.csv"))
  write.table(DECounts_real,file=paste0(FileName,"_DEG_table.csv"), sep="\t", quote=FALSE)
  print(paste0("Differentially expressed genes (DEGs) have been stored in: ", FileName, "_DEG_table.csv"))
  DEGs_display<-read.table(paste0(FileName, "_DEG_table.csv"))
  #}
```  

DESeq2 identification of DEGs finished. 
  
# **DESeq2 Results**

##### **PCA plot of normalized data**
```{r echo=F, results=F,out.width="95%",fig.width=6,fig.height=6, fig.align='center',cache=F}
Make_PCA(norm_data, paste0(" NormData ", Samp4compare, " vs ", Cont4compare), paste0("_NormData_", Samp4compare, "_vs_", Cont4compare), "bottomleft")
```

```{r echo=F, results=F,out.width="95%",fig.width=12,fig.height=6, fig.align='center',cache=F}
#Make_PCA_SampleID(norm_data, paste0(analysisID, " NormData ", Samp4compare, " vs ", Cont4compare), paste0(outputdir, analysisID, "_NormData_", Samp4compare, "_vs_", Cont4compare), "topleft")
```

**Heatmap of all genes**
```{r echo=F, out.width="95%", fig.width=8,fig.height=10 ,fig.align='center',cache=F}
heatmap(norm_data)
```

**Heatmap of DEGs**
```{r echo=F, out.width="95%", fig.width=8,fig.height=10 ,fig.align='center',cache=F}
DEG_counts<-norm_data[c(row.names(DEGs_display)),]
heatmap(DEG_counts)
```

## **Table of identified DEGs (FDR<0.01)**
```{r echo=F, results='asis'}
datatable(DEGs_display)
```

