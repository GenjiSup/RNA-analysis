# R script rsem.merge
predir <- "/share/analysis/Carlo/"
WORKDIR <- paste(predir,"3_rsem.star", sep="")
setwd(WORKDIR)

# Defining pattern
patt <- "*gene"
patt_isoform <- "*isoforms"
#GENES
# Listing all the files with the pattern
fichier <- list.files(pattern = patt)
# Show how many files there are ([1] 5, if there are 5 files)
taille <- length(fichier)
# Prepare merge_file_gene
# [1] = take the first value from the file
 temp1<- read.table(file = fichier[1], header = TRUE, sep = "\t")
 temp1 <- subset(temp1[,c(1,5)])
 # Extracting characters from 1-14 (letters/numbers) from the file name using substr
 colnames(temp1)<- c("Row.names",substr(fichier[1],8,20))
# [2] = take the second value from the file  
 temp2<- read.table(file =fichier[2], header = TRUE,  sep = "\t")
 temp2 <- subset(temp2[,c(1,5)])
 colnames(temp2)<- c("Row.names",substr(fichier[2],8,20))
  
 merge_data<- merge(temp1,temp2,by="Row.names", all=TRUE)

# Merge all gene.results files	
 for (m in 3:taille) {
            #m=3
	    temp1A <- read.table(fichier[m], header = TRUE, sep = "\t")
	    temp1A <- subset(temp1A[,c(1,5)])
	    print(paste("Processed files: ",fichier[m],sep=""))
	    colnames(temp1A)<- c("Row.names",substr(fichier[m],8,20))
	    
	    merge_data <- merge(merge_data,temp1A,by="Row.names",all=TRUE)
	}
	print("Mergefile.genes.results done")
  
  setwd(predir)

	write.table(merge_data, file =paste("gene", "_counts.txt", sep=""), row.names = FALSE, quote = FALSE, sep="\t") 
 
# ISOFORMS
setwd(WORK.DIR)

fichier2 = list.files(pattern = patt_isoform)
taille = length(fichier2)

# Prepare merge_data
temp1 = read.table(file = fichier2[1], header = TRUE, sep = "\t")
temp1 = subset(temp1[, c(1,5)])
colnames(temp1) = c("Row.names", substr(fichier2[1],8,20))

temp2 = read.table(file = fichier2[2], header = TRUE,  sep = "\t")
temp2 = subset(temp2[,c(1,5)])
colnames(temp2) = c("Row.names", substr(fichier2[2],8,20))
  
merge_iso = merge(temp1, temp2, by = "Row.names", all = TRUE)

# Merge all gene.results files	
for (m in 3:taille) {
	#m=3
	temp1A = read.table(fichier2[m], header = TRUE, sep = "\t")
	temp1A = subset(temp1A[,c(1,5)])	
	print(paste("Processed files: ", fichier2[m], sep=""))
	colnames(temp1A) = c("Row.names", substr(fichier2[m],8,20))
	    
	merge_iso = merge(merge_iso, temp1A, by = "Row.names", all = TRUE)
}

print("Mergefile isoforms.results done")

setwd(predir)

write.table(merge_iso, file = "isoform_counts_batch1.txt", row.names = FALSE, quote = FALSE, sep="\t") 
