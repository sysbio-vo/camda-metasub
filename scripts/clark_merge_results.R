library(plyr)

meta.data <- read.table("../otu_metadata/Boston_metadata.tsv", header = TRUE, stringsAsFactors = FALSE,
                        sep = "\t")

filenames <- paste(meta.data$ID, "_clark_OTU.csv", sep="")

tables <- c()
merged <- c()
for (f in filenames) {
  t <- read.table(paste("../boston_clark/", f , sep=""), fill = TRUE, 
             sep = ",", header = TRUE, stringsAsFactors = FALSE)
  colnames(t)[5] <- "Proportion_All"
  colnames(t)[6] <- "Proportion_Classified"
  t$Proportion_All[nrow(t)] <- t$Count[nrow(t)]
  t$Count[nrow(t)] <- t$Lineage[nrow(t)]
  t$Lineage[nrow(t)] <- "UNKNOWN"
  
  tables <- c(tables, list(t))
  merged <- merge(merged, t[,1:3], all=TRUE)
}



merged <- merge(tables[[1]][,1:3], tables[[2]][,1:3], all=TRUE)
