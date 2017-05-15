library(plyr)
library(metagenomeSeq)
library(stringr)

city = "boston"
meta.data <- read.table(paste("../metadata/", city, "_metadata.tsv", sep=""), header = TRUE,
                        stringsAsFactors = FALSE, sep = "\t")

filenames <- paste(meta.data$ID, "_", city, "_clark_OTU.csv", sep="")

tables <- c()
merged <- c()
for (f in filenames) {
  t <- read.table(paste("../otu/", city, "/", f , sep=""), fill = TRUE, 
             sep = ",", header = TRUE, stringsAsFactors = FALSE)
  colnames(t)[5] <- "All"
  colnames(t)[6] <- "Classified"
  t$All[nrow(t)] <- t$Count[nrow(t)]
  t$Count[nrow(t)] <- t$Lineage[nrow(t)]
  t$Lineage[nrow(t)] <- "UNKNOWN"
  #merged <- merge(merged, t[,1:3], all=TRUE, by.all="TaxID")
  merged <- merge(merged, t[,1:3], all=TRUE)
  t <- t[, c(-1, -3)]
  tables <- c(tables, list(t))
}

final <- merged
for (i in 1:length(tables)) {
  final <- merge(final, tables[[i]], all=TRUE, by="TaxID")
  colnames(final)[(i*3+1):(i*3+3)] = c(paste(c("Count", "All", "Classified"), ".", meta.data$ID[i], sep=""))
}

count <- final[,grepl("Count", names(final))]
colnames(count) <- meta.data$ID
count[is.na(count)] <- 0
rownames(count) <- final$TaxID

taxa <- str_split_fixed(final$Lineage, ";", 6)
taxa <- as.data.frame(taxa)
taxa <- cbind(final$Lineage, taxa)
taxa <- cbind(final$TaxID, taxa)
taxa <- cbind(taxa, final$Name)
taxa <- cbind(taxa, final$Name)
colnames(taxa) <- colnames(read.delim(file.path(dataDirectory, "CHK_otus.taxonomy.csv"),
                             stringsAsFactors = FALSE))

taxa$strain <- NA
taxa <- apply(taxa, 2, function(x) gsub("^$|^ $", NA, x))
taxa <- as.data.frame(taxa)

# Creating a MRexperiment object
meta <- meta.data
row.names(meta) <- meta.data$ID
meta <- meta[, -1]
phenotypeData = AnnotatedDataFrame(meta)
OTUdata = AnnotatedDataFrame(taxa)

count <- as.data.frame(sapply(count, as.numeric))

obj = newMRexperiment(count, phenoData=phenotypeData, featureData=OTUdata)
