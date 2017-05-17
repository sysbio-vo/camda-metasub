library(plyr)
library(stringr)

city = "sacramento"
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

write.table(final, paste("../otu/", city, "_merged_otu.tsv", sep=""), sep="\t",
            row.names = FALSE)
