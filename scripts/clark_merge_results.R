library(plyr)
library(metagenomeSeq)
library(stringr)

city = "ny"
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
colnames(taxa) <- colnames(read.delim(file.path(system.file("extdata", package = "metagenomeSeq"), "CHK_otus.taxonomy.csv"),
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

# Normalization
p = cumNormStatFast(obj)
obj = cumNorm(obj, p = p)

# Visualization

## Heatmap
trials = pData(obj)$surface_material
heatmapColColors = brewer.pal(12, "Set3")[as.integer(factor(trials))]
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)

pdf(paste("../plots/heatmap_", city, ".pdf", sep=""),width=9,height=8)
par(cex.main=0.8)
pl <- plotMRheatmap(obj = obj, n = 200, fun=mad, cexRow = 0.4, cexCol = 0.4, 
              trace = "none", col = heatmapCols, ColSideColors = heatmapColColors)

dev.off()

## CMDS

cl = factor(pData(obj)$surface_material)

pdf(paste("../plots/CMDS_", city, ".pdf", sep=""),width=9,height=8)
# plotOrd - can load vegan and set distfun = vegdist and use
# dist.method='bray'
plotOrd(obj, tran = TRUE, bg = cl, pch = 21)
# plotRare
res = plotRare(obj, cl = cl, pch = 21, bg = cl)
# Linear fits for plotRare / legend
tmp = lapply(levels(cl), function(lv) lm(res[, "ident"] ~ res[,"libSize"] - 1, subset = cl == lv))
for (i in 1:length(levels(cl))) {
  abline(tmp[[i]], col = i)
}
legend("topleft", c(levels(cl), "NA"), text.col = c(1:length(levels(cl))+1),
       box.col = NA)
dev.off()
