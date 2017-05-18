library(metagenomeSeq)
library(plyr)
library(stringr)

city = "boston"
meta.data <- read.table(paste("../metadata/", city, "_metadata.tsv", sep=""), header = TRUE,
                        stringsAsFactors = FALSE, sep = "\t")

otu.data <- read.table(paste("../otu/", city, "_merged_otu.tsv", sep=""), header = TRUE,
                       stringsAsFactors = FALSE, sep = "\t")

rownames(meta.data) <- meta.data$ID

count <- otu.data[,grepl("Count", names(otu.data))]
colnames(count) <- meta.data$ID
count[is.na(count)] <- 0
rownames(count) <- otu.data$TaxID


if (city=="ny") {
  sums <- colSums(t(count))
} else {
  sums <- colSums(t(count[,1:(ncol(count)-3)]))
}

most.abundant <- otu.data[, 1:2]
most.abundant$Count <- sums
most.abundant <- most.abundant[order(most.abundant[,3], decreasing=TRUE), ]
#most.abundant$Count <- round(most.abundant$Count/nrow(meta.data), 4)
most.abundant$Count <- round(most.abundant$Count/sum(most.abundant$Count)*100, 4)

write.table(most.abundant, paste("../otu/", city, "_sum_otu.tsv", sep=""), sep="\t", quote = FALSE,
            row.names = FALSE)

lev = levels(factor(meta.data$surface_material))
dna <- data.frame(surface_material=c(1:length(l)), mean=c(1:length(l)), sd=c(1:length(l)))
for (i in 1:length(lev)) {
  dna[i,1] <- lev[i]
  dna[i,2] <- mean(meta.data[which(meta.data$surface_material==lev[i]),]$hs_dna)
  dna[i,3] <- sd(meta.data[which(meta.data$surface_material==lev[i]),]$hs_dna)
}

write.table(dna, paste("../metadata/", city, "_hdna_surface.tsv", sep=""), sep="\t",
            row.names = FALSE)

taxa <- str_split_fixed(otu.data$Lineage, ";", 6)
taxa <- as.data.frame(taxa)
taxa <- cbind(otu.data$Lineage, taxa)
taxa <- cbind(otu.data$TaxID, taxa)
taxa <- cbind(taxa, otu.data$Name)
taxa <- cbind(taxa, otu.data$Name)
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
legend("topleft", levels(cl), text.col = 1:length(levels(cl)),
       box.col = NA)
dev.off()

### Differential abundances
data <- obj
rareFeatures = which(rowSums(MRcounts(data) > 0) < 10)
data = data[-rareFeatures, ]
p = cumNormStat(data, pFlag = TRUE, main = "Trimmed data")
data = cumNorm(data, p = p)

pData(data)$surface_material[24] <- "unknown"
material = pData(data)$surface_material
mod = model.matrix(~as.factor(material))
colnames(mod) = levels(factor(material))
# fitting the ZIG model
settings = zigControl(maxit = 10, verbose = TRUE)
res = fitZig(obj = data, mod = mod, control = settings)
# The output of fitZig contains a list of various useful
# items. hint: names(res). Probably the most useful is the
# limma 'MLArrayLM' object called fit.
zigFit = res$fit
finalMod = res$fit$design
contrast.matrix = makeContrasts(PVC - polyester, levels = finalMod)
fit2 = contrasts.fit(zigFit, contrast.matrix)
fit2 = eBayes(fit2)
t <- topTable(fit2, sort.by="logFC")
t$species <- taxa[rownames(t),]$species

write.table(t, paste("../otu/", city, "_diffabund_pvc_polyester.tsv", sep=""), sep="\t", quote = FALSE,
            row.names = FALSE)
