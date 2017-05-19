library(plyr)
source("plots_utils.R")

city = "ny"
meta.data <- read.table(paste("../metadata/", city, "_metadata.tsv", sep=""), header = TRUE,
                        stringsAsFactors = FALSE, sep = "\t")

otu.data <- read.table(paste("../otu/", city, "_merged_otu.tsv", sep=""), header = TRUE,
                       stringsAsFactors = FALSE, sep = "\t")

rownames(meta.data) <- meta.data$ID

count <- otu.data[,grepl("Count", names(otu.data))]
colnames(count) <- meta.data$ID
count[is.na(count)] <- 0
rownames(count) <- otu.data$TaxID

if (city=="boston") {
#  features = c("surface_material", "hs_dna", "GC", "surface_type",
#               "avg_seq_length", "unknown_species", "dups", "station",
#               "collection_date")
  features = c("surface_material", "surface_type")
  pca = prcomp(t(log(count+0.00001)))
  nc = 2; nr = 1
}
if (city=="ny") {
  features = c("surface_material", "borough", "GC", "surface_type",
               "avg_seq_length", "unknown_species", "Mseqs", "avg_dew_point",
               "sampling_place", "avg_air_temp_F", "ground_level", "avg_abs_humidity")
  pca = prcomp(t(log(count+0.00002)))
  nc = 3; nr = 4
}
if (city=="sacramento") {
  features = c("surface_material", "hs_dna", "GC", "surface_type",
               "avg_seq_length", "unknown_species", "dups", "station",
               "avg_abs_humidity")
  pca = prcomp(t(count))
  nc = 3; nr = 3
}

pl <- pcaPlots(pca, meta.data,
               features,
               ncol = nc)

save_plot(paste("../plots/PCA/", city, "_PCA.pdf", sep=""), ncol=nc, nrow = nr, 
          base_height=2.8, base_width= 6, pl[[1]])
