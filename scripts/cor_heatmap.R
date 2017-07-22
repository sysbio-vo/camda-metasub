library(heatmap3)
city = "ny"
cor.data <- read.table(paste("../metadata/", city, "_cor.tsv", sep=""), header = TRUE,
                        stringsAsFactors = FALSE, sep = "\t")
meta.data <- read.table(paste("../metadata/", city, "_metadata.tsv", sep=""), header = TRUE,
                        stringsAsFactors = FALSE, sep = "\t")


#Create a custom color scale
library(RColorBrewer)
materialCol <- brewer.pal(length(levels(factor(meta.data$surface_material))),"Set2")
typeCol <- colorRampPalette(brewer.pal(12, "Paired"))(length(levels(factor(meta.data$surface_type)))) 

surface_material <- factor(meta.data$surface_material, labels=materialCol)
surface_type <- factor(meta.data$surface_type, labels=typeCol)
group <- cbind(SurfaceMaterial=as.vector(surface_material), SurfaceType=as.vector(surface_type))

heatmap3(cor.data, sym = T, legendfun=function()showLegend(x="center",
                            legend=levels(factor(meta.data$surface_material)),
                            col=materialCol, y.intersp=0.7),
         ColSideColors = group, showRowDendro=FALSE)

