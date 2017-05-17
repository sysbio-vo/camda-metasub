require(ggplot2)
require(cowplot)
require(ggfortify)
library(grid)
library(reshape2)
library(RColorBrewer)

pcaPlots <- function(pca.data, pheno.data, meta.vars, title, ncol) {
  pheno.data[] <- lapply(pheno.data, as.character)
  plots <- c()
  ar <- -100
  for (i in meta.vars) {
    pl <- autoplot(pca.data, data = pheno.data, colour=i) +
      coord_fixed()
    newar <- getAspectRatio(pl)
    if (ar<newar) {
      ar <- newar
    }
    plots <- c(plots, list(pl))
  }
  if(missing(ncol)) {
    pl <- plot_grid(plotlist = plots, ncol=length(meta.vars), align="hv")
  } else {
    pl <- plot_grid(plotlist = plots, ncol=ncol, align="hv")
  }
  if (!missing(title)) {
    title <- ggdraw() + draw_label(title, fontface='bold')  
    pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  }
  pl <- pl + theme(plot.margin=margin(t=10, r=10, b=10, l=10))
  return(list(pl, ar))
}