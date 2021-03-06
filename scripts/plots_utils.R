require(ggplot2)
require(cowplot)
require(ggfortify)
library(grid)
library(reshape2)
library(RColorBrewer)

getAspectRatio <- function(p){
  gb <- ggplot_build(p)
  g <- ggplot_gtable(gb)
  
  nullw <- sapply(g$widths, attr, "unit")
  nullh <- sapply(g$heights, attr, "unit")
  
  # ar of plot
  if(any(nullw == "null"))
    ar <- unlist(g$widths[nullw == "null"]) / unlist(g$heights[nullh == "null"])
  
  # ar of plot + legend
  g$fullwidth <- convertWidth(sum(g$widths), "in", valueOnly=TRUE)
  g$fullheight <- convertHeight(sum(g$heights), "in", valueOnly=TRUE)
  ar <- g$fullwidth / g$fullheight
  
  return(ar)
}

pcaPlots <- function(pca.data, pheno.data, meta.vars, title, ncol) {
  plots <- c()
  ar <- -100
  for (i in meta.vars) {
    pl <- autoplot(pca.data, data = pheno.data, colour=i)
      #coord_fixed()
      #geom_text(aes(label=ID), size=2)
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
  pl <- pl + theme(plot.margin=margin(t=1, r=1, b=1, l=1))
  return(list(pl, ar))
}