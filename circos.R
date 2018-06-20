# 06/19/2018 Zhi Huang

library(circlize)
circlizeGenomics <- function(BED.data, factors, xlim, mySpecies, myTitle, circos_param_genelink, circos_param_genesymbol, font.scale, link.width, color.picker){
  # save(BED.data, file ="~/Desktop/BEDdata.Rdata")
  par(mar = c(1, 1, 1, 1))
  # reference: http://zuguang.de/circlize_book/book/initialize-genomic-plot.html#initialize-cytoband
  circos.clear()
  circos.par(cell.padding = c(0, 0, 0, 0))
  # circos.initializeWithIdeogram(plotType = NULL)
  if (mySpecies=="hg19") {
    circos.initializeWithIdeogram(species = mySpecies)
  }
  else if (mySpecies=="hg38"){
    circos.initializeWithIdeogram(species = mySpecies, chromosome.index = paste0("chr", c(1:22, "X", "Y")))
  }
  # text(0, 0, "Human Chromosomes", cex = 1)
  # Abbreviations of species. e.g. hg19 for human, mm10 for mouse.
  # If this value is specified, the function will download cytoBand.txt.gz
  # from UCSC website automatically. If there is no cytoband for user's species,
  # it will keep on trying to download chromInfo file.
  
  # we assume data is simply a data frame in BED format
  # (where the first column is the chromosome name, the
  # second and third column are start and end positions,
  # and the following columns are associated values)
  
  # circos.genomicRainfall(data.frame(hg38.ring[,c(4,6:7)]), pch = 16, cex = 0.4, col = "#FF000080")
  title(myTitle)
  circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA,
                         bg.col = rep("grey", 24), track.height = 0.05,
                         panel.fun = function(x, y) {
                           sector.name = get.cell.meta.data("sector.index")
                           xlim = get.cell.meta.data("xlim")
                           # circos.text(mean(xlim), 2.5,
                           #             facing = "outside", niceFacing = T,
                           #             sector.name, cex = 0.8, adj = c(0.5, 0))
                         })
  circos.genomicTrack(BED.data, track.height = 0.01, bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
                        i = getI(...)
                        # numeric.column is automatically passed to `circos.genomicPoints()`
                        circos.genomicPoints(region, value = 1, ...)
                      })
  if(circos_param_genelink){
    # circos.initializeWithIdeogram(plotType = NULL)
    for (i in 1:length(BED.data$chrom)){
      for (j in 1:length(BED.data$chrom)){
        if(j<i){
          circos.link(sector.index1=BED.data$chrom[i], point1=c(BED.data$txStart[i],BED.data$txEnd[i]),
                      sector.index2=BED.data$chrom[j], point2=c(BED.data$txStart[j],BED.data$txEnd[j]),
                      col = color.picker, lwd = link.width)
          # R color: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
        }
      }
    }
  }
  if(circos_param_genesymbol){
    circos.genomicLabels(BED.data, labels.column = 5, side = "inside",
                         cex = 0.8*font.scale,
                         col = "black",
                         line_col = "black"
                         # col = as.numeric(factor(BED.data[[1]])),
                         #line_col = as.numeric(factor(BED.data[[1]]))
                         )
  }
  
  
  circos.clear()
}