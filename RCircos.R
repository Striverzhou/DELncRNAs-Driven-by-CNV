#install.packages("RCircos")

setwd(".") # set your workspace
library("RCircos")

#Initialize RCircos core components
cytoBandIdeogram=read.table("refer.txt",sep="\t",header=T)
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

#Modifying RCircos core components
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size <- 0.5
#rcircos.params$heatmap.color <- "GreenWhiteRed"
RCircos.Reset.Plot.Parameters(rcircos.params)

pdf(file="lncRNA_RCircosOut.pdf", height=8, width=8)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#Scatter
RCircos.Scatter.Data=read.table("lncRNA_scatter.txt",sep="\t",header=T)
#RCircos.Scatter.Data[,1]=paste("chr",RCircos.Scatter.Data[,1],sep="")
data.col <- 4
track.num <- 1
side <- "in"
by.fold <- 1
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col,track.num, side, by.fold)

#Gene Labels
RCircos.Gene.Label.Data=read.table("lncRNA_Label.txt",sep="\t",header=T)
name.col <- 4
side <- "in"
track.num <- 2
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data,track.num, side)
track.num <- 3
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data,name.col,track.num, side)
dev.off()
