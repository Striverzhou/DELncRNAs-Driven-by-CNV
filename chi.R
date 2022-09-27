#!/usr/bin/env Rscript

setwd(".") # set your workspace
rt=read.table("chiInput.txt",sep="\t",header=T,row.names=1)

outTab=data.frame()
for(i in 1:nrow(rt)){
  x=matrix(c(rt[i,1],rt[i,2],rt[i,3],rt[i,4]), ncol = 2)
  chiTest=chisq.test(x)
  normalRatio=rt[i,1]/(rt[i,1]+rt[i,2])
  tumorRatio=rt[i,3]/(rt[i,3]+rt[i,4])
  Gene=row.names(rt[i,])
  Stat=chiTest$statistic
  Pvalue=chiTest$p.value
  outTab=rbind(outTab,cbind(Gene,normalRatio,tumorRatio,Stat,Pvalue))
}
pvalue=as.numeric(as.vector(outTab[,"Pvalue"]))
adjP=p.adjust(pvalue,method ="bonferroni")
outTab=cbind(outTab,adjPvalue=adjP)
write.table(outTab,file="chiResult.txt",sep="\t",quote=F,row.names=F)
diffTab=outTab[outTab[,"adjPvalue"]<0.05,]
write.table(diffTab,file="diffCNV_lncRNAs.txt",sep="\t",quote=F,row.names=F)
