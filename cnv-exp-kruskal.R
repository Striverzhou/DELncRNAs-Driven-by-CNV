setwd("your workspace") # set your workspace

geneFile="candidate_CNVdriven_DElncRNAs.txt"
cnvFile="cnvMatrix.txt"
expFile="difflncRNAExp.txt"
cnv = read.table(cnvFile, row.names=1 ,header=T,sep="\t",check.names=F)
RNA = read.table(expFile, row.names=1 ,header=T,sep="\t",check.names=F)
gene <- read.table(geneFile,header=T,row.names=1,sep="\t",check.names=F)

colnames(cnv)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(cnv))
colnames(RNA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(RNA))
sameSample=intersect(colnames(cnv),colnames(RNA))
sameGene=rownames(gene)

cnv=cnv[sameGene,sameSample]
RNA=RNA[sameGene,sameSample]
rownames(cnv)=paste(rownames(cnv),"cnv",sep="|")
rownames(RNA)=paste(rownames(RNA),"exp",sep="|")

rt=rbind(id=sameSample,cnv,RNA)
#write.table(rt,file="merge.txt",sep="\t",quote=F,col.names=F)

down_genes <- gene[gene$Regulation=="down",]
up_genes <- gene[gene$Regulation=="up",]

group=sapply(strsplit(colnames(rt),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
rt=rt[,group==0]
outTab=data.frame()

for (id in rownames(down_genes)) {
  geneName=id
  cnv=paste(geneName,"|cnv",sep="")
  exp=paste(geneName,"|exp",sep="")
  ct=rt[,rt[cnv,]<=0]
  ct[cnv,]=gsub("-2","-1",ct[cnv,])
  normal <- colnames(ct[cnv,which(ct[cnv,]==0)])
  delet <- colnames(ct[cnv,which(ct[cnv,]<0)])
  medn <- median(as.numeric(ct[exp,normal]))
  medd <- median(as.numeric(ct[exp,delet]))
  data=rbind(cnv=as.numeric(ct[cnv,]),exp=log2(as.numeric(ct[exp,])+1))
  
  data=t(data)
  ksTest<-kruskal.test(exp ~ cnv, data=data)
  ksPval=ksTest$p.value
  f=factor(data[,"cnv"])
  labels=levels(f)
  
  pval=0
  if(ksPval<0.001){
    pval=signif(ksPval,4)
    pval=format(pval, scientific = TRUE)
  }else{
    pval=round(ksPval,3)
  }
  outTab <- rbind(outTab,cbind(gene=geneName,regulation="del",
                               median_normal=medn,median_cnv=medd,
                               pvalue=pval))
}

for (id in rownames(up_genes)) {
  geneName=id
  cnv=paste(geneName,"|cnv",sep="")
  exp=paste(geneName,"|exp",sep="")
  ct=rt[,rt[cnv,]>=0]
  ct[cnv,]=gsub("2","1",ct[cnv,])
  normal <- colnames(ct[cnv,which(ct[cnv,]==0)])
  amp <- colnames(ct[cnv,which(ct[cnv,]>0)])
  medn <- median(as.numeric(ct[exp,normal]))
  meda <- median(as.numeric(ct[exp,amp]))
  data=rbind(cnv=as.numeric(ct[cnv,]),exp=log2(as.numeric(ct[exp,])+1))
  
  data=t(data)
  ksTest<-kruskal.test(exp ~ cnv, data=data)
  ksPval=ksTest$p.value
  f=factor(data[,"cnv"])
  labels=levels(f)
  pval=0
  if(ksPval<0.001){
    pval=signif(ksPval,4)
    pval=format(pval, scientific = TRUE)
  }else{
    pval=round(ksPval,3)
  }
  outTab <- rbind(outTab,cbind(gene=geneName,regulation="amp",
                               median_normal=medn,median_cnv=meda,
                               pvalue=pval))
}

write.table(outTab,file = "cnv-Exp-ks.txt",sep = "\t",quote = F,row.names = F)
