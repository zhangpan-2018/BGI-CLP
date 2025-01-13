.libPaths("/data/public/zhangpanyu/Tools/R4.1.1/library")

args = commandArgs(T)
if(length(args) != 3){
	print("this script need args")
	print("Rscript *R genelist(symbol) species(hsa/mmu) outdir")
	quit()
}
genefile=args[1]
sp=args[2]
outdir = args[3]
if(!dir.exists(outdir)){dir.create(outdir,recursive=T)}

if(sp == "hsa"){
	OrgDbsp="org.Hs.eg.db"
	organismsp="hsa"
}else{
	OrgDbsp="org.Mm.eg.db"
	organismsp="mmu"
}
## 参数：gene、物种、outdir
gene=read.table(genefile)
genelist=as.vector(unlist(gene))

## library
library(clusterProfiler)
library("org.Hs.eg.db")
library(org.Mm.eg.db)
library(topGO)

test = bitr(genelist,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb=OrgDbsp)
write.table(as.data.frame(test),paste(outdir,'id.xls',sep="/"),sep="\t",row.names = F,quote = F)

################## GO 富集分析
## go file
ego_ALL <- enrichGO(gene = test$ENTREZID, 
                    OrgDb = OrgDbsp, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 1,
                    readable = TRUE) #Gene ID 转成gene Symbol ，易读barplot(ego_ALL, showCategory = 10)

ifelse(!dir.exists(paste(outdir,"/GO",sep="")), dir.create(paste(outdir,"/GO",sep="")), FALSE)

gobar=barplot(ego_ALL, showCategory = 10,colorBy=pvalue,color = "pvalue") 
#gobar=barplot(ego_ALL, showCategory = 10) 
png(paste(outdir,"/GO/","/go.","barplot.png",sep=""))
gobar
dev.off()
pdf(paste(outdir,"/GO/","/go.","barplot.pdf",sep=""),height=5)
gobar
dev.off()

godot=dotplot(ego_ALL, showCategory = 10,color = "pvalue")
png(paste(outdir,"/GO/","/go.","dotplot.png",sep=""))
godot
dev.off()
pdf(paste(outdir,"/GO/","/go.","dotplot.pdf",sep=""),height=5)
godot
dev.off()

write.table(as.data.frame(ego_ALL),paste(outdir,"GO",'go.txt',sep="/"),sep="\t",row.names = F,quote = F)
for (i in c("BP","MF","CC")){
  #unlink(paste(outdir,"/GO/",i,sep=""), recursive = TRUE)
  #dir=paste(outdir,"/GO/",i,sep="")
  #if(file.exists(dir)){}else{dir.create(paste(outdir,"/GO/",i,sep=""),recursive = TRUE)}
  ifelse(!dir.exists(paste(outdir,"/GO/",i,sep="")), dir.create(paste(outdir,"/GO/",i,sep="")), FALSE)
  type = ego_ALL[ego_ALL$ONTOLOGY==i,]
  gotype <- enrichGO(test$ENTREZID, OrgDb = OrgDbsp, ont=i,pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1,keyType = 'ENTREZID')

  barplot=barplot(gotype, showCategory = 10,colorBy=pvalue,color = "pvalue")
  png(paste(outdir,"/GO/",i,"/go.","barplot.png",sep=""))
  print(barplot)
  dev.off()
  pdf(paste(outdir,"/GO/",i,"/go.","barplot.pdf",sep=""),height=5)
  print(barplot)
  dev.off()
  
  dotplot=dotplot(gotype, showCategory = 10,color = "pvalue")
  png(paste(outdir,"/GO/",i,"/go.","dotplot.png",sep=""))
  print(dotplot)
  dev.off()
  pdf(paste(outdir,"/GO/",i,"/go.","dotplot.pdf",sep=""),height=5)
  #pdf(paste(outdir,"/GO/",i,"/go.","dotplot.pdf",sep=""))
  print(dotplot)
  dev.off()
 
#  png(paste(outdir,"/GO/",i,"/go.","plotGOgraph.png",sep=""))
#  plotGOgraph(gotype)##有向无环
#  dev.off()
#  pdf(paste(outdir,"/GO/",i,"/go.","plotGOgrapht.pdf",sep="")) 
#  plotGOgraph(gotype)
#  dev.off()
  
  write.table(as.data.frame(type),paste(outdir,"/GO/",i,"/","go.",i,".txt",sep=""),sep="\t",row.names = F,quote = F)
}


##################### KEGG富集分析
kegg <- enrichKEGG(test$ENTREZID, organism = organismsp, keyType = 'kegg', pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,pvalueCutoff = 1, qvalueCutoff = 1,use_internal_data = FALSE)
keggdir=paste(outdir,"/KEGG",sep="")
#unlink(paste(outdir,"/KEGG",sep=""), recursive = TRUE)
ifelse(!dir.exists(paste(outdir,"/KEGG",sep="")), dir.create(paste(outdir,"/KEGG",sep="")), FALSE)
write.table(as.data.frame(kegg),paste(outdir,"KEGG",'kegg.path.xls',sep="/"),sep="\t",row.names = F,quote = F)

keggbar=barplot(kegg, showCategory = 10,,colorBy=pvalue,color = "pvalue") 
png(paste(outdir,"/KEGG/","/kegg.","barplot.png",sep=""))
keggbar
dev.off()
pdf(paste(outdir,"/KEGG/","/kegg.","barplot.pdf",sep=""),height=5)
keggbar
dev.off()

keggdot=dotplot(kegg, showCategory = 10,color = "pvalue")
png(paste(outdir,"/KEGG/","/kegg.","dotplot.png",sep=""))
keggdot
dev.off()
pdf(paste(outdir,"/KEGG/","/kegg.","dotplot.pdf",sep=""),height=5)
keggdot
dev.off()








