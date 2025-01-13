#.libPaths("/jdfsbjcas1/ST_BJ/P21H28400N0233/zhangpanyu/Rpkgs/")                                
#.libPaths("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/lib/Rlib")

library(tidyverse)
library(CellChat)
library(ggalluvial)
library(Seurat)
options(stringsAsFactors = FALSE)

indir="/jdfsbjcas1/ST_BJ/P21H28400N0233/zhangpanyu/Project/CLP/scRNA/TEC_DP_cellchat/0.data/"
outdir="/jdfsbjcas1/ST_BJ/P21H28400N0233/zhangpanyu/Project/CLP/scRNA/TEC_DP_cellchat/resultdir/Secreted Signaling"
if(!dir.exists(outdir)){dir.create(outdir,recursive=T)}

step1outdir=paste0(outdir,"/step1_SingSample/")
step2outdir=paste0(outdir,"/step2_Compare/")
if(!dir.exists(step1outdir)){dir.create(step1outdir,recursive=T)}
if(!dir.exists(step2outdir)){dir.create(step2outdir,recursive=T)}

compare=c("CLP_0w_TcellTEC","CLP_2w_TcellTEC")
subcelltype=c("DPrea1","DPrea2","DPrea3","DPrea4","DPsel1","DPsel2","DPsel3","DPsel4","DPsel5","mTEC-II-Aire+","mTEC-I-Aire-","cTEC","mTEC-III-Aire-")
#subcelltype=c("DPrea1","DPrea2","DPrea3","DPrea4","DPrea5","DPsel1","DPsel2","DPsel3","DPsel4","DPsel5","mTEC-II-Aire+","mTEC-I-Aire-","cTEC","mTEC-III-Aire-")

for (i in compare){
  scRNA=readRDS(paste0(indir,"/",i,".rds"))
  scRNA=scRNA[,scRNA$subcelltype %in% subcelltype]
  delgene=read.table("/jdfsbjcas1/ST_BJ/P21H28400N0233/zhangpanyu/Project/CLP/scRNA/TEC_DP_cellchat/delgene.list",header=F)
  keepgene=setdiff(rownames(scRNA),delgene$V1)
  scRNA=subset(scRNA, features = keepgene)
  
  ## create cellchat object
  data.input <- GetAssayData(scRNA, assay = "RNA", slot = "data")
  #meta=scRNA@meta.data
  meta=data.frame(rownames(scRNA@meta.data),scRNA$subcelltype)
  colnames(meta)=c("cellid","subcelltype")
  rownames(meta)=meta$cellid
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "subcelltype")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "subcelltype")
  
  ## database
  CellChatDB <- CellChatDB.mouse
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
  #CellChatDB.use = CellChatDB
  cellchat@DB <- CellChatDB.use 
  
  ## cellchat 
  cellchat <- subsetData(cellchat)
  future::plan("multicore", workers = 4)##设置多线程
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat = filterCommunication(cellchat,min.cells=10)

## net
  df.net=subsetCommunication(cellchat)
  write.table(df.net,paste0(step1outdir,"/",i,"_Net_LR.xls"),quote=FALSE,row.names=F,sep="\t")
  
  cellchat <- aggregateNet(cellchat) 
  groupSize <- as.numeric(table(cellchat@idents))
## netP
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.table(df.net,paste0(step1outdir,"/",i,"_NetP.xls"),quote=FALSE,row.names=F,sep="\t")

## pdf
  pdf(paste0(step1outdir,"/",i,"_NetCountWeight.pdf"),20,20)
  countfig=netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  weightfig=netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  print(countfig)
  print(weightfig)
  dev.off()
  
  saveRDS(cellchat,paste0(step1outdir,"/",i,"_cellchat.rds"))
}


############### step2 compare
prename=paste0(compare[1],"_",compare[2])
con=readRDS(paste0(step1outdir,"/",compare[1],"_cellchat.rds"))
case=readRDS(paste0(step1outdir,"/",compare[2],"_cellchat.rds"))

object.list <- list(con = con, case = case)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat, paste0(step2outdir,"/",prename,"_merge.rds"))

###  不同生物条件的细胞通信网络的相互作用数量和强度
pdf(paste0(step2outdir,"/",prename,"_Compare.pdf"))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()
### 不同细胞群之间的相互作用数量或相互作用强度不同
## 圆图
pdf(paste0(step2outdir,"/",prename,"Compare_circles.pdf"))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
## heatmap
pdf(paste0(step2outdir,"/",prename,"_Compare_heatmap.pdf"))
netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
dev.off()

# 比较配受体对的差异
pdf(paste0(step2outdir,"/",prename,"_Compare_LR.pdf"))
#compareInteractions(cellchat, show.legend = F, group = c(1,2))
#compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
#dev.off()
#netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in CLP2w", angle.x = 45, remove.isolate = T)
netVisual_bubble(cellchat,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in CLP2w", angle.x = 45, remove.isolate = T)
netVisual_bubble(cellchat,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in CLP2w", angle.x = 45, remove.isolate = T)
dev.off()

# 比较信号通路的差异
p1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
p2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

pdf(paste0(step2outdir,"/",prename,"_compareNet.pdf"),8,8)
p1
p2
dev.off()

# 差异配受体对表格
diff_interactions <- subsetCommunication(cellchat)
write.csv(diff_interactions, paste0(step2outdir,"/",prename,"_differential_interactions.csv"))

# 差异信号通路表格
pathway.union <- union(con@netP$pathways, case@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(con, signaling = pathway.union, title = "Control", width = 5, height = 6)
ht2 <- netAnalysis_signalingRole_heatmap(case, signaling = pathway.union, title = "Treatment", width = 5, height = 6)

plotGrandComparison(cellchat, show.legend = T, justified = T)
# 信号通路差异可视化
ht1 + ht2
# 网络图可视化
netVisual_diffInteraction(cellchat_merged, weight.scale = T)
netVisual_heatmap(cellchat_merged)
