.libPaths("/data/public/zhangpanyu/Tools/R4.1.1/library")

args=commandArgs(T)
if(length(args) != 4){ 
	print("this script need args")
	print("Rscript *R scrds sprds sctime routdir(outdir=routdir/sctime") 
	quit()
}

library(ImageGP)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(egg)
library(ggrepel)
library(Seurat)

rdsfile=args[1]# /data/public/zhangpanyu/Project/CLP/scRNA/DPsel_DEGEnrich_2w-vs-0w_4w-vs-0w/0.data/CLP_0d7d10d14d28d.DPsel.rds
dedir=args[2] # /data/public/zhangpanyu/Project/CLP/scRNA/DPsel_DEGEnrich_2w-vs-0w_4w-vs-0w/DegVolcano/DEG
group=args[3] # CLP_4w-vs-CLP_0w
outdir=args[4] # /data/public/zhangpanyu/Project/CLP/scRNA/DPsel_DEGEnrich_2w-vs-0w_4w-vs-0w/DegVolcano/VolcanoPlot
if(!dir.exists(outdir)){dir.create(outdir,recursive=T)}
casename=strsplit(group,"-vs-")[[1]][1]
contralname=strsplit(group,"-vs-")[[1]][2]

sel=readRDS(rdsfile)
celltype=unique(sel$celltype)
for (i in celltype){
  ## avg exp
  con=sel[,sel$time==contralname]
  case=sel[,sel$time==casename]
  concell=con[,con$celltype==i]
  casecell=case[,case$celltype==i]

  aveexpcon=AverageExpression(concell)
  aveexpcon=as.data.frame(aveexpcon$RNA)
  aveexpcase=AverageExpression(casecell)
  aveexpcase=as.data.frame(aveexpcase$RNA)


  ## de file
  defile=paste0(dedir,"/",group,"/DEG_",group,"_",i,".CSV")
  deg=read.csv(defile,header=T)
  #deg$updown="Nosig"
  #deg$updown[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0.2 & deg$pct.1>0.2] = "up"
  #deg$updown[deg$p_val_adj < 0.05 & deg$avg_log2FC < -0.2 & deg$pct.1>0.2] = "down"
  #deg <- subset(deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.2 & pct.1>0.2)

  expcon=aveexpcon[deg$X,]
  expcase=aveexpcase[deg$X,]
  
  ## plot
  plotfile=cbind(deg,expcon,expcase)
#  plotfile$expconlog=log2(plotfile$expcon+1)
#  plotfile$expcaselog=log2(plotfile$expcase+1)
  plotfile$expconlog=(plotfile$expcon)
  plotfile$expconlog[plotfile$expconlog>30]=30
  plotfile$expcaselog=(plotfile$expcase)
  plotfile$expcaselog[plotfile$expcaselog>30]=30
  pdf(paste0(outdir,"/",group,"_",i,"_Volcano.pdf"),4,3)
  p=sp_scatterplot(plotfile, xvariable = "expconlog", yvariable = "expcaselog",color_variable = "updown",title =paste0(i," ",group),color_variable_order = c("Nosig","up", "down"),manual_color_vector = c("grey","red","blue")) + coord_fixed(1)+ labs(x = contralname, y = casename)
  #p=sp_scatterplot(plotfile, xvariable = "expcaselog", yvariable = "expconlog",color_variable = "updown",title =paste0(i," ",group),color_variable_order = c("Nosig","up", "down"),manual_color_vector = c("grey","firebrick","dodgerblue")) + coord_fixed(1)+ labs(x = casename, y = contralname)
  print(p)
  dev.off()
  
}


#sel2w=sel[,sel$time=="CLP_14d"]
#sel12w=sel2w[,Idents(sel2w)=="DPsel1"]
#sel0w=sel[,sel$time=="CLP_0w"]
#sel10w=sel0w[,Idents(sel0w)=="DPsel1"]
#sel1deg=read.csv("/data/public/zhangpanyu/Project/CLP/scRNA/DPsel_DEGEnrich_2w-vs-0w_4w-vs-0w/DEG/CLP_10d-vs-CLP_0w/DEG_CLP_10d-vs-CLP_0w_DPsel1.CSV")
#
#aveexp0w=AverageExpression(sel0w)
#aveexp0w=as.data.frame(aveexp0w$RNA)
#aveexp2w=AverageExpression(sel2w)
#aveexp2w=as.data.frame(aveexp2w$RNA)
#
#exp0w=aveexp0w[sel1deg$X,]
#exp0w=exp0w$DPsel1
#exp2w=aveexp2w[sel1deg$X,]
#exp2w=exp2w$DPsel1
#
#### 
#plotfile=cbind(sel1deg,exp0w,exp2w)
#
#plotfile$exp0wlog=log2(plotfile$exp0w+1)
#plotfile$exp2wlog=log2(plotfile$exp2w+1)
#
##sp_scatterplot(plotfile, xvariable = "exp2w", yvariable = "exp0w", 
##                 color_variable = "updown",
##                title ="DPsel1 2w-vs-0w", 
##               color_variable_order = c("up", "down"),
##              manual_color_vector = c("firebrick","dodgerblue")) + 
##  coord_fixed(1)+ labs(x = "exp2w", y = "exp0w")
#
#pdf("/data/public/zhangpanyu/Project/CLP/scRNA/DPsel_DEGEnrich_2w-vs-0w_4w-vs-0w/test_voc.pdf")
#sp_scatterplot(plotfile, xvariable = "exp2wlog", yvariable = "exp0wlog", 
#                 color_variable = "updown",
#                title ="DPsel1 2w-vs-0w", 
#               color_variable_order = c("up", "down"),
#              manual_color_vector = c("firebrick","dodgerblue")) + 
#  coord_fixed(1)+ labs(x = "exp2wlog", y = "exp0wlog")
#dev.off()
#
#
