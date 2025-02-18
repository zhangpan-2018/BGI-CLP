.libPaths("zhangpanyu/Tools/R4.1.1/library")

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

rdsfile=args[1]# 
dedir=args[2] # 
group=args[3] # 
outdir=args[4] #
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


