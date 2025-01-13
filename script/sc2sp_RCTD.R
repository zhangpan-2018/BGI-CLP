.libPaths("/data/public/zhangpanyu/Tools/R4.1.1/library")

args=commandArgs(T)

#library(RCTD)
library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)
library("quadprog")
source("/data/public/zhangpanyu/Project/CLP/sc2sp/MultiROIRegionStatistics_BigCelltype_20240619/RCTDFunctionv2.R")
#source("/data/public/zhangpanyu/CLP/sc2sp/RCTD_allCluster/New1wROI/DPsel/new.R")

scrds=args[1]
sprds=args[2]
outdir=args[3]
sctime=args[4]
type=args[5]## clu/cell/sub
sc2spfi=args[6] ## roi经过allregion映射后，得到哪些是皮质哪些是髓质

sc2sp=readRDS(sc2spfi)
#dotsi=as.numeric(args[5])
roi=strsplit(outdir,"/")[[1]][length(strsplit(outdir,"/")[[1]])]
#roi=strsplit(dirname(outdir),"/")[[1]][length(strsplit(dirname(outdir),"/")[[1]])]
print(roi)
if(roi == "CLP_0w_ROI1"){dotsize = 0.7}else if(sctime=="CLP_0w"){dotsize = 0.6}else if(sctime == "CLP_1w"){dotsize = 1.5}else{dotsize=1}
#subcluster=c('CD8', 'DPbla', 'DN', 'Treg', 'CD4', 'DPrea', 'DPsel', 'ETP', 'NKT', 'Macro', 'B')
#subcelltype=c('B', 'CD4', 'CD8', 'DN', 'DPbla_CD2-', 'DPbla_CD2+',  'DPrea_IFN', 'DPrea_Ly6d_high1', 'DPrea_Ly6d_high2', 'DPrea_Ly6d_low','DPrea_Rag1', 'DPsel_Bcl2I11+', 'DPsel_Ccr7+', 'DPsel_Ikzf2+', 'DPsel_Item2a+', 'DPsel_Rag1_High', 'ETP', 'Macro', 'Treg','NKT')
keycelltype=c("DN","DPbla","DPrea","DPsel","CD8","CD4","Macro","Treg")

#subcelltype=c('DPbla_CD2+', 'DPbla_CD2-', 'DPrea_IFN', 'DPrea_Ly6d_high1', 'DPrea_Ly6d_high2', 'DPrea_Ly6d_low','DPrea_Rag1', 'DPsel_Bcl2I11+', 'DPsel_Ccr7+', 'DPsel_Ikzf2+', 'DPsel_Item2a+', 'DPsel_Rag1_High')
if(!dir.exists(outdir)){dir.create(outdir,recursive=T)}
setwd(outdir)
roi=strsplit(outdir,"/")[[1]][length(strsplit(outdir,"/")[[1]])]
############## scRNA
scRNA<-readRDS(scrds) 
sc <- subset(scRNA,subset = time == c(sctime))
#sc=sc[,sc$subcelltype %in% subcelltype]
sc=sc[,sc$celltype %in% keycelltype]
#### Load in/preprocess your data, this might vary based on your file type
counts <- as.matrix(sc@assays$RNA@counts )
meta_data <- sc@meta.data 
if(type == "clu"){
	cell_types <- meta_data$seurat_clusters; names(cell_types) <- rownames(meta_data)
}else if(type == "sub"){
	cell_types <- meta_data$subcelltype; names(cell_types) <- rownames(meta_data)
}else{
	cell_types <- meta_data$celltype; names(cell_types) <- rownames(meta_data)
}
cell_types <- as.factor(cell_types)
nUMI <- meta_data$nCount_RNA; names(nUMI) <- rownames(meta_data)
### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
saveRDS(reference,'ref.rds')

#reference=readRDS("/data/public/zhangpanyu/Tools/RCTD/ref.rds")
################## spRNA
spatial <- readRDS(sprds)
coords <- GetTissueCoordinates(spatial)
if(sctime == "CLP_0w" | sctime == "CLP_4w"){tmp=cbind(coords$imagecol,coords$imagerow);colnames(tmp)=c("imagerow","imagecol");rownames(tmp)=rownames(coords);coords=as.data.frame(tmp)}
spatialcounts <- as.matrix(spatial@assays$Spatial@counts)
spnUMI <- spatial@meta.data$nCount_Spatial
names(spnUMI)<-rownames(spatial@meta.data)
### Create SpatialRNA object
puck <- SpatialRNA(coords, spatialcounts, spnUMI)

################## predirct
myRCTD <- create.RCTD(puck, reference, max_cores = 1)
#myRCTD <- create.RCTD(puck, reference, max_cores = 20,test_mode = T)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') #doublet_mode = 'doublet'
save(myRCTD,file=paste0(sctime,"myRCTD.Rdata"))

############ result
results <- myRCTD@results
save(results,file=paste0(sctime,"RctdResult.Rdata"))
# normalize the cell type proportions to sum to 1.
#load(paste0(sctime,"myRCTD.Rdata"))
#load(paste0(sctime,"RctdResult.Rdata"))
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
cell_type_names <- myRCTD@cell_type_info$info[[2]] 
spatialRNA <- myRCTD@spatialRNA
resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)
# make the plots 
# Plots the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights,dotsize) 
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                      results$results_df) 

# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
############################################# 得到每个roi rctd与region结果  
regionration=cond_occur_ratio(cell_type_names,resultsdir,norm_weights,spatialRNA,sc2sp,roi,sctime)
write.table(regionration,file=paste0(resultsdir,"/",roi,"_regionratio.txt"),quote=F,row.names=F,sep="\t")
regioncount=cond_occur_count(cell_type_names,resultsdir,norm_weights,spatialRNA,sc2sp,roi)
write.table(regionration,file=paste0(resultsdir,"/",roi,"_regioncount.txt"),quote=F,row.names=F,sep="\t")
cortexplot=plot_cond_occur_region(cell_type_names,resultsdir,norm_weights,spatialRNA,sc2sp,"cortex")
medullaplot=plot_cond_occur_region(cell_type_names,resultsdir,norm_weights,spatialRNA,sc2sp,"medulla")

### 参考：https://www.jianshu.com/p/169cf7f394a1
### 其他可参考：https://www.jianshu.com/p/92a6df0dcb1b

