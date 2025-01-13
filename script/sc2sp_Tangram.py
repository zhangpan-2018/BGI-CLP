# -*- coding: utf-8 -*-                                                                               
'''
Description: 
Author: Panyu Zhang
Date: 2022-11
E-mail: zhangpanyu@genomics.cn
'''

import os,sys
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import torch
sys.path.append('./')  # uncomment for local import
import tangram as tg

#%load_ext autoreload
#%autoreload 2
#%matplotlib inline

tg.__version__

def getDefaultColors(n, type = 1):
    if type == 1:
        colors = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941","#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC","#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007","#809693", "#1B4400", "#4FC601", "#3B5DFF", "#FF2F80","#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9","#B903AA", "#DDEFFF", "#7B4F4B", "#A1C299", "#0AA6D8", "#F38400","#A1CAF1", "#C2B280", "#848482", "#E68FAC", "#0067A5","#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300","#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26","#191970", "#000080","#6495ED", "#1E90FF", "#00BFFF", "#00FFFF", "#FF1493","#FF00FF", "#A020F0", "#63B8FF", "#008B8B", "#54FF9F","#00FF00", "#76EE00", "#FFF68F"]
    return colors


args = sys.argv
scfile=args[1]
keygroup=args[2]
spfile=args[3]
defile=args[4]
outdir=args[5]
maskfile=args[6]

os.makedirs(outdir, exist_ok=True)
os.chdir(outdir)
def split_list(input_list, n):
    avg = len(input_list) // n
    remainder = len(input_list) % n
    result = [input_list[i * avg + min(i, remainder):(i + 1) * avg + min(i + 1, remainder)] for i in range(n)]
    return result
# 定义一个函数来替换最后一个'-'为'_low'，但只有当它是最后一个字符时
def replace_last_dash(s):
    if s.endswith('-'):
        return s[:-1] + '_low'
    return s
def replace_last_dash_2(s):
    if s.endswith('+'):
        return s[:-1] + '_high'
    return s
############### load sptail data 
sppre=os.path.basename(spfile)
sppre=sppre.split('.')[0]
ad_sp = sc.read_h5ad(spfile)
if "CLP_1w" in spfile:
    ad_sp.obs['x']=ad_sp.obs['imagerow']
    ad_sp.obs['y']=ad_sp.obs['imagecol']
elif "CLP_4w" in spfile:
    ad_sp.obs['x']=ad_sp.obs['imagerow']
    ad_sp.obs['y']=ad_sp.obs['imagecol']
xs = ad_sp.obs.x.values
ys = ad_sp.obs.y.values
plt.axis('off')
plt.scatter(xs, ys, s=.7);
plt.gca().invert_yaxis()

############### load scRNA data
scpre=os.path.basename(scfile)
scpre=scpre.split('.')[0]
ad_sc = sc.read_h5ad(scfile)
ad_sc.obs['subcelltype'] = ad_sc.obs['subcelltype'].replace('DPsel2+', 'DPsel2')
ad_sc.obs['subcelltype'] = ad_sc.obs['subcelltype'].replace('DPsel3+', 'DPsel3')
ad_sc.obs['subcelltype'] = ad_sc.obs['subcelltype'].replace('DPsel4+', 'DPsel4')
ad_sc.obs['subcelltype'] = ad_sc.obs['subcelltype'].replace('DPsel5+', 'DPsel5')

keytype=['DPrea1', 'CD4', 'DPbla_CD2-', 'DPsel4', 'DPsel2', 'DPrea3', 'DN','DPrea2', 'NKT', 'Treg', 'CD8', 'DPsel1', 'DPrea4', 'B', 'DPsel3', 'ETP', 'DPsel5', 'DPbla_CD2+', 'Macro','mTEC-II-Aire+','mTEC-I-Aire-', 'cTEC','mTEC-III-Aire-','DPbla']
ad_sc = ad_sc[ad_sc.obs.subcelltype.isin(keytype),:]
ad_sc.obs['subcelltype'] = ad_sc.obs['subcelltype'].apply(replace_last_dash_2)
ad_sc.obs['subcelltype'] = ad_sc.obs['subcelltype'].apply(replace_last_dash)
ad_sc
ad_sc.X=ad_sc.raw.X
np.unique(ad_sc.X.toarray()[0, :])
sc.pp.normalize_total(ad_sc)
ad_sc.obs.subcelltype.value_counts()

##############  Prepare to map
if defile != "No":
    df_genes = pd.read_csv(defile, index_col=0)
    markers = list(df_genes.index)
    tg.pp_adatas(ad_sc, ad_sp, genes=markers)
else:
    tg.pp_adatas(ad_sc, ad_sp)

outpre=scpre+sppre
################# Map
ad_map = tg.map_cells_to_space(
    adata_sc=ad_sc,
    adata_sp=ad_sp,
    mode='clusters',
    cluster_label='subcelltype',
    device='cpu',
    density_prior='rna_count_based',
    num_epochs=500,
)

################### plot
tg.project_cell_annotations(ad_map, ad_sp, annotation="subcelltype")
ad_sp.write(outpre+"_"+keygroup+"_ad_sp.h5ad")
annotation_list = list(pd.unique(ad_sc.obs['subcelltype']))
if "CLP_2w" in spfile:
    tg.plot_cell_annotation_sc(ad_sp, annotation_list,perc=0.02,spot_size= 3, scale_factor=1.0)
elif "CLP_4w" in spfile:
    tg.plot_cell_annotation_sc(ad_sp, annotation_list,perc=0.02,spot_size= 3, scale_factor=1.0)
elif "CLP_1w" in spfile:
    tg.plot_cell_annotation_sc(ad_sp, annotation_list,perc=0.02,spot_size= 3, scale_factor=1.0)
else:
    tg.plot_cell_annotation_sc(ad_sp, annotation_list,perc=0.02,spot_size= 80, scale_factor=1.0)

plt.savefig(outpre+"_"+keygroup+"_SeuratClusters_reflect.png",bbox_inches="tight")
#ad_map.write(outpre+"_"+keygroup+"_ad_map2.h5ad")
#ad_map=sc.read_h5ad(outpre+"_"+keygroup+"_ad_map2.h5ad")

## no normalize 
ad_sp.obs.drop(annotation_list, inplace=True, errors="ignore", axis=1)
## construct df_plot
ad_sp.obs['percent.mt']=1;ad_sp.obs['batch']=1;ad_sp.obs['uniform_density']=1;ad_sp.obs['rna_count_based_density']=1
ad_sp.obsm["tangram_ct_pred"]["imagerow"]=ad_sp.obs['x']
ad_sp.obsm["tangram_ct_pred"]["imagecol"]=ad_sp.obs['y']
df_plot = ad_sp.obsm["tangram_ct_pred"][annotation_list]
#pd.concat([ad_sp.obs, df_plot,ad_sp.obsm["tangram_ct_pred"]["imagerow"],ad_sp.obsm["tangram_ct_pred"]["imagecol"]], axis=1)
combined_df = pd.concat([ad_sp.obs, df_plot,ad_sp.obsm["tangram_ct_pred"]["imagerow"],ad_sp.obsm["tangram_ct_pred"]["imagecol"]], axis=1)
#combined_df.to_csv(outpre + "_" + keygroup + "_Tangram.txt",sep="\t",index=False)
combined_df.to_csv(outpre + "_" + keygroup + "_Tangram.txt",sep="\t",index=True)
perc=0
## normalize 
df_plot = df_plot.clip(df_plot.quantile(perc), df_plot.quantile(1 - perc), axis=1)## normalize
df_plot = (df_plot - df_plot.min()) / (df_plot.max() - df_plot.min())
df_plot_normalized= df_plot.div(df_plot.sum(axis=1), axis=0)
#df_plot_normalized.to_csv(outpre + "_" + keygroup + "_TangramNormalized.txt",sep="\t",index=False)
df_plot_normalized.to_csv(outpre + "_" + keygroup + "_TangramNormalized.txt",sep="\t",index=True)
# 遍历 annotation_list 保存pdf
# 计算每个子列表的长度
list_length = len(annotation_list)
sublist_length = list_length // 5  # 整除，保证切分为5等份

# 打印切分后的子列表
n=5
sublists = split_list(annotation_list, n)
for i, sublist in enumerate(sublists):
       tg.plot_cell_annotation_sc(ad_sp, sublist, perc=0.02, spot_size=50, scale_factor=1.0)
       #tg.plot_cell_annotation_sc(ad_sp, sublist, perc=0.02, spot_size=80, scale_factor=1.0)
       plt.savefig(outpre + "_" + keygroup + "_Num"+ str(i) + "_SeuratClusters_reflect.pdf", bbox_inches="tight")
       plt.clf()

