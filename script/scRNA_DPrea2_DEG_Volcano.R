library(ggplot2)
library(dplyr)
library(ggrepel)
#setwd("C:/Users/zhangpanyu/Desktop/CLP/火山图")
setwd("D:/BGI项目/CLP/火山图")
#DEG <- read.table("C:/Users/zhangpanyu/Desktop/CLP/火山图/CLP_0d7d10d14d28d.DPrea2.14dvs0d.fc02.txt", header=T, sep="\t")
#DEG <- read.table("C:/Users/zhangpanyu/Desktop/CLP/火山图/DEG_DPrea2_2w_0w.csv", header=T, sep=",")
DEG <- read.table("C:/Users/zhangpanyu/Desktop/CLP/火山图/DEG_DPrea2_2w_0w_new.csv", header=T, sep=",")
#DEG <- read.table("C:/Users/zhangpanyu/Desktop/CLP/火山图/DEG_DPrea2_2w_0w_noGm.csv", header=T, sep=",")
DEG <- read.table("D:/BGI项目/CLP/火山图/DEG_DPrea2_2w_0w_new.csv", header=T, sep=",")

## 按照pct.1过滤
DEG=DEG[DEG$pct.1>0.2,]
DEG$Gene=(DEG$X)

data <- 
  DEG %>% 
  mutate(change = as.factor(ifelse(p_val_adj < 0.05 & abs(avg_log2FC ) > 0.2,
                                   ifelse(p_val_adj < 0.05 &avg_log2FC  < -0.2,'Down','Up'),'No change')))# %>% 
  #rownames_to_column('gene')
table(data$change)

#up=data[data$change=='Up',]
#sortup=up[order(up$pct.1, decreasing = TRUE), ]
## 去掉GM，mt-，Rs-，Rp相关基因，按照pct.1从大到小排序，取上调和下调的top10

core_gene=c("Eef1a1","Tpt1","Tmsb4x","Tcf7","Ptma","Ccr9","Rack1","Eef2","Hsp90ab1","Ncl",
            "Trbc2","Cd8b1","H3f3a","Lck","Lat","Thy1","Ptprc","Rmnd5a","Trbc1","Fau")
logFC_t <- with(data,mean(abs(avg_log2FC)) + 2*sd(abs(avg_log2FC)))

## 这里使用的是动态阈值，也可以自定义例如logFC=1或2这种静态阈值
logFC_t <- round(logFC_t, 3) # 取前三位小数，这步也可以不运行
logFC_t #看一下动态阈值
this_tile="DPrea2_2w-vs-DPrea2_0w"

pdf("DPrea2_DEG_Volcano.pdf",7,6)
ggplot(
  # 数据、映射、颜色
  data, aes(x = avg_log2FC, y = -log10(p_val_adj), colour=change)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  scale_x_continuous(limits = c(-2, 2)) +
  # 辅助线
  geom_vline(xintercept=c(-0.2,0.2),linetype=2, color="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),linetype=2, color="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_classic(base_size = 14)+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+ggtitle(this_tile)+
  ## 关键基因标注
  geom_label_repel(data = filter(data, Gene %in% core_gene),
                 size = 3,box.padding = unit(0.5, "lines"),
                 segment.color = "black",
                 show.legend = FALSE,max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                 aes(label = Gene))+
  theme_bw()
dev.off()

