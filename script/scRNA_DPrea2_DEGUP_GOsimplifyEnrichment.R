library(simplifyEnrichment)
set.seed(888)
go_id = random_GO(500)
head(go_id)
setwd("C:/Users/zhangpanyu/Desktop/CLP/火山图")
go=read.csv("go.BP.txt",sep="\t",header = T)
siggo=go[go$pvalue < 0.05,]
siggoid=siggo$ID
mat = GO_similarity(siggoid)
simplifyGO(mat, word_cloud_grob_param = list(max_width = 80))

