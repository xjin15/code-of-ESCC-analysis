rm(list = ls())
setwd(dir = "I://JINXING/data/yhq/")
all_diff<- data.table::fread(file="outdata/菌群cluster_mrna_all_diff.csv")
library(clusterProfiler)
library(org.Hs.eg.db) 
library(tidyverse)
FCvalue <- 1
p <- 0.05
#将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
#将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因

all_diff$group <-  "not-significant"
all_diff$group[which((all_diff$P.Value< p) & (all_diff$logFC> FCvalue))] = "up-regulated"
all_diff$group[which((all_diff$P.Value< p) & (all_diff$logFC< -FCvalue))] = "down-regulated"
table(all_diff$group)
all_diff$group <- as.factor(all_diff$group)
all_diff$logP <- -log10(all_diff$P.Value)

degupdown <- all_diff[all_diff$group != "not-significant",]
#############################################功能富集
###1.GO
# GO三大注释——BP:生物学过程   CC:细胞学组分  MF：分子生物学功能

##三种基因ID类型
#ENTREZID:   4312           8318  
#SYMBOL:     MMP1           CDC45
#ENSEMBL:  ENSG00000196611  ENSG00000093009
library(clusterProfiler)
library(org.Hs.eg.db)
genename <- as.character(degupdown$ID)
#将基因名转换为ENTREZID格式
#基因ID类型为ENSEMBL的ID形式，选择BP功能组(ont="ALL"同时展示BP/CC/MF)
Go_resultBP <- enrichGO(genename, 'org.Hs.eg.db', 
                        keyType = "SYMBOL", 
                        ont="BP", 
                        pvalueCutoff=0.05) 
barplot(Go_resultBP, showCategory=10)
library(export)
graph2ppt(file="output/GObarplot.pptx", width = 8.5, aspectr = 1.5, append = T)
dotplot(Go_resultBP, showCategory=10)
graph2ppt(file="output/GOdotplot.pptx", width = 8.5, aspectr = 1.5, append = T)
###2.KEGG(需要用ENTREZID形式的数据)
library(org.Hs.eg.db)
library(clusterProfiler)
genename <- bitr(genename, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genename<-genename$ENTREZID
enrich_KEGG <- enrichKEGG(genename,
                          organism = "hsa",
                          pvalueCutoff = 0.1)
dotplot(enrich_KEGG, showCategory=10) 
graph2ppt(file="output/KEGGdotplot.pptx", width = 8.5, aspectr = 1.5, append = T)
