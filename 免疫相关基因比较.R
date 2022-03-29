setwd("I:\\JINXING/data/yhq/")
rm(list = ls())
sessionInfo()
library(tidyverse)
library(pacman)
p_load(tidyverse,ggplot2, limma, export)

library(pheatmap)
# 准备数据
load("outdata/step3-2mrna.Rdata")
rm(expr)
head(mrna)
mrna <- t(mrna)
######免疫基因数据读取：包括固有免疫功能、抗原提呈能力##########
# immugene <- c('IRF3','MYD88','TICAM1','TLR3','TLR5','TLR7',
#               'TLR8','DDX58','IFIH1','MAVS','CLEC7A','CLEC4E',
#               'CD209','CLEC10A','NLRP3','AIM2','PYCARD', #innate免疫共17个分子initiation of innate immunity
#               'HLA-A','HLA-B','HLA-C','TAP1','TAP2','B2M','HLA-DPA1',
#               'HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DQB2',
#               'HLA-DRB4','HLA-DRB6','HLA-E','HLA-F','HLA-J', #抗原提呈系列MHC-I/II antigen-presenting process的HLA共16个分子
#               'CTLA4','CD48','PDCD1','LAG3','HAVCR1','BTN2A2',
#               'LAIR1','BTN3A1','PDCD1LG2','BTN1A1','VTCN1','BTNL2', #免疫检查点12个immune co-inhibitors
#               'ICOS','TNFRSF9','CD70','CD80','TNFRSF13B','CD27',
#               'SLAMF1','TNFRSF14','TNFRSF8','ICOSLG','CD58','TNFSF13','TNFSF15',
#               'TNFRSF4','HAVCR1','TNFSF18','TNFSF4','CD86' ) #co-stimulators
# 数据读取
{
  immugeneMHC <- data.table::fread(
    "Antigen presentation
HLA-A
HLA-B
HLA-C
HLA-DPA1
HLA-DPB1
HLA-DQA1
HLA-DQA2
HLA-DQB1
HLA-DQB2
HLA-DRA
HLA-DRB1
HLA-DRB3
HLA-DRB4
HLA-DRB5
MICA
MICB
") %>% 
    pull()
  
} # 读取抗原呈递的MHC
{
  immugeneINHI <- data.table::fread(
    "inhibitory
  ADORA2A
ARG1
BTLA
CD274
CD276
CTLA4
EDNRB
HAVCR2
IDO1
IL10
IL13
IL4
KIR2DL1
KIR2DL2
KIR2DL3
LAG3
PDCD1
SLAMF7
TGFB1
TIGIT
VEGFA
VEGFB
C10orf54
VTCN1
"
  ) %>% pull()
} # 读取免疫检查点 抑制免疫作用的
{
  immugeneSTIM <- data.table::fread(
    "Stimulatory
  BTN3A1
BTN3A2
CCL5
CD27
CD28
CD40
CD40LG
CD70
CD80
CX3CL1
CXCL10
CXCL9
ENTPD1
GZMA
HMGB1
ICAM1
ICOS
ICOSLG
IFNA1
IFNA2
IFNG
IL12A
IL1A
IL1B
IL2
IL2RA
ITGB2
PRF1
SELP
TLR4
TNF
TNFRSF14
TNFRSF18
TNFRSF4
TNFRSF9
TNFSF4
TNFSF9
"
  ) %>% pull()} # 读取免疫检查点 促进免疫作用的


immugene <- c('IRF3','MYD88','TICAM1','TLR3','TLR5','TLR7',
              'TLR8','DDX58','IFIH1','MAVS','CLEC7A','CLEC4E',
              'CD209','CLEC10A','NLRP3','AIM2','PYCARD', #innate免疫共17个分子initiation of innate immunity
              'HLA-A','HLA-B','HLA-C','TAP1','TAP2','B2M','HLA-DPA1',
              'HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DQB2',
              'HLA-DRB4','HLA-DRB6','HLA-E','HLA-F','HLA-J', #抗原提呈系列MHC-I/II antigen-presenting process的HLA共16个分子
              "ADORA2A", "ARG1", "BTLA", "CD274", "CD276", "CTLA4", "EDNRB", 
              "HAVCR2", "IDO1", "IL10", "IL13", "IL4", "KIR2DL1", "KIR2DL2", 
              "KIR2DL3", "LAG3", "PDCD1", "SLAMF7", "TGFB1", "TIGIT", "VEGFA", 
              "VEGFB", "C10orf54", "VTCN1", # 免疫检查点抑制基因 24个
              "BTN3A1", "BTN3A2", "CCL5", "CD27", 
              "CD28", "CD40", "CD40LG", "CD70", "CD80", "CX3CL1", "CXCL10", 
              "CXCL9", "ENTPD1", "GZMA", "HMGB1", "ICAM1", "ICOS", "ICOSLG", 
              "IFNA1", "IFNA2", "IFNG", "IL12A", "IL1A", "IL1B", "IL2", "IL2RA", 
              "ITGB2", "PRF1", "SELP", "TLR4", "TNF", "TNFRSF14", "TNFRSF18", 
              "TNFRSF4", "TNFRSF9", "TNFSF4", "TNFSF9" ) #免疫检查点激活基因37个



# 1. 选择相应的基因并且读取-----------------------------------------------------------------------

immugene %>% duplicated()  %>%  table()
Innateimmune <- mrna[, match(immugene,colnames(mrna),nomatch = 0)] %>% as.data.frame()



is.na(Innateimmune) %>% table()

Innateimmune$group <- "GroupA"
Innateimmune$group[rownames(Innateimmune) %in% sample_group$sample[23:44]] <- "GroupB"
table(Innateimmune$group)

# Innateimmune[, 59] <- NULL
# 过滤掉第一组，保留第二、第三组内的数据
# Innateimmune<-Innateimmune[match(GSEpred$X,rownames(Innateimmune)),]
# Innateimmune$TMEcluster<-GSEpred$GSEpred
################ 2. 画直方图，比较两组的差异###############
bodat <- Innateimmune %>% 
  select(group, everything())
bodat2 <- bodat %>% 
  pivot_longer(cols = c(2:ncol(bodat)))

library(ggpubr)

x <- compare_means(value ~ group, group.by = "name",data = bodat2)
xnosig <- x[which(x$p > 0.05),]

######### 3.免疫检查点的基因########
Innateimmune <- mrna[, match(x = c(immugeneINHI,immugeneSTIM),colnames(mrna),nomatch = 0)] %>% as.data.frame()

Innateimmune$group <- "GroupA"
Innateimmune$group[rownames(Innateimmune) %in% sample_group$sample[23:44]] <- "GroupB"
is.na(Innateimmune) %>% table()
table(Innateimmune$group)

bodat <- Innateimmune %>% 
  select(group, everything())
bodat2 <- bodat %>% 
  pivot_longer(cols = c(2:ncol(Innateimmune)))
x <- compare_means(value ~ group, group.by = "name",data = bodat2)
xnosig <- x[which(x$p > 0.05),]
bodatsig <- bodat2 %>% 
  filter(!name %in% xnosig$name)
ggboxplot(bodatsig,x="group",y="value",fill = "group",size = 0.1,
          facet.by = "name",
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("InfiltrationScore")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))

######## 4.选择一些免疫相关基因做图比较两组差异########

bodatpd <- bodat2 %>% 
  filter(name %in% c("PDCD1","CTLA4","CD274","LAG3","TGFB1","VEGFA","VEGFB"))
ggboxplot(bodatpd,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("log2(FPKM+1)")+
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.6))
graph2ppt(file = "output/plots/6免疫基因")


############## 5.画热图 求的是分组均数到总中位数的距离###################
Innateimmune <- mrna[, match(immugene,colnames(mrna),nomatch = 0)] %>% as.data.frame()
is.na(Innateimmune) %>% table()
Innateimmune$group <- "GroupA"
Innateimmune$group[rownames(Innateimmune) %in% sample_group$sample[23:44]] <- "GroupB"
table(Innateimmune$group)
InnateimmuneA <- Innateimmune %>% filter(group == "GroupA")
InnateimmuneA$group <- NULL
InnateimmuneA[] <- lapply(InnateimmuneA, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
MeanA <- as.data.frame(apply(InnateimmuneA,2,mean))

InnateimmuneB <- Innateimmune%>%filter(group == "GroupB")
InnateimmuneB$group <- NULL
InnateimmuneB[] <- lapply(InnateimmuneB, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

MeanB <- as.data.frame(apply(InnateimmuneB,2,mean))

Mean <- cbind( MeanA, MeanB)
colnames(Mean)<-c("GroupA","GroupB")

Innateimmune$group <- NULL
Innateimmune[] <- lapply(Innateimmune, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
Median <- as.data.frame(apply(Innateimmune,2,median))
Mean$Mut <- Mean$Mut - Median$`apply(Innateimmune, 2, median)`
Mean$Wild <- Mean$Wild - Median$`apply(Innateimmune, 2, median)`


#热图绘制
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(Mean,cluster_rows = F, cluster_cols = F,
         show_colnames =T,show_rownames = T, 
         color = colorRampPalette(color.key)(50),
         cutree_cols = 0,
         #annotation_colors = color.annotation,
         #annotation_col = annotation_col,
         cellwidth = 50,
         cellheight = 18,
         gaps_row = c(17,32,67),
         filename = "output/mianyijinrun.pdf"
)
dev.off()


############### 6.相关性图############
aidat <- read.table(file = "cibersort/immucellAI/ImmuCellAI_icb_result.txt",header = T)
aidat <- aidat %>% 
  select(Th17,nTreg, InfiltrationScore)
TCGAtmb<-read.csv("outdata/TCGA_LUAD_tmv.csv", row.names = 1)

tmb <- TCGAtmb %>% 
  filter(Tumor_Sample_Barcode %in% allid) %>% 
  select(Tumor_Sample_Barcode, total) %>% 
  rename(SampleID = Tumor_Sample_Barcode,
         tmb =  total)
rownames(tmb) <- tmb[,1]  
tmb <- tmb[allid, -1] 
repla <- mean(tmb,na.rm = T)
tmb <- replace_na(data = tmb,replace = repla)

corgene <- c("PDCD1","CTLA4","CD274","TGFB1","VEGFA","VEGFB", 
             "HLA-DQB2","CD40","CX3CL1","ICAM1","ITGB2")
cordat <- mrna[, match(corgene,colnames(mrna),nomatch = 0)]
cordat <- cbind(aidat,tmb,cordat)

cormat <- round(cor(cordat), 2) 

ggcorrplot(cormat, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_bw(),
)
graph2ppt(file = "output/plots/6免疫基因.pptx",append = T)
### 样式复杂的相关性图
p_load(PerformanceAnalytics)
my_data <- cordat
chart.Correlation(my_data, histogram=TRUE, pch=19)


