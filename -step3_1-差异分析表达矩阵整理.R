### 加载数据和R包########
rm(list = ls())
load(file = "outdata/差异分析degroup.Rdata")
load(file = "outdata/expr_整理好.Rdata")
library(pacman)
p_load(limma,edgeR,ggplot2,ggthemes,ggpubr,ggrepel,export,data.table,tidyverse)
sample_group <- degroup
sample_group$group %>%  table()
levels(sample_group$group) <- c("clusterA","clusterB")
levels(sample_group$group) 
# A预后好，B预后差 比较的时候用B-A,一共44个人，AB各22人
# 设定阈值

rt <- as.matrix(expr)  
rt <- rt[,sample_group$sample]
rt[1:3,1:3]
dim(rt) # 58385    44
rt <- rt[rowSums(rt) > 0, ]    
dim(rt) # 55624    44
rt[1:3,1:3]


rownames(rt) %>% duplicated() %>% table()
rownames(expr) %>% duplicated() %>% table()

####### 基因ID转换及选择mrna #####
{
  load(file = "outdata/gtf注释文件.Rdata")
  ## 选择proteincoding的蛋白质
  gtf_gene_protein <- gtf_data %>% 
    dplyr::rename( Ensembl_ID = gene_id,
                   Symbol  = gene_name,
                   Biotype = gene_type  ) %>% 
    dplyr::filter(Biotype %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", 
                                 "IG_M_gene", "IG_V_gene", "IG_Z_gene", 
                                 "nonsense_mediated_decay", "nontranslating_CDS", 
                                 "non_stop_decay", "polymorphic_pseudogene", 
                                 "protein_coding", "TR_C_gene", "TR_D_gene", "TR_gene", 
                                 "TR_J_gene", "TR_V_gene")) %>%
    select(Ensembl_ID,Symbol,Biotype) %>% unique()
  index <- intersect(rownames(rt),gtf_gene_protein$Symbol)
  
  stopifnot(
    {
      duplicated(rownames(rt))== F
    }
  )
  
  # 查看行和列数
  dim(rt) # 19907    44
  mrna <- rt[index,]   ## 5万多选2万左右的基因为编码蛋白的mrna
  # 重复的基因名字进行合并且取均值
  dim(mrna)  #19907    44
  mrna <- limma::avereps(mrna) #19971 去重
  dim(mrna)  #19907    44
  
}

save(mrna,sample_group,file = "outdata/step3-2mrna.Rdata")
