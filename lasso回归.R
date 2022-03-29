setwd("G:\\stem celll\\Fig4\\monocle")
pbmcmeta<-read.csv("pbmcmeta.csv",row.names = 1)
stem<-AddMetaData(stem,metadata = pbmcmeta)


DimPlot(stem,group.by = "stage2")
DimPlot(stem)
pbmc<-stem

library(monocle)
library(Seurat)
library(tidyverse)
library(ggsci)
pbmc<-stem
Idents(pbmc)<-pbmc$stage2
############ befroe monocle analysis, you should selected some important genes for downstream analysis, such as marker genes####
logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_logFC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

######### sample cells for qucik analysis######
expm<-pbmc@assays$RNA@counts
#num = length(expm[1,])
#train = sample(1:num,(num/2)*1)
#exprMat<-expm[,train]  ###count data for downstream analysis########
exprMat<-expm 
pd<-pbmc@meta.data  ###metadata
pd$cell2<-pbmc@active.ident ###cell2 means celltype
pd<-pd[colnames(exprMat),]
data<-exprMat
pd<-as.data.frame(pd)
fd=data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd<-pd[colnames(data),]
colnames(pd)<-gsub("_","",colnames(pd))
fd=cbind(id=row.names(fd),fd)
pd<-new("AnnotatedDataFrame", data = pd)
fd<-new("AnnotatedDataFrame", data = fd)
cds<-newCellDataSet(data,phenoData = pd,featureData = fd)
names(pData(cds))[names(pData(cds))=="cell2"]="Cluster2"  ## you need change cell2 or cluster2 based on your own metadata####
pData(cds)[,"Cluster2"]=paste0("",pData(cds)[,"Cluster2"])
cds <- estimateSizeFactors(cds)
cds <- setOrderingFilter(cds, sig.markers$gene)  ###选择特征基因
cds <- reduceDimension(cds, max_components =2,reduction_method ="DDRTree")
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = 3)
plot_cell_trajectory(cds, color_by = "Cluster2",cell_size = 1.5)+scale_color_npg()
graph2pdf(file="monoclcluster.pdf",width=10,height=10)
######  you can run following code to change the root #####   
######cds <- orderCells(cds,root_state = xxx)########
plot_cell_trajectory(cds, color_by = "Pseudotime")+scale_colour_viridis_c(option = "inferno")
graph2pdf(file="monoclestemtumortime.pdf",width=10,height=10)
plot_cell_trajectory(cds, color_by = "stage2",cell_size = 1.5)+scale_color_npg()
graph2pdf(file="monoclestage2.pdf",width=10,height=10)

##########plot gene monocle#####
pic<-getwd()
setwd("G:/stem celll/Fig4/monocle/gene")
for (gene2 in top10$gene ) {
  to_be_tested<-gene2
  cds_subset <- cds[to_be_tested,]
  figure <-plot_genes_in_pseudotime(cds_subset,color_by ="Pseudotime",cell_size = 1.5)+scale_colour_viridis_c(option = "inferno")
  ggsave(paste(pic,"\\",gene2,".tiff",sep=""),width = 8, height=8,figure)
}
######################plot heatmap#############
phd<-plot_pseudotime_heatmap(cds[unique(top10$gene),],
                             num_clusters = 3,
                             cores = 1,
                             return_heatmap = TRUE,
                             show_rownames = TRUE)


phd$tree_row
clu<-cutree(phd$tree_row,k=3)
clusi<-data.frame(clu)
clusi[,1]<-as.character(clusi[,1])  #####Now you get the gene information of the heapmap and can use these genes for downstream analysis#####

##################plot fate heatmap############

cds<-estimateDispersions(cds, cores=4, relative_expr = TRUE)
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- cds[disp.genes,]  ####or you can use the marker genes####mycds_sub <- cds[top10$gene,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]

plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 3, show_rownames = T,return_heatmap = T)

graph2pdf(file="HEATMAP.pdf",width=10,height=15)
graph2ppt(file="HEATMAP.pptx",width=10,height=15)

#################### calculate the pathways and Pseudotime corrections#########
library(scCancer)
library(ggplot2)
gene<-getDefaultGeneSets(species = "human")
mai<-runGeneSets(pbmc,geneSets=gene, method = "average")
pbmc$Cluster<-pbmc@active.ident
meta<-pbmc@meta.data
meta2<-cbind(mai,meta)  
prefix = "GS__"
gs.name <- colnames(meta2)
gs.name <- grep(paste0("^", "GS__"), gs.name, value = TRUE)
pathway <- meta2[, gs.name]
colnames(pathway) <- substr(gs.name, 1 + nchar(prefix), 200)
pathway<-pathway[colnames(exprMat),]  ####selected the sample cell
limitData <- function(data, min = NULL, max = NULL){
  data2 <- data
  if(!is.null(min)){
    data2[data2 < min] <- min
  }
  if(!is.null(max)){
    data2[data2 > max] <- max
  }
  return(data2)
}   ####deal with the data
cds_exprs<-pathway
timedata<-cds@phenoData@data[rownames(cds_exprs),]
cds_exprs<-as.data.frame(cds_exprs)
cds_exprs$Pseudotime<-timedata$Pseudotime
pic<-getwd()
for (paths in colnames(pathway)) {
  Pseudotime<-"Pseudotime"
  Q <- quantile(cds_exprs[,paths], probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(cds_exprs[,paths])
  up <- Q[2]+1.5*iqr # Upper Range 
  low<- Q[1]-1.5*iqr # Lower Range 
  cds_exprs[,paths] <- limitData(cds_exprs[,paths], min = low, max = up)
  figure <- ggplot(cds_exprs, aes_string(Pseudotime,paths))+geom_smooth(method="loess",level=0.6,size=2.5,colour="#F8766D")+theme_classic()+theme(axis.line.x=element_line(linetype=1,color="black",size=1.5,lineend = 5),
                                                                                                                                                  axis.line.y=element_line(linetype=1,color="black",size=1.5,lineend = 5))   ###### use aes_string!!
  ggsave(paste(pic,"\\",paths,".tiff",sep=""), width = 8, height=8,figure)}  ####then you can choose the pathway for downstream analysis#######
ggsave(paste(pic,"\\",paths,".tiff",sep=""), width = 8, height=8,figure)
graph2ppt(x=figure,file=paste(pic,"\\",paths,".pptx",sep=""),width = 10, height=6)
graph2ppt(x=figure,file=paste(pic,"\\",paths,".pptx",sep=""),width = 10, height=6)
####################################################







