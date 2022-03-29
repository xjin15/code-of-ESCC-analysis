rm(list = ls())

getwd()
setwd(dir = "I://JINXING/data/yhq/")

group <- data.table::fread("data/yhq/TCGA肿瘤样本分组.csv")
survival <- read.csv("data/yhq/ESCAsurvival整理.csv")

test <- gsub(pattern = "\\.", replacement = "-", x=survival$id)

survival$id <- test

survival_1 <- subset(survival, survival$id == group$sample)
survival_1 <- subset(survival,  group$sample == id)
survival_1 <- survival[ survival$id %in% group$sample,]

surv_id <- survival$id
group_id <- group$sample
intersect(group_id,surv_id)

survival <- survival[survival$id %in% group$sample , ]
survival$X <- group$group[match(x= survival$id, table = group$sample, nomatch = 0)]
colnames(survival)[1] <- "group"


###整理生存数据
data_sur <- survival[,c(2,1,3,5)]
# OS时间/30 
data_sur$OS.time <- data_sur$OS.time / 30

library(survival)
library(survminer)


fit<-survfit(Surv(OS.time,OS)~group,data=data_sur)  
ggsurvplot(fit, data=data_sur, linetype = 1, 
           palette = c("jco"),size=1,
           surv.scale = c("percent"),
           pval = TRUE,
           legend.title = "Indentation",
           legend.labs = c("1","2"),
           break.time.by =20,xlim = c(0,60),
           risk.table = TRUE,
           risk.table.title = "Patients at risk",
           ylab = "Overall Survival, %",
           xlab = "Months",
           font.x = c(30,"plain","black"),
           font.y = c(30,"plain","black"),
           font.tickslab = c(30,"plain","black"),
           risk.table.fontsize = 8,font.legend =c(30,"plain","black"),font.main = c(30,"plain","black"),pval.size = 12)

export::graph2ppt(file = "data/yhq/output/surival_plot.pptx",
                  append = T, # 同一个文件新出的图加到PPT下一页中
                  width = 12, height = 8, )
write.csv(survival,file="data/yhq/data/43survivaldata.csv")

save(survival,file="data/yhq/outdata/step1_43survival_data.rds")
