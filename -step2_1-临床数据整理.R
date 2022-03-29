rm(list = ls())
###################  整理phenotype的数据 还原OS.time的计算过程  ###############
library(pacman)
p_load(ggplot2, ggpubr, dplyr, survival, survminer, stringr,export)
getwd()
setwd("I:\\JINXING/data/yhq/")
getwd()
######### 导入xena的phenotype数据
phenodata <- data.table::fread(file ="data/esca_characteristic.csv")
sur_data <- read.csv(file = "data/43survivaldata.csv")
sur_data <- sur_data[,-1]
sur_data$id <- gsub("\\-",".",sur_data$id,)
phe_data <- phenodata[match(sur_data$id,phenodata$id), ]
phe_data <- phe_data[,!8:9]
## 使用cbind直接列的连接。前提条件是id顺序一致

phe_data <- cbind(phe_data,sur_data[,c(1,3,5)])

#导入表达矩阵看这43个人是否都有RNA测序信息
library(tidyverse)
expr <- data.table::fread(file = "data/expr.csv")
expr[1:4,1:4]
class(expr)

expr <- column_to_rownames(expr,var="V1")
expr[1:4,1:4]
intersect(colnames(expr),sur_data$id) %>% length()

str(phe_data)

# ###group必须是因子factor，代表有无突变
str(phe_data)
phe_data$group <- as.factor(phe_data$group)
phe_data$group  %>%  table()
phe_data$OS <- as.integer(phe_data$OS)
phe_data$OS.time <- phe_data$OS.time / 30
boxplot(phe_data$OS.time)

library(dplyr)

data1 <- phe_data %>% filter(OS.time<=100)

fit<-survfit(Surv(OS.time,OS)~group, data = data1) 
pp_os <- ggsurvplot(fit, data=data1, linetype = 1, palette = c("jco"),size=1,
                    surv.scale = c("percent"),pval = TRUE,
                    legend.title = '', legend.labs=c("Mut","Wild"),
                    break.time.by =12,
                    xlim = c(0,60),risk.table = TRUE,
                    risk.table.title = "Patients at risk",
                    ylab = "Overall Survival, %",xlab = "Months",
                    font.x = c(20,"plain","black"),font.y = c(20,"plain","black"),
                    font.tickslab = c(20,"plain","black"),
                    risk.table.fontsize = 6,font.legend =c(20,"plain","black"),
                    font.main = c(20,"plain","black"),pval.size = 8)

######risk table 的修改
pp_os$table <- ggpar(
  pp_os$table,
  font.title    = c(13, "bold.italic", "green"),  ###risk table 标题的修改
  font.subtitle = c(15, "bold", "pink"),  ###risk table小标题的修改
  font.caption  = c(11, "plain", "darkgreen"), ####插入字的修改
  font.x        = c(18, "plain", "black"), ### risk table x的修改
  font.y        = c(18, "plain", "black"),### risk table y的修改
  font.xtickslab = c(18, "plain", "black"),### risk table x 坐标轴的修改
  font.ytickslab = c(18),  ### risk table y 坐标轴的修改
  legend=c(0.8,0.88),
  censor.size=3
)
print(pp_os) ##终于做出和survival一样的生存分析图了


# 整理数据 --------------------------------------------------------------------
colnames(phe_data) <- c("sample","age","ajcc_M","ajcc_N","ajcc_T",
                        "BMI","sex","smoke","stage","group","OS","OS.time")

# age
phe_data$age

#age_medi 根据中位年龄分组
median_age <- median(phe_data$age, na.rm = T) #中位年龄为66
phe_data$age_median <- ifelse( phe_data$age >= median_age, 'older', 'younger' )
table(phe_data$age_median,useNA = "ifany")

# smoke smokesta 分成吸烟和不吸烟
phe_data$smoke
table(phe_data$smoke, useNA = "ifany")
phe_data$smoke <- factor(phe_data$smoke,levels = c(0,1),labels = c("No","Yes"))
table(phe_data$smoke, useNA = "ifany")

# T,N,M T分期不进行改动
table(phe_data$ajcc_T, useNA = "ifany") 

table(phe_data$ajcc_N, useNA = "ifany")

table(phe_data$ajcc_M, useNA = "ifany")

#stage
table(phe_data$stage, useNA = "ifany")
fix(phe_data)
table(phe_data$stage, useNA = "ifany")
phe_data$stage2 <- phe_data$stage

phe_data$stage2 <- ifelse(phe_data$stage %in% c("stage i","stage ia","stage ib"),"I",
                          ifelse(phe_data$stage %in% c("stage iia","stage iib"),"II",
                                 ifelse(phe_data$stage %in% c("stage iii","stage iiia","stage iiib","stage iiic"),"III",
                                        ifelse(phe_data$stage %in% c("stage iv"),"IV",NA))))

phe_data <- phe_f
table(phe_data$stage2, useNA = "ifany")
table(phe_data$stage, useNA = "ifany")

## 字符向量改成factor形式
str(phe_data)
library(tidyverse)
phe_data1 <- phe_data %>% 
  mutate(sex = as.factor(sex),
         ajcc_M = as.factor(ajcc_M),
         ajcc_N = as.factor(ajcc_N),
         ajcc_T = as.factor(ajcc_T), 
         stage2 = factor(stage2,levels = c("I","II","III","IV")),
         age_median = factor(age_median),
                  ) %>% 
  select(sample,group,sex,age,age_median,BMI,smoke,ajcc_T,ajcc_N,ajcc_M,stage2,OS,OS.time,) %>% 
  column_to_rownames("sample")
# csv保存在data文件夹
# rds文件保存在outdata文件夹
phe_f <- phe_data1
write.csv(phe_f,file = "data/phe_data.csv")
save(phe_f,file = "outdata/step2_phedata.rds")

