install.packages("compareGroups")
library(compareGroups)
data("predimed")
str(predimed)
descrTable(~.,data = predimed)
descrTable(group ~.,data = predimed)
descrTable(famhist ~.,data = predimed)


rm(list = ls())
getwd()
setwd(dir = "I://JINXING/data/yhq/")
phenodata<-data.table::fread("data/phe_data.csv")
str(phenodata)
library(tidyverse)
phenodata[1:4,1:4]
class(phenodata)
phenodata<-column_to_rownames(phenodata,var="V1")
str(phenodata)
phenodata$group <- as.factor(phenodata$group)
phenodata$group  %>%  table()
phenodata$OS <- as.integer(phenodata$OS)
boxplot(phenodata$OS.time)
str(phenodata)
phenodata$sex <- as.factor(phenodata$sex)
phenodata$ajcc_M = as.factor(phenodata$ajcc_M)
phenodata$ajcc_M  %>%  table()

phenodata$ajcc_N = as.factor(phenodata$ajcc_N)
phenodata$ajcc_N  %>%  table()

phenodata$ajcc_T = as.factor(phenodata$ajcc_T) 
phenodata$ajcc_T  %>%  table()

phenodata$stage2 = as.factor(phenodata$stage2)
phenodata$stage2 %>% table()
phenodata$age_median = as.factor(phenodata$age_median) 
str(phenodata)
descrTable(~.,data =phenodata)

shapiro.test(phenodata[phenodata$group==1,"BMI"])
# W = 0.88522, p-value = 0.0265 p值＜0.05 説明是非正態分佈 

shapiro.test(phenodata[phenodata$group==1,"age"])
# W = 0.97685, p-value = 0.874 p值＞0.05 説明是正太分佈


phetable2 <- descrTable(group ~.,data =phenodata,method = NA)
export2word(phetable2,file = "I:/JINXING/data/yhq/output/Table1.docx",) #导出时注意 修改word名字
