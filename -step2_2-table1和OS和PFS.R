p_load(do,tidyverse,ggplot2)
rm_all()
# 载入数据 --------------------------------------------------------------------
load("outdata/step2_phedata.rds")
# load("outdata/step4_grphe_f <- phe_data1
write.csv(phe_f,file = "data/phe_data.csv")
save(phe_f,file = "outdata/step2_phedata.rds")


# table one 制作 ------------------------------------------------------------------
library(compareGroups)
str(phe_f)
dt_factor <- phe_f

phetable2 <- descrTable(formula = group~. , # 公式，左边为分组，右边为变量
           data = dt_factor,  # 数据集
           method = 4, # 根据数据实际情况，自动选择统计方法
           alpha = 0.05,  # 显著性水平
           Q1 = 0.25, Q3 = 0.75,  # 默认输出p25和p75的分位数结果
           show.n = T, # 显示样本量
           show.ci = F, # 显示置信区间，默认是F
           conf.level = 0.95, # 置信区间范围
           type = 2, # 分类变量会显示频数和百分比
           show.p.overall = TRUE, # 显示P值
           digits.p = 3, # p值小数点位数
           sd.type = 2, # 1位mean（sd），2位mean±sd
)#加载包
phetable1 <- descrTable(group ~ ., data = dt_factor, method = 2,type = 2 )
print(phetable1)
print(phetable2)

export2word(phetable2,file = "I:/JINXING/data/yhq/output/Table1.docx",) #导出时注意 修改word名字

