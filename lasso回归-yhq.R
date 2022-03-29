rm(list = ls())
# LASSO挑选合适的基因作为signature ------------------------------------------------
options(digits = 4) # 设置小数点后4位
p_load(ggplot2,ggthemes,ggpubr,ggrepel,export,tidyverse,glmnet,corrplot)
 # 批量加载要用的包
# 数据准备 --------------------------------------------------------------------
load("outdata/expr_整理好.Rdata") #读取已经保存好的表达矩阵
load("outdata/step3-2mrna.Rdata") # 读取mrna，只有编码蛋白的基因
load("outdata/差异分析degroup.Rdata") # 读取分组信息
library(corrplot) #相关系数分析
library(glmnet) # 岭回归、 LASSO 、弹性网络模型
library(MASS) 

# epxr和mrna转成 行为样本，列为基因及其表达量的matrix。

mrna <- na.omit(mrna) # 去除有na值的行
mrna <- t(mrna) # 转置
mrna[1:4,1:4]

set.seed(123) #random number generator 
x <- as.matrix(mrna)  # 预测变量为多个因素，放在x中
# 如果只用DEG的133个基因做预测，那就
load("outdata/step3.2_差异mrna.rds")
X2<-intersect(DEG$genes,colnames(mrna))
x3 <- x[ ,X2]
y <- as.numeric(sample_group$group) # y就是因变量

y

######## 10折交叉验证挑选合适的λ值 #######
# 交叉验证

set.seed(8)
cvlasso <- cv.glmnet(x = x,y = y,
                     family = "binomial", # 因变量是二分类
                     alpha = 1, # lasso回归
                     type.measure = "class", # 也可以不选
                     nfolds = 5,)
plot(cvlasso, xvar = "lambda", label = T) # 回归系数随着logλ的变化曲线
graph2ppt(file = "output/LASSO.pptx", append = T)

lasso_deg <- glmnet(x = x,y = y,family = "binomial",alpha = 1,standardize = T,)
plot(lasso_deg,xvar = "lambda")
graph2ppt(file = "output/LASSO.pptx", append = T)

model_lasso_min <- glmnet(x = x, y= y,
                          alpha = 1,
                          lambda=cvlasso$lambda.min,
                          family = "binomial")
model_lasso_1se <- glmnet(x = x, y= y,
                          alpha = 1,
                          lambda=cvlasso$lambda.1se,
                          family = "binomial")
head(model_lasso_min$beta)
head(model_lasso_1se$beta)
choose_gene_min <- rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se <- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]

betamin <- model_lasso_min$beta[as.numeric(model_lasso_min$beta)!=0]
beta1se <- model_lasso_1se$beta[as.numeric(model_lasso_1se$beta)!=0]

length(choose_gene_1se);choose_gene_1se;beta1se
length(choose_gene_min);choose_gene_min;betamin

lasso.prob <- predict(cvlasso, newx=x , s=c(cvlasso$lambda.min,cvlasso$lambda.1se) )
lasso.prob
library(pROC)
plot.roc(roc(y,lasso.prob[,1]),
         xlab = "特异度",
         ylab = "灵敏度",
         reuse.auc=T,
         print.thres = T,
         print.auc = T)

plot.roc(roc(y,lasso.prob[,2]),
         xlab = "特异度",
         ylab = "灵敏度",
         reuse.auc=T,
         print.thres = T,
         print.auc = T,
         # add =T
         )
plot(smooth(roc(y,lasso.prob[,2])),
     add = T)

graph2ppt(file = "output/LASSO.pptx", append = T)

########### 用筛选到的特征基因做logistic回归分析#########
length(choose_gene_1se);choose_gene_1se;beta1se
length(choose_gene_min);choose_gene_min;betamin

# 最终选择 min2 的19个基因表达量为突变分组条件
choose_gene_min
options(digits = 3)
betamin <- model_lasso_min$beta[as.numeric(model_lasso_min$beta)!=0] 
betamin
x[1:4,1:4]
## 写公式，算分数 #####
formular_call = paste(round(betamin,digits = 4), choose_gene_min,sep = " * ",collapse= ") + (");formular_call
# 便于复制到word中
paste(round(betamin,digits = 4), choose_gene_min,sep = "×",collapse = "）+（")
paste(choose_gene_min,collapse = "、")

## 计算fpkm中的score
x[1:4,1:4]

logi_x <- x %>% 
  as.data.frame() %>% 
  # rownames_to_column("ID") %>% 
  transmute(score = (0.1093 * SNX3) + 
  (0.2403 * AKIRIN2) + (-0.0616 * TMEM87B) +   (-0.4857 * STEAP3) + 
  (0.1149 * PPME1) + (0.0511 * LGALS7B) + (-0.168 * ARFRP1) + 
  (-0.1601 * KIAA0556) +  (-0.0292 * CLN3) + (0.0682 * STX11) + 
  (-0.2265 * `RP11-295P9.3`) + (-0.1583 * `RP11-434D12.1`)) 

logi_x$group <- as.factor(sample_group$group)
# 建logistic模型
library(pROC)
logiroc <- roc(logi_x$group,logi_x$score) 

plot.roc(logiroc,print.thres=T,print.auc =T,main="ROC曲线",col="#008600")
plot(logiroc,
     print.thres=TRUE,
     main="ROC曲线",
     col="#008600")
graph2ppt(file = "output/plots/LASSO专利", append = T)


