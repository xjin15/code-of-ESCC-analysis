library(rtracklayer)
gtf_data = import('..//gencode.v22.annotation.gtf.gz') #gtf的路径
#这里使用import导入gtf文件， 生成一个GRangs对象
options(stringsAsFactors = F)
gtf_data = as.data.frame(gtf_data)
gtf_data 
save(gtf_data,file = "outdata/gtf注释文件.Rdata")