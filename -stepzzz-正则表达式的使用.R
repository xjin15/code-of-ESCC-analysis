#####正式进入正则表达式######
#元字符
#(1)表达自身
grep("ab",x)
#(2)任意字符：.
grep("a.b",x,value = TRUE)
#(3)表示数量:*表示0及以上，+表示1及以上，?表示0或1个，{1,3}表示1-3个字符,{2}表示仅2个字符。这些都要放在后面，例如".*",表示任意数量的任意字符
grep("a.*b",x,value = TRUE)
grep("a.+b",x,value = TRUE)
grep("a.?b",x,value = TRUE)
grep("a.{1,3}b",x,value = TRUE)
grep("a.{1}b",x,value = TRUE)
#(4)表示^开头和结尾$
grep("^1",x,value = TRUE)
grep("4$",x,value = TRUE)
#(5)表示包含任意以下字符[2c]:2或c
grep("a[2c]b",c("a2b","a1cb","ab","a111b","a1b","acb"),value = TRUE)
grep("a[2c]*b",c("a2b","a1cb","ab","a111b","a1b","acb"),value = TRUE)#几个都行
grep("a[0-9]*b",c("a2b","a1cb","ab","a111b","a1b","acb"),value = TRUE)#中间有1-9的数字
grep("a[a-z]b",c("a2b","a1cb","ab","a111b","a1b","acb"),value = TRUE)#中间有a-z的字幕
grep("a[^c]b",c("a2b","a1cb","ab","a111b","a1b","acb"),value = TRUE)#^表示否定，除了c以外的都行
#(6)表示字符组合
grep("(13|13ab).4",x,value = TRUE)
gsub("(13|13ab).4","\\1",x)

#转义字符
#(1)\\$,\\^，分别表示$,^本身
grep("a\\$b",c("acb","a?b","a??b","a$b"),value = TRUE)
#(2)\\表示反斜杠本身，\\?表示?本身
grep("a\\?b",c("acb","a?b","a??b","a$b"),value = TRUE)
grep("a\\?+b",c("acb","a?b","a??b","a$b"),value = TRUE)
#(3)其他
grep("^\\w+$",x,value = TRUE)#从头到尾都是文字/数字的字符，所以不能有空格
grep("\\W{1}",x,value = TRUE)#包含空格的字符
grep("\\d",x,value = TRUE)#包含数字的字符
grep("^\\d+$",x,value = TRUE)#从头到尾都是数字的字符
grep("^\\D+$",x,value = TRUE)#从头到尾没有数字的字符
grep("^\\D+$",grep("^\\w+$",x,value = TRUE),value = TRUE)#从头到尾只有文字没有数字和空格

#用过的案例
gsub("_.*$","",x)#切掉_以后的部分
#去除括号中的内容
gsub("\\(.*\\)","",colnames(Achilles))
#删去.（点）后的内容
genename<-gsub("\\..*","",genename)
#删去/后的内容
GPL570$gene<-gsub("/.*$","",GPL570$gene)#删去/后的内容
#去除引号
sb$target<-gsub("\"","",sb$target)#去除引号

regexpr#给出位置
str_extract_all