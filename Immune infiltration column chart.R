library(tidyr)
library(reshape2)

source('Cibersort.R')
result1 <- CIBERSORT('LM22.txt','DATA.txt', perm = 1000, QN = T)  #perm置换次数=0，QN分位数归一化=TRUE

a <- as.data.frame(result1)
a<-a[,c(1:22)]
head(a)
a$X<-rownames(a)
fen<-read.csv("fen.csv")
a<-merge(fen,a,"X")

mydata1<-melt(
  a,
  id.vars=c("X","fen"),
  variable.name="immunecell",
  value.name="tpm"
)

library(ggpubr)
library(ggplot2)

ylabname <- paste("immunecell", "expression")
colnames(mydata1) <- c("Sample", "Groups", "immunecell","tpm")
# 计算p value
pvalues <- sapply(mydata1$immunecell, function(x) {
  res <- wilcox.test(as.numeric(tpm) ~ Groups, data = subset(mydata1, immunecell == x)) #两组，wilcox.test或t.test；多组，kruskal.test或aov(one-way ANOVA test)
  res$p.value
})
pv <- data.frame(gene = mydata1$immunecell, pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                  labels=c('****','***', '**', '*', 'ns'))
mydata1<-mydata1[,-1]
# 画box plot
p.box <- ggplot(mydata1, aes(x=immunecell, y=tpm, color=Groups, fill=Groups)) +
  geom_boxplot(alpha = .5) + #半透明
  theme_classic() + #或theme_bw()
  scale_fill_brewer(palette = "Set1") + #按类填充颜色
  scale_color_brewer(palette = "Set1") + #按类给边框着色
  
  theme(axis.text.x = element_text(colour="black", size = 15,
                                   #名太挤，旋转45度
                                   angle = 45, hjust =1, vjust =1)) +
  geom_text(aes(x=gene, y=max(mydata1$tpm) * 1.1,
                label = pv$sigcode),
            data=pv, 
            inherit.aes=F) +
  ylab(ylabname)
p.box
ggsave("immunecellbox.pdf", width = 14, height = 5)
# 画带散点的box plot
p.box.dot <- p.box + geom_point(shape = 21, size=.5, # 点的形状和大小
                                position = position_jitterdodge(), # 让点散开
                                alpha = .5) #半透明
p.box.dot
ggsave("immunecellsanbox.pdf", width = 14, height = 5)


library(IOBR)

data<-read.table("DATA.txt",row.names = 1,header=T,sep = '\t')

im_mcpcounter <- deconvo_tme(eset = data,
                             method = "mcpcounter"
)

