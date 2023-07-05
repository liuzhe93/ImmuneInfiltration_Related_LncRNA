setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset")
rm(list=ls()) 
library("survival")
library("survminer")
rt=read.csv("TCGA_exp_surv.csv")
rownames(rt)<-rt$id
rt$id<-NULL
rt$futime<-rt$futime/365
fit <- coxph(Surv(futime, fustat) ~ CYB561D2 + LINC00638 + DANCR, data = rt)
sum.surv <- summary(fit)
c_index <- sum.surv$concordance
c_index
#         C      se(C) 
#0.59206244 0.03194791 

table = data.frame(
  tcga = c("TCGA-COAD", 0.59206244),
  gse1 = c("GSE39582", 0.59942799)
)
print(table) # 查看 table 数据
pdf("C-index.pdf")
barplot(height = c(0.59206244, 0.59942799),  
        names.arg = c('TCGA', 'GSE39582'),  
        col = 'steelblue',  
        border = '#ffffff',   
        xlab = 'dataset', 
        ylab = 'C-index',  
        main = 'C-index',  
        horiz = FALSE,  
        ylim = c(0, 0.7), 
)
dev.off()


