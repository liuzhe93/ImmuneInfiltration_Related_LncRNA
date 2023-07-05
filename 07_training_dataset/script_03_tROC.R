setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/07_training_dataset")

rm(list=ls())

library(timeROC)
library(survival)
library(survivalROC)

risk<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",
                 sep="\t",header=T)
dim(risk)
#[1] 562   8
head(risk)
time_roc_res <- timeROC(
  T = risk$futime,
  delta = risk$fustat,
  marker = risk$riskScore,
  cause = 1,
  weighting="marginal",
  times = c(3, 5, 7),
  #times = c(1,2,3,4,5,6,7,8,9,10),
  ROC = TRUE,
  iid = TRUE
)

time_roc_res$AUC
#      t=3       t=5       t=7 
#0.7001561 0.7018321 0.6507481

confint(time_roc_res, level = 0.95)$CI_AUC
#     2.5% 97.5%
#t=3 55.64 84.39
#t=5 57.01 83.36
#t=7 52.64 77.51

time_ROC_df <- data.frame(
  TP_3year = time_roc_res$TP[, 1],
  FP_3year = time_roc_res$FP[, 1],
  TP_5year = time_roc_res$TP[, 2],
  FP_5year = time_roc_res$FP[, 2],
  TP_7year = time_roc_res$TP[, 3],
  FP_7year = time_roc_res$FP[, 3]
)
library(ggplot2)
pdf("ROC_curve_year357_TCGA.pdf")
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_7year, y = TP_7year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 7 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()


