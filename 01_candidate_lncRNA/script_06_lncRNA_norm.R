setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA")
rm(list=ls())

mydata<-read.csv("exp_with_celltype.csv", header = T, row.names = 1)
mydata<-mydata[-1,]
colnames(mydata)<-mydata[1,]
mydata<-mydata[-1,]
write.csv(mydata, "Matrix_rowLncRNA_colCellType.csv", quote = F)

mydata_num<-as.data.frame(lapply(mydata, as.numeric))
rownames(mydata_num)<-rownames(mydata)
#######1. cell type: Myeloid dendritic cells
data_exp<-subset(mydata_num, select = c("Myeloid.dendritic.cells_1", "Myeloid.dendritic.cells_2",
                                        "Myeloid.dendritic.cells_3", "Myeloid.dendritic.cells_4"))
data_exp$Myeloid_dendritic_cells<-rowMeans(data_exp)
data_Myeloid_dendritic_cells<-subset(data_exp, select = Myeloid_dendritic_cells)

#######2. cell type: CD8 T cell activated
data_exp<-subset(mydata_num, select = c("CD8.T.cell.activated_1", "CD8.T.cell.activated_2",
                                        "CD8.T.cell.activated_3", "CD8.T.cell.activated_4",
                                        "CD8.T.cell.activated_5", "CD8.T.cell.activated_6"))
data_exp$CD8_T_cell_activated<-rowMeans(data_exp)
data_CD8_T_cell_activated<-subset(data_exp, select = CD8_T_cell_activated)

#######3. cell type: CD8 T cell resting
data_exp<-subset(mydata_num, select = c("CD8.T.cell.resting_1", "CD8.T.cell.resting_2",
                                        "CD8.T.cell.resting_3", "CD8.T.cell.resting_4"))
data_exp$CD8_T_cell_resting<-rowMeans(data_exp)
data_CD8_T_cell_resting<-subset(data_exp, select = CD8_T_cell_resting)

#######4. cell type: Monocytes
data_exp<-subset(mydata_num, select = c("Monocytes_1", "Monocytes_2",
                                        "Monocytes_3", "Monocytes_4",
                                        "Monocytes_5", "Monocytes_6"))
data_exp$Monocytes<-rowMeans(data_exp)
data_Monocytes<-subset(data_exp, select = Monocytes)

#######5. cell type: T helper 17
data_exp<-subset(mydata_num, select = c("T.helper.17_1", "T.helper.17_2",
                                        "T.helper.17_3", "T.helper.17_4",
                                        "T.helper.17_5", "T.helper.17_6",
                                        "T.helper.17_7", "T.helper.17_8",
                                        "T.helper.17_9"))
data_exp$T_helper_17<-rowMeans(data_exp)
data_T_helper_17<-subset(data_exp, select = T_helper_17)

#######6. cell type: Dendritic cells resting
data_exp<-subset(mydata_num, select = c("Dendritic.cells.resting_1", "Dendritic.cells.resting_2",
                                        "Dendritic.cells.resting_3", "Dendritic.cells.resting_4",
                                        "Dendritic.cells.resting_5", "Dendritic.cells.resting_6"))
data_exp$Dendritic_cells_resting<-rowMeans(data_exp)
data_Dendritic_cells_resting<-subset(data_exp, select = Dendritic_cells_resting)

#######7. cell type: Dendritic cells activated
data_exp<-subset(mydata_num, select = c("Dendritic.cells.activated_1", "Dendritic.cells.activated_2",
                                        "Dendritic.cells.activated_3", "Dendritic.cells.activated_4"))
data_exp$Dendritic_cells_activated<-rowMeans(data_exp)
data_Dendritic_cells_activated<-subset(data_exp, select = Dendritic_cells_activated)

#######8. cell type: Immature dendritic cells
data_exp<-subset(mydata_num, select = c("Immature.dendritic.cells_1", "Immature.dendritic.cells_2",
                                        "Immature.dendritic.cells_3", "Immature.dendritic.cells_4",
                                        "Immature.dendritic.cells_5", "Immature.dendritic.cells_6"))
data_exp$Immature_dendritic_cells<-rowMeans(data_exp)
data_Immature_dendritic_cells<-subset(data_exp, select = Immature_dendritic_cells)

#######9. cell type: NK resting
data_exp<-subset(mydata_num, select = c("NK.resting_1"))
data_exp$NK_resting<-rowMeans(data_exp)
data_NK_resting<-subset(data_exp, select = NK_resting)

#######10. cell type: NK activated
data_exp<-subset(mydata_num, select = c("NK.activated_1", "NK.activated_2",
                                        "NK.activated_3", "NK.activated_4",
                                        "NK.activated_5", "NK.activated_6",
                                        "NK.activated_7", "NK.activated_8",
                                        "NK.activated_9",  "NK.activated_10",
                                        "NK.activated_11"))
data_exp$NK_activated<-rowMeans(data_exp)
data_NK_activated<-subset(data_exp, select = NK_activated)

#######11. cell type: T gamma delta
data_exp<-subset(mydata_num, select = c("T.gamma.delta_1", "T.gamma.delta_2",
                                        "T.gamma.delta_3", "T.gamma.delta_4",
                                        "T.gamma.delta_5", "T.gamma.delta_6",
                                        "T.gamma.delta_7", "T.gamma.delta_8",
                                        "T.gamma.delta_9", "T.gamma.delta_10"))
data_exp$T_gamma_delta<-rowMeans(data_exp)
data_T_gamma_delta<-subset(data_exp, select = T_gamma_delta)

#######12. cell type: Mast cells activated
data_exp<-subset(mydata_num, select = c("Mast.cells.activated_1", "Mast.cells.activated_2",
                                        "Mast.cells.activated_3", "Mast.cells.activated_4"))
data_exp$Mast_cells_activated<-rowMeans(data_exp)
data_Mast_cells_activated<-subset(data_exp, select = Mast_cells_activated)

#######13. cell type: B cell activated
data_exp<-subset(mydata_num, select = c("B.cell.activated_1", "B.cell.activated_2",
                                        "B.cell.activated_3", "B.cell.activated_4",
                                        "B.cell.activated_5", "B.cell.activated_6",
                                        "B.cell.activated_7", "B.cell.activated_8",
                                        "B.cell.activated_9"))
data_exp$B_cell_activated<-rowMeans(data_exp)
data_B_cell_activated<-subset(data_exp, select = B_cell_activated)

#######14. cell type: Eosinophils
data_exp<-subset(mydata_num, select = c("Eosinophils_1", "Eosinophils_2",
                                        "Eosinophils_3"))
data_exp$Eosinophils<-rowMeans(data_exp)
data_Eosinophils<-subset(data_exp, select = Eosinophils)

#######15. cell type: CD4 T cell resting
data_exp<-subset(mydata_num, select = c("CD4.T.cell.resting_1", "CD4.T.cell.resting_2",
                                        "CD4.T.cell.resting_3", "CD4.T.cell.resting_4"))
data_exp$CD4_T_cell_resting<-rowMeans(data_exp)
data_CD4_T_cell_resting<-subset(data_exp, select = CD4_T_cell_resting)

#######16. cell type: CD4 T cell activated
data_exp<-subset(mydata_num, select = c("CD4.T.cell.activated_1", "CD4.T.cell.activated_2",
                                        "CD4.T.cell.activated_3", "CD4.T.cell.activated_4",
                                        "CD4.T.cell.activated_5", "CD4.T.cell.activated_6",
                                        "CD4.T.cell.activated_7"))
data_exp$CD4_T_cell_activated<-rowMeans(data_exp)
data_CD4_T_cell_activated<-subset(data_exp, select = CD4_T_cell_activated)

#######17. cell type: NKT activated
data_exp<-subset(mydata_num, select = c("NKT.activated_1", "NKT.activated_2",
                                        "NKT.activated_3", "NKT.activated_4",
                                        "NKT.activated_5", "NKT.activated_6"))
data_exp$NKT_activated<-rowMeans(data_exp)
data_NKT_activated<-subset(data_exp, select = NKT_activated)

#######18. cell type: Plasmacytoid dendritic cells
data_exp<-subset(mydata_num, select = c("Plasmacytoid.dendritic.cells_1", "Plasmacytoid.dendritic.cells_2",
                                        "Plasmacytoid.dendritic.cells_3", "Plasmacytoid.dendritic.cells_4",
                                        "Plasmacytoid.dendritic.cells_5", "Plasmacytoid.dendritic.cells_6",
                                        "Plasmacytoid.dendritic.cells_7", "Plasmacytoid.dendritic.cells_8"))
data_exp$Plasmacytoid_dendritic_cells<-rowMeans(data_exp)
data_Plasmacytoid_dendritic_cells<-subset(data_exp, select = Plasmacytoid_dendritic_cells)

#######19. cell type: Neutrophils
data_exp<-subset(mydata_num, select = c("Neutrophils_1", "Neutrophils_2",
                                        "Neutrophils_3", "Neutrophils_4",
                                        "Neutrophils_5", "Neutrophils_6",
                                        "Neutrophils_7"))
data_exp$Neutrophils<-rowMeans(data_exp)
data_Neutrophils<-subset(data_exp, select = Neutrophils)


merged_data<-cbind(data_Myeloid_dendritic_cells, data_CD8_T_cell_activated,
                   data_CD8_T_cell_resting, data_Monocytes,
                   data_T_helper_17, data_Dendritic_cells_resting,
                   data_Dendritic_cells_activated, data_Immature_dendritic_cells,
                   data_NK_resting, data_NK_activated,
                   data_T_gamma_delta, data_Mast_cells_activated,
                   data_B_cell_activated, data_Eosinophils,
                   data_CD4_T_cell_resting, data_CD4_T_cell_activated,
                   data_NKT_activated, data_Plasmacytoid_dendritic_cells,
                   data_Neutrophils)
rmax <- vector()
rmin <- vector()
for (i in 1:nrow(merged_data)){
  rmax = c(rmax,max(merged_data[i,]))
  rmin = c(rmin,min(merged_data[i,]))
}
temp <- cbind(merged_data,rmax,rmin)
data_norm<-merged_data
for(i in 1:nrow(data_norm)){
  for(j in 1:ncol(data_norm)){
    data_norm[i,j] = (data_norm[i,j]-temp[i,21])/(temp[i,20]-temp[i,21])
  }
}
write.csv(data_norm, "data_norm.csv", quote = F)



#1. Myeloid dendritic cells_1
#2. CD8 T cell activated_1
#3. CD8 T cell resting_1
#4. Monocytes_1
#5. T helper 17_1
#6. Dendritic cells resting_1
#7. Dendritic cells activated_1
#8. Immature dendritic cells_1
#9. NK resting_1
#10. NK activated_1
#11. T gamma delta_1
#12. Mast cells activated_1
#13. B cell activated_1
#14. Eosinophils
#15. CD4 T cell resting_1
#16. CD4 T cell activated_1
#17. NKT activated_1
#18. Plasmacytoid dendritic cells_1
#19. Neutrophils_1








