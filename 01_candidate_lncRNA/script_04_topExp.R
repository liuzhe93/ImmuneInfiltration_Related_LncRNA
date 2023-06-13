setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA")
rm(list=ls())

mydata<-read.csv("exp_with_celltype.csv", header = T, row.names = 1)
mydata<-mydata[-1,]
colnames(mydata)<-mydata[1,]
mydata<-mydata[-1,]
write.csv(mydata, "Matrix_rowLncRNA_colCellType.csv", quote = F)

mydata_num<-as.data.frame(lapply(mydata, as.numeric))
rownames(mydata_num)<-rownames(mydata)
#n_candi: the number of candidate immune-related lncRNAs
n_candi<-round(dim(data_Myeloid_dendritic_cells_sort)[1]*0.3)
n_candi
#[1] 427

#######1. cell type: Myeloid dendritic cells
data_exp<-subset(mydata_num, select = c("Myeloid.dendritic.cells_1", "Myeloid.dendritic.cells_2",
                                        "Myeloid.dendritic.cells_3", "Myeloid.dendritic.cells_4"))
data_exp$Myeloid_dendritic_cells<-rowMeans(data_exp)
data_Myeloid_dendritic_cells<-subset(data_exp, select = Myeloid_dendritic_cells)
data_Myeloid_dendritic_cells_sort<-data_Myeloid_dendritic_cells
data_Myeloid_dendritic_cells_sort$LncRNA<-row.names(data_Myeloid_dendritic_cells_sort)
data_Myeloid_dendritic_cells_sort<-data_Myeloid_dendritic_cells_sort[order(-data_Myeloid_dendritic_cells_sort$Myeloid_dendritic_cells),]
dim(data_Myeloid_dendritic_cells_sort)
#[1] 1422    2
top30per_Myeloid_dendritic_cells<-rownames(data_Myeloid_dendritic_cells_sort)[1:n_candi]
  
#######2. cell type: CD8 T cell activated
data_exp<-subset(mydata_num, select = c("CD8.T.cell.activated_1", "CD8.T.cell.activated_2",
                                        "CD8.T.cell.activated_3", "CD8.T.cell.activated_4",
                                        "CD8.T.cell.activated_5", "CD8.T.cell.activated_6"))
data_exp$CD8_T_cell_activated<-rowMeans(data_exp)
data_CD8_T_cell_activated<-subset(data_exp, select = CD8_T_cell_activated)
data_CD8_T_cell_activated_sort<-data_CD8_T_cell_activated
data_CD8_T_cell_activated_sort$LncRNA<-row.names(data_CD8_T_cell_activated_sort)
data_CD8_T_cell_activated_sort<-data_CD8_T_cell_activated_sort[order(-data_CD8_T_cell_activated_sort$CD8_T_cell_activated),]
dim(data_CD8_T_cell_activated_sort)
#[1] 1422    2
top30per_CD8_T_cell_activated<-rownames(data_CD8_T_cell_activated_sort)[1:n_candi]

#######3. cell type: CD8 T cell resting
data_exp<-subset(mydata_num, select = c("CD8.T.cell.resting_1", "CD8.T.cell.resting_2",
                                        "CD8.T.cell.resting_3", "CD8.T.cell.resting_4"))
data_exp$CD8_T_cell_resting<-rowMeans(data_exp)
data_CD8_T_cell_resting<-subset(data_exp, select = CD8_T_cell_resting)
data_CD8_T_cell_resting_sort<-data_CD8_T_cell_resting
data_CD8_T_cell_resting_sort$LncRNA<-row.names(data_CD8_T_cell_resting_sort)
data_CD8_T_cell_resting_sort<-data_CD8_T_cell_resting_sort[order(-data_CD8_T_cell_resting_sort$CD8_T_cell_resting),]
dim(data_CD8_T_cell_resting_sort)
#[1] 1422    2
top30per_CD8_T_cell_resting<-rownames(data_CD8_T_cell_resting_sort)[1:n_candi]

#######4. cell type: Monocytes
data_exp<-subset(mydata_num, select = c("Monocytes_1", "Monocytes_2",
                                        "Monocytes_3", "Monocytes_4",
                                        "Monocytes_5", "Monocytes_6"))
data_exp$Monocytes<-rowMeans(data_exp)
data_Monocytes<-subset(data_exp, select = Monocytes)
data_Monocytes_sort<-data_Monocytes
data_Monocytes_sort$LncRNA<-row.names(data_Monocytes_sort)
data_Monocytes_sort<-data_Monocytes_sort[order(-data_Monocytes_sort$Monocytes),]
dim(data_Monocytes_sort)
#[1] 1422    2
top30per_Monocytes<-rownames(data_Monocytes_sort)[1:n_candi]

#######5. cell type: T helper 17
data_exp<-subset(mydata_num, select = c("T.helper.17_1", "T.helper.17_2",
                                        "T.helper.17_3", "T.helper.17_4",
                                        "T.helper.17_5", "T.helper.17_6",
                                        "T.helper.17_7", "T.helper.17_8",
                                        "T.helper.17_9"))
data_exp$T_helper_17<-rowMeans(data_exp)
data_T_helper_17<-subset(data_exp, select = T_helper_17)
data_T_helper_17_sort<-data_T_helper_17
data_T_helper_17_sort$LncRNA<-row.names(data_T_helper_17_sort)
data_T_helper_17_sort<-data_T_helper_17_sort[order(-data_T_helper_17_sort$T_helper_17),]
dim(data_T_helper_17_sort)
#[1] 1422    2
top30per_T_helper_17<-rownames(data_T_helper_17_sort)[1:n_candi]

#######6. cell type: Dendritic cells resting
data_exp<-subset(mydata_num, select = c("Dendritic.cells.resting_1", "Dendritic.cells.resting_2",
                                        "Dendritic.cells.resting_3", "Dendritic.cells.resting_4",
                                        "Dendritic.cells.resting_5", "Dendritic.cells.resting_6"))
data_exp$Dendritic_cells_resting<-rowMeans(data_exp)
data_Dendritic_cells_resting<-subset(data_exp, select = Dendritic_cells_resting)
data_Dendritic_cells_resting_sort<-data_Dendritic_cells_resting
data_Dendritic_cells_resting_sort$LncRNA<-row.names(data_Dendritic_cells_resting_sort)
data_Dendritic_cells_resting_sort<-data_Dendritic_cells_resting_sort[order(-data_Dendritic_cells_resting_sort$Dendritic_cells_resting),]
dim(data_Dendritic_cells_resting_sort)
#[1] 1422    2
top30per_Dendritic_cells_resting<-rownames(data_Dendritic_cells_resting_sort)[1:n_candi]

#######7. cell type: Dendritic cells activated
data_exp<-subset(mydata_num, select = c("Dendritic.cells.activated_1", "Dendritic.cells.activated_2",
                                        "Dendritic.cells.activated_3", "Dendritic.cells.activated_4"))
data_exp$Dendritic_cells_activated<-rowMeans(data_exp)
data_Dendritic_cells_activated<-subset(data_exp, select = Dendritic_cells_activated)
data_Dendritic_cells_activated_sort<-data_Dendritic_cells_activated
data_Dendritic_cells_activated_sort$LncRNA<-row.names(data_Dendritic_cells_activated_sort)
data_Dendritic_cells_activated_sort<-data_Dendritic_cells_activated_sort[order(-data_Dendritic_cells_activated_sort$Dendritic_cells_activated),]
dim(data_Dendritic_cells_activated_sort)
#[1] 1422    2
top30per_Dendritic_cells_activated<-rownames(data_Dendritic_cells_activated_sort)[1:n_candi]

#######8. cell type: Immature dendritic cells
data_exp<-subset(mydata_num, select = c("Immature.dendritic.cells_1", "Immature.dendritic.cells_2",
                                        "Immature.dendritic.cells_3", "Immature.dendritic.cells_4",
                                        "Immature.dendritic.cells_5", "Immature.dendritic.cells_6"))
data_exp$Immature_dendritic_cells<-rowMeans(data_exp)
data_Immature_dendritic_cells<-subset(data_exp, select = Immature_dendritic_cells)
data_Immature_dendritic_cells_sort<-data_Immature_dendritic_cells
data_Immature_dendritic_cells_sort$LncRNA<-row.names(data_Immature_dendritic_cells_sort)
data_Immature_dendritic_cells_sort<-data_Immature_dendritic_cells_sort[order(-data_Immature_dendritic_cells_sort$Immature_dendritic_cells),]
dim(data_Immature_dendritic_cells_sort)
#[1] 1422    2
top30per_Immature_dendritic_cells<-rownames(data_Immature_dendritic_cells_sort)[1:n_candi]

#######9. cell type: NK resting
data_exp<-subset(mydata_num, select = c("NK.resting_1"))
data_exp$NK_resting<-rowMeans(data_exp)
data_NK_resting<-subset(data_exp, select = NK_resting)
data_NK_resting_sort<-data_NK_resting
data_NK_resting_sort$LncRNA<-row.names(data_NK_resting_sort)
data_NK_resting_sort<-data_NK_resting_sort[order(-data_NK_resting_sort$NK_resting),]
dim(data_NK_resting_sort)
#[1] 1422    2
top30per_NK_resting<-rownames(data_NK_resting_sort)[1:n_candi]

#######10. cell type: NK activated
data_exp<-subset(mydata_num, select = c("NK.activated_1", "NK.activated_2",
                                        "NK.activated_3", "NK.activated_4",
                                        "NK.activated_5", "NK.activated_6",
                                        "NK.activated_7", "NK.activated_8",
                                        "NK.activated_9",  "NK.activated_10",
                                        "NK.activated_11"))
data_exp$NK_activated<-rowMeans(data_exp)
data_NK_activated<-subset(data_exp, select = NK_activated)
data_NK_activated_sort<-data_NK_activated
data_NK_activated_sort$LncRNA<-row.names(data_NK_activated_sort)
data_NK_activated_sort<-data_NK_activated_sort[order(-data_NK_activated_sort$NK_activated),]
dim(data_NK_activated_sort)
#[1] 1422    2
top30per_NK_activated<-rownames(data_NK_activated_sort)[1:n_candi]

#######11. cell type: T gamma delta
data_exp<-subset(mydata_num, select = c("T.gamma.delta_1", "T.gamma.delta_2",
                                        "T.gamma.delta_3", "T.gamma.delta_4",
                                        "T.gamma.delta_5", "T.gamma.delta_6",
                                        "T.gamma.delta_7", "T.gamma.delta_8",
                                        "T.gamma.delta_9", "T.gamma.delta_10"))
data_exp$T_gamma_delta<-rowMeans(data_exp)
data_T_gamma_delta<-subset(data_exp, select = T_gamma_delta)
data_T_gamma_delta_sort<-data_T_gamma_delta
data_T_gamma_delta_sort$LncRNA<-row.names(data_T_gamma_delta_sort)
data_T_gamma_delta_sort<-data_T_gamma_delta_sort[order(-data_T_gamma_delta_sort$T_gamma_delta),]
dim(data_T_gamma_delta_sort)
#[1] 1422    2
top30per_T_gamma_delta<-rownames(data_T_gamma_delta_sort)[1:n_candi]

#######12. cell type: Mast cells activated
data_exp<-subset(mydata_num, select = c("Mast.cells.activated_1", "Mast.cells.activated_2",
                                        "Mast.cells.activated_3", "Mast.cells.activated_4"))
data_exp$Mast_cells_activated<-rowMeans(data_exp)
data_Mast_cells_activated<-subset(data_exp, select = Mast_cells_activated)
data_Mast_cells_activated_sort<-data_Mast_cells_activated
data_Mast_cells_activated_sort$LncRNA<-row.names(data_Mast_cells_activated_sort)
data_Mast_cells_activated_sort<-data_Mast_cells_activated_sort[order(-data_Mast_cells_activated_sort$Mast_cells_activated),]
dim(data_Mast_cells_activated_sort)
#[1] 1422    2
top30per_Mast_cells_activated<-rownames(data_Mast_cells_activated_sort)[1:n_candi]

#######13. cell type: B cell activated
data_exp<-subset(mydata_num, select = c("B.cell.activated_1", "B.cell.activated_2",
                                        "B.cell.activated_3", "B.cell.activated_4",
                                        "B.cell.activated_5", "B.cell.activated_6",
                                        "B.cell.activated_7", "B.cell.activated_8",
                                        "B.cell.activated_9"))
data_exp$B_cell_activated<-rowMeans(data_exp)
data_B_cell_activated<-subset(data_exp, select = B_cell_activated)
data_B_cell_activated_sort<-data_B_cell_activated
data_B_cell_activated_sort$LncRNA<-row.names(data_B_cell_activated_sort)
data_B_cell_activated_sort<-data_B_cell_activated_sort[order(-data_B_cell_activated_sort$B_cell_activated),]
dim(data_B_cell_activated_sort)
#[1] 1422    2
top30per_B_cell_activated<-rownames(data_B_cell_activated_sort)[1:n_candi]

#######14. cell type: Eosinophils
data_exp<-subset(mydata_num, select = c("Eosinophils_1", "Eosinophils_2",
                                        "Eosinophils_3"))
data_exp$Eosinophils<-rowMeans(data_exp)
data_Eosinophils<-subset(data_exp, select = Eosinophils)
data_Eosinophils_sort<-data_Eosinophils
data_Eosinophils_sort$LncRNA<-row.names(data_Eosinophils_sort)
data_Eosinophils_sort<-data_Eosinophils_sort[order(-data_Eosinophils_sort$Eosinophils),]
dim(data_Eosinophils_sort)
#[1] 1422    2
top30per_Eosinophils<-rownames(data_Eosinophils_sort)[1:n_candi]

#######15. cell type: CD4 T cell resting
data_exp<-subset(mydata_num, select = c("CD4.T.cell.resting_1", "CD4.T.cell.resting_2",
                                        "CD4.T.cell.resting_3", "CD4.T.cell.resting_4"))
data_exp$CD4_T_cell_resting<-rowMeans(data_exp)
data_CD4_T_cell_resting<-subset(data_exp, select = CD4_T_cell_resting)
data_CD4_T_cell_resting_sort<-data_CD4_T_cell_resting
data_CD4_T_cell_resting_sort$LncRNA<-row.names(data_CD4_T_cell_resting_sort)
data_CD4_T_cell_resting_sort<-data_CD4_T_cell_resting_sort[order(-data_CD4_T_cell_resting_sort$CD4_T_cell_resting),]
dim(data_CD4_T_cell_resting_sort)
#[1] 1422    2
top30per_CD4_T_cell_resting<-rownames(data_CD4_T_cell_resting_sort)[1:n_candi]

#######16. cell type: CD4 T cell activated
data_exp<-subset(mydata_num, select = c("CD4.T.cell.activated_1", "CD4.T.cell.activated_2",
                                        "CD4.T.cell.activated_3", "CD4.T.cell.activated_4",
                                        "CD4.T.cell.activated_5", "CD4.T.cell.activated_6",
                                        "CD4.T.cell.activated_7"))
data_exp$CD4_T_cell_activated<-rowMeans(data_exp)
data_CD4_T_cell_activated<-subset(data_exp, select = CD4_T_cell_activated)
data_CD4_T_cell_activated_sort<-data_CD4_T_cell_activated
data_CD4_T_cell_activated_sort$LncRNA<-row.names(data_CD4_T_cell_activated_sort)
data_CD4_T_cell_activated_sort<-data_CD4_T_cell_activated_sort[order(-data_CD4_T_cell_activated_sort$CD4_T_cell_activated),]
dim(data_CD4_T_cell_activated_sort)
#[1] 1422    2
top30per_CD4_T_cell_activated<-rownames(data_CD4_T_cell_activated_sort)[1:n_candi]

#######17. cell type: NKT activated
data_exp<-subset(mydata_num, select = c("NKT.activated_1", "NKT.activated_2",
                                        "NKT.activated_3", "NKT.activated_4",
                                        "NKT.activated_5", "NKT.activated_6"))
data_exp$NKT_activated<-rowMeans(data_exp)
data_NKT_activated<-subset(data_exp, select = NKT_activated)
data_NKT_activated_sort<-data_NKT_activated
data_NKT_activated_sort$LncRNA<-row.names(data_NKT_activated_sort)
data_NKT_activated_sort<-data_NKT_activated_sort[order(-data_NKT_activated_sort$NKT_activated),]
dim(data_NKT_activated_sort)
#[1] 1422    2
top30per_NKT_activated<-rownames(data_NKT_activated_sort)[1:n_candi]

#######18. cell type: Plasmacytoid dendritic cells
data_exp<-subset(mydata_num, select = c("Plasmacytoid.dendritic.cells_1", "Plasmacytoid.dendritic.cells_2",
                                        "Plasmacytoid.dendritic.cells_3", "Plasmacytoid.dendritic.cells_4",
                                        "Plasmacytoid.dendritic.cells_5", "Plasmacytoid.dendritic.cells_6",
                                        "Plasmacytoid.dendritic.cells_7", "Plasmacytoid.dendritic.cells_8"))
data_exp$Plasmacytoid_dendritic_cells<-rowMeans(data_exp)
data_Plasmacytoid_dendritic_cells<-subset(data_exp, select = Plasmacytoid_dendritic_cells)
data_Plasmacytoid_dendritic_cells_sort<-data_Plasmacytoid_dendritic_cells
data_Plasmacytoid_dendritic_cells_sort$LncRNA<-row.names(data_Plasmacytoid_dendritic_cells_sort)
data_Plasmacytoid_dendritic_cells_sort<-data_Plasmacytoid_dendritic_cells_sort[order(-data_Plasmacytoid_dendritic_cells_sort$Plasmacytoid_dendritic_cells),]
dim(data_Plasmacytoid_dendritic_cells_sort)
#[1] 1422    2
top30per_Plasmacytoid_dendritic_cells<-rownames(data_Plasmacytoid_dendritic_cells_sort)[1:n_candi]

#######19. cell type: Neutrophils
data_exp<-subset(mydata_num, select = c("Neutrophils_1", "Neutrophils_2",
                                        "Neutrophils_3", "Neutrophils_4",
                                        "Neutrophils_5", "Neutrophils_6",
                                        "Neutrophils_7"))
data_exp$Neutrophils<-rowMeans(data_exp)
data_Neutrophils<-subset(data_exp, select = Neutrophils)
data_Neutrophils_sort<-data_Neutrophils
data_Neutrophils_sort$LncRNA<-row.names(data_Neutrophils_sort)
data_Neutrophils_sort<-data_Neutrophils_sort[order(-data_Neutrophils_sort$Neutrophils),]
dim(data_Neutrophils_sort)
#[1] 1422    2
top30per_Neutrophils<-rownames(data_Neutrophils_sort)[1:n_candi]

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
write.csv(merged_data, "ExpPro_rowGene_colImmuneCell.csv", quote = F)

top30_merge_uniq<-unique(c(top30per_Myeloid_dendritic_cells, top30per_CD8_T_cell_activated, 
                         top30per_CD8_T_cell_resting, top30per_Monocytes,
                         top30per_T_helper_17, top30per_Dendritic_cells_resting,
                         top30per_Dendritic_cells_activated, top30per_Immature_dendritic_cells,
                         top30per_NK_resting, top30per_NK_activated,
                         top30per_T_gamma_delta, top30per_Mast_cells_activated,
                         top30per_B_cell_activated, top30per_Eosinophils,
                         top30per_CD4_T_cell_resting, top30per_CD4_T_cell_activated,
                         top30per_NKT_activated, top30per_Plasmacytoid_dendritic_cells,
                         top30per_Neutrophils))
top30_merge_uniq<-as.data.frame(top30_merge_uniq)
write.csv(top30_merge_uniq, "top30_expressed_lncRNA_in_each_ImmuCell.csv", quote = F)

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








