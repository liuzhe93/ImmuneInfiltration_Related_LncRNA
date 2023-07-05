setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA")
rm(list=ls())

library("pheatmap")
hklncRNA<-read.csv("hklncRNA_inf.csv")
hklncRNA_list<-unique(hklncRNA$X)
speclncRNA<-read.csv("speclncRNA_inf.csv")
speclncRNA_list<-unique(speclncRNA$X)

###############housekeeping and celltype-specific lncRNAs' expression levels in 19 types of immune cells
exp<-read.csv("Matrix_rowLncRNA_colCellType.csv", row.names = 1)
exp_hklncRNA<-exp[hklncRNA_list,]
exp_spelncRNA<-exp[speclncRNA_list,]
mydata<-rbind(exp_hklncRNA, exp_spelncRNA)
data_sort<-subset(mydata, select = c("Myeloid.dendritic.cells_1","Myeloid.dendritic.cells_2","Myeloid.dendritic.cells_3",
                                     "Myeloid.dendritic.cells_4","CD8.T.cell.activated_1","CD8.T.cell.activated_2",
                                     "CD8.T.cell.activated_3","CD8.T.cell.activated_4","CD8.T.cell.activated_5",
                                     "CD8.T.cell.activated_6","CD8.T.cell.resting_1","CD8.T.cell.resting_2",
                                     "CD8.T.cell.resting_3","CD8.T.cell.resting_4","Monocytes_1","Monocytes_2","Monocytes_3",
                                     "Monocytes_4","Monocytes_5","Monocytes_6","T.helper.17_1","T.helper.17_2","T.helper.17_3",
                                     "T.helper.17_4","T.helper.17_5","T.helper.17_6","T.helper.17_7","T.helper.17_8",
                                     "T.helper.17_9","Dendritic.cells.resting_1","Dendritic.cells.resting_2",
                                     "Dendritic.cells.resting_3","Dendritic.cells.resting_4","Dendritic.cells.resting_5",
                                     "Dendritic.cells.resting_6","Dendritic.cells.activated_1","Dendritic.cells.activated_2",
                                     "Dendritic.cells.activated_3",'Dendritic.cells.activated_4',"Immature.dendritic.cells_1",
                                     "Immature.dendritic.cells_2","Immature.dendritic.cells_3","Immature.dendritic.cells_4",
                                     "Immature.dendritic.cells_5","Immature.dendritic.cells_6","NK.resting_1","NK.activated_1",
                                     "NK.activated_2","NK.activated_3","NK.activated_4","NK.activated_5","NK.activated_6",
                                     "NK.activated_7","NK.activated_8","NK.activated_9","NK.activated_10","NK.activated_11",
                                     "T.gamma.delta_1","T.gamma.delta_2","T.gamma.delta_3","T.gamma.delta_4","T.gamma.delta_5",
                                     "T.gamma.delta_6","T.gamma.delta_7","T.gamma.delta_8","T.gamma.delta_9","T.gamma.delta_10",
                                     "Mast.cells.activated_1","Mast.cells.activated_2","Mast.cells.activated_3",
                                     "Mast.cells.activated_4","B.cell.activated_1","B.cell.activated_2","B.cell.activated_3",
                                     "B.cell.activated_4","B.cell.activated_5","B.cell.activated_6","B.cell.activated_7",
                                     "B.cell.activated_8","B.cell.activated_9","Eosinophils_1","Eosinophils_2","Eosinophils_3",
                                     "CD4.T.cell.resting_1","CD4.T.cell.resting_2","CD4.T.cell.resting_3","CD4.T.cell.resting_4",
                                     "CD4.T.cell.activated_1","CD4.T.cell.activated_2","CD4.T.cell.activated_3",
                                     "CD4.T.cell.activated_4","CD4.T.cell.activated_5","CD4.T.cell.activated_6",
                                     "CD4.T.cell.activated_7","NKT.activated_1","NKT.activated_2","NKT.activated_3",
                                     "NKT.activated_4","NKT.activated_5","NKT.activated_6","Plasmacytoid.dendritic.cells_1",
                                     "Plasmacytoid.dendritic.cells_2","Plasmacytoid.dendritic.cells_3",
                                     "Plasmacytoid.dendritic.cells_4","Plasmacytoid.dendritic.cells_5",
                                     "Plasmacytoid.dendritic.cells_6","Plasmacytoid.dendritic.cells_7",
                                     "Plasmacytoid.dendritic.cells_8","Neutrophils_1","Neutrophils_2","Neutrophils_3",
                                     "Neutrophils_4","Neutrophils_5","Neutrophils_6","Neutrophils_7"))
colnam<-colnames(data_sort)
colnam_cell<-colnam
for(i in 1:length(colnam)){
  colnam_cell[i]<-strsplit(colnam[i],"_")[[1]][1]
}
annotation_col = data.frame(Type=colnam_cell)
annotation_col$Type<-as.factor(annotation_col$Type)
rownames(annotation_col) = colnam
rownam<-c(rep("house_keeping_lncRNAs",221),rep("cell_type_specificity_lncRNAs",221))
annotation_row = matrix(rownam,nrow = 442,ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(data_sort)
colnames(annotation_row)<-"gene_type"
annotation_row<-as.data.frame(annotation_row)
ann_colors = list(Type = c(Myeloid.dendritic.cells = "#ed1299", CD8.T.cell.activated = "#09f9f5", CD8.T.cell.resting = "#246b93",
                           Monocytes = "#cc8e12", T.helper.17 = "#d561dd", Dendritic.cells.resting = "#c93f00",
                           Dendritic.cells.activated = "#ddd53e", Immature.dendritic.cells = "#4aef7b",
                           NK.resting =  "#e86502", NK.activated ="#9ed84e", T.gamma.delta = "#39ba30", 
                           Mast.cells.activated = "#6ad157", B.cell.activated = "#8249aa", Eosinophils = "#99db27",
                           CD4.T.cell.resting = "#ce2523", CD4.T.cell.activated = "#f7aa5d", NKT.activated = "#cebb10",
                           Plasmacytoid.dendritic.cells = "#03827f", Neutrophils = "#931635"),
                  gene_type = c(house_keeping_lncRNAs = "red", cell_type_specificity_lncRNAs  = "blue"))
pdf(file="heatmap_immune.pdf")
pheatmap(data_sort, clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 221, 
         angle_col = "45", annotation_colors = ann_colors,  
         show_rownames = F, show_colnames = F, main = "Title", fontsize = 6)
dev.off()




