###############################microarry data organization############################################################
setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA")
rm(list=ls())
##############################lncRNA annnotation information###########################################################
# download from https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/HG-U133_Plus_2-na36-annot-csv.zip
anno_file<-read.csv("LncRNAannotation/NetAffx/HG-U133_Plus_2-na36-annot-csv/HG-U133_Plus_2.na36.annot.csv", skip = 25)
anno_file<-subset(anno_file, select = c("Probe.Set.ID", "Gene.Symbol"))
write.csv(anno_file, "Probe_anno.csv", row.names = F, quote = F)
#human reference genome 38 downloaded from https://www.gencodegenes.org/human/release_21.html
lncRNA<-read.table("LncRNAannotation/GENCODE/gencode.v21.long_noncoding_RNAs.gtf", header = F, sep = "\t")
colnames(lncRNA)<-c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
lncRNA_sel<-as.data.frame(lncRNA$attributes)
colnames(lncRNA_sel)<-"attributes"
library("dplyr")
library("tidyr")
lncRNA_sel<-lncRNA_sel %>%
  separate(col = attributes,into = c("gene_id","transcript_id", "gene_type", "gene_status", "gene_name", "transcript_type",
                                     "transcript_status", "transcript_name", "level", "havana_gene"), sep = ";")
table(lncRNA_sel$gene_type)
#gene_type 3prime_overlapping_ncrna                 gene_type antisense               gene_type known_ncrna 
#                               118                               43254                                   6 
#gene_type lincRNA                gene_type non_coding      gene_type processed_transcript 
#            55426                                   3                               12319 
#gene_type sense_intronic         gene_type sense_overlapping                       gene_type TEC 
#                    3664                                1336                                3214 
lncRNA_sel$gene_name_cor<-lncRNA_sel$gene_name
for(i in 1:dim(lncRNA_sel)[1]){
  lncRNA_sel$gene_name_cor[i]<-gsub(" gene_name ","",lncRNA_sel$gene_name[i])
}
head(lncRNA_sel$gene_name_cor)
genelist<-unique(lncRNA_sel$gene_name_cor)
write.csv(genelist, "LncRNA_genename.csv", quote = F, row.names = F)

## match probe-gene
## python script_02_format_probe_gene.py Probe_anno.csv Probe_anno_format.csv
## read files
anno_file<-read.csv("Probe_anno_format.csv", header = T )
genelist<-read.csv("LncRNA_genename.csv", header = T)
merged<-merge(anno_file, genelist, by.x = "gene", by.y = "x")
write.csv(merged, "LncRNA_probe.csv", quote = F, row.names = F)



