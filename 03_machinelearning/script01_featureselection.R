setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning")
rm(list=ls())

###########lncRNA expression level information
allfeatures<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/02_upregulated_in_immunecells/Upreg_in_immune.csv", header = T)
allfeatures<-subset(allfeatures, select = "Gene.Name")
dim(allfeatures)
#[1] 87   1
genelist<-allfeatures$Gene.Name

cancer_exp<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/02_upregulated_in_immunecells/cancer_exp.csv", row.names = 1)
cancer_exp$genename<-rownames(cancer_exp)
data_can<-cancer_exp[cancer_exp$genename %in% genelist,]
data_can_fil<-na.omit(data_can)
data_can_fil$genename<-NULL
mydata<-t(data_can_fil)
mydata<-as.data.frame(mydata)
mydata$ID<-row.names(mydata)
for(i in 1:dim(mydata)[1]){
  temp<-strsplit(mydata$ID[i], "_")
  id<-temp[[1]][1]
  mydata$ID[i]<-id
}
dim(mydata)
#[1] 566  72
mydata<-mydata[,c(72,1:71)]
rownames(mydata)<-NULL

###########clinical factors information
sampleinfo<-read.table("SampleInfor.txt", header = T, sep = "\t")
sam_inf<-t(sampleinfo)
colnames(sam_inf)<-sam_inf[1,]
sam_inf<-sam_inf[-1,]
#RFS（无复发生存）、OS（总生存期）
colnames(sam_inf)<-c("Sample_geo_accession","Sample_status","Sample_submission_date","Sample_last_update_date",
                     "Sample_type","Sample_channel_count","Sample_source_name_ch1","Sample_organism_ch1",
                     "organism","dataset","gender","age","tnm.stage","tnm.t","tnm.n","tnm.m","tumor.location",
                     "chemotherapy.adjuvant","chemotherapy.adjuvant.type","rfs.event","rfs.delay","os.event","os.delay (months)",
                     "mmr.status","cimp.status","cin.status","tp53.mutation","tp53.mutation.dna","tp53.mutation.exon.number",
                     "tp53.mutation.protein","kras.mutation","kras.mutation.dna","kras.mutation.exon.number","kras.mutation.protein",
                     "braf.mutation","braf.mutation.dna","braf.mutation.exon.number","braf.mutation.protein","cit.molecularsubtype",
                     "dependancy sample")
sam_inf<-as.data.frame(sam_inf)
info_sel<-subset(sam_inf, select = c("Sample_geo_accession","gender", "age", "tnm.stage", "tnm.t","tnm.n","tnm.m","os.event","os.delay (months)"))
info_sel$gender<-gsub("Sex: ","",info_sel$gender)
info_sel$age<-gsub("^age.* ","",info_sel$age)
info_sel$tnm.stage<-gsub("^tnm.* ","",info_sel$tnm.stage)
info_sel$tnm.t<-gsub("^tnm.* ","",info_sel$tnm.t)
info_sel$tnm.n<-gsub("^tnm.* ","",info_sel$tnm.n)
info_sel$tnm.m<-gsub("^tnm.* ","",info_sel$tnm.m)
info_sel$os.event<-gsub("^os.* ","",info_sel$os.event)
info_sel$`os.delay (months)`<-gsub("^os.* ","",info_sel$`os.delay (months)`)


data_merged<-merge(info_sel, mydata, by.x = "Sample_geo_accession", by.y = "ID")
write.csv(data_merged, "exp_cli.csv", quote = F, row.names = F)

input<-data_merged[,c(1,9,8,10:80)]
colnames(input)[1:3]<-c("id",	"futime",	"fustat")
write.table(input, "uniSigExp.txt", quote = F, row.names = F, sep = "\t")
# two machine-learning based algorithms were used for feature selection





