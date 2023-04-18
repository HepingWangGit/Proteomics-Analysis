
setwd("/Users/hwang/Desktop/Somalogic Data/Somanic")

#Load Library
#Read Data
meta_data<-read.csv("output_file/Nash_EDTA_Metadata_clean_HW.csv")
File<-"output_file/somaLogic.EDTA.limmaNorm.txt"
data <- read.table(File)
data<-as.data.frame(t(data))
index<-match(rownames(data),meta_data$UniqueID)
meta_data<-meta_data[index,]
Annotation<-read.table("output_file/Annotation.txt",sep="\t",quote = "",comment.char = "",header = T)

#Select through top500 most variance data
vari_data<-apply(data,2,function(x){var(x)})
vari_data<-vari_data[order(vari_data,decreasing = T)]
seq_var<-names(vari_data[1:500])
index<-match(seq_var,colnames(data))
var_data<-data[,index]

#fdr of the coefficients of proteomics
D1ID_AXA<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="axa1125 24g bid"),]
D1ID_AXA_s<-D1ID_AXA$SubjectID
index1<-match(D1ID_AXA$UniqueID,rownames(var_data))
D1_data<-var_data[index1,]

PDFF_W16<-D1ID_AXA$classifier_pdff_W16

pro_data<-D1_data
df<-as.data.frame(cbind(pro_data,PDFF_W16))
coef_pvalue<-c()
coef<-c()
for(i in 1:ncol(pro_data)){
  glm.fit <- glm(factor(PDFF_W16) ~pro_data[,i] , data=df, family = binomial)
  coef_pvalue[i]<-coef(summary(glm.fit))[2,4]
  coef[i]<-coef(glm.fit)[2]
}
Seqtag<-colnames(pro_data)
Seqtag<-gsub("^.+\\.","",Seqtag,perl=T)

coef_fdrvalue<-p.adjust(coef_pvalue,method="BH")
coefficient<-data.frame(SeqId=Seqtag,coefficient=coef,coef_pvalue=coef_pvalue,coef_fdrvalue=coef_fdrvalue)
coefficient_Annoted<-merge(coefficient,
                           Annotation,by="SeqId",
                           all.x=T)
write.csv(coefficient_Annoted,"output_file/D1_VarianceFilter_cPDFF_W16_coefficient.csv")

#Cor.test for Point-biserial correlation
cor_pvalue<-c()
cor_value<-c()
for(i in 1:500){
  cor_pvalue[i]<-cor.test(pro_data[,i],PDFF_W16)$p.value
  cor_value[i]<-cor.test(pro_data[,i],PDFF_W16)$estimate
}

Seqtag<-colnames(pro_data)
Seqtag<-gsub("^.+\\.","",Seqtag,perl=T)

cor_fdrvalue<-p.adjust(cor_pvalue,method="BH")
correlation<-data.frame(SeqId=Seqtag,correlation=cor_value,cor_pvalue=cor_pvalue,cor_fdrvalue=cor_fdrvalue)
correlation_Annoted<-merge(correlation,
                           Annotation,by="SeqId",
                           all.x=T)

