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
pkdd_data<-read.csv("output_file/AXA1125_003_PKPD_Clean_HW.csv")

#Feature selection

#Baseline Features
pkdd_data<-pkdd_data[,c(colnames(pkdd_data)[grep("B",colnames(pkdd_data))])]
index<-match(meta_data$SubjectID,pkdd_data$SUBJID)
pkdd_data<-pkdd_data[index,]
meta_data<-cbind(meta_data,pkdd_data[,c(6:57)])
meta_data<-meta_data[,-which(colnames(meta_data)=="BPDFF")]
meta_data<-meta_data[,-which(colnames(meta_data)=="BCT1")]
meta_data<-meta_data[,-which(colnames(meta_data)=="BALT")]

##AXA1125 Timeline index pair matching D1/W8

D1ID_AXA<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="axa1125 24g bid"),]
D1ID_AXA_s<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="axa1125 24g bid"),]$SubjectID
index1<-which(rownames(data)%in%D1ID_AXA$UniqueID)
D1_data<-data[index1,]

W8_AXA<-meta_data[which(meta_data$Time_point=="W8"&meta_data$Sample_Group=="axa1125 24g bid"),]
W8_AXA<-W8_AXA[match(D1ID_AXA_s,W8_AXA$SubjectID),]
index2<-which(rownames(data)%in%W8_AXA$UniqueID)
W8_data<-data[index2,]

Diff_exp<-W8_data-D1_data #Prot_data

PDFF<-W8_AXA$PDFF-D1ID_AXA$PDFF #Responce
c_PDFF<-D1ID_AXA$classifier_pdff_W8#Classifier

##Step1. Correlation test between delta(W8-D1) proteomics and delta pdff 
###Proteomics Features Preselection
library(caret)
#p.value<-c()
#corr.stat<-c()
#for(i in 1:dim(Diff_exp)[2]){
#  p.value[i]<-cor.test(Diff_exp[,i],PDFF,method = "spearman")$p.value
#  corr.stat[i]<-cor.test(Diff_exp[,i],PDFF,method = "spearman")$estimate
#}
#corrTest_result<-data.frame(SeqId=colnames(Diff_exp),p.value=p.value,corr.stat=corr.stat)

#index<-which(corrTest_result$p.value<=0.05)
#Diff_exp<-Diff_exp[,index]

#corMatrix<-cor(Diff_exp) #check inter correlation of proteins
#hc = findCorrelation(corMatrix, cutoff=0.8) # putt any value as a "cutoff" 
#hc = sort(hc)
#reduced_Data = Diff_exp[,-c(hc)]

##Select clinical features
cli_data<-D1ID_AXA[,c(9:11,24:33,36:40,43:70)] #A total of 46 Clinical features at baseline
cli_data<-cli_data[,1:7]

df<-as.data.frame(cbind(cli_data,c_PDFF))
coef_pvalue<-c()
coef<-c()
for(i in 1:ncol(cli_data)){
  glm.fit <- glm(factor(c_PDFF) ~cli_data[,i] , data=df, family = binomial)
  coef_pvalue[i]<-coef(summary(glm.fit))[2,4]
  coef[i]<-coef(glm.fit)[2]
}
colnames(cli_data)[which(coef_pvalue<=0.1)]
cor(cli_data[,which(coef_pvalue<=0.1)])#indicating BMI and BWC is highly correlated Thus we drop BMI select BWC

##Select Protein features
df<-as.data.frame(cbind(Diff_exp,c_PDFF))
coef_pvalue<-c()
coef<-c()
for(i in 1:ncol(Diff_exp)){
  glm.fit <- glm(factor(c_PDFF) ~Diff_exp[,i] , data=df, family = binomial)
  coef_pvalue[i]<-coef(summary(glm.fit))[2,4]
  coef[i]<-coef(glm.fit)[2]
}
Seqtag<-colnames(Diff_exp)
Seqtag<-gsub("^.+\\.","",Seqtag,perl=T)
coefficient<-data.frame(SeqId=Seqtag,coefficient=coef,coef_pvalue=coef_pvalue)
coefficient_Annoted<-merge(coefficient,
                           Annotation,by="SeqId",
                           all.x=T)
write.csv(coefficient_Annoted,"output_file/D1W8_cPDFF_coefficient.csv")

coefficient_Annoted[coefficient_Annoted$coef_pvalue<=0.05,]
#Feature Selection

ctrl <- trainControl(method = "LOOCV")
#WaistC
model1 <- train(as.factor(c_PDFF) ~ WaistC, data = df, method = "multinom", trControl = ctrl)
summary(model1)

fit_wa<-glm(factor(c_PDFF) ~WaistC +ACE2 , data=df, family = binomial)
fit_w<-glm(factor(c_PDFF) ~WaistC  , data=df, family = binomial)
lrtest(fit_wa,fit_w) #fit_wa fit significantly better than the waist only model
