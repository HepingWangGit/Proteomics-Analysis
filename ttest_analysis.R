library(readat)
#Read Data

setwd("/Users/hwang/Desktop/SomaLogic Data/Somanic")

File<-"output_file/somaLogic.EDTA.limmaNorm.txt"
data <- read.table(File)
data<-as.data.frame(t(data))
meta_data<-read.csv("output_file/Nash_EDTA_Metadata_clean_HW.csv",row.names = 1)

index<-match(rownames(data),meta_data$UniqueID)
meta_data<-meta_data[index,]

Annotation<-read.table("output_file/Annotation.txt",sep="\t",quote = "",comment.char = "",header = T)

#TTest

#Treatment (AXA1125 vs Placebo)
summary(as.factor(meta_data$Sample_Group))

Ttest_p<-function(df,ind1,ind2){
  group1<-df[ind1]
  group2<-df[ind2]
  sd1<-sd(group1)
  sd2<-sd(group2)
  sd3<-sd(c(group1,group2))
  if(min(sd1,sd2,sd3)<0.01){
    return(1.)
  }
  else{
  p_value<-t.test(group1,group2)$p.value
  return(p_value)
  }
}

Ttest_Stat<-function(df,ind1,ind2){
  group1<-df[ind1]
  group2<-df[ind2]
  sd1<-sd(group1)
  sd2<-sd(group2)
  sd3<-sd(c(group1,group2))
  if(min(sd1,sd2,sd3)<0.01){
    return(0)
  }
  else{
  statistic<-t.test(group1,group2)$statistic
  return(statistic)
  }
}

Pairedttest_p<-function(df,ind1,ind2){
  group1<-df[ind1]
  group2<-df[ind2]
  sd1<-sd(group1)
  sd2<-sd(group2)
  sd3<-sd(c(group1,group2))
  if(min(sd1,sd2,sd3)<0.01){
    return(1.)
  }
  else{
    p_value<-t.test(group1,group2,paired = TRUE, 
                    alternative = "two.sided")$p.value
    return(p_value)
  }
}

Pairedttest_Stat<-function(df,ind1,ind2){
  group1<-df[ind1]
  group2<-df[ind2]
  sd1<-sd(group1)
  sd2<-sd(group2)
  sd3<-sd(c(group1,group2))
  if(min(sd1,sd2,sd3)<0.01){
    return(0)
  }
  else{
    statistic<-t.test(group1,group2, paired = TRUE,
                      alternative = "two.sided")$statistic
    return(statistic)
  }
}

Log2Ratio<-function(df,ind1,ind2){
  group1<-df[ind1]
  group2<-df[ind2]
  log2ratio<-mean(group1)-mean(group2)
  return(log2ratio)
}

SeqID<-colnames(data)
SeqID<-gsub("^.+\\.","",SeqID,perl=T)

AXAID<-meta_data[which(meta_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
PlaceboID<-meta_data[which(meta_data$Sample_Group=="placebo 24g bid"),]$UniqueID

index1<-which(rownames(data)%in%AXAID)
index2<-which(rownames(data)%in%PlaceboID)
P_value<-apply(data,2, function(x)Ttest_p(x,index1,index2))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index1,index2))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index1,index2))
fdr <- p.adjust(P_value,method="BH")

ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)

write.table(ttest_Result,"output_file/Treatment_ttest_all.txt",quote=F,sep="\t",row.names=F)

ttest_Sig_Result<-ttest_Result[which(ttest_Result$fdr.value<=0.05),]
write.table(ttest_Sig_Result,"output_file/Treatment_ttest_Significant.txt",quote=F,sep="\t",row.names=F)

#Annotation Attachment

ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_all_Annotated.txt",quote=F,sep="\t",row.names=F)

ttest_Result_Sig_Annoted<-merge(ttest_Sig_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Sig_Annoted,"output_file/Treatment_ttest_Sig_Annotated.txt",quote=F,sep="\t",row.names=F)

#Treatment Timeline*Treatment
summary(as.factor(meta_data$Time_point))
D1ID_AXA<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
index1<-which(rownames(data)%in%D1ID_AXA)

D1ID_PLA<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="placebo 24g bid"),]$UniqueID
index2<-which(rownames(data)%in%D1ID_PLA)

W8ID_AXA<-meta_data[which(meta_data$Time_point=="W8"&meta_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
index3<-which(rownames(data)%in%W8ID_AXA)

W8ID_PLA<-meta_data[which(meta_data$Time_point=="W8"&meta_data$Sample_Group=="placebo 24g bid"),]$UniqueID
index4<-which(rownames(data)%in%W8ID_PLA)

W16ID_AXA<-meta_data[which(meta_data$Time_point=="W16"&meta_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
index5<-which(rownames(data)%in%W16ID_AXA)

W16ID_PLA<-meta_data[which(meta_data$Time_point=="W16"&meta_data$Sample_Group=="placebo 24g bid"),]$UniqueID
index6<-which(rownames(data)%in%W16ID_PLA)

##D1
P_value<-apply(data,2, function(x)Ttest_p(x,index1,index2))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index1,index2))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index1,index2))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_D1.txt",quote=F,sep="\t",row.names=F)

##W8
P_value<-apply(data,2, function(x)Ttest_p(x,index3,index4))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index3,index4))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index3,index4))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_W8.txt",quote=F,sep="\t",row.names=F)

##W16
P_value<-apply(data,2, function(x)Ttest_p(x,index5,index6))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index5,index6))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index5,index6))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_W16.txt",quote=F,sep="\t",row.names=F)

#Timeline

#AXA1125
#Timeline index pair matching
#D1/D8
D1ID_AXA<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
D1ID_AXA_s<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="axa1125 24g bid"),]$SubjectID
index1<-which(rownames(data)%in%D1ID_AXA)

D1ID_PLA<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="placebo 24g bid"),]$UniqueID
D1ID_PLA_s<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="placebo 24g bid"),]$SubjectID
index2<-which(rownames(data)%in%D1ID_PLA)

W8_AXA<-meta_data[which(meta_data$Time_point=="W8"&meta_data$Sample_Group=="axa1125 24g bid"),]
W8_AXA<-W8_AXA[match(D1ID_AXA_s,W8_AXA$SubjectID),]
index3<-which(rownames(data)%in%W8_AXA$UniqueID)

W8_PLA<-meta_data[which(meta_data$Time_point=="W8"&meta_data$Sample_Group=="placebo 24g bid"),]
W8_PLA<-W8_PLA[match(D1ID_PLA_s,W8_PLA$SubjectID),]
index4<-which(rownames(data)%in%W8_PLA$UniqueID)

##D1/W8 ttest
P_value<-apply(data,2, function(x)Pairedttest_p(x,index1,index3))
t_Stat<-apply(data,2, function(x)Pairedttest_Stat(x,index1,index3))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index1,index3))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_AXA_D1W8.txt",quote=F,sep="\t",row.names=F)

P_value<-apply(data,2, function(x)Pairedttest_p(x,index2,index4))
t_Stat<-apply(data,2, function(x)Pairedttest_Stat(x,index2,index4))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index2,index4))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_PLA_D1W8.txt",quote=F,sep="\t",row.names=F)

#D1/W16
#Paired ttest index matching
D1ID_AXA<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="axa1125 24g bid"),]
W16_AXA<-meta_data[which(meta_data$Time_point=="W16"&meta_data$Sample_Group=="axa1125 24g bid"),]
W16_AXA<-W16_AXA[na.omit(match(D1ID_AXA$SubjectID,W16_AXA$SubjectID)),]
D1ID_AXA<-D1ID_AXA[na.omit(match(W16_AXA$SubjectID,D1ID_AXA$SubjectID)),]
index1<-which(rownames(data)%in%D1ID_AXA$UniqueID)
index3<-which(rownames(data)%in%W16_AXA$UniqueID)

D1ID_PLA<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="placebo 24g bid"),]
W16_PLA<-meta_data[which(meta_data$Time_point=="W16"&meta_data$Sample_Group=="placebo 24g bid"),]
W16_PLA<-W16_PLA[na.omit(match(D1ID_PLA$SubjectID,W16_PLA$SubjectID)),]
D1ID_PLA<-D1ID_PLA[na.omit(match(W16_PLA$SubjectID,D1ID_PLA$SubjectID)),]
index2<-which(rownames(data)%in%D1ID_PLA$UniqueID)
index4<-which(rownames(data)%in%W16_PLA$UniqueID)

#D1/W16 ttest
P_value<-apply(data,2, function(x)Pairedttest_p(x,index1,index3))
t_Stat<-apply(data,2, function(x)Pairedttest_Stat(x,index1,index3))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index1,index3))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_AXA_D1W16.txt",quote=F,sep="\t",row.names=F)

P_value<-apply(data,2, function(x)Pairedttest_p(x,index2,index4))
t_Stat<-apply(data,2, function(x)Pairedttest_Stat(x,index2,index4))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index2,index6))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_PLA_D1W16.txt",quote=F,sep="\t",row.names=F)

#W8/W16
#Index Matching for W8/W16
W8_AXA<-meta_data[which(meta_data$Time_point=="W8"&meta_data$Sample_Group=="axa1125 24g bid"),]
W16_AXA<-meta_data[which(meta_data$Time_point=="W16"&meta_data$Sample_Group=="axa1125 24g bid"),]
W16_AXA<-W16_AXA[na.omit(match(W8_AXA$SubjectID,W16_AXA$SubjectID)),]
W8_AXA<-W8_AXA[na.omit(match(W16_AXA$SubjectID,W8_AXA$SubjectID)),]
index1<-which(rownames(data)%in%W8_AXA$UniqueID)
index3<-which(rownames(data)%in%W16_AXA$UniqueID)

W8_PLA<-meta_data[which(meta_data$Time_point=="W8"&meta_data$Sample_Group=="placebo 24g bid"),]
W16_PLA<-meta_data[which(meta_data$Time_point=="W16"&meta_data$Sample_Group=="placebo 24g bid"),]
W16_PLA<-W16_PLA[na.omit(match(W8_PLA$SubjectID,W16_PLA$SubjectID)),]
W8_PLA<-W8_PLA[na.omit(match(W16_PLA$SubjectID,W8_PLA$SubjectID)),]
index2<-which(rownames(data)%in%W8_PLA$UniqueID)
index4<-which(rownames(data)%in%W16_PLA$UniqueID)

#ttest
P_value<-apply(data,2, function(x)Pairedttest_p(x,index1,index3))
t_Stat<-apply(data,2, function(x)Pairedttest_Stat(x,index3,index3))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index1,index3))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_AXA_W8W16.txt",quote=F,sep="\t",row.names=F)

P_value<-apply(data,2, function(x)Pairedttest_p(x,index2,index4))
t_Stat<-apply(data,2, function(x)Pairedttest_Stat(x,index2,index4))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index2,index4))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_PLA_W8W16.txt",quote=F,sep="\t",row.names=F)

#AXA vs PLA within Diabetes group

Dia_data<-meta_data[meta_data$T2DIABFL=="Y",]
NonDia_data<-meta_data[meta_data$T2DIABFL=="N",]

#Diabetes
#Treatment Timeline*Treatment
D1ID_AXA<-Dia_data[which(Dia_data$Time_point=="D1"&Dia_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
index1<-which(rownames(data)%in%D1ID_AXA)

D1ID_PLA<-Dia_data[which(Dia_data$Time_point=="D1"&Dia_data$Sample_Group=="placebo 24g bid"),]$UniqueID
index2<-which(rownames(data)%in%D1ID_PLA)

W8ID_AXA<-Dia_data[which(Dia_data$Time_point=="W8"&Dia_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
index3<-which(rownames(data)%in%W8ID_AXA)

W8ID_PLA<-Dia_data[which(Dia_data$Time_point=="W8"&Dia_data$Sample_Group=="placebo 24g bid"),]$UniqueID
index4<-which(rownames(data)%in%W8ID_PLA)

W16ID_AXA<-Dia_data[which(Dia_data$Time_point=="W16"&Dia_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
index5<-which(rownames(data)%in%W16ID_AXA)

W16ID_PLA<-Dia_data[which(Dia_data$Time_point=="W16"&Dia_data$Sample_Group=="placebo 24g bid"),]$UniqueID
index6<-which(rownames(data)%in%W16ID_PLA)

##D1
P_value<-apply(data,2, function(x)Ttest_p(x,index1,index2))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index1,index2))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index1,index2))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_D1_Dia.txt",quote=F,sep="\t",row.names=F)

##W8
P_value<-apply(data,2, function(x)Ttest_p(x,index3,index4))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index3,index4))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index3,index4))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_W8_Dia.txt",quote=F,sep="\t",row.names=F)

##W16
P_value<-apply(data,2, function(x)Ttest_p(x,index5,index6))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index5,index6))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index5,index6))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_W16_Dia.txt",quote=F,sep="\t",row.names=F)

#Non-Diabetes group

#Treatment Timeline*Treatment
D1ID_AXA<-NonDia_data[which(NonDia_data$Time_point=="D1"&NonDia_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
index1<-which(rownames(data)%in%D1ID_AXA)

D1ID_PLA<-NonDia_data[which(NonDia_data$Time_point=="D1"&NonDia_data$Sample_Group=="placebo 24g bid"),]$UniqueID
index2<-which(rownames(data)%in%D1ID_PLA)

W8ID_AXA<-NonDia_data[which(NonDia_data$Time_point=="W8"&NonDia_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
index3<-which(rownames(data)%in%W8ID_AXA)

W8ID_PLA<-NonDia_data[which(NonDia_data$Time_point=="W8"&NonDia_data$Sample_Group=="placebo 24g bid"),]$UniqueID
index4<-which(rownames(data)%in%W8ID_PLA)

W16ID_AXA<-NonDia_data[which(NonDia_data$Time_point=="W16"&NonDia_data$Sample_Group=="axa1125 24g bid"),]$UniqueID
index5<-which(rownames(data)%in%W16ID_AXA)

W16ID_PLA<-NonDia_data[which(NonDia_data$Time_point=="W16"&NonDia_data$Sample_Group=="placebo 24g bid"),]$UniqueID
index6<-which(rownames(data)%in%W16ID_PLA)

##D1
P_value<-apply(data,2, function(x)Ttest_p(x,index1,index2))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index1,index2))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index1,index2))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_D1_NonDia.txt",quote=F,sep="\t",row.names=F)

##W8
P_value<-apply(data,2, function(x)Ttest_p(x,index3,index4))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index3,index4))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index3,index4))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_W8_NonDia.txt",quote=F,sep="\t",row.names=F)

##W16
P_value<-apply(data,2, function(x)Ttest_p(x,index5,index6))
t_Stat<-apply(data,2, function(x)Ttest_Stat(x,index5,index6))
log2ratio<-apply(data,2, function(x)Log2Ratio(x,index5,index6))
fdr <- p.adjust(P_value,method="BH")
ttest_Result<-data.frame(SeqId=SeqID,p_value=P_value,t.statistics=t_Stat,fdr.value=fdr,log2ratio=log2ratio)
ttest_Result_Annoted<-merge(ttest_Result,Annotation,by="SeqId",all.x=T)
write.table(ttest_Result_Annoted,"output_file/Treatment_ttest_W16_NonDia.txt",quote=F,sep="\t",row.names=F)


