#Data Summary

File<-"output_file/somaLogic.EDTA.limmaNorm.txt"
data <- read.table(File)
data<-as.data.frame(t(data))

meta_data<-read.csv("output_file/Nash_EDTA_Metadata_clean_HW.csv",row.names = 1)
index<-match(rownames(data),meta_data$UniqueID)
meta_data<-meta_data[index,]

table(meta_data$Sample_Group,meta_data$Time_point)

#Heatmap for correlated proteins

#TOP30 Heatmap of PDFF selected proteins and how they correlated with other features
pdff_axa_w8<-read.csv("output_file/PDFF_correlations_AXA_W8.csv",row.names = 1)
alt_axa_w8<-read.csv("output_file/ALT_correlations_AXA_W8.csv",row.names = 1)
ct1_axa_w8<-read.csv("output_file/CT1_correlations_AXA_W8.csv",row.names = 1)
HOMAIR_axa_w8<-read.csv("output_file/HOMAIR_correlations_AXA_W8.csv",row.names = 1)
HBA1C_axa_w8<-read.csv("output_file/HBA1C_correlations_AXA_W8.csv",row.names = 1)
PROC3_axa_w8<-read.csv("output_file/PROC3_correlations_AXA_W8.csv",row.names = 1)


pdff_pla_w8<-read.csv("output_file/PDFF_correlations_PLA_W8.csv",row.names = 1)
alt_pla_w8<-read.csv("output_file/ALT_correlations_PLA_W8.csv",row.names = 1)
ct1_pla_w8<-read.csv("output_file/CT1_correlations_PLA_W8.csv",row.names = 1)
HOMAIR_pla_w8<-read.csv("output_file/HOMAIR_correlations_PLA_W8.csv",row.names = 1)
HBA1C_pla_w8<-read.csv("output_file/HBA1C_correlations_PLA_W8.csv",row.names = 1)
PROC3_pla_w8<-read.csv("output_file/PROC3_correlations_PLA_W8.csv",row.names = 1)

#pdff
diff_seq<-setdiff(pdff_axa_w8[which(pdff_axa_w8$p.value<=0.05),]$SeqId,pdff_pla_w8[which(pdff_pla_w8$p.value<=0.05),]$SeqId)
diff_axa<-pdff_axa_w8[which(pdff_axa_w8$SeqId%in%diff_seq),]

top30<-diff_axa[1:30,]$SeqId
tag<-diff_axa[1:30,]$Target

pdff=pdff_axa_w8[which(pdff_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(pdff)<-c("pdff","Target")
alt=alt_axa_w8[which(alt_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(alt)<-c("alt","Target")
ct1=ct1_axa_w8[which(ct1_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(ct1)<-c("ct1","Target")
homair=HOMAIR_axa_w8[which(HOMAIR_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(homair)<-c("homair","Target")
hba1c=HBA1C_axa_w8[which(HBA1C_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(hba1c)<-c("hba1c","Target")
proc3=PROC3_axa_w8[which(PROC3_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(proc3)<-c("proc3","Target")


df<-merge(pdff,alt,by="Target",all.x=T)  
df<-merge(df,ct1,by="Target",all.x=T)  
df<-merge(df,homair,by="Target",all.x=T)  
df<-merge(df,hba1c,by="Target",all.x=T)  
df<-merge(df,proc3,by="Target",all.x=T)  

rownames(df)<-df$Target
df<-df[order(df$pdff),]
df<-as.matrix(df[,-1])
heatmap(df, Colv = NA, Rowv = NA, scale="column")
#legend(x="bottomright", legend=c("min", "med", "max"),fill=heat.colors(3))
#write.csv(df,"output_file/pdff_top30.csv")

#alt
diff_seq<-setdiff(alt_axa_w8[which(alt_axa_w8$p.value<=0.05),]$SeqId,alt_pla_w8[which(alt_pla_w8$p.value<=0.05),]$SeqId)
diff_axa<-alt_axa_w8[which(alt_axa_w8$SeqId%in%diff_seq),]

top30<-diff_axa[1:30,]$SeqId
tag<-diff_axa[1:30,]$Target

pdff=pdff_axa_w8[which(pdff_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(pdff)<-c("pdff","Target")
alt=alt_axa_w8[which(alt_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(alt)<-c("alt","Target")
ct1=ct1_axa_w8[which(ct1_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(ct1)<-c("ct1","Target")
homair=HOMAIR_axa_w8[which(HOMAIR_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(homair)<-c("homair","Target")
hba1c=HBA1C_axa_w8[which(HBA1C_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(hba1c)<-c("hba1c","Target")
proc3=PROC3_axa_w8[which(PROC3_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(proc3)<-c("proc3","Target")


df<-merge(pdff,alt,by="Target",all.x=T)  
df<-merge(df,ct1,by="Target",all.x=T)  
df<-merge(df,homair,by="Target",all.x=T)  
df<-merge(df,hba1c,by="Target",all.x=T)  
df<-merge(df,proc3,by="Target",all.x=T)  

rownames(df)<-df$Target
df<-df[order(df$alt),]
df<-as.matrix(df[,-1])
heatmap(df, Colv = NA, Rowv = NA, scale="column")

#CT1
diff_seq<-setdiff(ct1_axa_w8[which(ct1_axa_w8$p.value<=0.05),]$SeqId,ct1_pla_w8[which(ct1_pla_w8$p.value<=0.05),]$SeqId)
diff_axa<-ct1_axa_w8[which(ct1_axa_w8$SeqId%in%diff_seq),]

top30<-ct1_axa_w8[1:30,]$SeqId
tag<-ct1_axa_w8[1:30,]$Target

pdff=pdff_axa_w8[which(pdff_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(pdff)<-c("pdff","Target")
alt=alt_axa_w8[which(alt_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(alt)<-c("alt","Target")
ct1=ct1_axa_w8[which(ct1_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(ct1)<-c("ct1","Target")
homair=HOMAIR_axa_w8[which(HOMAIR_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(homair)<-c("homair","Target")
hba1c=HBA1C_axa_w8[which(HBA1C_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(hba1c)<-c("hba1c","Target")
proc3=PROC3_axa_w8[which(PROC3_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(proc3)<-c("proc3","Target")


df<-merge(pdff,alt,by="Target",all.x=T)  
df<-merge(df,ct1,by="Target",all.x=T)  
df<-merge(df,homair,by="Target",all.x=T)  
df<-merge(df,hba1c,by="Target",all.x=T)  
df<-merge(df,proc3,by="Target",all.x=T)  

rownames(df)<-df$Target
df<-df[order(df$ct1),]
df<-as.matrix(df[,-1])
heatmap(df, Colv = NA, Rowv = NA, scale="column")

#HOMAIR
diff_seq<-setdiff(HOMAIR_axa_w8[which(HOMAIR_axa_w8$p.value<=0.05),]$SeqId,HOMAIR_pla_w8[which(HOMAIR_pla_w8$p.value<=0.05),]$SeqId)
diff_axa<-HOMAIR_axa_w8[which(HOMAIR_axa_w8$SeqId%in%diff_seq),]

top30<-HOMAIR_axa_w8[1:30,]$SeqId
tag<-HOMAIR_axa_w8[1:30,]$Target

pdff=pdff_axa_w8[which(pdff_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(pdff)<-c("pdff","Target")
alt=alt_axa_w8[which(alt_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(alt)<-c("alt","Target")
ct1=ct1_axa_w8[which(ct1_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(ct1)<-c("ct1","Target")
homair=HOMAIR_axa_w8[which(HOMAIR_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(homair)<-c("homair","Target")
hba1c=HBA1C_axa_w8[which(HBA1C_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(hba1c)<-c("hba1c","Target")
proc3=PROC3_axa_w8[which(PROC3_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(proc3)<-c("proc3","Target")


df<-merge(pdff,alt,by="Target",all.x=T)  
df<-merge(df,ct1,by="Target",all.x=T)  
df<-merge(df,homair,by="Target",all.x=T)  
df<-merge(df,hba1c,by="Target",all.x=T)  
df<-merge(df,proc3,by="Target",all.x=T)  

rownames(df)<-df$Target
df<-df[order(df$homair),]
df<-as.matrix(df[,-1])
heatmap(df, Colv = NA, Rowv = NA, scale="column")

#HBA1C
diff_seq<-setdiff(HBA1C_axa_w8[which(HBA1C_axa_w8$p.value<=0.05),]$SeqId,HBA1C_pla_w8[which(HBA1C_pla_w8$p.value<=0.05),]$SeqId)
diff_axa<-HBA1C_axa_w8[which(HBA1C_axa_w8$SeqId%in%diff_seq),]

pdff=pdff_axa_w8[which(pdff_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(pdff)<-c("pdff","Target")
alt=alt_axa_w8[which(alt_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(alt)<-c("alt","Target")
ct1=ct1_axa_w8[which(ct1_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(ct1)<-c("ct1","Target")
homair=HOMAIR_axa_w8[which(HOMAIR_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(homair)<-c("homair","Target")
hba1c=HBA1C_axa_w8[which(HBA1C_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(hba1c)<-c("hba1c","Target")
proc3=PROC3_axa_w8[which(PROC3_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(proc3)<-c("proc3","Target")


df<-merge(pdff,alt,by="Target",all.x=T)  
df<-merge(df,ct1,by="Target",all.x=T)  
df<-merge(df,homair,by="Target",all.x=T)  
df<-merge(df,hba1c,by="Target",all.x=T)  
df<-merge(df,proc3,by="Target",all.x=T)  

rownames(df)<-df$Target
df<-df[order(df$hba1c),]
df<-as.matrix(df[,-1])
heatmap(df, Colv = NA, Rowv = NA, scale="column")

#prco3
diff_seq<-setdiff(PROC3_axa_w8[which(PROC3_axa_w8$p.value<=0.05),]$SeqId,PROC3_pla_w8[which(PROC3_pla_w8$p.value<=0.05),]$SeqId)
diff_axa<-PROC3_axa_w8[which(PROC3_axa_w8$SeqId%in%diff_seq),]

pdff=pdff_axa_w8[which(pdff_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(pdff)<-c("pdff","Target")
alt=alt_axa_w8[which(alt_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(alt)<-c("alt","Target")
ct1=ct1_axa_w8[which(ct1_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(ct1)<-c("ct1","Target")
homair=HOMAIR_axa_w8[which(HOMAIR_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(homair)<-c("homair","Target")
hba1c=HBA1C_axa_w8[which(HBA1C_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(hba1c)<-c("hba1c","Target")
proc3=PROC3_axa_w8[which(PROC3_axa_w8$SeqId%in%top30),c("correlation","Target")]
colnames(proc3)<-c("proc3","Target")


df<-merge(pdff,alt,by="Target",all.x=T)  
df<-merge(df,ct1,by="Target",all.x=T)  
df<-merge(df,homair,by="Target",all.x=T)  
df<-merge(df,hba1c,by="Target",all.x=T)  
df<-merge(df,proc3,by="Target",all.x=T)  

rownames(df)<-df$Target
df<-df[order(df$proc3),]
df<-as.matrix(df[,-1])
heatmap(df, Colv = NA, Rowv = NA, scale="column")
#Scatterplot for certain proteins

#Read Data

File<-"output_file/somaLogic.EDTA.limmaNorm.txt"
data <- read.table(File)
data<-as.data.frame(t(data))

meta_data<-read.csv("output_file/Nash_EDTA_Metadata_clean_HW.csv",row.names = 1)
index<-match(rownames(data),meta_data$UniqueID)
meta_data<-meta_data[index,]

Annotation<-read.table("output_file/Annotation.txt",sep="\t",quote = "",comment.char = "",header = T)

#AXA1125 Timeline index pair matching D1/W8

D1ID_AXA<-meta_data[which(meta_data$Time_point=="D1"&meta_data$Sample_Group=="axa1125 24g bid"),]
D1ID_AXA_s<-D1ID_AXA$SubjectID
index1<-which(rownames(data)%in%D1ID_AXA$UniqueID)
D1_data<-data[index1,]

W8_AXA<-meta_data[which(meta_data$Time_point=="W8"&meta_data$Sample_Group=="axa1125 24g bid"),]
W8_AXA<-W8_AXA[match(D1ID_AXA_s,W8_AXA$SubjectID),]
index2<-which(rownames(data)%in%W8_AXA$UniqueID)
W8_data<-data[index2,]

Diff_exp<-W8_data-D1_data

PDFF<-W8_AXA$PDFF-D1ID_AXA$PDFF
HBA1C<-W8_AXA$HBA1C-D1ID_AXA$HBA1C
ALT<-W8_AXA$ALT-D1ID_AXA$ALT
CT1<-W8_AXA$CT1-D1ID_AXA$CT1
HOMAIR<-W8_AXA$HOMAIR-D1ID_AXA$HOMAIR

#(1) SEM5B
x<-Diff_exp$`SeqId.20527-47`
y<-PDFF
plot(x, y, main="Scatterplot of SEM5B vs PDFF", 
     xlab="log2(delta_SEM5B) ", ylab=" delta_PDFF ", pch=19)
abline(v=0, col=c("red"), lty=c(2), lwd=c(3))

#(2) 
x<-Diff_exp$`SeqId.6440-31`
y<-PDFF
plot(x, y, main="Scatterplot of MFAP5 vs PDFF", 
     xlab="log2(delta_MFAP5) ", ylab=" delta_PDFF ", pch=19)
abline(v=0, col=c("red"), lty=c(2), lwd=c(3))









