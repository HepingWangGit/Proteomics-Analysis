#Library

library(ggplot2)
library(ggfortify)

#Read Data

setwd("/Users/hwang/Desktop/SomaLogic Data/Somanic")
File<-"output_file/somaLogic.EDTA.limmaNorm.txt"
data <- read.table(File)
data<-as.data.frame(t(data))
meta_data<-read.csv("SOMAscan_NASH_Data/Sample Submission Form_Axcella_Russell_NASH_EDTA.csv")

index<-match(rownames(data),meta_data$UniqueID)
meta_data<-meta_data[index,]

sample_name<-paste0(meta_data$Sample_Group,meta_data$Time_point)
dt<-data
df <- as.data.frame(cbind(dt,sample_name))
pca_res <- prcomp(dt, scale. = TRUE)
g <- autoplot(pca_res,data=df,colour='sample_name',size=2)+ggtitle("sample_name")+
  theme(text=element_text(size=14))
ggsave("Figure/PCA.sample.name.jpg",g,device="jpeg")


index<-which(meta_data$Sample_Group=="axa1125 24g bid")
dt_AXA<-as.matrix(dt[index,])
Time_point<-meta_data$Time_point[index]
df_AXA<-as.data.frame(cbind(dt_AXA,Time_point))
pca_res <- prcomp(dt_AXA, scale. = TRUE)
g <- autoplot(pca_res,data=df_AXA,colour='Time_point',size=2)+ggtitle("time")+
  theme(text=element_text(size=14))
ggsave("Figure/PCA.time_AXA.jpg",g,device="jpeg")

index<-which(meta_data$Sample_Group=="placebo 24g bid")
dt_PLA<-as.matrix(dt[index,])
Time_point<-meta_data$Time_point[index]
df_PLA<-as.data.frame(cbind(dt_PLA,Time_point))
pca_res <- prcomp(dt_PLA, scale. = TRUE)
g <- autoplot(pca_res,data=df_PLA,colour='Time_point',size=2)+ggtitle("time")+
  theme(text=element_text(size=14))
ggsave("Figure/PCA.time_PLA.jpg",g,device="jpeg")

combination <- meta_data$Sample_Group
df <- as.data.frame(cbind(dt,combination))
pca_res <- prcomp(dt, scale. = TRUE)
g <- autoplot(pca_res,data=df,colour='combination',size=2)+ggtitle("combination")+
  theme(text=element_text(size=14))
ggsave("Figure/PCA.combination.jpg",g,device="jpeg")

# PCA Group Distance Analysis

Group_Distance<-function(disMatrix,GroupID1,GroupID2){
  DisGroup<-c()
  
  D<-0
  for(i in 1:length(GroupID1)){
    for(j in 1:length(GroupID2)){
      index_i<-which(rownames(disMatrix)==GroupID1[i])
      index_j<-which(colnames(disMatrix)==GroupID2[j])
      D<-D+disMatrix[index_i,index_j]
    }
  }
  D<-D/(length(GroupID1)*length(GroupID2))
  
  
  D_Axa<-0
  Comb_pair<-combn(length(GroupID1),2)
  for (j in 1:ncol(Comb_pair)) {
    index_a<-c(Comb_pair[,j])[1]
    index_b<-c(Comb_pair[,j])[2]
    index_i<-which(rownames(disMatrix)==GroupID1[index_a])
    index_j<-which(colnames(disMatrix)==GroupID1[index_b])
    D_Axa<-D_Axa+disMatrix[index_i,index_j]
  }
  
  D_Axa<-2*D_Axa/(length(GroupID1)*(length(GroupID1)-1))
  
  D_Placebo<-0
  Comb_pair<-combn(length(GroupID2),2)
  for (j in 1:ncol(Comb_pair)) {
    index_a<-c(Comb_pair[,j])[1]
    index_b<-c(Comb_pair[,j])[2]
    index_i<-which(rownames(disMatrix)==GroupID2[index_a])
    index_j<-which(colnames(disMatrix)==GroupID2[index_b])
    D_Placebo<-D_Placebo+disMatrix[index_i,index_j]
  }
  
  D_Placebo<-2*D_Placebo/(length(GroupID2)*(length(GroupID2)-1))
  
  DisGroup<-D/(D_Axa+D_Placebo)
  
  return(DisGroup)
}

DisPerm<-function(data,GroupID1,GroupID2,PermTime){
  DisGroup<-c()
  data<-data[which(rownames(data)%in%GroupID1|rownames(data)%in%GroupID2),]
  for (p in 1:PermTime){
  index<-sample(c(1:length(rownames(data))),replace = F)  
  Radom_data<-data[index,]
  rownames(Radom_data)<-rownames(data)
  pcorMatrix<-cor(t(Radom_data),method = "pearson")
  disMatrix<-1-pcorMatrix
  
  DisGroup[p]<-Group_Distance(disMatrix = disMatrix,GroupID1=GroupID1,GroupID2=GroupID2)
}
return(DisGroup)
}

#AXA1125 vs Placebo Group distance

Axa1125ID<-meta_data[meta_data$Sample_Group=="axa1125 24g bid",]$UniqueID
PlaceboID<-meta_data[meta_data$Sample_Group=="placebo 24g bid",]$UniqueID
dt<-data

pcorMatrix<-cor(t(dt),method = "pearson")
disMatrix<-1-pcorMatrix

Dis_Group<-Group_Distance(disMatrix = disMatrix,GroupID1=Axa1125ID,GroupID2=PlaceboID)
NullDis<-DisPerm(data = dt,GroupID1=Axa1125ID,GroupID2 = PlaceboID, PermTime = 10000)
P_value<-sum(NullDis>=Dis_Group)/length(NullDis) #0.0199

#Time point PCA group distance
D1ID<-meta_data[meta_data$Time_point=="D1",]$UniqueID
W8ID<-meta_data[meta_data$Time_point=="W8",]$UniqueID
W16ID<-meta_data[meta_data$Time_point=="W16",]$UniqueID
Dis_Group_D1W8<-Group_Distance(disMatrix = disMatrix,GroupID1=D1ID,GroupID2=W8ID)
Dis_Group_D1W16<-Group_Distance(disMatrix = disMatrix,GroupID1=D1ID,GroupID2=W16ID)
Dis_Group_W8W16<-Group_Distance(disMatrix = disMatrix,GroupID1=W8ID,GroupID2=W16ID)

NullDis_D1W8<-DisPerm(data = dt,GroupID1=D1ID,GroupID2 = W8ID, PermTime = 1000)
NullDis_D1W16<-DisPerm(data = dt,GroupID1=D1ID,GroupID2 = W16ID, PermTime = 1000)
NullDis_W8W16<-DisPerm(data = dt,GroupID1=W8ID,GroupID2 = W16ID, PermTime = 1000)

P_value_D1W8<-sum(NullDis_D1W8>=Dis_Group_D1W8)/length(NullDis_D1W8) 
P_value_D1W16<-sum(NullDis_D1W16>=Dis_Group_D1W16)/length(NullDis_D1W16) 
P_value_W8W16<-sum(NullDis_W8W16>=Dis_Group_W8W16)/length(NullDis_W8W16) 

P_value=data.frame(Group=c("D1 vs W8","D1 vs W16","W8 vs W16"),
                   Group.Distance=c(Dis_Group_D1W8,Dis_Group_D1W16,Dis_Group_W8W16),
                   p.value=c(P_value_D1W8,P_value_D1W16,P_value_W8W16))
write.table(P_value,file="output_file/Timepoint_PCA_GroupDistance.txt",sep="\t")

# Sample PCA group seperation
AXA_D1_ID<-meta_data[meta_data$Time_point=="D1"&meta_data$Sample_Group=="axa1125 24g bid",]$UniqueID
AXA_W8_ID<-meta_data[meta_data$Time_point=="W8"&meta_data$Sample_Group=="axa1125 24g bid",]$UniqueID
AXA_W16_ID<-meta_data[meta_data$Time_point=="W16"&meta_data$Sample_Group=="axa1125 24g bid",]$UniqueID

PLA_D1_ID<-meta_data[meta_data$Time_point=="D1"&meta_data$Sample_Group=="placebo 24g bid",]$UniqueID
PLA_W8_ID<-meta_data[meta_data$Time_point=="W8"&meta_data$Sample_Group=="placebo 24g bid",]$UniqueID
PLA_W16_ID<-meta_data[meta_data$Time_point=="W16"&meta_data$Sample_Group=="placebo 24g bid",]$UniqueID

Dis_Group_1<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_D1_ID,GroupID2=AXA_W8_ID)
Dis_Group_2<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_D1_ID,GroupID2=AXA_W16_ID)
Dis_Group_3<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_D1_ID,GroupID2=PLA_D1_ID)
Dis_Group_4<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_D1_ID,GroupID2=PLA_W8_ID)
Dis_Group_5<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_D1_ID,GroupID2=PLA_W16_ID)
Dis_Group_6<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_W8_ID,GroupID2=AXA_W16_ID)
Dis_Group_7<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_W8_ID,GroupID2=PLA_D1_ID)
Dis_Group_8<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_W8_ID,GroupID2=PLA_W8_ID)
Dis_Group_9<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_W8_ID,GroupID2=PLA_W16_ID)
Dis_Group_10<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_W16_ID,GroupID2=PLA_D1_ID)
Dis_Group_11<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_W16_ID,GroupID2=PLA_W8_ID)
Dis_Group_12<-Group_Distance(disMatrix = disMatrix,GroupID1=AXA_W16_ID,GroupID2=PLA_W16_ID)
Dis_Group_13<-Group_Distance(disMatrix = disMatrix,GroupID1=PLA_D1_ID,GroupID2=PLA_W8_ID)
Dis_Group_14<-Group_Distance(disMatrix = disMatrix,GroupID1=PLA_D1_ID,GroupID2=PLA_W16_ID)
Dis_Group_15<-Group_Distance(disMatrix = disMatrix,GroupID1=PLA_W8_ID,GroupID2=PLA_W16_ID)

NullDis_1<-DisPerm(data = dt,GroupID1=AXA_D1_ID,GroupID2=AXA_W8_ID, PermTime = 1000)
NullDis_2<-DisPerm(data = dt,GroupID1=AXA_D1_ID,GroupID2=AXA_W16_ID, PermTime = 1000)
NullDis_3<-DisPerm(data = dt,GroupID1=AXA_D1_ID,GroupID2=PLA_D1_ID, PermTime = 1000)
NullDis_4<-DisPerm(data = dt,GroupID1=AXA_D1_ID,GroupID2=PLA_W8_ID, PermTime = 1000)
NullDis_5<-DisPerm(data = dt,GroupID1=AXA_D1_ID,GroupID2=PLA_W16_ID, PermTime = 1000)
NullDis_6<-DisPerm(data = dt,GroupID1=AXA_W8_ID,GroupID2=AXA_W16_ID, PermTime = 1000)
NullDis_7<-DisPerm(data = dt,GroupID1=AXA_W8_ID,GroupID2=PLA_D1_ID, PermTime = 1000)
NullDis_8<-DisPerm(data = dt,GroupID1=AXA_W8_ID,GroupID2=PLA_W8_ID, PermTime = 1000)
NullDis_9<-DisPerm(data = dt,GroupID1=AXA_W8_ID,GroupID2=PLA_W16_ID, PermTime = 1000)
NullDis_10<-DisPerm(data = dt,GroupID1=AXA_W16_ID,GroupID2=PLA_D1_ID, PermTime = 1000)
NullDis_11<-DisPerm(data = dt,GroupID1=AXA_W16_ID,GroupID2=PLA_W8_ID, PermTime = 1000)
NullDis_12<-DisPerm(data = dt,GroupID1=AXA_W16_ID,GroupID2=PLA_W16_ID, PermTime = 1000)
NullDis_13<-DisPerm(data = dt,GroupID1=PLA_D1_ID,GroupID2=PLA_W8_ID, PermTime = 1000)
NullDis_14<-DisPerm(data = dt,GroupID1=PLA_D1_ID,GroupID2=PLA_W16_ID, PermTime = 1000)
NullDis_15<-DisPerm(data = dt,GroupID1=PLA_W8_ID,GroupID2=PLA_W16_ID, PermTime = 1000)

P_value_1<-sum(NullDis_1>=Dis_Group_1)/length(NullDis_1) 
P_value_2<-sum(NullDis_2>=Dis_Group_2)/length(NullDis_2) 
P_value_3<-sum(NullDis_3>=Dis_Group_3)/length(NullDis_3) 
P_value_4<-sum(NullDis_4>=Dis_Group_4)/length(NullDis_4) 
P_value_5<-sum(NullDis_5>=Dis_Group_5)/length(NullDis_5) 
P_value_6<-sum(NullDis_6>=Dis_Group_6)/length(NullDis_6) 
P_value_7<-sum(NullDis_7>=Dis_Group_7)/length(NullDis_7) 
P_value_8<-sum(NullDis_8>=Dis_Group_8)/length(NullDis_8) 
P_value_9<-sum(NullDis_9>=Dis_Group_9)/length(NullDis_9) 
P_value_10<-sum(NullDis_10>=Dis_Group_10)/length(NullDis_10) 
P_value_11<-sum(NullDis_11>=Dis_Group_11)/length(NullDis_11) 
P_value_12<-sum(NullDis_12>=Dis_Group_12)/length(NullDis_12) 
P_value_13<-sum(NullDis_13>=Dis_Group_13)/length(NullDis_13) 
P_value_14<-sum(NullDis_14>=Dis_Group_14)/length(NullDis_14) 
P_value_15<-sum(NullDis_15>=Dis_Group_15)/length(NullDis_15) 

P_value=data.frame(Group=c("AXA D1 vs W8","AXA D1 vs W16","AXA D1 vs PLA D1",
                           "AXA D1 vs PLA W8","AXA D1 vs PLA W16","AXA W8 vs AXA W16",
                           "AXA W8 vs PLA D1","AXA W8 vs PLA W8","AXA W8 vs PLA W16",
                           "AXA W16 vs PLA D1","AXA W16 vs PLA W8","AXA W16 vs PLA W16",
                           "PLA D1 vs W8","PLA D1 vs W16","PLA W8 vs W16"),
                   Group.Distance=c(Dis_Group_1,Dis_Group_2,Dis_Group_3,Dis_Group_4,
                                    Dis_Group_5,Dis_Group_6,Dis_Group_7,
                                    Dis_Group_8,Dis_Group_9,Dis_Group_10,
                                    Dis_Group_11,Dis_Group_12,Dis_Group_13,Dis_Group_14,
                                    Dis_Group_15),
                   p.value=c(P_value_1,P_value_2,P_value_3,P_value_4,P_value_5,
                             P_value_6,P_value_7,P_value_8,P_value_9,P_value_10,
                             P_value_11,P_value_12,P_value_13,P_value_14,P_value_15))
write.table(P_value,file="output_file/Allsample_PCA_GroupDistance.txt",sep="\t")

