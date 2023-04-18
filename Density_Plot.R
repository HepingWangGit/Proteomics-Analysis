# Prep Library

library(readat)
library(reshape2)
library(limma)
library(dplyr)
library(tidyr)
library(stringr)

#Read Data

setwd("/Users/hwang/Desktop/Somanic")
#adatFile1<-"SOMAscan_NASH_Data/SS-217128_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.adat"
adatFile1<-"SOMAscan_Nash_Data/SS-217128_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
plasma1 <- readAdat(adatFile1)
meta_data<-read.csv("SOMAscan_NASH_Data/Sample Submission Form_Axcella_Russell_NASH_EDTA.csv")

Annotation<-attr(plasma1,"SequenceData")
write.table(Annotation,file="output_file/Annotation.txt",quote=F,sep="\t",row.names=F)
#Data transformation

longPlasma1 <- melt(plasma1)
(interestingSeqs <-
    getSequencesWithLargestBetweenGroupVariation(
      longPlasma1, n = 2))

#Orig Data Density Plot

df<-data.frame(SeqID=longPlasma1$SeqId,SampleID<-longPlasma1$SampleId,Intensity=log2(1+longPlasma1$Intensity))
write.table(df,"output_file/EDTA_log.txt",row.names = F)

g <- ggplot(df, aes(x=Intensity, fill=SampleID))+
  geom_density(alpha=.01)+
  theme(legend.position = "none")+
  labs(title="SomaLogic Normalized data density plot",x="log2 transformed values",y="density",subtitle="one density plot per sample")+
  theme(text = element_text(size=14))
ggsave("Figure/somaLogic.EDTA.log.density.jpg",g,device="jpeg")


#Limma quantile normalization

df<-log2(1+plasma1[,-c(1:33)])
df<-t(df)
colnames(df)<-plasma1$SampleId
m <- as.matrix(df)
m <- normalizeQuantiles(m)
write.table(m,"output_file/somaLogic.EDTA.limmaNorm.txt", quote=F, sep="\t", col.names = T,row.names = T)

#LIMMANORM Data Density Plot
m<-read.table("output_file/somaLogic.EDTA.limmaNorm.txt")
df<-as.data.frame(m)
dt <- melt(df)
g <- ggplot(dt, aes(x=value, fill=variable)) +
  geom_density(alpha=.01)+
  theme(legend.position = "none")+
  labs(title="After log transformation and quantile normalization",x="normalized log2 transformed values",y="density",subtitle="one density plot per sample")+
  theme(text = element_text(size=14))
ggsave("Figure/somaLogic.EDTAnorm.log.density.jpg",g,device="jpeg")
