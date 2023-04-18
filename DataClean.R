# Packages Install

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install("readat")
BiocManager::install("limma")
library(devtools)
#install_bitbucket("graumannlabtools/readat")
library(readat)
library(reshape2)
library(dplyr)
library(tidyr)

#Read Data

setwd("/Users/hwang/Desktop/Somalogic Data/Somanic")
adatFile1<-"SOMAscan_EDTA_Data/SS-217128_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
plasma1 <- readAdat(adatFile1)
meta_data<-read.csv("SOMAscan_NASH_Data/Sample Submission Form_Axcella_Russell_NASH_EDTA.csv")
# Removing 11 sequences that failed QC.

#Statistical Sample Summary

table(meta_data$Time_point,meta_data$Sample_Group)


longPlasma1 <- melt(plasma1)
(interestingSeqs <-
    getSequencesWithLargestBetweenGroupVariation(
      longPlasma1, n = 2)[, .(SeqId, Target)])

(interestingSeqs <-
    getSequencesWithLargestBetweenGroupVariation(
      longPlasma2, n = 2)[, .(SeqId, Target)])

# PKPD Dataset Preview Adding for Metadata
pkdd_data<-read.csv("SOMAscan_Nash_Data/AXA1125_003 PKPD data.csv")
index<-match(meta_data$SubjectID,pkdd_data$SUBJID)
pkdd_data<-pkdd_data[index,]
meta_data$T2DIABFL<-pkdd_data$T2DIABFL
write.csv(meta_data,"SOMAscan_NASH_Data/Sample Submission Form_Axcella_Russell_NASH_EDTA.csv")

pkdd_data<-read.csv("SOMAscan_Nash_Data/AXA1125_003 PKPD data.csv")
index<-which(pkdd_data$SUBJID%in%meta_data$SubjectID)
pkdd_data<-pkdd_data[index,]
pkdd_data_wide<-pkdd_data[,c("SUBJID",colnames(pkdd_data)[grep("ALT|HOMAIR|PDFF|CT1|HBA1C|PROC3",colnames(pkdd_data))])]

sample<-unique(pkdd_data$SUBJID)
pkdd_data_long<-lapply(c("B","WK8","WK16"), function(x){
  tmp<-data.frame(SUBJID=sample,TIMEPOINT=x,pkdd_data_wide[match(sample,pkdd_data_wide$SUBJID),match(paste0(x,c("ALT","HOMAIR","PDFF","CT1","HBA1C","PROC3")),colnames(pkdd_data_wide))])
  colnames(tmp)[c(-1:-2)]<-c("ALT","HOMAIR","PDFF","CT1","HBA1C","PROC3")
  return(tmp)
})
pkdd_data_long<-rbind(pkdd_data_long[[1]],pkdd_data_long[[2]],pkdd_data_long[[3]])
pkdd_data_long[which(pkdd_data_long$TIMEPOINT=="B"),]$TIMEPOINT<-"D1"
pkdd_data_long[which(pkdd_data_long$TIMEPOINT=="WK8"),]$TIMEPOINT<-"W8"
pkdd_data_long[which(pkdd_data_long$TIMEPOINT=="WK16"),]$TIMEPOINT<-"W16"

pkdd_data_long<-rename(pkdd_data_long,c("Time_point"="TIMEPOINT"))
pkdd_data_long<-rename(pkdd_data_long,c("SubjectID"="SUBJID"))

meta_full<-full_join(meta_data, pkdd_data_long, by=c("Time_point","SubjectID"))

write.csv(meta_full,"output_file/Nash_EDTA_Metadata_clean_HW.csv")

#Absolute Change
pkdd_data<-read.csv("SOMAscan_Nash_Data/AXA1125_003 PKPD data.csv")

pkdd_data$delta_PDFF_W8<-pkdd_data$WK8PDFF-pkdd_data$BPDFF
pkdd_data$delta_PDFF_W16<-pkdd_data$WK16PDFF-pkdd_data$BPDFF

pkdd_data$delta_ALT_W8<-pkdd_data$WK8ALT-pkdd_data$BALT
pkdd_data$delta_ALT_W16<-pkdd_data$WK16ALT-pkdd_data$BALT

pkdd_data$delta_CT1_W8<-pkdd_data$WK8CT1-pkdd_data$BCT1
pkdd_data$delta_CT1_W16<-pkdd_data$WK16CT1-pkdd_data$BCT1

pkdd_data$delta_HOMAIR_W8<-pkdd_data$WK8HOMAIR-pkdd_data$BHOMAIR
pkdd_data$delta_HOMAIR_W16<-pkdd_data$WK16HOMAIR-pkdd_data$BHOMAIR

pkdd_data$delta_HBA1C_W8<-pkdd_data$WK8HBA1C-pkdd_data$BHBA1C
pkdd_data$delta_HBA1C_W16<-pkdd_data$WK16HBA1C-pkdd_data$BHBA1C

pkdd_data$delta_PROC3_W8<-pkdd_data$WK8PROC3-pkdd_data$BPROC3
pkdd_data$delta_PROC3_W16<-pkdd_data$WK16PROC3-pkdd_data$BPROC3

write.csv(pkdd_data,"output_file/AXA1125_003_PKPD_Clean_HW.csv")

#
meta_data<-read.csv("output_file/Nash_EDTA_Metadata_clean_HW.csv")

pkdd_data$PDFF_W8<-ifelse(pkdd_data$delta_PDFF_W8/pkdd_data$BPDFF<=-0.3,1,0)
pkdd_data$PDFF_W16<-ifelse(pkdd_data$delta_PDFF_W16/pkdd_data$BPDFF<=-0.3,1,0)

pkdd_data$ALT_W8<-ifelse(pkdd_data$delta_ALT_W8<=-17,1,0)
pkdd_data$ALT_W16<-ifelse(pkdd_data$delta_ALT_W16<=-17,1,0)

pkdd_data$CT1_W8<-ifelse(pkdd_data$delta_CT1_W8<=-80,1,0)
pkdd_data$CT1_W16<-ifelse(pkdd_data$delta_CT1_W16<=-80,1,0)

pkdd_data<-pkdd_data[match(meta$SubjectID,pkdd_data$SUBJID),]

meta$classifier_pdff_W8<-pkdd_data$PDFF_W8
meta$classifier_pdff_W16<-pkdd_data$PDFF_W16

meta$classifier_alt_W8<-pkdd_data$ALT_W8
meta$classifier_alt_W16<-pkdd_data$ALT_W16

meta$classifier_ct1_W8<-pkdd_data$CT1_W8
meta$classifier_ct1_W16<-pkdd_data$CT1_W16

write.csv(meta,"output_file/Nash_EDTA_Metadata_clean_HW.csv")
