library(caret)
library(ggplot2)

getBLModel <- function(filepath,meth,K,reps){
  
  dat <- read.csv(file=filepath,row.names = 1)
  train_control <- NULL
  if (meth == "LOOCV"){
    train_control <- trainControl(method = meth, 
                                  number = 1,
                                  savePredictions = "final")
  }else{
    if (meth == "boot" | meth == "boot632"){
      train_control <- trainControl(method = meth, 
                                    number = K,
                                    savePredictions = "final")
    }else{
      train_control <- trainControl(method = meth, 
                                    number = K, repeats=reps,
                                    savePredictions = "final")
    }
  }
  model <- train(form = as.factor(PDFF_W16) ~ BMI+ +SPA4L+CDCP1+`DC.SIGN`,
                 data = dat,
                 trControl = train_control,
                 method = "glm",
                 family = "binomial",
                 metric = "Kappa")
  model
}

saveStatsAndModel <- function(s, method){
  
  file <- paste0("output_file/W16_W8_Valid/",method,".coef.res")
  nms <- names(s$coef)
  vals <- unlist(s$coef)
  mat <- matrix(nrow=1,ncol=length(nms))
  for (i in 1:length(nms)){
    mat[1,i] <- vals[i]
  }
  df <- as.data.frame(mat)
  colnames(df) <- nms
  write.table(df,file,quote=F,sep="\t",row.names=F)
  df <- s$df
  file <- paste0("output_file/W16_W8_Valid/",method,".stat.res")
  write.table(df,file,quote=F,sep="\t",row.names=F)
}

saveAllStats <- function(nsim){
  
  mod <- getBLModel("output_file/PreData_cPDFF_W16.csv","LOOCV",1,1)
  s <- getModelStats(mod)
  saveStatsAndModel(s,"LOOCV")
  
  mod <- getBLModel("output_file/PreData_cPDFF_W16.csv","boot",nsim,1)
  s <- getModelStats(mod)
  saveStatsAndModel(s,"boot")
  
  mod <- getBLModel("output_file/PreData_cPDFF_W16.csv","boot632",nsim,1)
  s <- getModelStats(mod)
  saveStatsAndModel(s,"boot632")
  
  mod <- getBLModel("output_file/PreData_cPDFF_W16.csv","repeatedcv",3,nsim)
  s <- getModelStats(mod)
  saveStatsAndModel(s,"repcvk3")
  
  mod <- getBLModel("output_file/PreData_cPDFF_W16.csv","repeatedcv",4,nsim)
  s <- getModelStats(mod)
  saveStatsAndModel(s,"repcvk4")
}


getStatsFromCM <- function(tb){
  
  TP <- tb[1,1]
  FP <- tb[1,2]
  FN <- tb[2,1]
  TN <- tb[2,2]
  P <- TP+FN
  N <- TN+FP
  PP <- TP+FP
  PN <- TN+FN
  n <- P+N
  SN <- TP/P
  SP <- TN/N
  OA <- (TP+TN)/n
  EA <- (N/n)*(PN/n)+(P/n)*(PP/n)
  KP <- (OA-EA)/(1-EA)
  PR <- TP/(TP+FP)
  NPV <- TN/(TN+FN)
  
  df <- data.frame(TP,FP,FN,TN,KP,OA,EA,SN,SP,PR,NPV)
  stat <- colnames(df)
  mn <- apply(df,2,mean)
  mn <- signif(mn,digits = 2)
  stdv <- apply(df,2,sd)
  df <- data.frame(stat,mn,stdv)
  colnames(df) <- c("stat","mean","sd")
  rownames(df) <- NULL
  df
}

getStatsFromCMs <- function(t){
  
  TP <- t[,1]
  FP <- t[,3]
  FN <- t[,2]
  TN <- t[,4]
  P <- TP+FN
  N <- TN+FP
  PP <- TP+FP
  PN <- TN+FN
  n <- P+N
  SN <- TP/P
  SP <- TN/N
  OA <- (TP+TN)/n
  EA <- (N/n)*(PN/n)+(P/n)*(PP/n)
  KP <- (OA-EA)/(1-EA)
  PR <- TP/(TP+FP)
  NPV <- TN/(TN+FN)
  
  df <- data.frame(TP,FP,FN,TN,KP,OA,EA,SN,SP,PR,NPV)
  stat <- colnames(df)
  mn <- apply(df,2,mean,na.rm=TRUE)
  mn <- signif(mn,digits=2)
  stdv <- apply(df,2,sd,na.rm=TRUE)
  stdv <- signif(stdv,digits=2)
  df <- data.frame(stat,mn,stdv)
  colnames(df) <- c("stat","mean","sd")
  rownames(df) <- NULL
  df
}

getModelStats <- function(mod){
  
  res <- list()
  res[["coef"]] <- mod$finalModel$coefficients
  if (mod$control$method == "LOOCV"){
    pred <- mod$pred$pred
    obs <- mod$pred$obs
    tb <- table(pred,obs)
    df <- getStatsFromCM(tb)
    res[["df"]] <- df
  }else{
    df <- getStatsFromCMs(mod$resampledCM)
    res[["df"]] <- df
  }
  res
}
