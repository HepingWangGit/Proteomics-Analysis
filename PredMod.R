library(caret)
library(cluster)

#   Functions 

getFormula <- function(prots){
  s <- "as.factor(PDFF_W16) ~ `cli_data$BBMI`+ "
  matrix<-t(combn(c(1:length(prots)),2))
  result<-apply(matrix,1,function(aa){
  prot1<-prots[aa[1]]
  prot2<-prots[aa[2]]
  s2<-paste0(prot1," + ",prot2, collapse = " ")
})
 s <- paste0(s,result, collapse = " ")
  as.formula(s)
}


getModel <- function(dat,prots,meth,K,reps){
  
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
  fm <- getFormula(prots)
  model <- train(form = fm,
                 data = dat,
                 trControl = train_control,
                 method = "glm",
                 family = "binomial",
                 metric = "Kappa")
  model
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

getSimpleLogisticRegression <- function(dat,prot){
  
  form <- getFormula(prot)
  mod <- glm(formula=form,data=dat,family = binomial)
  mod
}

saveSimpleLogisticRegressions <- function(){
  
  df <- getFullDataDF()
  prots <- colnames(df)
  prots <- prots[-1]
  dat <- loadDFForProtsAtStage(prots,"AC")
  pv.AC <- NULL
  sd.AC <- NULL
  mean.AC <- NULL
  cv.AC <- NULL
  for (prot in prots){
    mod <- getSimpleLogisticRegression(dat,prot)
    s <- summary(mod)
    p <- s$coefficients[2,4]
    pv.AC <- c(pv.AC,signif(p,digits=2))
    stdev <- sd(dat[,prot])
    mn <- mean(dat[,prot])
    cv <- stdev/mn
    sd.AC <- c(sd.AC,signif(stdev,digits=2))
    mean.AC <- c(mean.AC,signif(mn,digits=2))
    cv.AC <- c(cv.AC,signif(cv,digits=2))
  }
  fdr.AC <- p.adjust(pv.AC,method="BH")
  dat <- loadDFForProtsAtStage(prots,"BL")
  pv.BL <- NULL
  for (prot in prots){
    mod <- getSimpleLogisticRegression(dat,prot)
    s <- summary(mod)
    p <- s$coefficients[2,4]
    pv.BL <- c(pv.BL,signif(p,digits=2))
  }
  fdr.BL <- p.adjust(pv.BL,method="BH")
  pv.min <- pmin(pv.AC,pv.BL)
  pv.max <- pmax(pv.AC,pv.BL)
  df <- data.frame(prots,pv.AC,pv.BL,pv.min,pv.max,fdr.AC,fdr.BL)
  df <- df[order(df$pv.max),]
  write.table(df,"res/SimpleLogisticRegression.res",sep="\t",quote=F,row.names=F)
}


getTPV <- function(dat,prot){
  
  y <- dat[,prot]
  fat <- dat[,"fatigue"]
  y1 <- y[which(fat==1)]
  y0 <- y[which(fat==0)]
  tt <- t.test(y1,y0,var.equal=T)
  pv <- signif(tt$p.value,digits=2)
  pv
}

saveTTests <- function(){
  
  df <- getFullDataDF()
  prots <- colnames(df)
  prots <- prots[-1]
  dat <- loadDFForProtsAtStage(prots,"AC")
  pv.AC <- NULL
  sd.AC <- NULL
  for (prot in prots){
    p <- getTPV(dat,prot)
    v <- dat[,prot]
    stdev <- sd(v)
    pv.AC <- c(pv.AC,signif(p,digits=2))
    sd.AC <- c(sd.AC,signif(stdev,digits=2))
  }
  fdr.AC <- p.adjust(pv.AC,method="BH")
  dat <- loadDFForProtsAtStage(prots,"BL")
  pv.BL <- NULL
  for (prot in prots){
    p <- getTPV(dat,prot)
    pv.BL <- c(pv.BL,signif(p,digits=2))
  }
  fdr.BL <- p.adjust(pv.BL,method="BH")
  dat <- loadDFForProtsAtStage(prots,"CV")
  pv.CV <- NULL
  for (prot in prots){
    p <- getTPV(dat,prot)
    pv.CV <- c(pv.CV,signif(p,digits=2))
  }
  fdr.CV <- p.adjust(pv.CV,method="BH")
  pv.max <- pmax(pv.AC,pv.BL)
  pv.max <- pmax(pv.max,pv.CV)
  df <- data.frame(prots,pv.max,pv.AC,pv.BL,pv.CV,fdr.AC,fdr.BL,fdr.CV)
  df <- df[order(df$pv.max),]
  write.table(df,"res/TTest.res",sep="\t",quote=F,row.names=F)
}

stepwiseRegression <- function(stage,thresh,outF){
  
  df <- read.table("res/SimpleLogisticRegression.res",head=T,sep="\t",quote="")
  cnm <- paste0("pv.",stage)
  pv <- df[,cnm]
  tokeep <- which(pv < thresh)
  prots <- df$prots[tokeep]
  dat <- loadDFForProtsAtStage(prots,stage)
  meth <- "repeatedcv"
  K <- 5
  reps <- 100
  bestProt <- NULL
  selected <- NULL
  scores <- NULL
  for (j in 1:10){
    bestKP <- -10
    for (i in 1:length(prots)){
      prot <- prots[i]
      selected <- c(selected,prot)
      mod <- getModel(dat,selected,meth,K,reps)
      s <- getModelStats(mod)
      KP = s$df[s$df$stat=="KP","mean"]
      if (KP > bestKP){
        bestKP <- KP
        bestProt <- prot
      }
      print(bestKP)
      selected <- selected[-which(selected==prot)]
    }
    selected <- c(selected,bestProt)
    prots <- prots[-which(prots==bestProt)]
    scores <- c(scores,bestKP)
  }
  df <- data.frame(selected,scores)
  write.table(df,outF,sep="\t",quote=F,row.names=F)
}

regressionPairs <- function(stage,thresh,outF){
  
  df <- read.table("res/TTest.res",head=T,sep="\t",quote="")
  cnm <- paste0("pv.",stage)
  pv <- df[,cnm]
  tokeep <- which(pv < thresh)
  prots <- df[,"prots"]
  prots <- prots[tokeep]
  dat <- loadDFForProtsAtStage(prots,stage)
  meth <- "repeatedcv"
  K <- 5
  reps <- 10
  print(length(prots))
  n <- length(prots)*(length(prots)-1)/2
  c <- 0
  prot1 <- NULL
  prot2 <- NULL
  kappa <- NULL
  bestKP <- -1
  for (i in 1:(length(prots)-1)){
    proti <- prots[i]
    for (j in (i+1):length(prots)){
      c <- c+1
      protj <- prots[j]
      pair <- c(proti,protj)
      mod <- getModel(dat,pair,meth,K,reps)
      s <- getModelStats(mod)
      KP = s$df[s$df$stat=="KP","mean"]
      if (KP > bestKP){
        bestKP <- KP
      }
      print(paste0(c," / ",n,", ",bestKP))
      prot1 <- c(prot1,proti)
      prot2 <- c(prot2,protj)
      kappa <- c(kappa,KP)
    }
  }
  df <- data.frame(prot1,prot2,kappa)
  write.table(df,outF,quote=F,sep="\t",row.names=F)
}

regressionTriplets <- function(stage,thresh,outF){
  
  df <- read.table("res/TTest.res",head=T,sep="\t",quote="")
  cnm <- paste0("pv.",stage)
  pv <- df[,cnm]
  tokeep <- which(pv < thresh)
  prots <- df[,"prots"]
  prots <- prots[tokeep]
  dat <- loadDFForProtsAtStage(prots,stage)
  meth <- "repeatedcv"
  K <- 5
  reps <- 10
  n <- length(prots)*(length(prots)-1)*(length(prots)-2)/6
  c <- 0
  prot1 <- NULL
  prot2 <- NULL
  prot3 <- NULL
  kappa <- NULL
  bestKP <- -1
  for (i in 1:(length(prots)-2)){
    proti <- prots[i]
    for (j in (i+1):(length(prots)-1)){
      protj <- prots[j]
      for (k in (j+1):length(prots)){
        protk <- prots[k]
        triplet <- c(proti,protj,protk)
        c <- c+1
        mod <- getModel(dat,triplet,meth,K,reps)
        s <- getModelStats(mod)
        KP = s$df[s$df$stat=="KP","mean"]
        if (KP > bestKP){
          bestKP <- KP
        }
        print(paste0(c," / ",n,", ",bestKP))
        prot1 <- c(prot1,proti)
        prot2 <- c(prot2,protj)
        prot3 <- c(prot3,protk)
        kappa <- c(kappa,KP)
      }
    }
  }
  df <- data.frame(prot1,prot2,prot3,kappa)
  write.table(df,outF,quote=F,sep="\t",row.names=F)
}

regressionQuadruplets <- function(stage,thresh,outF){
  
  df <- read.table("res/TTest.res",head=T,sep="\t",quote="")
  cnm <- paste0("pv.",stage)
  pv <- df[,cnm]
  tokeep <- which(pv < thresh)
  prots <- df[,"prots"]
  prots <- prots[tokeep]
  dat <- loadDFForProtsAtStage(prots,stage)
  meth <- "repeatedcv"
  K <- 5
  reps <- 4
  n <- length(prots)*(length(prots)-1)*(length(prots)-2)*(length(prots)-3)/24
  c <- 0
  prot1 <- NULL
  prot2 <- NULL
  prot3 <- NULL
  prot4 <- NULL
  kappa <- NULL
  bestKP <- -1
  for (i in 1:(length(prots)-3)){
    proti <- prots[i]
    for (j in (i+1):(length(prots)-2)){
      protj <- prots[j]
      for (k in (j+1):(length(prots)-1)){
        protk <- prots[k]
        for (l in (k+1):length(prots)){
          protl <- prots[l]
          quadruplet <- c(proti,protj,protk,protl)
          c <- c+1
          mod <- getModel(dat,quadruplet,meth,K,reps)
          s <- getModelStats(mod)
          KP = s$df[s$df$stat=="KP","mean"]
          if (KP > bestKP){
            bestKP <- KP
          }
          prot1 <- c(prot1,proti)
          prot2 <- c(prot2,protj)
          prot3 <- c(prot3,protk)
          prot4 <- c(prot4,protl)
          kappa <- c(kappa,KP)
        }
      }
      print(paste0(c," / ",n,", ",bestKP))
    }
  }
  df <- data.frame(prot1,prot2,prot3,kappa)
  write.table(df,outF,quote=F,sep="\t",row.names=F)
}

medoidSummary <- function(cm, n){
  
  pm <- pam(cm,n,diss=TRUE)
  res <- list()
  meds <- pm$medoids
  res[["meds"]] <- meds
  prot2cluster <- as.list(pm$clustering)
  names(prot2cluster) <- names(pm$clustering)
  res[["prot2cluster"]] <- prot2cluster
  cluster <- pm$clustering
  prot <- names(pm$clustering)
  df <- data.frame(cluster,prot)
  cluster2prots <- split(df$prot,df$cluster)
  names(cluster2prots) <- meds
  res[["cluster2prots"]] <- cluster2prots
  prot2med <- lapply(prot2cluster,function(x){meds[x]})
  res[["prot2med"]] <- prot2med
  df <- as.data.frame(pm$clusinfo)
  id <- c(1:length(meds))
  df <- cbind(id,meds,df)
  res[["info"]] <- df
  res
}

goDown <- function(stage,prots){
  
  dat <- getFullDFatStage(stage)
  mod <- getModel(dat,prots,"repeatedcv",5,1000)
  s <- getModelStats(mod)
  KP = s$df[s$df$stat=="KP","mean"]
  print(KP)
  for (i in 1:length(prots)){
    sel <- prots
    torm <- which(sel == prots[i])
    sel <- sel[-torm]
    mod <- getModel(dat,sel,"repeatedcv",5,1000)
    s <- getModelStats(mod)
    KP = s$df[s$df$stat=="KP","mean"]
    print(KP)
  }
}

temp1 <- function(){
  dat <- getFullDFatStage("BL")
  prots <- c("SIT1_IRE", "PDGFB_CVD2", "PRKCQ_IRE", "NPTXR_MET",
             "NCR1_IRE", "DDX58_IRE", "CAPG_ODA", "NTF4_IRE",
             "IL8_INF", "BIRC2_IRE")
  prots <- prots[-10]
  prots <- prots[-3]
  prots <- prots[-7]
  # mod <- getModel(dat,prots,"repeatedcv",10,1000)
  # s <- getModelStats(mod)
  goDown("BL",prots)
}


orthogonalSelection <- function(stage,thresh,n){
  
  dat <- getFullDFatStage("BL")
  df <- read.table("res/TTest.res",head=T,sep="\t",quote="")
  t <- paste0("pv.",stage)
  pv <- df[,t]
  tokeep <- which(pv < thresh)
  df <- df[tokeep,]
  prots <- df[,"prots"]
  fatigue <- dat[,"fatigue"]
  dat <- dat[,-c(1,2,3)]
  tokeep <- which(colnames(dat) %in% prots)
  dat <- dat[,tokeep]
  prots <- colnames(dat)
  cfr <- NULL
  for (i in 1:length(prots)){
    sct <- cor.test(fatigue,dat[,prots[i]],method="spearman")
    cfr <- c(cfr,abs(sct$estimate))
  }
  cm <- cor(dat,method="spearman")
  cm <- abs(cm)
  m <- min(dim(cm))
  
  getAvFFCor <- function(x){
    
    u <- combn(x,2)
    res <- 0
    for (j in 1:length(u[1,])){
      k <- u[1,j]
      l <- u[2,j]
      res <- res + cm[k,l]
    }
    res <- res / choose(length(x),2)
    res
  }
  
  getAvFRCor <- function(x){
    
    res <- mean(cfr[x])
    res
  }
  
  cbn <- combn(c(1:m),n)
  N <- length(cbn[1,])
  besti <- 0
  bestf <- 0
  for (i in 1:N){
    v <- cbn[,i]
    num <- getAvFRCor(v)
    denom <- getAvFFCor(v)*(n+n*(n-1))
    denom <- sqrt(denom)
    f <- num/denom
    if (f > bestf){
      bestf <- f
      besti <- i
    }
  }
  res <- NULL
  v <- cbn[,besti]
  for (j in v){
    res <- c(res,colnames(cm)[j])
  }
  print(bestf)
  res
}

univariatePredictor <- function(stage,thresh){
  
  dat <- getFullDFatStage(stage)
  dftt <- read.table("res/TTest.res",head=T,sep="\t",quote="")
  clnm <- paste0("pv.",stage)
  pv <- dftt[,clnm]
  tokeep <- which(pv < thresh)
  dftt <- dftt[tokeep,]
  prots <- dftt[,"prots"]
  tt.p <- dftt[,clnm]
  slr.p <- NULL
  KP <- NULL
  OA <- NULL
  EA <- NULL
  SN <- NULL
  SP <- NULL
  PR <- NULL
  NPV <- NULL
  meth <- "repeatedcv"
  K <- 5
  reps <- 1000
  for (i in 1:length(prots)){
    prot <- prots[i]
    print(prot)
    slr <- getSimpleLogisticRegression(dat,prot)
    s <- summary(slr)
    slr.p <- c(slr.p,signif(s$coefficients[2,4]))
    mod <- getModel(dat,prot,meth,K,reps)
    s <- getModelStats(mod)
    d <- s$df
    KP <- c(KP,signif(d[d$stat=="KP","mean"],digits=2))
    OA <- c(OA,signif(d[d$stat=="OA","mean"],digits=2))
    EA <- c(EA,signif(d[d$stat=="EA","mean"],digits=2))
    SN <- c(SN,signif(d[d$stat=="SN","mean"],digits=2))
    SP <- c(SP,signif(d[d$stat=="SP","mean"],digits=2))
    PR <- c(PR,signif(d[d$stat=="PR","mean"],digits=2))
    NPV <- c(NPV,signif(d[d$stat=="NPV","mean"],digits=2))
  }
  dfo <- data.frame(prots,slr.p,KP,OA,EA,SN,SP,PR,NPV)
  f <- paste0("res/UnivariatePred.",stage,".res")
  write.table(dfo,f,quote=F,sep="\t",row.names=F)
}

univariatePredictorS <- function(stage,thresh){
  
  dat <- getFullDFatStage(stage)
  dftt <- read.table("res/TTest.res",head=T,sep="\t",quote="")
  clnm <- paste0("pv.",stage)
  pv <- dftt[,clnm]
  tokeep <- which(pv < thresh)
  dftt <- dftt[tokeep,]
  prots <- dftt[,"prots"]
  tt.p <- dftt[,clnm]
  slr.p <- NULL
  KP <- NULL
  OA <- NULL
  EA <- NULL
  SN <- NULL
  SP <- NULL
  PR <- NULL
  NPV <- NULL
  meth <- "repeatedcv"
  K <- 5
  reps <- 100
  for (i in 1:length(prots)){
    prot <- prots[i]
    print(prot)
    slr <- getSimpleLogisticRegression(dat,prot)
    s <- summary(slr)
    slr.p <- c(slr.p,signif(s$coefficients[2,4]))
    mod <- getModel(dat,prot,meth,K,reps)
    s <- getModelStats(mod)
    d <- s$df
    KP <- c(KP,signif(d[d$stat=="KP","mean"],digits=2))
    OA <- c(OA,signif(d[d$stat=="OA","mean"],digits=2))
    EA <- c(EA,signif(d[d$stat=="EA","mean"],digits=2))
    SN <- c(SN,signif(d[d$stat=="SN","mean"],digits=2))
    SP <- c(SP,signif(d[d$stat=="SP","mean"],digits=2))
    PR <- c(PR,signif(d[d$stat=="PR","mean"],digits=2))
    NPV <- c(NPV,signif(d[d$stat=="NPV","mean"],digits=2))
  }
  dfo <- data.frame(prots,slr.p,KP,OA,EA,SN,SP,PR,NPV)
  f <- paste0("res/UnivariatePred.S.",stage,".res")
  write.table(dfo,f,quote=F,sep="\t",row.names=F)
}

refinePairs <- function(n){
  
  df <- read.table("res/Pairs.K5.T.AC.res",head=T,sep="\t",quote="")
  df <- df[order(df$kappa,decreasing=T),]
  t <- 0
  for (i in 1:(2*n)){
    s <- c(df$prot1[1:i],df$prot2[1:i])
    x <- length(unique(s))
    if (x >= n){
      t <- i
      break
    }
  }
  df <- df[1:t,]
  KP <- NULL
  OA <- NULL
  EA <- NULL
  SN <- NULL
  SP <- NULL
  PR <- NULL
  NPV <- NULL
  dat <- getFullDFatStage("AC")
  meth <- "repeatedcv"
  K <- 5
  reps <- 1000
  for (i in 1:t){
    prots <- c(df$prot1[i],df$prot2[i])
    print(prots)
    mod <- getModel(dat,prots,meth,K,reps)
    s <- getModelStats(mod)
    u <- s$df
    KP <- c(KP,signif(u$mean[5],digits=2))
    OA <- c(OA,signif(u$mean[6],digits=2))
    EA <- c(EA,signif(u$mean[7],digits=2))
    SN <- c(SN,signif(u$mean[8],digits=2))
    SP <- c(SP,signif(u$mean[9],digits=2))
    PR <- c(PR,signif(u$mean[10],digits=2))
    NPV <- c(NPV,signif(u$mean[11],digits=2))
  }
  prot1 <- df[,"prot1"]
  prot2 <- df[,"prot2"]
  df <- data.frame(prot1,prot2,KP,OA,EA,SN,SP,PR,NPV)
  write.table(df,"res/RefinedPairs.AC.res",sep="\t",quote=F,row.names=F)
}


refineTriplets <- function(n){
  
  df <- read.table("res/Triplets.K5.T01.AC.res",head=T,sep="\t",quote="")
  df <- df[order(df$kappa,decreasing=T),]
  t <- 0
  for (i in 1:(2*n)){
    s <- c(df$prot1[1:i],df$prot2[1:i],df$prot3[1:i])
    x <- length(unique(s))
    if (x >= n){
      t <- i
      break
    }
  }
  df <- df[1:t,]
  KP <- NULL
  OA <- NULL
  EA <- NULL
  SN <- NULL
  SP <- NULL
  PR <- NULL
  NPV <- NULL
  dat <- getFullDFatStage("AC")
  meth <- "repeatedcv"
  K <- 5
  reps <- 1000
  for (i in 1:t){
    prots <- c(df$prot1[i],df$prot2[i],df$prot3[i])
    print(prots)
    mod <- getModel(dat,prots,meth,K,reps)
    s <- getModelStats(mod)
    u <- s$df
    KP <- c(KP,signif(u$mean[5],digits=2))
    OA <- c(OA,signif(u$mean[6],digits=2))
    EA <- c(EA,signif(u$mean[7],digits=2))
    SN <- c(SN,signif(u$mean[8],digits=2))
    SP <- c(SP,signif(u$mean[9],digits=2))
    PR <- c(PR,signif(u$mean[10],digits=2))
    NPV <- c(NPV,signif(u$mean[11],digits=2))
  }
  prot1 <- df[,"prot1"]
  prot2 <- df[,"prot2"]
  prot3 <- df[,"prot3"]
  df <- data.frame(prot1,prot2,prot3,KP,OA,EA,SN,SP,PR,NPV)
  df <- df[order(df$KP,decreasing=T),]
  write.table(df,"res/RefinedTriplets.AC.res",sep="\t",quote=F,row.names=F)
}

# refineTriplets(20)

# refinePairsS(20)

# dat <- getFullDFatStage("AC")
# prots <- c("CDH2_MET","ANGPTL7_MET")
# mod <- getModel(dat,prots,"repeatedcv",5,100)
# s <- getModelStats(mod)
# print(s$df)

# univariatePredictor("AC",0.1)
# univariatePredictorS("AC",0.05)

#regressionPairs("AC",0.02,"res/Pairs.K5.AC.res")
# stepwiseRegression("BL",0.1,"res/SRK5.BL.res")