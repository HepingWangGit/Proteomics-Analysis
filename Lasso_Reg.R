library(glmnet)


lasso_regression<-function(x,
                           y,
                           nfold = 5,
                           alphalist = seq(0,1,by=0.1),
                           Annotation
){
  result_list<-list()
  alphalist<-alphalist #optimal alpha
  foldid<-sample(1:nfold,size=length(y),replace=TRUE) #pre-defined foldid for reproducibility
  elasticnet <- lapply(alphalist, function(a){cv.glmnet(x, y, alpha=a,lambda.min.ratio=.001,nfolds = nfold,foldid = foldid)})
  min_cvm<-c() 
  for (i in 1:11) {
    #print(min(elasticnet[[i]]$cvm))
    min_cvm<-c(min_cvm,min(elasticnet[[i]]$cvm))
  }
  result_list[["min_cvm"]]<-min(min_cvm)
  
  optimal_alpha<-alphalist[which.min(min_cvm)]
  result_list[["optimal_alpha"]]<-optimal_alpha
  
  cv_model <- cv.glmnet(x, y, alpha = optimal_alpha, nfolds = 5,foldid   = foldid)
  best_lambda <- cv_model$lambda.min
  result_list[["best_lambda"]] <-best_lambda 
  g1<-plot(cv_model)#produce plot of test MSE by lambda value
  result_list[["MSE_CV"]] <-g1
  best_model <- glmnet(Diff_exp, PDFF, alpha = optimal_alpha, 
                       lambda = best_lambda)
  result_list[["best_model"]]<-best_model
  coefficient<-coef(best_model)#find coefficients of best model
  coefficient<-coefficient[coefficient[,1]!=0,]
  SeqID=names(coefficient)
  SeqID<-gsub("^.+\\.","",SeqID,perl=T)
  coefficient<-data.frame(SeqId=SeqID,Coefficient=coefficient)
  coefficient_Annoted<-merge(coefficient,
                             Annotation,by="SeqId",
                             all.x=T)
  result_list[["coefficient"]]<-coefficient_Annoted
  return(result_list)
}
