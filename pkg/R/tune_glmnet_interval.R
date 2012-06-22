tune.glmnet.interval<-function(parms, x, y,
                                         weights, 
                                         offset = NULL, 
                                         lambda = NULL, 
                                         type.measure = c("mse", "deviance", "class", "auc", "mae"),
                                         seed=12345, 
                                         nfolds = 10, 
                                         foldid=NULL, 
                                         grouped = TRUE, 
                                         type.min=c("lambda.min", "lambda.1se"),
                                         verbose=FALSE,
                                         ...){
  
  # 1. decode the parameters ############################################################
  
  print("parms.coding")
  print(parms.coding)
  
  alpha<-parms[1]
  names(alpha)<- NULL
  
  if (verbose) print(paste("alpha=",alpha))
  
  #  2. find optimal lambda for given alpha #######################################################################
  # for EPSGO, Problem = cv.glmnet , output mean misclassification error
  
   
  # find optimal lambda for given alpha
  set.seed(seed)
  cv<-cv.glmnet(x=x,y=y,family=family,  
                alpha=alpha,
                offset = NULL,
                lambda = NULL, 
                type.measure =type.measure,
                nfolds = nfolds, 
                foldid = foldid,
                grouped = grouped)
  
  
  opt.lambda<-ifelse(type.min=="lambda.min", cv$lambda.min, cv$lambda.1se )
  
  # q.val= mean cross-validated error over the folds
  q.val<-cv$cvm[which(cv$lambda == opt.lambda) ]
  
  #  3. fit  the model for given alpha and opt.lambda ########################################################
  # fit the model for given alpha and opt.lambda
  fit<-glmnet(x=x,y=y,family=family,
              alpha=alpha, 
              lambda=opt.lambda)
  
  #OR extract the coefs from the cv.glmnet object
  #cof <- coef(cv, s=opt.lambda)
  #names.cof <- rownames(cof)
  #cofn <- cof[which(cof!=0)]
  #names(cofn) <- names.cof[which(cof!=0)]
  
  #summary(coef(fit))
  
  #ret<-list(q.val=q.val, model=list(fit,cv)) 
  #ret<-list(q.val=q.val, model=cv) 
  ret<-list(q.val=q.val, model=fit)
  
  return(ret)
}

#stolen from MCREstimate package
# used for stratified(balanced) classification
my.balanced.folds <- function(class.column.factor, cross.outer)
{
  # get balanced folds from pamr
  sampleOfFolds  <- get("balanced.folds",en=asNamespace("pamr"))(class.column.factor, nfolds=cross.outer)
  permutated.cut <- rep(0,length(class.column.factor))
  for (sample in 1:cross.outer)
  {
    cat(sample,"\n")
    permutated.cut[sampleOfFolds[[sample]]] <- sample
  }
  return(permutated.cut)
}

