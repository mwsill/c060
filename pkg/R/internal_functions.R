
###########################################################################################################
.get.cofn<-function(model){
  # get coef for a model from fit object after running interval search 
  f1 <- model$cvreg
  res<-f1$glmnet.fit
  cof <- coef(res, s=model$lambda)
  names.cof <- rownames(cof)
  cofn <- cof[which(cof!=0)]
  names(cofn) <- names.cof[which(cof!=0)] 
  bet <- res$beta[match(names(cofn), rownames(res$beta)),]
  
  return(cofn)
}  

