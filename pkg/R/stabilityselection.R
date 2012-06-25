#Change log:
# 13-06-2012 (MZ):
# - changes to plot.stabpath() in a fashion similar to glmnet:::plotCoef; allow users to choose colours

require(glmnet)
require(parallel)

stability.path <- function(y,x,mc.cores=getOption("mc.cores", 2L),size=0.632,steps=100,weakness=1,...){
	fit <- glmnet(x,y,...)
	p <- ncol(x)
	#draw subsets
  	subsets <- sapply(1:steps,function(v){sample(1:nrow(x),nrow(x)*size)})
  
	# parallel computing depending on OS
	# UNIX/Mac
	if (.Platform$OS.type!="windows") {
	  res <- mclapply(1:steps,mc.cores=mc.cores,glmnet.subset,subsets,x,y,lambda=fit$lambda,weakness,p,...)
	} else {
	  # Windows  
	  cl  <- makePSOCKcluster(mc.cores)
	  clusterExport(cl,c("glmnet","drop0"))
	  res <- parLapply(cl, 1:steps,glmnet.subset,subsets,x,y,lambda=fit$lambda,weakness,p,...)
	  stopCluster(cl)
	}
  
  #merging
	stabpath <- res[[1]]
	qmat <- matrix(ncol=ncol(res[[1]]),nrow=steps)
	qmat[1,] <- colSums(res[[1]])
	for(i in 2:length(res)){
  		qmat[i,] <- colSums(res[[i]])
		stabpath <- stabpath + res[[i]]
	}
	stabpath <- stabpath/length(res)
	qs <- colMeans(qmat)
	out <- list(fit=fit,stabpath=stabpath,qs=qs)	
	class(out) <- "stabpath" 
	return(out)
}

#internal function used by lapply 
glmnet.subset <- function(index,subsets,x,y,lambda,weakness,p,...){
  if(length(dim(y))==2|class(y)=="Surv"){
    glmnet(x[subsets[,index],],y[subsets[,index],],lambda=lambda
           ,penalty.factor= 1/runif(p,weakness,1),...)$beta!=0
  }else
    glmnet(x[subsets[,index],],y[subsets[,index]],lambda=lambda
           ,penalty.factor= 1/runif(p,weakness,1),...)$beta!=0
}

#performs error control and returns estimated set of stable variables and corresponding lambda
stability.selection <- function(stabpath,fwer,pi_thr=0.8){
	p <- dim(stabpath$fit$beta)[1]
	qv <- sqrt(pi_thr*fwer*p)
	lpos <- which(stabpath$qs>qv)[1]
	stable <- which(stabpath$stabpath[,lpos]>=pi_thr)
	out <- list(stable=stable,lambda=stabpath$fit$lambda[lpos],lpos=lpos)
	return(out)
}

#plot penalization and stability path 
plot.stabpath <- function(stabpath, fwer, pi_thr=0.8, xvar=c("lambda", "norm", "dev"), col.all="black", col.sel="red",...){
  sel <- stability.selection(stabpath,fwer,pi_thr)
  
  beta = as.matrix(stabpath$fit$beta)
  p <- dim(beta)[1]
  which = nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = apply(abs(beta), 2, sum)
    iname = "L1 Norm"
  }, lambda = {
    index = log(stabpath$fit$lambda)
    iname = expression(paste("log ",lambda))
  }, dev = {
    index = stabpath$fit$dev
    iname = "Fraction Deviance Explained"
  })
  
  #stability path
  cols <- rep(col.all,p)
  cols[sel$stable] <- col.sel
  lwds <- rep(1,p)
  lwds[sel$stable] <- 2
  par(mfrow=c(2,1))
  matplot(y=t(beta), x=index
          ,type="l",col=cols,lwd=lwds,lty=1,ylab=expression(paste(beta[i]))
          ,xlab=iname,main="Penalization Path",cex.lab=1,cex.axis=1,...)
  matplot(y=as.matrix(t(stabpath$stabpath)), x=index
          ,type="l",col=cols,lwd=lwds,lty=1,ylab=expression(paste(hat(Pi)))
          ,xlab=iname,main="Stability Path",ylim=c(0,1),cex.lab=1,cex.axis=1,...)
  abline(h=pi_thr,col="darkred",lwd=1,lty=1)
  abline(v=index[sel$lpos],col="darkred",lwd=1,lty=1)
  #text(x=20,y=0.9,paste(expression(paste(lambda)),"=",paste(round(sel[[2]],digits=3)),sep=""),cex=0.75)
  
  return(sel)
}

