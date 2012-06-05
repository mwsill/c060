library(glmnet)
library(parallel)

stability.path <- function(y,x,mc.cores=getOption("mc.cores", 2L),size=0.632,steps=100,weakness=1,...){
	fit <- glmnet(x,y,...)
	p <- ncol(x)
	#draw subsets
  	subsets <- sapply(1:steps,function(v){sample(1:nrow(x),nrow(x)*size)})
	res <- mclapply(1:steps,mc.cores=mc.cores,glmnet.subset,subsets,x,y
                  ,lambda=fit$lambda,weakness,p,...)
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
plot.stabpath <- function(stabpath,fwer,pi_thr=0.8,...){
	sel <- stability.selection(stabpath,fwer,pi_thr)
	p <- dim(stabpath$fit$beta)[1]
	#stability path
	cols <- rep("black",p)
	cols[sel$stable] <- "red"
	par(mfrow=c(2,1))
	matplot(t(stabpath$fit$beta),xaxt="n",type="l",col=cols,lty=1,ylab=expression(paste(beta[i]))
        ,xlab=expression(paste(lambda)),main="Penalization Path",cex.lab=1,cex.axis=1,...)
	matplot(t(stabpath$stabpath),xaxt="n",type="l",col=cols,lty=1,ylab=expression(paste(hat(Pi)))
        ,xlab=expression(paste(lambda)),main="Stability Path",ylim=c(0,1),cex.lab=1,cex.axis=1,...)
	abline(h=pi_thr,col="darkred",lwd=1,lty=1)
	abline(v=sel$lpos,col="darkred",lwd=1,lty=1)
	#text(x=20,y=0.9,paste(expression(paste(lambda)),"=",paste(round(sel[[2]],digits=3)),sep=""),cex=0.75)
	return(sel)
}

#gaussian
x=matrix(rnorm(100*10000),100,10000)
y=rnorm(100)
#1 core
system.time(res <- stability.path(y,x,weakness=1,mc.cores=1))
#2 cores
system.time(res <- stability.path(y,x,weakness=1,mc.cores=2))
# additional arguments for glmnet nlambda alpha 
system.time(res <- stability.path(y,x,weakness=1,mc.cores=2,nlambda=50,alpha=0.5))
stability.selection(res,fwer=0.2)
plot.stabpath(res,fwer=0.2)

#binomial
y=sample(1:2,100,replace=TRUE)
res <- stability.path(y,x,weakness=1,mc.cores=2,family="binomial")
stability.selection(res,fwer=0.2)
plot.stabpath(res,fwer=0.2)
    
#multinomial
#y=sample(1:4,100,replace=TRUE)
#res <- stability.path(y,x,weakness=0.2,family="multinomial")
#stability.selection(res,fwer=0.2)
#plot.stabpath(res,fwer=0.2)
    
#poisson
N=100; p=10000
nzc=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
f = x[,seq(nzc)]%*%beta
mu=exp(f)
y=rpois(N,mu)
res <- stability.path(y,x,weakness=1,mc.cores=2
                      ,family="poisson")
stability.selection(res,fwer=0.2)
plot.stabpath(res,fwer=0.2)

#Cox
library(survival)
set.seed(10101)
N=100;p=10000
nzc=p/3
x=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
fx=x[,seq(nzc)]%*%beta/3
hx=exp(fx)
ty=rexp(N,hx)
tcens=rbinom(n=N,prob=.3,size=1)# censoring indicator
y=cbind(time=ty,status=1-tcens) # y=Surv(ty,1-tcens) with library(survival)

res <- stability.path(y,x,weakness=1,mc.cores=2,family="cox")
stability.selection(res,fwer=0.2)
plot.stabpath(res,fwer=0.2)










