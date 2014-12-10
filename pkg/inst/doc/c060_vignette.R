### R code from vignette source 'c060_vignette.Rnw'

###################################################
### code chunk number 1: setup
###################################################
rm(list=ls())
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
#setwd("./vignettes")
# load libraries
library(cacheSweave)
library(Biobase)
library(genefilter)
library(GEOquery)
library(limma)
library(xtable)
library(survival)

library(glmnet)
require(penalizedSVM)
library(parallel) #for stabilityselection.R
library(peperr) #for peperr_glmnet.R
#library(tgp) #for EPSGO.R
#library(mlegp) #for EPSGO.R
library(pamr) #for tune.glmnet.interval.R
library(lattice) # for plotting 

library(c060)

dir.create("figures",showWarnings = FALSE)

# set cache folder; needs different driver, to run as:
# Sweave("c060_vignette.Rnw",driver=cacheSweaveDriver)
setCacheDir("cache")


###################################################
### code chunk number 2: GEOdataGSE12417
###################################################
if (file.exists(".Rdata/expressionSet.Rdata")) {
    load(".Rdata/expressionSet.Rdata") 
} else {
    library(GEOquery)
    
    dir.create(".GEOdata",showWarnings = FALSE)
    dir.create(".Rdata",showWarnings = FALSE)
    
    # download processed data files from url or if available get data from locally saved copies
    geo_exprs_data  <- getGEO("GSE12417", destdir=".GEOdata",AnnotGPL=T)[[1]]
    annotation(geo_exprs_data) <- "hgu133plus2"
    
    #clinical characteristics including survival data are in "characteristics_ch1" in the phenoData 
    clin.df <- as.data.frame(strsplit2(as.character(pData(geo_exprs_data)$characteristics_ch1), split=";"))
    clin.df[,1] <- strsplit2(clin.df[,1]," ")[,7]
    clin.df[,2] <- as.numeric(sub("=", "", strsplit2(clin.df[,2]," ")[,3]))
    clin.df[,3] <- as.numeric(strsplit2(clin.df[,3]," ")[,4])
    clin.df[,4] <- as.numeric(strsplit2(clin.df[,4]," ")[,4])
    colnames(clin.df) <- c("FAB", "age", "os", "os_status")
    rownames(clin.df) <- rownames(pData(geo_exprs_data))
    pData(geo_exprs_data) <- clin.df

    save(geo_exprs_data, file=".Rdata/expressionSet.Rdata")  
    unlink(".GEOdata", recursive=TRUE) #delete GEOdata
} 


###################################################
### code chunk number 3: preprocessData
###################################################
# unspecific filtering: select most varying probe sets
nfeatures <- 10000
rowvars    <- rowVars(exprs(geo_exprs_data))
topprobes   <- which(rowvars>=sort(rowvars,decreasing=T)[nfeatures])
eset       <- geo_exprs_data[topprobes,]     


###################################################
### code chunk number 4: glmnet0
###################################################
fit <- glmnet(y=Surv(pData(eset)$os, pData(eset)$os_status), 
              x=t(exprs(eset)), family="cox")


###################################################
### code chunk number 5: glmnet
###################################################
set.seed(1234)
cvres <- cv.glmnet(y=Surv(pData(eset)$os, pData(eset)$os_status), 
                   x=t(exprs(eset)), family="cox", nfolds=10)
res <- cvres$glmnet.fit
plot(cvres)


###################################################
### code chunk number 6: glmnet0
###################################################
cof <- coef(res, s=cvres$lambda.min)
names.cof <- rownames(cof)
cofn <- cof[which(cof!=0)]
names(cofn) <- names.cof[which(cof!=0)]


###################################################
### code chunk number 7: glmnet2
###################################################
print(cofn, digits=3)
cofn.lasso <- cofn
#cofn.ix <- predict(res, type="nonzero", s=cvres$lambda.min)


###################################################
### code chunk number 8: glmnet3 (eval = FALSE)
###################################################
## cof <- coef(res, s=cvres$lambda.min)
## Plot.coef.glmnet(cvfit=cvres, betas=rownames(cof)[which(cof!=0)])


###################################################
### code chunk number 9: glmnet4
###################################################
Plot.coef.glmnet(cvfit=cvres, betas=names(cofn))


###################################################
### code chunk number 10: peperrMandatoryParallel
###################################################
obj  <- peperr(response=Surv(pData(eset)$os, pData(eset)$os_status),
               x=data.frame(eset$age,t(exprs(eset))),
               fit.fun=fit.glmnet, args.fit=list(standardize=FALSE, family="cox",
               penalty.factor=rep(0:1, times=c(1,dim(eset)[1]))),
               complexity=complexity.glmnet, 
               args.complexity=list(standardize=FALSE, nfolds=10, family="cox", 
               penalty.factor=rep(0:1, times=c(1,dim(eset)[1]))),
               RNG="fixed", seed=0815, cpus=3, parallel=TRUE, clustertype="SOCK", 
               load.list=list(packages=c("c060")), 
               indices=resample.indices(n=dim(eset)[2],
               sample.n=1000, method="sub632"))


###################################################
### code chunk number 11: peperrPlot1 (eval = FALSE)
###################################################
## Plot.peperr.curves(obj, at.risk=TRUE, allErrors=FALSE, bootRuns=FALSE,
##                    bootQuants=TRUE, bootQuants.level=0.95, leg.cex=0.7)


###################################################
### code chunk number 12: peperrPlot2
###################################################
Plot.peperr.curves(obj, at.risk=TRUE, allErrors=FALSE, bootRuns=FALSE,
                   bootQuants=TRUE, bootQuants.level=0.95, leg.cex=0.7)


###################################################
### code chunk number 13: stabilitySelection1
###################################################
y <- cbind(time=pData(eset)$os, status=pData(eset)$os_status)
set.seed(1234)
spath <- stabpath(y=y, x=t(exprs(eset)), mc.cores=2,
                        family="cox", weakness=.8)


###################################################
### code chunk number 14: stabilitySelection2
###################################################
stabsel(spath, error=1, type="pfer", pi_thr=0.6)$stable


###################################################
### code chunk number 15: stabilitySelection3
###################################################
plot(spath, error=1, type="pfer", pi_thr=0.6, xvar="lambda", col.all="gray")


###################################################
### code chunk number 16: interval_search_setup
###################################################
  x <- t(exprs(eset))
  y <- cbind(time=pData(eset)$os,status=pData(eset)$os_status)

  bounds <- t(data.frame(alpha=c(0, 1)))
  colnames(bounds)<-c("lower","upper")
  
  nfolds <- 10
  set.seed(1234)
  foldid <- balancedFolds(class.column.factor=y[,2], cross.outer=nfolds)


###################################################
### code chunk number 17: interval_search_cox_show (eval = FALSE)
###################################################
##   fit <- epsgo(Q.func="tune.glmnet.interval", 
##              bounds=bounds, 
##              parms.coding="none", 
##              seed = 1234, 
##              fminlower = -100,
##              x = x, y = y, family = "cox", 
##              foldid = foldid,
##              type.min = "lambda.1se",
##              type.measure = "deviance")


###################################################
### code chunk number 18: interval_search_cox
###################################################
if (file.exists(".Rdata/fit_interval_search.RData")) {
  load(".Rdata/fit_interval_search.RData") 
} else {
  print("start interval search")
  fit <- epsgo(Q.func="tune.glmnet.interval",  
             bounds=bounds, 
             parms.coding="none", 
             show = "none", 
             N = NULL,   
             seed = 1234, 
             fminlower = -100,
             # Q.func arguments
             x = x, y = y, family = "cox", 
             foldid = foldid,
             type.min = "lambda.1se",
             type.measure = "deviance",
             verbose = TRUE)
  #names(fit)
  save(fit,  file=".Rdata/fit_interval_search.RData") 
}  


###################################################
### code chunk number 19: interval_search_cox_extract_show
###################################################
sumint <- summary(fit, verbose=TRUE)


###################################################
### code chunk number 20: interval_search_cox_fit_opt
###################################################
#select the model with optimal parameters from the object fit
opt.model <- sumint$opt.models[[1]]
cofn <- coef(opt.model$fit)


###################################################
### code chunk number 21: interval_search_cox_check_overlap_stability
###################################################
names(cofn.lasso) %in% rownames(cofn)
'206932_at' %in% rownames(cofn)


###################################################
### code chunk number 22: interval_search_cox_out_plot1
###################################################
plot(sumint)


###################################################
### code chunk number 23: interval_search_cox_out_plot2
###################################################
plot(sumint,type="points") 


