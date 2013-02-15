\name{plot.peperr.curves}
\alias{plot.peperr.curves}
\title{Plot method for prediction error curves of a peperr object}
\description{Plots individual and aggregated prediction error estimates based on bootstrap samples.}
\usage{
plot.peperr.curves(x, at.risk=TRUE, allErrors=FALSE,...)
}
\arguments{
\item{x}{\code{peperr} object.}
\item{at.risk}{number at risk to be display. default is TRUE.}
\item{allErrors}{Display .632, no information and average out-of-bag error in addition. default is FALSE.}
\item{\dots}{additional arguments, not used.}

}
\details{
This function is literally taken from \code{plot.peperr} in the \code{peperr} package.
The display of prediction error curves is adapted to allow for numbers at risk.
}
\examples{
\dontrun{

# example from glmnet package
set.seed(10101)
library(glmnet)
library(survival)
library(peperr)

N=1000;p=30
nzc=p/3
x=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
fx=x[,seq(nzc)]%*%beta/3
hx=exp(fx)
ty=rexp(N,hx)
tcens=rbinom(n=N,prob=.3,size=1)# censoring indicator
y=Surv(ty,1-tcens)

peperr.object <- peperr(response=y, x=x, 
                        fit.fun=fit.glmnet, args.fit=list(family="cox"), 
                        complexity=complexity.glmnet,  
                        args.complexity=list(family="cox",nfolds=10),
                        indices=resample.indices(n=N, method="sub632", sample.n=10))

plot.peperr.curves(peperr.object)
# peperr
plot(peperr.object)
}
}
\keyword{models} \keyword{regression} \keyword{survival}

\author{
Thomas Hielscher 
\email{t.hielscher@dkfz.de}
}

\seealso{\code{\link[peperr]{peperr}}}
