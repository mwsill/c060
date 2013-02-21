\name{summary.int.search}
\alias{summary.int.search}
\title{
Summary method for interval search models
}
\description{
Produces a summary of a fitted interval search model
}
\usage{
summary.int.search(fit,digits = max(3, getOption("digits") - 3), verbose=TRUE,first.n=5)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
 an object as returned by the function \code{EPSGO}.
  }
  \item{digits}{
 digits after the comma
  }
   \item{verbose}{
 default set to TRUE.
  }
   \item{first.n}{
 show first.n entries , default 5.
  }
}
\value{
A List of following elements 
  \item{info}{
  a data frame of four objects for optimal models\cr
    \describe{  
      \item{alpha}{
      a vector of alphas
      }
      \item{lambda}{
      a vector of  penalization parameter lambda
      }
      \item{deviances}{
      a vector of deviances
      }
      \item{n.features}{
      a vector of number of features selected in each optimal model 
      }
    }
  }  
\item{opt.alpha}{
  an optimal value for alpha
}
\item{opt.lambda}{
  an optimal value for lambda
}
\item{opt.error}{
  an optimal value for alpha
}
\item{opt.alpha}{
  an optimal value for error, hier minimal diviance
}
\item{opt.models}{
  a list of optimal models with the same optimal error 
}  
}
\author{
Natalia Becker  \
\email{natalia.becker@dkfz.de}
}

\seealso{ \code{\link{EPSGO}}
}
\examples{
\dontrun{
#gaussian
#set.seed(1234)
#x=matrix(rnorm(100*1000,0,1),100,1000)
#y <- x[1:100,1:1000] \%*\% matrix(c(rep(2,5),rep(-2,5),rep(.1,990)))
#y<-ifelse(y>0, "1","-1")

#bounds <- t(data.frame(alpha=c(0, 1)))
#colnames(bounds) <- c("lower", "upper") 

#Q.func<- "tune.glmnet.interval"
#fit <- EPSGO(Q.func,bounds,parms.coding="none",x=x,y=y)
#summary.int.search(fit)
}}
\keyword{summary}
