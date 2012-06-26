\name{stability.path}
\alias{stability.path}
\title{
Stability path for glmnet models
}
\description{
The function calculates the stability path for glmnet models
}
\usage{
stability.path(y,x,size=0.632,steps=100,weakness=1,mc.cores=getOption("mc.cores", 2L),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response variable.
}
  \item{x}{
matrix of predictor variables.
}
  \item{size}{
proportion of samples drawn in every subsample.
}
  \item{steps}{
number of subsamples.
}
  \item{weakness}{
weakness parameter for the randomised lasso.
}
  \item{mc.cores}{
number of cores used for the parallelization
}
  \item{...}{
further arguments that are passed to glmnet. 
}
}
\description{
performs stability selection...
}
\value{
returns a list...
  \item{}{}
  \item{}{}
  \item{}{}
}
\author{
Martin Sill \
\email{m.sill@dkfz.de}
}
\references{
Meinshausen and B\"uhlmann, Stability Selection 
}
\seealso{ \code{\link{glmnet}}
}

\examples{
}

\keyword{stability selection}
