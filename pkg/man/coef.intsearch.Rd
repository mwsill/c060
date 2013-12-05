\name{coef.opt.intsearch}
\alias{coef.opt.intsearch}

\title{
Get coefficients for a model 
}
\description{
Get coefficients for a model after applying interval search for tuning parameters 
}
\usage{
coef.opt.intsearch(object,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{ an object as returned by the function \code{summary.intsearch}.}
\item{\dots}{additional argument(s)}
}
\value{
 named vector of non-zero coeficients for the optimal lambda}
\author{
Natalia Becker  \
\email{natalia.becker@dkfz.de}
}

\seealso{ \code{\link{EPSGO}}, \code{\link{summary.intsearch}}
}

\keyword{system}
