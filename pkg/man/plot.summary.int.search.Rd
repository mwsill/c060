\name{plot.summary.int.search}
\alias{plot.summary.int.search}
\alias{plot.points.int.search}
\title{
Plot Summary object for interval search models
}
\description{
Produces a plot for summary object of a fitted interval search model.
Plot 'visited' points against iteration steps. start.N points are initial points selected before interval search starts.
}
\usage{
plot.summary.int.search(summary.int)
plot.points.int.search(summary.int, start.N=21)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
 an object as returned by the function \code{summary.int.search}.
  }
  \item{srat.N}{
 number of initial points, default 21.
  }
}

\author{
Natalia Becker  \
\email{natalia.becker@dkfz.de}
}

\seealso{ \code{\link{EPSGO}}, \code{\link{summary.int.search}}
}

\keyword{plot}
