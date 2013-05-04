\name{incrementalMultivariateRegression}
\alias{incrementalMultivariateRegression}
\alias{incrementalMultivariateRegression-methods}
\alias{incrementalMultivariateRegression,data.frame-method}
\alias{incrementalMultivariateRegression,matrix-method}
\alias{incrementalMultivariateRegression,numeric-method}
\alias{incrementalMultivariateRegression,vector-method}

\title{
incremental Multivariate Regression
}
\description{
processed the incremental Multivariate Regression
}
\usage{
incrementalMultivariateRegression(dataSet, yName, alphaValue, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
input data set
}
  \item{yName}{
Name of the column which should be predicted
}
  \item{alphaValue}{
significance level
}
  \item{\dots}{
not used
}
}

\author{
Alexander Entzian \email{a.entzian@alexander-entzian.de}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{incReg}}
\code{\link{car}}
\code{\link{stats}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }
\keyword{ multivariate }
