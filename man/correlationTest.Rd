\name{correlationTest}
\alias{correlationTest}
\alias{correlationTest-methods}
\alias{correlationTest,data.frame-method}
\alias{correlationTest,matrix-method}
\alias{correlationTest,numeric-method}
\alias{correlationTest,vector-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
correlation Test
}
\description{
Check if a correlation exits between two descriptors
}
\usage{
correlationTest(dataSet, yName, cutlineCorr, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Input data set (should be numeric values)
}
  \item{yName}{
Name of the column which should be predicted
}
  \item{cutlineCorr}{
The column, from a correlation pair (value more than the cutline), which has less correlation agains the y column is removed
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
\code{\link{car}}, 
\code{\link{incReg}}
}

\keyword{ manip }

