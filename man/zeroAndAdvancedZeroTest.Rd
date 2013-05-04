\name{zeroAndAdvancedZeroTest}
\alias{zeroAndAdvancedZeroTest}
\alias{zeroAndAdvancedZeroTest-methods}
\alias{zeroAndAdvancedZeroTest,data.frame-method}
\alias{zeroAndAdvancedZeroTest,matrix-method}
\alias{zeroAndAdvancedZeroTest,numeric-method}
\alias{zeroAndAdvancedZeroTest,vector-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
zero and advanced zero test
}
\description{
Check if the single columns have enough variations
}
\usage{
zeroAndAdvancedZeroTest(dataSet, yName, cutlineZeroTest, cutlineAdvZeroTest, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
input data set
}
  \item{yName}{
Name of the column which should be predicted
}
  \item{cutlineZeroTest}{
not more than the given cutline should be zero
}
  \item{cutlineAdvZeroTest}{
there have to be no descriptor which have more than the given cutline on the same cases
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

}

\keyword{ manip }