\name{multicollinearityTest}
\alias{multicollinearityTest}
\alias{multicollinearityTest-methods}
\alias{multicollinearityTest,data.frame-method}
\alias{multicollinearityTest,matrix-method}
\alias{multicollinearityTest,numeric-method}
\alias{multicollinearityTest,vector-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
multicollinearity Test
}
\description{
Check if a multicollinearity exits between more than two descriptors. 
If more than two descriptors have multicol. then the descriptor with the lowest correlation to "y" removed. 
After that multicollinearity is recalculated.
}
\usage{
multicollinearityTest(dataSet, yName, cutlineMulticoll, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
input data set
}
  \item{yName}{
Name of the column which should be predicted
}
  \item{cutlineMulticoll}{
if a descriptor have a higher value than this cutline, multicollinearity is there
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

