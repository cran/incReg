\name{deletRowsAndColumsAndOrCalculatMean}
\alias{deletRowsAndColumsAndOrCalculatMean}
\alias{deletRowsAndColumsAndOrCalculatMean-methods}
\alias{deletRowsAndColumsAndOrCalculatMean,data.frame-method}
\alias{deletRowsAndColumsAndOrCalculatMean,matrix-method}
\alias{deletRowsAndColumsAndOrCalculatMean,numeric-method}
\alias{deletRowsAndColumsAndOrCalculatMean,vector-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
delet rows and colums and/or calculat mean
}
\description{
Delete all rows and cols with NA or calculate the mean for each cell (special handling for the y column)
}
\usage{
deletRowsAndColumsAndOrCalculatMean(dataSet, yName, 
			missingValueTyp, missingValueIsCol, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
input data set
}
  \item{yName}{
Name of the column which should be predicted
}
  \item{missingValueTyp}{
remove ("rm") or average ("avg")
}
  \item{missingValueIsCol}{
if missingValueTyp is "rm" then this value is important: TRUE -> column will be removed, FALSE -> the row
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

