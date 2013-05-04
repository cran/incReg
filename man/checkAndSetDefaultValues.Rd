\name{checkAndSetDefaultValues}
\alias{checkAndSetDefaultValues}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
check and set default values
}
\description{
checked some parameters if they are in their specified limits. Otherwise, the default value is used.
}
\usage{
checkAndSetDefaultValues(listOfData)
}

\arguments{
  \item{listOfData}{
		should be a list of the parameter or a object of the class \code{call}
	}

}
\details{
The values to be checked are currently: csvFile, missingValueTyp, cutlineZeroTest, cutlineAdvZeroTest , cutlineCorr, cutlineMulticoll, alphaValue
}
\value{
the returned value is a list of the adapted values
}

\author{
Alexander Entzian \email{a.entzian@alexander-entzian.de}
}

\seealso{
	\code{\link{incReg}}
}

\keyword{ manip }
