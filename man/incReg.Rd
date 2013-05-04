\name{incReg}
\alias{incReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
incReg constructor
}
\description{
Constructor for users. Creates the incReg-object and start the process. 
}
\usage{
incReg(yName, dataSet = data.frame(), csvFile = "None", doAnalysis = TRUE, 
	   shouldScaling = FALSE, shouldCentering = FALSE, missingValueTyp = "rm", 
	   missingValueIsCol = TRUE, cutlineZeroTest = 0.8, cutlineAdvZeroTest = 0.25, 
	   cutlineCorr = 0.49, cutlineMulticoll = 10, alphaValue = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yName}{
column name of the column which is to be predicted 
}
  \item{dataSet}{
input data set, alternative see \code{csvFile}
}
  \item{csvFile}{
input csv File, alternative see \code{dataSet}
}
  \item{doAnalysis}{
should the process be performed directly (TRUE or FALSE) 
}
  \item{shouldScaling}{
should the data set be scaled (TRUE or FALSE)
}
  \item{shouldCentering}{
should the data set be centered (TRUE or FALSE)
}
  \item{missingValueTyp}{
what should the method do with missing values: remove ("rm") or calculate the average ("avg")
}
  \item{missingValueIsCol}{
if \code{missingValueTyp} is "rm" then this value is important: TRUE -> column will be removed, FALSE -> the row
}
  \item{cutlineZeroTest}{
not more than the given cutline should be zero else the descriptor have to be removed
}
  \item{cutlineAdvZeroTest}{
there have to be no descriptor which have more than the given cutline on the same cases
}
  \item{cutlineCorr}{
the column, from a correlation pair (value more than the cutline), which has less correlation against the y column is removed
}
  \item{cutlineMulticoll}{
if a descriptor have a higher value than this cutline, multicollinearity is there
}
  \item{alphaValue}{
significance level
}
}

\author{
Alexander Entzian \email{a.entzian@alexander-entzian.de}
}

\seealso{
\code{\link{incReg-package}}
\code{\link{IncReg-class}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

	data(exampleIncRegData)

	result <- incReg(yName = "exp.pkb", dataSet = exampleIncRegData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip}

