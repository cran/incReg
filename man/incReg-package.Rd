\name{incReg-package}
\alias{incReg-package}
\docType{package}
\title{
incremental multivariate Regression
}
\description{
The incReg package calculate the incremental multivariate Regression based on a given dataset or a .csv file, based on the multiple linear regression
}
\details{
\tabular{ll}{
Package: \tab incReg\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-05\cr
License: \tab LGPL (>= 3)\cr
Depends: \tab R (>= 2.14.0), car, methods\cr
}

}
\author{
Alexander Entzian \email{a.entzian@alexander-entzian.de}
}

\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ multivariate }

\seealso{
\code{\link{incReg}}
\code{\link{IncReg-class}}
}
\examples{

	data(exampleIncRegData)
	
	result <- incReg(yName = "exp.pkb", dataSet = exampleIncRegData)
	summary(result)
	finalModel <- getFinalModel(result)
	
}
