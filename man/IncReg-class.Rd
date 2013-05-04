\name{IncReg-class}
\Rdversion{1.1}
\docType{class}
\alias{IncReg-class}
\alias{[,IncReg-method}
\alias{[<-,IncReg-method}
\alias{getAlphaValue,IncReg-method}
\alias{getCollResult,IncReg-method}
\alias{getCsvFile,IncReg-method}
\alias{getCutlineAdvZeroTest,IncReg-method}
\alias{getCutlineCorr,IncReg-method}
\alias{getCutlineMulticoll,IncReg-method}
\alias{getCutlineZeroTest,IncReg-method}
\alias{getDataSet,IncReg-method}
\alias{getFinalModel,IncReg-method}
\alias{getMissingValueIsCol,IncReg-method}
\alias{getMissingValueTyp,IncReg-method}
\alias{getModelNumberX,IncReg-method}
\alias{getModelResult,IncReg-method}
\alias{getMultiCollResult,IncReg-method}
\alias{getPreModelResults,IncReg-method}
\alias{getShouldCentering,IncReg-method}
\alias{getShouldScaling,IncReg-method}
\alias{getSquareBracket,IncReg-method}
\alias{getYname,IncReg-method}
\alias{print,IncReg-method}
\alias{processing,IncReg-method}
\alias{setAlphaValue<-,IncReg-method}
\alias{setCollResult<-,IncReg-method}
\alias{setCsvFile<-,IncReg-method}
\alias{setCutlineAdvZeroTest<-,IncReg-method}
\alias{setCutlineCorr<-,IncReg-method}
\alias{setCutlineMulticoll<-,IncReg-method}
\alias{setCutlineZeroTest<-,IncReg-method}
\alias{setDataSet<-,IncReg-method}
\alias{setMissingValueIsCol<-,IncReg-method}
\alias{setMissingValueTyp<-,IncReg-method}
\alias{setModelResult<-,IncReg-method}
\alias{setMultiCollResult<-,IncReg-method}
\alias{setShouldCentering<-,IncReg-method}
\alias{setShouldScaling<-,IncReg-method}
\alias{setSquareBracket,IncReg-method}
\alias{setYname<-,IncReg-method}
\alias{show,IncReg-method}
\alias{summary,IncReg-method}

\title{Class \code{"IncReg"}}
\description{
The incReg class calculate the incremental multivariate Regression based on a given dataset or a .csv file, based on the multiple linear regression
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IncReg", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.dataSet}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{.csvFile}:}{Object of class \code{"character"} ~~ }
    \item{\code{.yName}:}{Object of class \code{"character"} ~~ }
    \item{\code{.shouldScaling}:}{Object of class \code{"logical"} ~~ }
    \item{\code{.shouldCentering}:}{Object of class \code{"logical"} ~~ }
    \item{\code{.missingValueTyp}:}{Object of class \code{"character"} ~~ }
    \item{\code{.missingValueIsCol}:}{Object of class \code{"logical"} ~~ }
    \item{\code{.cutlineZeroTest}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.cutlineAdvZeroTest}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.cutlineCorr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.cutlineMulticoll}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.alphaValue}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.multiCollResult}:}{Object of class \code{"list"} ~~ }
    \item{\code{.collResult}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{.modelResult}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{IncRegOption}"}, directly.
Class \code{"\linkS4class{IncRegFinalmodelResult}"}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "IncReg")}: ... }
    \item{[<-}{\code{signature(x = "IncReg")}: ... }
    \item{getAlphaValue}{\code{signature(object = "IncReg")}: ... }
    \item{getCollResult}{\code{signature(object = "IncReg")}: ... }
    \item{getCsvFile}{\code{signature(object = "IncReg")}: ... }
    \item{getCutlineAdvZeroTest}{\code{signature(object = "IncReg")}: ... }
    \item{getCutlineCorr}{\code{signature(object = "IncReg")}: ... }
    \item{getCutlineMulticoll}{\code{signature(object = "IncReg")}: ... }
    \item{getCutlineZeroTest}{\code{signature(object = "IncReg")}: ... }
    \item{getDataSet}{\code{signature(object = "IncReg")}: ... }
    \item{getFinalModel}{\code{signature(object = "IncReg")}: ... }
    \item{getMissingValueIsCol}{\code{signature(object = "IncReg")}: ... }
    \item{getMissingValueTyp}{\code{signature(object = "IncReg")}: ... }
    \item{getModelNumberX}{\code{signature(object = "IncReg")}: ... }
    \item{getModelResult}{\code{signature(object = "IncReg")}: ... }
    \item{getMultiCollResult}{\code{signature(object = "IncReg")}: ... }
    \item{getPreModelResults}{\code{signature(object = "IncReg")}: ... }
    \item{getShouldCentering}{\code{signature(object = "IncReg")}: ... }
    \item{getShouldScaling}{\code{signature(object = "IncReg")}: ... }
    \item{getSquareBracket}{\code{signature(x = "IncReg")}: ... }
    \item{getYname}{\code{signature(object = "IncReg")}: ... }
    \item{print}{\code{signature(x = "IncReg")}: ... }
    \item{processing}{\code{signature(IncReg = "IncReg")}: ... }
    \item{setAlphaValue<-}{\code{signature(object = "IncReg")}: ... }
    \item{setCollResult<-}{\code{signature(object = "IncReg")}: ... }
    \item{setCsvFile<-}{\code{signature(object = "IncReg")}: ... }
    \item{setCutlineAdvZeroTest<-}{\code{signature(object = "IncReg")}: ... }
    \item{setCutlineCorr<-}{\code{signature(object = "IncReg")}: ... }
    \item{setCutlineMulticoll<-}{\code{signature(object = "IncReg")}: ... }
    \item{setCutlineZeroTest<-}{\code{signature(object = "IncReg")}: ... }
    \item{setDataSet<-}{\code{signature(object = "IncReg")}: ... }
    \item{setMissingValueIsCol<-}{\code{signature(object = "IncReg")}: ... }
    \item{setMissingValueTyp<-}{\code{signature(object = "IncReg")}: ... }
    \item{setModelResult<-}{\code{signature(object = "IncReg")}: ... }
    \item{setMultiCollResult<-}{\code{signature(object = "IncReg")}: ... }
    \item{setShouldCentering<-}{\code{signature(object = "IncReg")}: ... }
    \item{setShouldScaling<-}{\code{signature(object = "IncReg")}: ... }
    \item{setSquareBracket}{\code{signature(x = "IncReg")}: ... }
    \item{setYname<-}{\code{signature(object = "IncReg")}: ... }
    \item{show}{\code{signature(object = "IncReg")}: ... }
    \item{summary}{\code{signature(object = "IncReg")}: ... }
	 }
}

\author{
Alexander Entzian \email{a.entzian@alexander-entzian.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{incReg}}
\code{\link{incReg-package}}
}
\examples{
showClass("IncReg")
}
\keyword{classes}
