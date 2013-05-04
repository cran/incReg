\name{IncRegOption-class}
\Rdversion{1.1}
\docType{class}
\alias{IncRegOption-class}

\title{Class \code{"IncRegOption"}}
\description{
save the options of the incReg class
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IncRegOption", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.shouldScaling}:}{Object of class \code{"logical"} ~~ }
    \item{\code{.shouldCentering}:}{Object of class \code{"logical"} ~~ }
    \item{\code{.missingValueTyp}:}{Object of class \code{"character"} ~~ }
    \item{\code{.missingValueIsCol}:}{Object of class \code{"logical"} ~~ }
    \item{\code{.cutlineZeroTest}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.cutlineAdvZeroTest}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.cutlineCorr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.cutlineMulticoll}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.alphaValue}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "IncRegOption" in the signature.
}

\author{
Alexander Entzian \email{a.entzian@alexander-entzian.de}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{incReg}}
}

\examples{
showClass("IncRegOption")
}
\keyword{classes}
