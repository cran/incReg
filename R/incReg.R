
###############################################################################
##
## Author: Alexander Entzian
## email: a.entzian@alexander-entzian.de
##
###############################################################################
#
#
# The incReg package is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
#
# The incReg package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You can find a copy of the GNU Lesser General Public License at <http://www.r-project.org/Licenses/> or <http://www.gnu.org/licenses/>.
#
#################################################
#
# Incremental Multivariate Regression
# the package calculate the incremental multivariate Regression based on a given dataset or a .csv file
#
# Test Methods:
# * Check NA-Values: delete all rows and cols with NA or calculate the mean for each cell (special handling for the y column)
# * Zero- and advanced Zero Test: check if the single columns have enough variations (Zero-Test: not more than the given cutline should be zero, ...)
# * Multicollinearity Test: check if a multicollinearity exits between more than two descriptors
# * correlation Test: check if a correlation exits between two descriptors
#
# NoTest Methods:
# * Scaling and centering of the dataset
# * Incremental Multivariate Regression
# * get the model
# * summary
#
# Package Use
# * result <- incReg(yName = "exp.pkb", csvFile = "incReg.csv")
# * summary(result)
# * finalModel <- getFinalModel(result)
#
#################################################


#################################################
##
## Debug
##
#################################################

DEBUG <- FALSE

if(DEBUG) {
	rm(list=ls(all=TRUE))
}


########################################
##
## Class IncRegFinalmodelResult
## in this class, the results will be saved
##
#########################################

setClass(Class = "IncRegFinalmodelResult",	
		representation = representation(
				.multiCollResult = "list",
				.collResult = "matrix",
				.modelResult = "list"
		)
)


########################################
##
## Class IncRegOption
## contains all the values ​​that can be set as an option
##
#########################################

setClass(Class = "IncRegOption",	
		representation = representation(
				.shouldScaling = "logical",
				.shouldCentering = "logical",
				.missingValueTyp = "character", # "rm" or "avg"
				.missingValueIsCol = "logical", #if missingValueTyp is "rm" then this value is important: TRUE -> column will be removed, FALSE -> the row
				.cutlineZeroTest = "numeric",
				.cutlineAdvZeroTest = "numeric",
				.cutlineCorr = "numeric",
				.cutlineMulticoll = "numeric",
				.alphaValue = "numeric"
		
		
		),
		prototyp = prototype(
				.shouldScaling = FALSE,
				.shouldCentering = FALSE,
				.missingValueTyp = "rm", # "rm" or "avg"
				.missingValueIsCol = TRUE, # TRUE column will be deleted , FALSE the row
				.cutlineZeroTest = 0.8,
				.cutlineAdvZeroTest = 0.25,
				.cutlineCorr = 0.49,
				.cutlineMulticoll = 10,
				.alphaValue = 0.05
		),
		
		validity = function(object) {
				
			if(!(object@.missingValueTyp == "avg" || object@.missingValueTyp == "rm")) {
				errorPrintOfTheClasses("missingValueTyp", FALSE)
				stop()
			}
			
			if(object@.cutlineZeroTest > 1.0 || object@.cutlineZeroTest < 0.0) {
				errorPrintOfTheClasses("cutlineZeroTest", FALSE)
				stop()
			}
			
			if(object@.cutlineAdvZeroTest > 1.0 || object@.cutlineAdvZeroTest < 0.0) {
				errorPrintOfTheClasses("cutlineAdvZeroTest", FALSE)
				stop()
			}
			
			if(object@.cutlineCorr > 1.0 || object@.cutlineCorr < 0.0) {
				errorPrintOfTheClasses("cuttlineCorr", FALSE)
				stop()
			}
			
			if(object@.cutlineMulticoll < 0.0) {
				errorPrintOfTheClasses("cutlineMulticoll", FALSE)
				stop()
			}
			
			if(object@.alphaValue > 1.0 || object@.alphaValue < 0.0) {
				errorPrintOfTheClasses("alphaValue", FALSE)
				stop()
			}
			return(TRUE)
		}
)


########################################
##
## Class IncReg
## The main class
##
#########################################

setClass(Class = "IncReg",	
		representation = representation(
				.dataSet = "data.frame",
				.csvFile = "character", # if this value is set, then the dataSet variable will be ignored
				.yName = "character"
		),
		
		contains = c("IncRegOption", "IncRegFinalmodelResult"),
		
		prototyp = prototype(
				.csvFile = "None"		
		),
		
		validity = function(object) {
			
			if(nchar(object@.csvFile) > 0 && object@.csvFile != "None") {
				if(!file.exists(object@.csvFile)) {
					stop("... csv-file dont exists!")
				}
			}
				
			if(!(nrow(object@.dataSet) > 1)) {
				stop("... too few rows!")
			}
				
			if(!(ncol(object@.dataSet) > 1)) {
				stop("... too few cols!")
			}
			return(TRUE)
		}
)


########################################
##
## Constructor for users
##
#########################################

incReg <- function(yName,
				   dataSet = data.frame(),
				   csvFile = "None",
				   doAnalysis = TRUE,
				   shouldScaling = FALSE, 
				   shouldCentering = FALSE, 
				   missingValueTyp = "rm", 
				   missingValueIsCol = TRUE, 
				   cutlineZeroTest = 0.8, 
		   		   cutlineAdvZeroTest = 0.25, 
				   cutlineCorr = 0.49, 
				   cutlineMulticoll = 10, 
				   alphaValue = 0.05) {
	
	mm <- match.call(expand.dots = FALSE)
	mm <- checkAndSetDefaultValues(mm)
	
	if (csvFile != "None") {
		mm$dataSet <- mm$csvFile
		mm$csvFile <- csvFile
	}

	inRegObject <- new(Class = "IncReg", 
					.dataSet = ownIfElse(is.null(mm$dataSet), dataSet, mm$dataSet), 
					.csvFile = ownIfElse(is.null(mm$csvFile), csvFile, mm$csvFile), 
					.yName = ownIfElse(is.null(mm$yName), yName, mm$yName), 
					.shouldScaling = ownIfElse(is.null(mm$shouldScaling), shouldScaling, mm$shouldScaling), 
					.shouldCentering = ownIfElse(is.null(mm$shouldCentering), shouldCentering, mm$shouldCentering), 
					.missingValueTyp = ownIfElse(is.null(mm$missingValueTyp), missingValueTyp, mm$missingValueTyp),
					.missingValueIsCol = ownIfElse(is.null(mm$missingValueIsCol), missingValueIsCol, mm$missingValueIsCol),
					.cutlineZeroTest = ownIfElse(is.null(mm$cutlineZeroTest), cutlineZeroTest, mm$cutlineZeroTest),
					.cutlineAdvZeroTest = ownIfElse(is.null(mm$cutlineAdvZeroTest), cutlineAdvZeroTest, mm$cutlineAdvZeroTest),
					.cutlineCorr = ownIfElse(is.null(mm$cutlineCorr), cutlineCorr, mm$cutlineCorr),
					.cutlineMulticoll = ownIfElse(is.null(mm$cutlineMulticoll), cutlineMulticoll, mm$cutlineMulticoll),
					.alphaValue = ownIfElse(is.null(mm$alphaValue), alphaValue, mm$alphaValue))
	if(doAnalysis) {
		return(processing(inRegObject))
	} else {
		return(inRegObject)
	}
			
}

reanalyse <- function(object) {
	if(class(object) == "IncReg") {
		return(processing(object))
	} else {
		ownCat("Only a object from type \"IncReg\" can be reanalysed")
	}
}

#########################################
##
## Check and set the default values
##
## change some values at the validity part is not possible, so it has to be done before the object is created
## the validity check themselve should be there, for example when some one use the class directly 
## and not over the "incReg"-function way
##
#########################################

checkAndSetDefaultValues <- function(listOfData) {
	mmAdapted <- mapply(function(xx, yy) {
		switch(xx,
			
			dataSet = 	return(NULL),	
				
			csvFile =	if(nchar(yy) > 0) {
							error <- try(read.csv(yy, sep=";", dec=","))
							if(isTRUE(all.equal(class(error), "try-error"))) {
								stop("... problem with the csv-file!")
							} else {			
								return(as.data.frame(error))
							}
						},	
			
			missingValueTyp = if(!(yy == "avg" || yy == "rm")) {
									default <- "rm"
									errorPrintOfTheClasses("missingValueTyp", TRUE, default)
									return(default)
								},
						
			cutlineZeroTest = if(as.numeric(yy) > 1.0 || as.numeric(yy) < 0.0) {
									default <- 0.8
									errorPrintOfTheClasses("cutlineZeroTest", TRUE, default)
									return(default)
								},
						
			cutlineAdvZeroTest = if(as.numeric(yy) > 1.0 || as.numeric(yy) < 0.0) {
									default <- 0.25
									errorPrintOfTheClasses("cutlineAdvZeroTest", TRUE, default)
									return(default)
								 },
						
		 	cutlineCorr = if(as.numeric(yy) > 1.0 || as.numeric(yy) < 0.0) {
							  default <- 0.49
							  errorPrintOfTheClasses("cutlineCorr", TRUE, default)
							  return(default)
						  },
						
			cutlineMulticoll = if(as.numeric(yy) < 0.0) {
									default <- 10
									errorPrintOfTheClasses("cutlineMulticoll", TRUE, default)
									return(default)
								},
						
			alphaValue = if(as.numeric(yy) > 1.0 || as.numeric(yy) < 0.0) {
							default <- 0.05
							errorPrintOfTheClasses("alphaValue", TRUE, default)
							return(default)
						},
						
			suppressWarnings(ifelse(is.na(as.numeric(yy)), yy, as.numeric(yy)))
		)
	}, names(listOfData), as.character(as.list(listOfData)))

	return(mmAdapted)
}

#########################################
##
## Error Print
##
#########################################

errorPrintOfTheClasses <- function(descriptor, isStandard = FALSE, standardValue = NULL) {
	
	switch(descriptor,
			missingValueTyp = ownCat("\"missingValueTyp\": only \"avg\" or \"rm\" is possible", FALSE),
			cutlineZeroTest = ownCat("\"cutlineZeroTest\": only values between 0.0 and 1.0 are possible", FALSE),
			cutlineAdvZeroTest = ownCat("\"cutlineAdvZeroTest\": only values between 0.0 and 1.0 are possible", FALSE),
			cutlineCorr = ownCat("\"cuttlineCorr\": only values between 0.0 and 1.0 are possible", FALSE),
			cutlineMulticoll = ownCat("\"cutlineMulticoll\": only values greater than 0.0 are possible", FALSE),
			alphaValue = ownCat("\"alphaValue\": only values between 0.0 and 1.0 are possible", FALSE)
	)
	
	if(isStandard) {
		ownCat(paste(", standard value is now used (", standardValue, ")", sep=""), FALSE)
	} 
	ownCat("!")
}


#########################################
##
## Print and Show methods
##
#########################################

setMethod(
		f = "print",
		signature = "IncReg",
		function(x, ...) {
			cat(paste("[1] dataSet = ", x@.dataSet, "\n", sep="\""))
			cat(paste("[2] csvFile = ", "\"", x@.csvFile, "\"", "\n", sep=""))
			cat(paste("[3] yName = ", "\"",x@.yName, "\"", "\n", sep=""))
			cat(paste("[4] shouldScaling = ", x@.shouldScaling, "\n", sep=""))
			cat(paste("[5] shouldCentering = ", x@.shouldCentering, "\n", sep=""))
			cat(paste("[6] missingValueTyp = ", "\"", x@.missingValueTyp, "\"", "\n", sep=""))
			cat(paste("[7] missingValueIsCol = ", x@.missingValueIsCol, "\n", sep=""))
			cat(paste("[8] cutlineZeroTest = ", x@.cutlineZeroTest, "\n", sep=""))
			cat(paste("[9] cutlineAdvZeroTest = ", x@.cutlineAdvZeroTest, "\n", sep=""))
			cat(paste("[10] cutlineCorr = ", x@.cutlineCorr, "\n", sep=""))
			cat(paste("[11] cutlineMulticoll = ", x@.cutlineMulticoll, "\n", sep=""))
			cat(paste("[12] alphaValue = ", x@.alphaValue, "\n", sep=""))
			cat(paste("[13] multiCollResult = ", x@.multiCollResult, "\n", sep=""))
			cat(paste("[13] collResult = ", x@.collResult, "\n", sep=""))
			cat(paste("[14] modelResult = ", x@.modelResult, "\n", sep=""))
#			
		}	
)

setMethod(
		f = "show",
		signature = "IncReg",
		function(object) {
			print(object)
		}	
)

#########################################
##
## Getter methods (accessor-methods)
##
#########################################

setGeneric(name = "getDataSet",	def = function(object) {standardGeneric("getDataSet")})
setMethod(
		f = "getDataSet",
		signature = "IncReg",
		definition = function(object) {
			return(object@.dataSet)
		}
)

setGeneric(name = "getCsvFile",	def = function(object) {standardGeneric("getCsvFile")})
setMethod(
		f = "getCsvFile",
		signature = "IncReg",
		definition = function(object) {
			return(object@.csvFile)
		}
)

setGeneric(name = "getYname",	def = function(object) {standardGeneric("getYname")})
setMethod(
		f = "getYname",
		signature = "IncReg",
		definition = function(object) {
			return(object@.yName)
		}
)

setGeneric(name = "getShouldScaling",	def = function(object) {standardGeneric("getShouldScaling")})
setMethod(
		f = "getShouldScaling",
		signature = "IncReg",
		definition = function(object) {
			return(object@.shouldScaling)
		}
)

setGeneric(name = "getShouldCentering",	def = function(object) {standardGeneric("getShouldCentering")})
setMethod(
		f = "getShouldCentering",
		signature = "IncReg",
		definition = function(object) {
			return(object@.shouldCentering)
		}
)

setGeneric(name = "getMissingValueTyp",	def = function(object) {standardGeneric("getMissingValueTyp")})
setMethod(
		f = "getMissingValueTyp",
		signature = "IncReg",
		definition = function(object) {
			return(object@.missingValueTyp)
		}
)

setGeneric(name = "getMissingValueIsCol",	def = function(object) {standardGeneric("getMissingValueIsCol")})
setMethod(
		f = "getMissingValueIsCol",
		signature = "IncReg",
		definition = function(object) {
			return(object@.missingValueIsCol)
		}
)

setGeneric(name = "getCutlineZeroTest",	def = function(object) {standardGeneric("getCutlineZeroTest")})
setMethod(
		f = "getCutlineZeroTest",
		signature = "IncReg",
		definition = function(object) {
			return(object@.cutlineZeroTest)
		}
)

setGeneric(name = "getCutlineAdvZeroTest",	def = function(object) {standardGeneric("getCutlineAdvZeroTest")})
setMethod(
		f = "getCutlineAdvZeroTest",
		signature = "IncReg",
		definition = function(object) {
			return(object@.cutlineAdvZeroTest)
		}
)

setGeneric(name = "getCutlineCorr",	def = function(object) {standardGeneric("getCutlineCorr")})
setMethod(
		f = "getCutlineCorr",
		signature = "IncReg",
		definition = function(object) {
			return(object@.cutlineCorr)
		}
)

setGeneric(name = "getCutlineMulticoll",	def = function(object) {standardGeneric("getCutlineMulticoll")})
setMethod(
		f = "getCutlineMulticoll",
		signature = "IncReg",
		definition = function(object) {
			return(object@.cutlineMulticoll)
		}
)

setGeneric(name = "getAlphaValue",	def = function(object) {standardGeneric("getAlphaValue")})
setMethod(
		f = "getAlphaValue",
		signature = "IncReg",
		definition = function(object) {
			return(object@.alphaValue)
		}
)

setGeneric(name = "getMultiCollResult",	def = function(object) {standardGeneric("getMultiCollResult")})
setMethod(
		f = "getMultiCollResult",
		signature = "IncReg",
		definition = function(object) {
			return(object@.multiCollResult)
		}
)

setGeneric(name = "getCollResult",	def = function(object) {standardGeneric("getCollResult")})
setMethod(
		f = "getCollResult",
		signature = "IncReg",
		definition = function(object) {
			return(object@.collResult)
		}
)

setGeneric(name = "getModelResult",	def = function(object) {standardGeneric("getModelResult")})
setMethod(
		f = "getModelResult",
		signature = "IncReg",
		definition = function(object) {
			return(object@.modelResult)
		}
)

setGeneric(name = "getSquareBracket",	def = function(x,i, j, drop) {standardGeneric("getSquareBracket")})
setMethod(
		f = "getSquareBracket",
		signature = "IncReg",
		definition = function(x,i, j, drop) {
			if(i == "dataSet" || i == 1) {return(getDataSet(x))}
			if(i == "csvFile" || i == 2) {return(getCsvFile(x))}
			if(i == "yName" || i == 3) {return(getYname(x))}
			if(i == "shouldScaling" || i == 4) {return(getShouldScaling(x))}
			if(i == "shouldCentering" || i == 5) {return(getShouldCentering(x))}
			if(i == "missingValueTyp" || i == 6) {return(getMissingValueTyp(x))}
			if(i == "missingValueIsCol" || i == 7) {return(getMissingValueIsCol(x))}
			if(i == "cutlineZeroTest" || i == 8) {return(getCutlineZeroTest(x))}
			if(i == "cutlineAdvZeroTest" || i == 9) {return(getCutlineAdvZeroTest(x))}
			if(i == "cutlineCorr" || i == 10) {return(getCutlineCorr(x))}
			if(i == "cutlineMulticoll" || i == 11) {return(getCutlineMulticoll(x))}
			if(i == "alphaValue" || i == 12) {return(getAlphaValue(x))}
			if(i == "multiCollResult" || i == 13) {return(getMultiCollResult(x))}
			if(i == "collResult" || i == 14) {return(getCollResult(x))}
			if(i == "modelResult" || i == 15) {return(getModelResult(x))}
		}
)

setMethod(
		f = "[",
		signature = "IncReg",
		definition = function(x, i, j, drop) {
			return(getSquareBracket(x, i, j, drop))
		}
)


#########################################
##
## Setter methods (replacement-methods)
##
#########################################

setGeneric(name = "setDataSet<-",	def = function(object, value) {standardGeneric("setDataSet<-")})
setReplaceMethod(
		f = "setDataSet",
		signature = "IncReg",
		definition = function(object, value) {
			object@.dataSet <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setCsvFile<-",	def = function(object, value) {standardGeneric("setCsvFile<-")})
setReplaceMethod(
		f = "setCsvFile",
		signature = "IncReg",
		definition = function(object, value) {
			object@.csvFile <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setYname<-",	def = function(object, value) {standardGeneric("setYname<-")})
setReplaceMethod(
		f = "setYname",
		signature = "IncReg",
		definition = function(object, value) {
			object@.yName <- value
			validObject(object)
			return(object)
		}
)



setGeneric(name = "setShouldScaling<-",	def = function(object, value) {standardGeneric("setShouldScaling<-")})
setReplaceMethod(
		f = "setShouldScaling",
		signature = "IncReg",
		definition = function(object, value) {
			object@.shouldScaling <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setShouldCentering<-",	def = function(object, value) {standardGeneric("setShouldCentering<-")})
setReplaceMethod(
		f = "setShouldCentering",
		signature = "IncReg",
		definition = function(object, value) {
			object@.shouldCentering <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setMissingValueTyp<-",	def = function(object, value) {standardGeneric("setMissingValueTyp<-")})
setReplaceMethod(
		f = "setMissingValueTyp",
		signature = "IncReg",
		definition = function(object, value) {
			object@.missingValueTyp <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setMissingValueIsCol<-",	def = function(object, value) {standardGeneric("setMissingValueIsCol<-")})
setReplaceMethod(
		f = "setMissingValueIsCol",
		signature = "IncReg",
		definition = function(object, value) {
			object@.missingValueIsCol <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setCutlineZeroTest<-",	def = function(object, value) {standardGeneric("setCutlineZeroTest<-")})
setReplaceMethod(
		f = "setCutlineZeroTest",
		signature = "IncReg",
		definition = function(object, value) {
			object@.cutlineZeroTest <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setCutlineAdvZeroTest<-",	def = function(object, value) {standardGeneric("setCutlineAdvZeroTest<-")})
setReplaceMethod(
		f = "setCutlineAdvZeroTest",
		signature = "IncReg",
		definition = function(object, value) {
			object@.cutlineAdvZeroTest <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setCutlineCorr<-",	def = function(object, value) {standardGeneric("setCutlineCorr<-")})
setReplaceMethod(
		f = "setCutlineCorr",
		signature = "IncReg",
		definition = function(object, value) {
			object@.cutlineCorr <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setCutlineMulticoll<-",	def = function(object, value) {standardGeneric("setCutlineMulticoll<-")})
setReplaceMethod(
		f = "setCutlineMulticoll",
		signature = "IncReg",
		definition = function(object, value) {
			object@.cutlineMulticoll <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setAlphaValue<-",	def = function(object, value) {standardGeneric("setAlphaValue<-")})
setReplaceMethod(
		f = "setAlphaValue",
		signature = "IncReg",
		definition = function(object, value) {
			object@.alphaValue <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setMultiCollResult<-",	def = function(object, value) {standardGeneric("setMultiCollResult<-")})
setReplaceMethod(
		f = "setMultiCollResult",
		signature = "IncReg",
		definition = function(object, value) {
			object@.multiCollResult <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setCollResult<-",	def = function(object, value) {standardGeneric("setCollResult<-")})
setReplaceMethod(
		f = "setCollResult",
		signature = "IncReg",
		definition = function(object, value) {
			object@.collResult <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setModelResult<-",	def = function(object, value) {standardGeneric("setModelResult<-")})
setReplaceMethod(
		f = "setModelResult",
		signature = "IncReg",
		definition = function(object, value) {
			object@.modelResult <- value
			validObject(object)
			return(object)
		}
)

setGeneric(name = "setSquareBracket",	def = function(x,i, j, value) {standardGeneric("setSquareBracket")})
setMethod(
		f = "setSquareBracket",
		signature = "IncReg",
		definition = function(x,i, j, value) {			
			if(i == "dataSet" || i == 1) {setDataSet(x) <- value}
			if(i == "csvFile" || i == 2) {setCsvFile(x) <- value}
			if(i == "yName" || i == 3) {setYname(x) <- value}
			if(i == "shouldScaling" || i == 4) {setShouldScaling(x) <- value}
			if(i == "shouldCentering" || i == 5) {setShouldCentering(x) <- value}
			if(i == "missingValueTyp" || i == 6) {setMissingValueTyp(x) <- value}
			if(i == "missingValueIsCol" || i == 7) {setMissingValueIsCol(x) <- value}
			if(i == "cutlineZeroTest" || i == 8) {setCutlineZeroTest(x) <- value}
			if(i == "cutlineAdvZeroTest" || i == 9) {setCutlineAdvZeroTest(x) <- value}
			if(i == "cutlineCorr" || i == 10) {setCutlineCorr(x) <- value}
			if(i == "cutlineMulticoll" || i == 11) {setCutlineMulticoll(x) <- value}
			if(i == "alphaValue" || i == 12) {setAlphaValue(x) <- value}
			if(i == "multiCollResult" || i == 13) {setMultiCollResult(x) <- value}
			if(i == "collResult" || i == 14) {setCollResult(x) <- value}
			if(i == "modelResult" || i == 15) {setModelResult(x) <- value}
			return(x)
		}
)

setReplaceMethod(
		f = "[",
		signature = "IncReg",
		definition = function(x, i ,j , value) {
			return(setSquareBracket(x, i, j, value))
		}
)


#########################################
##
## Processing
##
## The whole processing
## First the variables will be checked
## After that the model will be calculated
##
#########################################

setGeneric(name = "processing", def = function(IncReg, ...) {standardGeneric("processing")})
setMethod(
		f = "processing",
		signature = "IncReg",
		definition = function(IncReg, ...) {	
			mm <- match.call(expand.dots = FALSE)
			mm$IncReg <- eval(mm$IncReg, parent.frame())
	
			# 1. Step: covert all values in the dataframe to numeric values, all other values convered to NAs
			mm$IncReg["dataSet"] <- convertAllToNumeric(mm$IncReg["dataSet"])

			# 2. Step: delete alle rows and cols with NA or calculate the mean for each cell (special handling for the y column)
			mm$IncReg["dataSet"] <- deletRowsAndColumsAndOrCalculatMean(mm$IncReg["dataSet"], mm$IncReg["yName"], mm$IncReg["missingValueTyp"], mm$IncReg["missingValueIsCol"])
			
			# 3. Step: Zero- und advanced Zero Test, check if the single columns have enough variations(Zero-Test: not more than the given cutline should be zero, ...)
			mm$IncReg["dataSet"] <- zeroAndAdvancedZeroTest(mm$IncReg["dataSet"], mm$IncReg["yName"], mm$IncReg["cutlineZeroTest"], mm$IncReg["cutlineAdvZeroTest"])

			# 4. Step: multicollinearity, check if exits between more than two descriptors dependencies
			multiResult <- multicollinearityTest(mm$IncReg["dataSet"], mm$IncReg["yName"], mm$IncReg["cutlineMulticoll"])
			mm$IncReg["dataSet"] <- multiResult$dataSet	
			if(!is.null(multiResult$result)) {
				mm$IncReg[multiResult$resultTyp] <- multiResult$result
			} 
			
			# 5. Step: correlation, check if exits between two descriptors dependencies
			multiResult <- correlationTest(mm$IncReg["dataSet"], mm$IncReg["yName"], mm$IncReg["cutlineCorr"])
			mm$IncReg["dataSet"] <- multiResult$dataSet	
			if(!is.null(multiResult$result)) {
				mm$IncReg[multiResult$resultTyp] <- multiResult$result
			} 

			# 6. Step: scaling and centering the dataset
			mm$IncReg["dataSet"] <- scalingAndCenteringData(mm$IncReg["dataSet"], mm$IncReg["shouldScaling"], mm$IncReg["shouldCentering"])
			
			# 7. Step: Prework is done now do the incremental Multivariate Regression
			mm$IncReg["modelResult"] <- incrementalMultivariateRegression(mm$IncReg["dataSet"], mm$IncReg["yName"], mm$IncReg["alphaValue"])
			
			return(mm$IncReg)
		}
)


#########################################
##
## convertAllToNumeric
##
## convert all not numeric values to numeric or NA
## thereby character convert to NAs
## also empty values like '' will be convert to NAs
##
#########################################

setGeneric(name = "convertAllToNumeric", def = function(dataSet, ...) {standardGeneric("convertAllToNumeric")})
setMethod(
		f = "convertAllToNumeric",
		signature = "data.frame",
		definition = function(dataSet, ...) {	
			for(nn in seq(along=dataSet)) {
				if(class(dataSet[,nn]) == "factor") {
					suppressWarnings(dataSet[,nn] <- as.numeric(levels(dataSet[,nn]))[dataSet[,nn]])
				} else if (class(dataSet[,nn]) != "numeric") {
					suppressWarnings(dataSet[,nn] <- as.numeric(dataSet[,nn]))
				}
			}
			return(dataSet)
		}
)

setMethod(
		f = "convertAllToNumeric",
		signature = "matrix",
		definition = function(dataSet, ...) {
			return(convertAllToNumeric(as.data.frame(dataSet, ...)))
		}
)

setMethod(
		f = "convertAllToNumeric",
		signature = "vector",
		definition = function(dataSet, ...) {
			return(convertAllToNumeric(t(t(dataSet)), ...))
		}
)

setMethod(
		f = "convertAllToNumeric",
		signature = "numeric",
		definition = function(dataSet, ...) {
			return(convertAllToNumeric(t(t(dataSet)), ...))
		}
)


#########################################
##
## deletRowsAndColumsAndOrCalculatMean
##
## Delete all rows and cols with NA or calculate the mean for each cell (special handling for the y column)
## 
## important options
## missingValueTyp: remove ("rm") or average ("avg")
## missingValueIsCol: if missingValueTyp is "rm" then this value is important: TRUE -> column will be removed, FALSE -> the row
##
#########################################

setGeneric(name = "deletRowsAndColumsAndOrCalculatMean", def = function(dataSet, yName, missingValueTyp, missingValueIsCol, ...) {standardGeneric("deletRowsAndColumsAndOrCalculatMean")})
setMethod(
		f = "deletRowsAndColumsAndOrCalculatMean",
		signature(dataSet = "data.frame"),
		definition = function(dataSet, yName, missingValueTyp, missingValueIsCol, ...) {	
				
			# TRUE when all cells in one row are NA (without the y column)
			boolVrow <- apply(is.na(getValuesWithoutYAsDataFrame(dataSet, yName)), 1, "all")
			names(boolVrow) <- rownames(dataSet)
			
			# TRUE when all NAs in the y column in combination with the line above
			boolVrow <- boolVrow | is.na(dataSet[,yName])
			
			# TRUE when all cells in one column are NA 
			boolVcol <- apply(is.na(dataSet),2,"all")
			
			if(all(boolVrow) || all(boolVcol[!(colnames(dataSet) %in% yName)]) || boolVcol[yName]) {
				stop("... all rows/cols are NA!")
			} else {
				if(any(boolVrow)) {
					printReduction("row", names(boolVrow)[boolVrow], "NA")
				}
				
				if(any(boolVcol)) {
					printReduction("col", names(boolVcol)[boolVcol], "NA")
				}
				
				dataSet <- dataSet[!boolVrow,!boolVcol]
				if(!missingValueIsCol && missingValueTyp == "rm") { # when at least one value missed in a column than delete the whole column
					boolV <- apply(is.na(dataSet),1,"any")
					if(all(boolV)) {
						stop("... all rows are NA, based on the specification that, when one NA is in a row the whole row schould removed")
					} else {
						if(any(boolV)){
							printReduction("row", names(boolV)[boolV], "NA")
						}
						dataSet <- dataSet[!boolV,]
					}
				} else if(missingValueIsCol && missingValueTyp == "rm") { #when at least one value missed in a row than delete the whole row
					boolV <- apply(is.na(dataSet),2,"any")
					if(all(boolV[!(colnames(dataSet) %in% yName)])) {
						stop("... all cols are NA, based on the specification that, when one NA is in a col the whole col schould removed")
					} else {
						if(any(boolV)){
							printReduction("col", names(boolV)[boolV], "NA")
						}
						dataSet <- dataSet[,!boolV]
					}
				} else if(missingValueTyp == "avg") { #calculate the mean of all given values in one column and use this value for the NAs in this column (only for the x-values (descriptors))
					dataSet <- as.data.frame(apply(dataSet, 2, function(x) {
										x[is.na(x)] <- mean(x, na.rm=TRUE)
										return(x)
									}))
				}
				return(dataSet)	
			}		
		}
)

setMethod(
		f = "deletRowsAndColumsAndOrCalculatMean",
		signature(dataSet = "matrix"),
		definition = function(dataSet, yName, missingValueTyp, missingValueIsCol, ...) {
			return(deletRowsAndColumsAndOrCalculatMean(as.data.frame(dataSet), yName, missingValueTyp, missingValueIsCol, ...))
		}
)

setMethod(
		f = "deletRowsAndColumsAndOrCalculatMean",
		signature(dataSet = "vector"),
		definition = function(dataSet, yName, missingValueTyp, missingValueIsCol, ...) {
			return(deletRowsAndColumsAndOrCalculatMean(t(t(dataSet)), yName, missingValueTyp, missingValueIsCol, ...))
		}
)

setMethod(
		f = "deletRowsAndColumsAndOrCalculatMean",
		signature(dataSet = "numeric"),
		definition = function(dataSet, yName, missingValueTyp, missingValueIsCol, ...) {
			return(deletRowsAndColumsAndOrCalculatMean(t(t(dataSet)), yName, missingValueTyp, missingValueIsCol, ...))
		}
)


#########################################
##
## zeroAndAdvancedZeroTest
##
## Check if the single columns have enough variations
## 
## important options
## cutlineZeroTest: not more than the given cutline should be zero
## cutlineAdvZeroTest: there have to be no descriptor which have more than the given cutline on the same cases
##
#########################################

setGeneric(name = "zeroAndAdvancedZeroTest", def = function(dataSet, yName, cutlineZeroTest, cutlineAdvZeroTest, ...) {standardGeneric("zeroAndAdvancedZeroTest")})
setMethod(
		f = "zeroAndAdvancedZeroTest",
		signature(dataSet = "data.frame"),
		definition = function(dataSet, yName, cutlineZeroTest, cutlineAdvZeroTest, ...) {	
						
			boolV <- apply(getValuesWithoutYAsDataFrame(dataSet, yName), 2, function(xx, cutZer, cutAdvZer) {
						
							divisionTable <- table(xx);
							divisionTable <- divisionTable / sum(divisionTable)			
							
							if(!is.na(divisionTable["0"])) {
								if(divisionTable["0"] > cutZer) {
									return(FALSE)
								} else {
									divisionTable <- divisionTable[!(names(divisionTable) %in% "1")]
								}
							} 
							
							if(sum(divisionTable > cutAdvZer, na.rm=TRUE) > 0) {
								return(FALSE)
							} else {
								return(TRUE)
							}
						}, cutlineZeroTest, cutlineAdvZeroTest)
			
			if(all(!boolV)) {
				stop("... all columns dont fulfill the zero or zero advanced test!")
			} else {
				if(any(!boolV)) {
					printReduction("col", names(boolV)[!boolV], "zero")
				}
				
				boolV <- (colnames(dataSet) %in% yName) | (colnames(dataSet) %in% names(boolV)[boolV]) 
				dataSet <- dataSet[,boolV]	
					
				return(dataSet)
			}
		}
)

setMethod(
		f = "zeroAndAdvancedZeroTest",
		signature(dataSet = "matrix"),
		definition = function(dataSet, yName, cutlineZeroTest, cutlineAdvZeroTest, ...) {	
			return(zeroAndAdvancedZeroTest(as.data.frame(dataSet), yName, cutlineZeroTest, cutlineAdvZeroTest, ...))
		}
)

setMethod(
		f = "zeroAndAdvancedZeroTest",
		signature(dataSet = "vector"),
		definition = function(dataSet, yName, cutlineZeroTest, cutlineAdvZeroTest, ...) {	
			return(zeroAndAdvancedZeroTest(t(t(dataSet)), yName, cutlineZeroTest, cutlineAdvZeroTest, ...))
		}
)

setMethod(
		f = "zeroAndAdvancedZeroTest",
		signature(dataSet = "numeric"),
		definition = function(dataSet, yName, cutlineZeroTest, cutlineAdvZeroTest, ...) {	
			return(zeroAndAdvancedZeroTest(t(t(dataSet)), yName, cutlineZeroTest, cutlineAdvZeroTest, ...))
		}
)


#########################################
##
## multicollinearityTest
##
## check if a multicollinearity exits between more than two descriptors. 
## If more than two descriptors have multicol. then the descriptor with the lowest correlation to "y" removed. 
## After that multicollinearity is recalculated.
## 
## important option
## cutlineMulticoll: if a descriptor have a higher value than this cutline, multicollinearity is there
##
#########################################

setGeneric(name = "multicollinearityTest", def = function(dataSet, yName, cutlineMulticoll, ...) {standardGeneric("multicollinearityTest")})
setMethod(
		f = "multicollinearityTest",
		signature(dataSet = "data.frame"),
		definition = function(dataSet, yName, cutlineMulticoll, ...) {	
			
			if(dim(getValuesWithoutYAsDataFrame(dataSet, yName))[2] > 1) {
				resultList <- list()
				rr <- 0
				
#				library("car")
				
				repeat {
					rr <- rr + 1
					fmla <- getFormula(dataSet, yName)
					result <- vif(lm(fmla, data=dataSet))			
					boolV <- (result > cutlineMulticoll)
					resultList[[rr]] <- result
					if(!any(boolV)) {
						break
					} else {
						if(sum(boolV) > 1) {
							correlationWithY <- cor(dataSet[,yName], getValuesWithoutYAsDataFrame(dataSet, yName), method="pearson")
							orderedValues <- correlationWithY[,boolV][order(correlationWithY[,boolV], na.last = TRUE,decreasing = TRUE)]
							boolV <- names(boolV) %in% names(orderedValues)[1]
							names(boolV) <- names(result)
						}
 						printReduction("col", names(boolV)[boolV], because = "multiCol")
						
						boolV <- (colnames(dataSet) %in% yName) | (colnames(dataSet) %in% names(boolV)[!boolV]) 
						dataSet <- dataSet[,boolV]	
						
						if(dim(getValuesWithoutYAsDataFrame(dataSet, yName))[2] <= 1) {
							break	
						}
					}	
				}
				
				return(list(dataSet = dataSet, result = resultList, resultTyp = "multiCollResult"))
				
			} else {
				return(list(dataSet = dataSet, result = NULL, resultTyp = NULL))
			}
		}
)

setMethod(
		f = "multicollinearityTest",
		signature(dataSet = "matrix"),
		definition = function(dataSet, yName, cutlineMulticoll, ...) {	
			return(multicollinearityTest(as.data.frame(dataSet), yName, cutlineMulticoll, ...))
		}
)

setMethod(
		f = "multicollinearityTest",
		signature(dataSet = "vector"),
		definition = function(dataSet, yName, cutlineMulticoll, ...) {	
			return(multicollinearityTest(t(t(dataSet)), yName, cutlineMulticoll, ...))
		}
)

setMethod(
		f = "multicollinearityTest",
		signature(dataSet = "numeric"),
		definition = function(dataSet, yName, cutlineMulticoll, ...) {	
			return(multicollinearityTest(t(t(dataSet)), yName, cutlineMulticoll, ...))
		}
)


#########################################
##
## correlationTest
##
## check if a correlation exits between two descriptors
## 
## important option
## cutlineCorr: the column, from a correlation pair (value more than the cutline), which has less correlation against the y column is removed
##
#########################################

setGeneric(name = "correlationTest", def = function(dataSet, yName, cutlineCorr, ...) {standardGeneric("correlationTest")})
setMethod(
		f = "correlationTest",
		signature(dataSet = "data.frame"),
		definition = function(dataSet, yName, cutlineCorr, ...) {	
			
			corrAll <- cor(dataSet, method="pearson")
			correlationWithY <- cor(dataSet[,yName], getValuesWithoutYAsDataFrame(dataSet, yName), method="pearson")
			corMatrixBoolean <- abs(cor(getValuesWithoutYAsDataFrame(dataSet, yName), method="pearson")) > cutlineCorr
			corMatrixValues <- cor(getValuesWithoutYAsDataFrame(dataSet, yName), method="pearson")
			corDataFrame <- data.frame(xx=character(0), yy=character(0), cor=numeric(0))
			
			for(xx in seq(along=corMatrixBoolean[1,])) {
				for(yy in seq(along=corMatrixBoolean[,1])) {
					if(corMatrixBoolean[xx,yy] && 
							rownames(corMatrixBoolean)[xx] != colnames(corMatrixBoolean)[yy] &&
							!(xx %in% corDataFrame[,"yy"] && yy %in% corDataFrame[,"xx"])) {
						corDataFrame <- rbind(corDataFrame, data.frame(xx = xx, yy = yy, cor = corMatrixValues[xx,yy]))
					}
				}
			}
						
			if(length(corDataFrame[,1]) > 0) {
				corDataFrame <- corDataFrame[order(corDataFrame$cor),]
				removeBoolean <- rep.int(FALSE, ncol(dataSet)) 
				
				repeat {
					if(length(corDataFrame$xx) == 0) {
						break
					}	
					
					if(abs(correlationWithY[corDataFrame$xx[1]]) >= abs(correlationWithY[corDataFrame$yy[1]])) {
						## remove the column from yy
						removeBoolean[(corDataFrame$yy[1]+1)] <- TRUE
						corDataFrame <- reduceTheDataFrame(corDataFrame, corDataFrame$yy[1])
					} else {
						## remove the column from xx
						removeBoolean[(corDataFrame$xx[1]+1)] <- TRUE
						corDataFrame <- reduceTheDataFrame(corDataFrame, corDataFrame$xx[1])
					}
					
				}

				ownCat("... following columns removed based on the correlation alignment")
				ownCat(paste("... first all correlated columns larger than ", cutlineCorr, " are sorted", sep=""))
				ownCat("... after that, the column, from a correlation pair, which has less correlation agains the y column is removed")
				ownCat(paste("... removed column(s): ", colnames(dataSet)[removeBoolean], sep=""))
				
				dataSet <- dataSet[, !removeBoolean]
			}
			
			if(is.null(corMatrixValues)) {
				return(list(dataSet = dataSet, result = NULL, resultTyp = NULL))
			} else {
				return(list(dataSet = dataSet, result = corrAll, resultTyp = "collResult"))
			}
		}
)

setMethod(
		f = "correlationTest",
		signature(dataSet = "matrix"),
		definition = function(dataSet, yName, cutlineCorr, ...) {	
			return(correlationTest(as.data.frame(dataSet), yName, cutlineCorr, ...))
		}
)

setMethod(
		f = "correlationTest",
		signature(dataSet = "vector"),
		definition = function(dataSet, yName, cutlineCorr, ...) {	
			return(correlationTest(t(t(dataSet)), yName, cutlineCorr, ...))
		}
)

setMethod(
		f = "correlationTest",
		signature(dataSet = "numeric"),
		definition = function(dataSet, yName, cutlineCorr, ...) {	
			return(correlationTest(t(t(dataSet)), yName, cutlineCorr, ...))
		}
)


#########################################
##
## scalingAndCenteringData
##
## Scaling and centering of the dataset
## 
## important option
## shouldScaling: TRUE or FALSE
## shouldCentering: TRUE or FALSE
##
#########################################

setGeneric(name = "scalingAndCenteringData", def = function(dataSet, shouldScaling, shouldCentering, ...) {standardGeneric("scalingAndCenteringData")})
setMethod(
		f = "scalingAndCenteringData",
		signature(dataSet = "data.frame"),
		definition = function(dataSet, shouldScaling, shouldCentering, ...) {	
			
			if(shouldScaling || shouldCentering) {
				dataSet <- as.data.frame(scale(dataSet, center = shouldCentering, scale = shouldScaling))
			}
			
			printScaleAndCentering("scale", shouldScaling)
			printScaleAndCentering("center", shouldCentering)
			
			return(dataSet)
		}
)

setMethod(
		f = "scalingAndCenteringData",
		signature(dataSet = "matrix"),
		definition = function(dataSet, shouldScaling, shouldCentering, ...) {	
			return(scalingAndCenteringData(as.data.frame(dataSet), shouldScaling, shouldCentering, ...))
		}
)

setMethod(
		f = "scalingAndCenteringData",
		signature(dataSet = "vector"),
		definition = function(dataSet, shouldScaling, shouldCentering, ...) {	
			return(scalingAndCenteringData(t(t(dataSet)), shouldScaling, shouldCentering, ...))
		}
)

setMethod(
		f = "scalingAndCenteringData",
		signature(dataSet = "numeric"),
		definition = function(dataSet, shouldScaling, shouldCentering, ...) {	
			return(scalingAndCenteringData(t(t(dataSet)), shouldScaling, shouldCentering, ...))
		}
)


#########################################
##
## incrementalMultivariateRegression
##
## Search the best multivariate regression modell of the given dataset. For this, first a modell with all descriptors calculate, after that the descriptor with the bigest F-value removed 
## and a new modell will be calculate. Finished when all F-values are less than the given alpha value or all descriptors removed
## 
## important option
## alphaValue: significance level
##
#########################################


setGeneric(name = "incrementalMultivariateRegression", def = function(dataSet, yName, alphaValue, ...) {standardGeneric("incrementalMultivariateRegression")})
setMethod(
		f = "incrementalMultivariateRegression",
		signature(dataSet = "data.frame"),
		definition = function(dataSet, yName, alphaValue, ...) {	
			
			options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))
#			library("car")
			resultList <- list()
			rr <- 0
			
			repeat {
				rr <- rr + 1
				ownCat(paste("Model:",rr))
				repeat {
					fmla <- getFormula(dataSet, yName)
					model <- lm(formula = fmla, data = dataSet) # singular.ok=TRUE
					if(sum(is.na(model$coefficients)) > 0) {
						ownCat("... singularities found, following columns are removed:")
						ownCat(colnames(dataSet)[is.na(model$coefficients)])
						
						dataSet <- dataSet[,!is.na(model$coefficients)]
						
					} else {
						break
					}
				}
				anovaResult <- Anova(model, typ="III", singular.ok=TRUE)
				
				booleanV <- sapply(rownames(anovaResult), function(x,y) {
							if(x  %in% y) {
								return(FALSE)
							} else {
								return(TRUE)
							}				
						}, c("(Intercept)", "Residuals"), USE.NAMES=FALSE)
				
				anovaResultSmall <- anovaResult[booleanV,]
				anovaResultSmall <- anovaResultSmall[order(anovaResultSmall$`Pr(>F)`, na.last=TRUE,decreasing = TRUE),]
				
				if(anovaResultSmall$`Pr(>F)`[1] > alphaValue) {
					ownCat(paste("... no signification! worst p-value = ", round(anovaResultSmall$`Pr(>F)`[1], 3)," but alpha = ", alphaValue, ", following descriptor removed from model: ", rownames(anovaResultSmall)[1], sep=""))
					resultList[[rr]] <- list(formula = fmla, model = model, anova = anovaResult, final = FALSE)
					
					if(length(anovaResultSmall$`Pr(>F)`) > 1) {
						dataSet <- dataSet[, !(colnames(dataSet) %in% rownames(anovaResultSmall)[1])]
					} else {
						ownCat(paste("... found no model with significantly values", sep=""))
						return(resultList)
						break;
					}
				} else {
					resultList[[rr]] <- list(formula = fmla, model = model, anova = anovaResult, final = TRUE)
					ownCat(paste("... the used descriptors for the model are significantly (p-value = ", round(anovaResultSmall$`Pr(>F)`[1], 3), " alpha = ", alphaValue, ")", sep=""))
					ownCat("... the new model is: ")
					print(model$coefficients)
					ownCat("for more information use getFinalModel(object); getModelNumberX(object, number); summary(object); getPreModelResults(object)!")
					return(resultList)
					break
				}
			}
		}
)

setMethod(
		f = "incrementalMultivariateRegression",
		signature(dataSet = "matrix"),
		definition = function(dataSet, yName, alphaValue, ...) {	
			return(incrementalMultivariateRegression(as.data.frame(dataSet), yName, alphaValue, ...))
		}
)

setMethod(
		f = "incrementalMultivariateRegression",
		signature(dataSet = "vector"),
		definition = function(dataSet, yName, alphaValue, ...) {	
			return(incrementalMultivariateRegression(t(t(dataSet)), yName, alphaValue, ...))
		}
)

setMethod(
		f = "incrementalMultivariateRegression",
		signature(dataSet = "numeric"),
		definition = function(dataSet, yName, alphaValue, ...) {	
			return(incrementalMultivariateRegression(t(t(dataSet)), yName, alphaValue, ...))
		}
)


#########################################
##
## methods for print/get results
##
#########################################

setGeneric(name = "getPreModelResults", def = function(object) {standardGeneric("getPreModelResults")})
setMethod(
		f = "getPreModelResults",
		signature = "IncReg",
		definition = function(object) {		
			return(list(multiCollResult = getMultiCollResult(object), collResult = getCollResult(object)))
		}
)

setGeneric(name = "getModelNumberX", def = function(object, number) {standardGeneric("getModelNumberX")})
setMethod(
		f = "getModelNumberX",
		signature = "IncReg",
		definition = function(object, number) {
			modelResult <- getModelResult(object)
			if(length(modelResult) > 0) {
				if(number <= length(modelResult) && number > 0) {
					return(modelResult[[number]])
				} else {
					ownCat(paste("Found no model with number ", number, sep=""))
				}
			} else {
				ownCat("Found no model")
			}
		}
)

setGeneric(name = "getFinalModel", def = function(object) {standardGeneric("getFinalModel")})
setMethod(
		f = "getFinalModel",
		signature = "IncReg",
		definition = function(object) {
			modelResult <- getModelResult(object)
			if(length(modelResult) > 0) {
				if(modelResult[[length(modelResult)]]$final) {
					return(modelResult[[length(modelResult)]])
				} else {
					ownCat("Found no final model")
				}
			} else {
				ownCat("Found no model")
			}
		}
)

setGeneric(name = "summary", def = function(object) {standardGeneric("summary")})
setMethod(
		f = "summary",
		signature = "IncReg",
		definition = function(object) {
			value <- getFinalModel(object)
			if(!is.null(value)) {
				return(summary(value))
			} else {
				ownCat("Summary not possible")
			}
		}
)

setMethod(
		f = "summary",
		signature(object = "list"),
		definition = function(object) {		
			summary.lm(object$model)
		}
)

#########################################
##
## miscellaneous
##
#########################################

ownCat <- function(text, endline=TRUE, printAsDataFrame=FALSE){
	if(printAsDataFrame) {
		print(text)
	} else {
		cat(unlist(text))	
	}
	
	if (endline) {
		cat("\n")
	}
}

ownIfElse <- function(test, yes, no) {
	if(test) {
		return(yes)
	} else {
		return(no)
	}
}

printReduction <- function(typ, description, because = "NA") {
	if(because == "zero") {
		because <- "zero or zero advanced test"
	} else if(because == "NA") {
		because <- "NA values"
	} else if(because == "multiCol") {
		because <- "multicollinearity test"
	}

	ownCat(paste("... following ", typ, "(s) are removed (because of ", because,"): ", paste(description, collapse=", "), sep=""))
}

printScaleAndCentering <- function(typ, boolean) {
	text <- ""
	
	if(!boolean) {
		text <- "not "
	}
	
	if(typ == "scale") {
		text <- paste(text, "scaled", sep = "")
	} else {
		text <- paste(text, "centered", sep = "")
	}
		
	ownCat(paste("... the dataset has been ", text, sep=""))
}



"%checkAndStop%" <- function(toVerify, checkValue) {
	if(toVerify == checkValue) {
		stop(paste("... the termination condition is statisfied! checkValue: ", checkValue, sep=""))
	}
}

getValuesWithoutYAsDataFrame <- function(dataSet, yName) {

	if(sum(!(colnames(dataSet) %in% yName)) > 1) {
		return(dataSet[,!(colnames(dataSet) %in% yName)])
	} else {
		newDataSet <- t(t(dataSet[,!(colnames(dataSet) %in% yName)]))
		rownames(newDataSet) <- rownames(dataSet)
		colnames(newDataSet) <- colnames(dataSet)[!(colnames(dataSet) %in% yName)]
		return(newDataSet)
	}
}

getFormula <- function(dataSet, yName) {
	return(as.formula(paste(yName, " ~ ", paste(colnames(getValuesWithoutYAsDataFrame(dataSet, yName)), collapse= "+"))))
}


reduceTheDataFrame <- function(dataSet, removeColumnIndex) {
	
	if(removeColumnIndex %in% dataSet$yy) {
		dataSet <- dataSet[!(dataSet$yy %in% removeColumnIndex),]
	}		
	
	if(removeColumnIndex %in% dataSet$xx) {
		dataSet <- dataSet[!(dataSet$xx %in% removeColumnIndex),]
	}
	
	return(dataSet)
}
