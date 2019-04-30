## This file contains checks of input function parameters that are
## used by several other functions.

check_method <- function(method) {
    if (! method %in% c("pros", "retros")) {
        stop("Unkown method specified; the method parameter ",
             "should be one of: 'pros', ",
             ", or 'retros'" )
    }
}


check_copula <- function(copfit) {
  if (! copfit %in% c("Gaussian", "Student", "Gumbel", "Frank", "Clayton")) {
    stop("Unkown copula model specified; the copfit parameter ",
         "should be one of: 'Gaussian', 'Student', 'Gumbel', 'Frank'",
         ", or 'Clayton'" )
  }
}


check_link <- function(link, Design) {
  if ((Design=="CC") | (Design=="MT")){ 
  if (! link %in% c("probit","logit","cloglog")) {
    stop("Unkown link function specified; the link parameter ",
         "should be one of: 'probit', 'logit'",
         ", or 'cloglog'" )
  }
  }
}


check_pheno <- function(y) {
    if (is.null(y)) {
        stop("Argument y is empty, you did not specify a phenotype and/or genotype")
    }
}

check_covariates <- function(covariates, y, method) {
    if (!is.null(covariates)) {
        if (!class(covariates) == "matrix") {
            stop("Argument covariates is not a matrix")
        }
    }

    if (nrow(covariates) != length(y)) {
        msg <- paste0("The number of rows of the covariate matrix (",
                      nrow(covariates),
                      ") is not equal to the number of elements in the",
                      " phenotype (",
                      length(y),
                      ")")
        stop(msg)
    }
  if ((!is.null(covariates)) && (method=="retros")) {
    message("Attention: The Retrospective likelihood, method=", method,", does not handle covariates, 
            they will be ignored in the analysis. We recommend to use the copula-based 
            prospective method in presence of covariates")  
    covariates =NULL
  }
  return(covariates)
}


check_prev <- function(prev,Design,method) {
  if (is.null(prev) && ((Design=="CC") | (Design=="MT"))) {
    if(method=="pros") {
  stop("The primary phenotype prevalance should be specified for the copula-based prospective approach, 
       method =", method, "for Case-Control or the Maultiple-Trait designs")
    }  
  }  
}


check_cutoffs <- function(cutoffs,Design) {
    if (!is.null(cutoffs) && (Design=="CC")) {
            message("Case-control desgin does not need cutoffs")
      cutoffs = NULL
    }

  if (is.null(cutoffs) && ((Design=="ET") | (Design=="MT"))) {
      stop("The cutoffs parameter shoul be specified for the data analysis 
      if the Extreme-Trait (ET) or Multiple-Trait (MT) designs are considered.")
    }
  
  if (!is.null(cutoffs) && (Design=="ET")) {
    if(!is.vector(cutoffs) | length(cutoffs) != 2) {
      stop("The 'cutoffs' parameter should be a numeric vector of legnth 2, c(ylb, yub):
           ylb and yub are the lower and the upper primary trait thresholds")
    }
  }
  
  if (!is.null(cutoffs) && (Design=="MT")) {
    if(length(cutoffs) != 1) {
      stop("The 'cutoffs' parameter should be a numeric scalar, yub:
           yub is the upper secondary trait threshold")
    }
  }
}


