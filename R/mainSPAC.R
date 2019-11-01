##' Main function of the SPAC package
##'
##' The SPAC function is the main function of the SPAC
##' used package.
##'
##' @title SPAC main function
##' @param y1 vector designates the primary phenotype (one row-entry per individual), of
##'     length \eqn{n}.
##' @param y2 vector designates the secondary phenotype (one row-entry per individual), of
##'     length \eqn{n}.
##' @param G vector deignates the SNP of interest for 
##'     which the association test specified by the \code{method}
##'     parameter should be run. The vector has length \eqn{n}.
##' @param covariates matrix of covariates NOT-including intercept (dimension:
##'     \eqn{n \times p}, with \eqn{p} the number of covariates). 
##' @param link character, specifies the link function to be used for modelling the marginal distribution of the binary primary phenotype for the CC and MT designs. 
##' The available link functions are link=c("probit","logit","cloglog"). The defaut is link = "probit", which th liabiltiy latent model.
##' @param method character, selects the method to use for correcting the sampling mechanism bias. Can be one of the following:
##' \itemize{
##' \item \code{"pros"} (default), copula-based prospective method
##' \item \code{"retros"}, copula-based retrospective method. 
##' Since the copula-based restrspective method does not take into account for the covariates, if the user specifies method="retros", 
##' in presence of cavariates, the latter will simply be ignored from the analysis  
##' }
##' @param Design character, specifies the sampling design of the data. 
##' Can be one of the following:
##' \itemize{
##' \item \code{"CC"} (default), for the Case-Control sampling mechanism
##' \item \code{"ET"}, for the Extreme-Trait sampling mechanism
##' \item \code{"MT"}, for the Multiple-Trait sampling mechanism 
##' }
##' @param prev scalar between 0 and 1, specifies the primary phenotype prevalance. 
##' It is needed for the \code{method} = "prosp" and \code{Design} equals "CC" or "MT".  
##' Default is \code{prev} = NULL. If it is not specified, it will be estimated form the data.
##' @param cutoffs vector or scalar, depending on the sampling mechanism design \code{Design}.
##'  Can be one of the following:
##' \itemize{
##' \item \code{c(ylb,yub)}, a vector of length \eqn{2} for the extreme-trait sampling design, if \code{Design = "ET"}. 
##' ylb, yub:  the lower (ylb) and the upper (yub) primary trait thresholds
##' \item \code{yub}, a scalar for the muliple-trait sampling design, if \code{Design = "MT"}.
##' yub: the upper (yub) secondary trait threshold.
##' }
##' @param copfit character, selects the copula to use for modelling priamry-secondary phenptypes dependence. 
##' Can be one of the following:
##' \itemize{
##' \item \code{"Gaussian"} (default)
##' \item \code{"Student"}, Student copula with degree of freedom equals 10
##' \item \code{"Clayton"}, Clayton copula 
##' \item \code{"Gumbel"}, Gumbel copula 
##' \item \code{"Frank"}, Frank copula
##' }
##' @return A list containing results of the association test
##'     specified by the \code{method} parameter using the copula model for modelling primary-secondary dependency 
##'     by the \code{copfit} parameter for the sampling mechanism \code{Design}. The
##'     output list contains the following results:
##'     \itemize{
##'     \item \code{intercept.SNP.SecP}: the intercept estimate from the marginal model of the secondary phenotype, and its standard error
##'     \item \code{SNP.SecP}: the effect size estimate of the SNP on the secondary phenotype, and its standard error
##'     \item \code{P.value.SecP}: the p-value of the SNP/secondary phenotype association test
##'     \item \code{intercept.SNP.PrP}: the intercept estimate from the marginal mode of the primamry phenotype, and its standard error
##'     \item \code{SNP.PrP}: the effect size estimate of the SNP on the primamry phenotype, and its standard error
##'     \item \code{P.value.PrP}: the p-value of the SNP/primary phenotype association test
##'     \item \code{alpha}: the estimated primary-secondary dependence parameter, as a parameter of the used copula
##'     \item \code{tau}: Kendall's tau, which is a function of the estimated primary-secondary dependence parameter alpha
##'     \item \code{df2}: the estimated degree of freedom of the secondary-phenotype marginal distribution, Student-t
##'     \item \code{AIC}: the AIC criterion calculated based on the prospective or rterospective copula-model 
##'     }
##' @author Karim Oualkacha
##' @export
SPAC <- function(y1=NULL, y2=NULL,
                   G=NULL,
                   covariates = NULL,
                   prev = NULL,
                   prev2 = NULL,
                   p.lower = NULL,
                   p.upper = NULL,
                   cutoffs = NULL,
                   link = "probit",
                   copfit = "Gaussian",
                   method = "pros",
                   Design = "CC"
)
{
  ## Parameter checks
  check_pheno(y1) 
  check_pheno(y2)
  check_pheno(G)
  covariates <- check_covariates(covariates, y1, method)
  check_copula(copfit)
  check_link(link,Design)
  check_method(method)
  check_prev(prev,Design,method)
  check_p.lower_uppers(p.lower,p.upper,Design,method)
  check_cutoffs(cutoffs,Design)
  
  #----- calculation of p-values
  message("Starting association analysis of the SNP...")
  
    switch(Design,
         MT={
           outfit = gwas_cop_mt_snps(y1=y1, y2=y2, marker=G, covar=covariates, y2ub=cutoffs, link=link, cop=copfit,
                    lik=method, prev = prev, prev2 =prev2) # I have to add prev to this function
           
           results <- list(intercept.SNP.SecP = c(outfit$b02, outfit$se02),
                           SNP.SecP = c(outfit$b12,outfit$se2), 
                           P.value.SecP = outfit$Pvalue2,
                           intercept.SNP.PrP = c(outfit$b01, outfit$se01), 
                           SNP.PrP = c(outfit$b11,outfit$se1),
                           P.value.PrP = outfit$Pvalue1,
                           alpha = outfit$dep[1],
                           tau = outfit$dep[2],
                           df2 = outfit$df2,
                           AIC = outfit$AIC)
           
           results
         },
         ET={
           outfit = gwas_cop_et_snps(y1=y1, y2=y2, marker=G, covar=covariates, yub=cutoffs[2], ylb=cutoffs[1], cop=copfit, 
                      lik=method, p.lower =p.lower, p.upper =p.upper)
           
           results <- list(intercept.SNP.SecP = c(outfit$b02, outfit$se02),
                           SNP.SecP = c(outfit$b12,outfit$se2), 
                           P.value.SecP = outfit$Pvalue2,
                           intercept.SNP.PrP = c(outfit$b01, outfit$se01), 
                           SNP.PrP = c(outfit$b11,outfit$se1),
                           P.value.PrP = outfit$Pvalue1,
                           alpha = outfit$dep[1],
                           tau = outfit$dep[2],
                           df2 = outfit$df2,
                           AIC = outfit$AIC)
           
           results
         },
         CC={
           outfit = gwas_cop_cc_snps(y1=y1, y2=y2, marker=G, covar=covariates, link=link, cop=copfit, 
                    lik=method, prev = prev)

                      results <- list(intercept.SNP.SecP = c(outfit$b02, outfit$se02),
                           SNP.SecP = c(outfit$b12,outfit$se2), 
                           P.value.SecP = outfit$Pvalue2,
                           intercept.SNP.PrP = c(outfit$b01, outfit$se01), 
                           SNP.PrP = c(outfit$b11,outfit$se1),
                           P.value.PrP = outfit$Pvalue1,
                           alpha = outfit$dep[1],
                           tau = outfit$dep[2],
                           df2 = outfit$df2,
                           AIC = outfit$AIC)
           
           results
           }
  )
  
  return(results)
}




