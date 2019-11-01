
## Update: 2019 04 07


#############################################################################
### Note : One can use BiCopHfunc function in the vineCopula R package##
#############################################################################
#rm(list=ls())
#library(VineCopula)

## The cumulative distribution function
copCDF<-function(u1,u2,alpha,cop="Gaussian",par2=0){

  if(cop=="Gaussian") cdf<-BiCopCDF(u1=u1, u2=u2, family=1, par=alpha)
  else if(cop=="Student")  cdf<-BiCopCDF(u1=u1, u2=u2, family=2, par=alpha, par2=round(par2,0))
  else if(cop=="Clayton")  cdf<-BiCopCDF(u1=u1, u2=u2, family=3, par=alpha)
  else if(cop=="Gumbel")   cdf<-BiCopCDF(u1=u1, u2=u2, family=4, par=alpha+1)
  else if(cop=="Frank")    cdf<-BiCopCDF(u1=u1, u2=u2, family=5, par=alpha)
  else if(cop=="Joe")      cdf<-BiCopCDF(u1=u1, u2=u2, family=6, par=alpha+1)
  else stop("Unrecognized copula family")

  return(cdf)
}
#copCDF(u1=0.1,u2=0.4,alpha=0.4,cop="Gaussia",par2=0)
## The PDF function
copPDF<-function(u1,u2,alpha,cop="Gaussian",par2=0){

  if(cop=="Gaussian") pdf<-BiCopPDF(u1=u1, u2=u2, family=1, par=alpha)
  else if(cop=="Student") pdf<-BiCopPDF(u1=u1, u2=u2, family=2, par=alpha, par2=round(par2,0))
  else if(cop=="Clayton") pdf<-BiCopPDF(u1=u1, u2=u2, family=3, par=alpha)
  else if(cop=="Gumbel")  pdf<-BiCopPDF(u1=u1, u2=u2, family=4, par=alpha+1)
  else if(cop=="Frank")   pdf<-BiCopPDF(u1=u1, u2=u2, family=5, par=alpha)
  else if(cop=="Joe")     pdf<-BiCopPDF(u1=u1, u2=u2, family=6, par=alpha+1)
  else stop("Unrecognized copula family")

  return(pdf)
}

#conditional distribution of C(u,v) given U=u   deriv C(u,v)/deriv v
Hfunc<-function(u1,u2,alpha,cop="Gaussian",par2=0){

  if(cop=="Gaussian")     out<-BiCopHfunc(u1=u1, u2=u2, family=1, par=alpha)$hfunc2
  else if(cop=="Student") out<-BiCopHfunc(u1=u1, u2=u2, family=2, par=alpha, par2=round(par2,0))$hfunc2
  else if(cop=="Clayton") out<-BiCopHfunc(u1=u1, u2=u2, family=3, par=alpha)$hfunc2
  else if(cop=="Gumbel")  out<-BiCopHfunc(u1=u1, u2=u2, family=4, par=alpha+1)$hfunc2
  else if(cop=="Frank")   out<-BiCopHfunc(u1=u1, u2=u2, family=5, par=alpha)$hfunc2
  else if(cop=="Joe")     out<-BiCopHfunc(u1=u1, u2=u2, family=6, par=alpha+1)$hfunc2
  else stop("Unrecognized copula family")

  return(out)
}
#Hfunc(u1=0.1,u2=0.5,alpha=0.1,cop="Student",par2=4)

# for simulation
simCop<-function(N,alpha,cop,par2=0){

  if(cop=="Gaussian"){
    p.cop = normalCopula(param=alpha, dim=2)
    u.alpha = rCopula(N,p.cop)
  }
  else if(cop=="Student"){
    p.cop = tCopula(param=alpha, dim=2, df=par2)
    u.alpha = rCopula(N,p.cop)
  }
  else if(cop=="Clayton"){
    p.cop =archmCopula(family="clayton", dim=2, param=alpha)
    u.alpha = rCopula(N,p.cop)
  }
  else if(cop=="Frank"){
    p.cop = archmCopula(family="frank", dim=2, param=alpha)
    u.alpha = rCopula(N,p.cop)
  }
  else if(cop=="Gumbel"){
    p.cop = archmCopula(family="gumbel", dim=2, param=alpha+1)
    u.alpha = rCopula(N,p.cop)
  }
  else if(cop=="Joe"){
    p.cop = archmCopula(family="joe", dim=2, param=alpha+1)
    u.alpha = rCopula(N,p.cop)
  }
  u.alpha
}



#######################################################################################################################
####  copPar2KT: given dependence parameter of some copula, this function evaluate the Kendall's tau in term of alpha
####  using functions of Vinecopula package
########################################################################################################################

copPar2KT<-function(alpha,cop="Gaussian"){

  if(cop=="Gaussian")     KT<-BiCopPar2Tau(1,alpha)
  else if(cop=="Student") KT<-BiCopPar2Tau(2,alpha)
  else if(cop=="Clayton") KT<-copClayton@tau(alpha)
  else if(cop=="Gumbel")  KT<-copGumbel@tau(alpha+1)
  else if(cop=="Frank")   KT<-copFrank@tau(alpha)
  else if(cop=="Joe")     KT<-copJoe@tau(alpha+1)
  else stop("Unrecognized copula family")

  return(KT)
}

#######################################################################################################################
####  given kendall's tau of some copula, this function return the correponding dependence parameter
####  using functions of Vinecopula package
########################################################################################################################
ftau<-function(tau,alpha,cop){tau-copPar2KT(alpha=alpha,cop=cop)}

copKT2Par<-function(tau,cop="Gaussian"){

  if(cop=="Gaussian")     par<-BiCopTau2Par(1,tau)
  else if(cop=="Student") par<-BiCopTau2Par(2,tau)
  else if(cop=="Clayton") par<-copClayton@iTau(tau)
  else if(cop=="Gumbel")  par<-copGumbel@iTau(2-1/tau)
  else if(cop=="Frank") par<-copFrank@iTau(tau)
  else if(cop=="Joe") par<-uniroot(ftau,lower=1e-9,upper=700,extendInt="yes",tau=tau,cop="Joe")$root
  else stop("Unrecognized copula family")
  par
}


#----------------
# last update 19/03/2019 : change compared to the version "SPA_copST_et_samples-KO-5.R" is
# output of the gwas_cop_et_snps() is now a list instead of a vector
#----------------
# last update 11/02/2019 : change compared to the version "SPA_copST_et_samples-KO-4.R" is
# 1) the output which now provides beta01 and beat02; reporting bias for b02 was required by Reviewer 1
# 2) the covariates can be now adjusted for
# 3) the prospectrive likelihood is now fitted with fixed prevalance at ist true value
#----------------
#----------------
# KO: no changes are made for ET design codes of FT, thus this file remains the same as FT send it, only the name is changed to match names of CC/MT designs
#----------------
#----------------
### Copyright 20016 - 2017 Fode Tounkara
### Created: December 2016
### Last modified:
###
### Perform copula-based method for secondary phenotype association analysis using Extreme trait samples
###


### This function implement the log-likelihood for copula

cop.negloglik.et.spl<-function(pars,y1,y2,marker,covar, yub, ylb,lik="pros", cop="Gaussian", p.lu){

  ### *** inputs ***
  ###
  ###  pars : vector of model parameters
  ###  y1,y2 : vector of primary (y1) and  secondary (y2) traits values
  ###  marker : vector of genotype value  : {0,1,2,NA}
  ###  covar: matrix, or data.frame, of covariables values: one covariable per column
  ### yub, ylb: the upper (yub) and the lower (ylb) primary trait thresholds
  ### pct.lower, pct.upper: the upper (yub) and the lower (ylb) percentages for seletive samples
  ### lik : vector of three ascertainment corrections : naive, retrospective and prospective
  ###
  ### *** output ***
  ###
  ### the negative log-likelihood

  covar.mat<-as.matrix(cbind(marker,covar))
  m<-ncol(covar.mat)

  alpha=pars[1];

  log.sigma.y1 <- pars[2]; sigma.y1 <- exp(log.sigma.y1);
  log.sigma.y2 <- pars[3]; sigma.y2 <- exp(log.sigma.y2);

  beta.10<-pars[4]; beta.11.vec<-pars[(4+1):(4+m)];
  beta.20<-pars[4+m+1];  beta.21.vec<-pars[(4+m+2):(4+m+1+m)];

  df2<- pars[4+m+1+m+1];

  par2=0
  if(cop=="Student"){
    par2<- pars[4+m+1+m+2] }


  if(lik=="retros" || lik=="pros"){
    lglgc.maf <- pars[4+m+1+m+2]
    maf <- exp( lglgc.maf )/(1+exp( lglgc.maf ))
    par2=0
    if(cop=="Student")
      par2<- pars[4+m+1+m+2+1];
  }


  if ( abs(alpha)> 1 &  (cop=="Gaussian" | cop=="Student") )
    return(Inf)

  if ( ( alpha == 0 | abs(alpha)> 100 ) & (cop == "Frank") )
    return(Inf)

  if( (alpha < 0.0000001 |  alpha>99 ) &   (cop=="Clayton" | cop=="Gumbel"  | cop=="Joe") )
    return(Inf)

  if( (df2!=0) & (df2 <=2 ) )
    return(Inf)

  if( (par2!=0) & (par2 <=2 ) )
    return(Inf)


  mu1<-beta.10+covar.mat%*%beta.11.vec
  mu2<-beta.20+covar.mat%*%beta.21.vec

  sigma22<-sqrt( sigma.y2^2*(df2-2)/df2 )
  if (sigma22 <= 0) sigma22 = 1e-4

  F.y1.given.gc<-pnorm(y1,mu1,sigma.y1);
  f.y1.given.gc<-dnorm(y1,mu1,sigma.y1);

  F.y2.given.gc<-pst(y2, mu2, sigma=sigma22, nu=df2);
  f.y2.given.gc<-dst(y2, mu2, sigma=sigma22, nu=df2);

  ## naive likelihood
  joint_density<-f.y1.given.gc*f.y2.given.gc*copPDF(u1=F.y1.given.gc,u2=F.y2.given.gc,
                                                    alpha=alpha,cop=cop,par2=par2 )


  #### sampling ascertainment corrections

  Likelihood <- rep(NA, length(y1))

  if(lik=="naive"){
    Likelihood<-joint_density
  } else if(lik=="retros") {

    ## Genotype frequency
    prob.marker<-(marker==0)*(1-maf)^2  + (marker==1)*2*maf*(1-maf) + (marker==2)*maf^2

    ## estimate of the unconditional percentages of individuals with  upper and lower y1 extreme values
    g<-0:2
    Pg<-c((1-maf)^2,2*maf*(1-maf),maf^2)
    mu1g<-beta.10+g*beta.11.vec[1]
    zlbg=(ylb-mu1g)/sigma.y1
    zubg=(yub-mu1g)/sigma.y1

    prob.y1.lower <- sum( pnorm(zlbg)*Pg );
    prob.y1.upper <-1-sum( (pnorm(zubg) )*Pg );
    #prob.y1.lower <- 0.1
    #prob.y1.upper <-1-0.1

    prob.y1 <- (y1<ylb)*prob.y1.lower + (y1>yub)*prob.y1.upper

    Likelihood<-joint_density*prob.marker/prob.y1

  } else if(lik=="pros") {

    nlb<-length(y1[y1<ylb])
    nub<-length(y1[y1>yub])
    net<-nlb+nub

    ### 11/02/2019: set the prev to its true value, and omit the lines 132-->137
    ## estimate of the unconditional percentages of individuals with  upper and lower y1 extreme values
    #g<-0:2
    #Pg<-c((1-maf)^2,2*maf*(1-maf),maf^2)
    #mu1g<-beta.10+g*beta.11.vec[1]

    #prob.y1.lower <- sum( pnorm(ylb,mu1g,sigma.y1)*Pg );
    #prob.y1.upper <-1-sum( (pnorm(yub,mu1g,sigma.y1) )*Pg );
    prob.y1.lower <- p.lu[1]
    prob.y1.upper <-p.lu[2]

    # spl probabilities
    net<-nlb+nub
    prob.spl.given.y1.lower<-nlb/(prob.y1.lower)
    prob.spl.given.y1.upper<-nub/(prob.y1.upper)

    prob.spl.given.y1<-(y1<ylb)*prob.spl.given.y1.lower + (y1>yub)*prob.spl.given.y1.upper

    # conditional probabilities : percentages of individuals with  upper and lower y1 extreme values gievn marker
    prob.y1.lower.given.marker <-pnorm(ylb,mu1,sigma.y1)
    prob.y1.upper.given.marker <-1-pnorm(yub,mu1,sigma.y1)

    num<-prob.spl.given.y1*joint_density
    den<-prob.spl.given.y1.lower*prob.y1.lower.given.marker  + prob.spl.given.y1.upper*prob.y1.upper.given.marker

    Likelihood<-num/den

  } else { stop(" Unrecognized ascertainment correction") }

  return( -sum(log(Likelihood)) )
}
#cop.negloglik.et.spl(pars.init,y1,y2,marker,covar, yub, ylb, lik="retros", cop="Student")

### This function estimates model parmeters and perform the secondary phenotype test
cop.fit.et.spl<-function(y1,y2,marker,covar, yub, ylb,
                         cop="Gaussian", lik="pros", p.lu){
  ### *** inputs ***
  ###  y1,y2 : vector of primary (y1) and  secondary (y2) traits values
  ###  marker : vector of genotype value  : {0,1,2,NA}
  ###  covar: matrix, or data.frame, of covariables values: one covariable per column
  ### yub, ylb: the upper (yub) and the lower (ylb) primary trait thresholds
  ### pct.lower, pct.upper: the upper (yub) and the lower (ylb) percentages for seletive samples
  ### lik : vector of three ascertainment corrections : naive, retrospective and prospective
  ### cop : vector of five copula model : Gaussian, Student, Clayton, Gumbel, Frank, Joe
  ###
  ### *** output ***
  ### the list of the estimates of regression parameters for markers/y1, markers/y2 and covariable/y2,
  ###, the estimates of the variance of the primary and seconadary trait distribution
  ### the estimates of correlation between y1 and y2, and the estimates of the estimates of AIV values

  y1=as.numeric(y1)
  y2=as.numeric(y2)

  covar.mat<-as.matrix(cbind(marker))
  m<-ncol(covar.mat)

  ## Get marker names
  # colnames(covar.mat) <- paste("SNP",1:m,sep="")

  if(!is.null(covar)){
    covar.mat<-as.matrix(cbind(marker,covar))
    m<-ncol(covar.mat)
    ## colnames(covar.mat) <- paste("SNP",1:(m-ncol(covar)),sep="") ### KO 12/02/2019: omit this line
  }

  inds=which(  is.na(y1) | is.na(y2) )

  if(length(inds)>0) {
    warning(length(inds), " samples are missing primary or secondary or SNP information; they are ignored.\n")
    y1 <- y1[-inds]
    y2 <- y2[-inds]
    covar.mat<-covar.mat[-inds,]
  }

  # fit multivariate t

  ## initial values
  # start value for dependence parameter
  #tau.init=cor(cbind(y1,y2), method="kendall")[1,2]
  #alpha.init=copKT2Par(tau=tau.init,cop=cop)

  ## initial values
  alpha.init=cor( cbind(y1,y2) )[1,2]


  if(cop=="Student"){
    t.cop <- tCopula(dim=2)
    m1 <- pobs(as.matrix(cbind(y1,y2)))
    tau.init=cor(m1, method="kendall")[1,2]
    start.vals = c(tau.init, 3)
    names(start.vals) = c("rho","df")
    fit <- fitCopula(t.cop,m1,
                     method="mpl",
                     start=start.vals, optim.method="L-BFGS-B",
                     lower=c(-0.99, 2),
                     upper=c(0.99, 2000))
    alpha.init=coef(fit)[1]
    par2.init= coef(fit)[2]

  }

  log.sigma.y1.init=0;
  log.sigma.y2.init=0;

  beta1.init<-lm(y1~covar.mat+1)$coefficients
  beta2.init<-lm(y2~covar.mat+1)$coefficients

  nloglik=function(pars1){
    df2=pars1[1]
    beta.20<-pars1[2];  beta.21.vec<-pars1[(2+1):(2+m)];
    mu2<-beta.20+covar.mat%*%beta.21.vec
    sigma22<-sqrt( (df2-2)/df2 )
    -sum( log( dst(y2, mu2, sigma=sigma22, nu=df2) ) )
  }

  #pary2.init=try( optim(c(10,beta2.init), nloglik )$par, silent=TRUE)
  pary2.init=try( optim(c(10,beta2.init), nloglik, # KO added constrained bonderies here
                        method = "L-BFGS-B",lower=c(3, rep(-Inf, m+1)), upper = c(Inf, rep(Inf, m+1)))$par, silent=TRUE)


  if(class(pary2.init)=="try-error"){
    fit.t = fitdistr(y2, densfun="t")
    df2.init = coef(fit.t)["df"]
    pary2.init=c(coef(fit.t)["df"],beta2.init)
  }

  df2.init =pary2.init[1]
  beta2.init=pary2.init[-1]



  lglgc.maf.init=0;

  pars.init<- c(alpha.init, log.sigma.y1.init, log.sigma.y2.init, beta1.init, beta2.init,df2.init);

  if(cop=="Student"){
    pars.init<- c(alpha.init,log.sigma.y1.init,log.sigma.y2.init, beta1.init, beta2.init,df2.init, par2.init );
  }

  if(lik=="retros" || lik=="pros"){
    maf.init <- mean(covar.mat[,1])/2 # KO: added this initial value for MAF
    lglgc.maf.init <- log(maf.init/(1-maf.init))
    pars.init<- c(alpha.init,log.sigma.y1.init,log.sigma.y2.init, beta1.init, beta2.init,df2.init,lglgc.maf.init);

    if(cop=="Student"){
      pars.init<- c(alpha.init,log.sigma.y1.init,log.sigma.y2.init, beta1.init, beta2.init, df2.init,lglgc.maf.init, par2.init );
    }
  }


  outoptim<-optim(pars.init, cop.negloglik.et.spl,
                  y1=y1,y2=y2,marker=marker,covar=covar, yub=yub, ylb=ylb,
                  cop=cop, lik=lik,p.lu = p.lu,
                  method="BFGS",
                  hessian=T);

  maxloglik<--outoptim$value
  pars.est=outoptim$par

  #mvar<-solve(outoptim$hes)
  mvar<-ginv(outoptim$hes)
  SE<-sqrt(diag(mvar))
  statist=pars.est[1:(2*m+5)]/SE[1:(2*m+5)]
  p.value <- 1-pchisq(statist^2,df=1)

  alpha.est<- pars.est[1];

  tau.est=copPar2KT(alpha=alpha.est,cop=cop)

  log.sigma.y1.est <- pars.est[2]; sigma.y1.est <- exp(log.sigma.y1.est);
  log.sigma.y2.est <- pars.est[3]; sigma.y2.est <- exp(log.sigma.y2.est);

  ## estimates of beta10 and beta11
  beta1<-matrix(NA,ncol=4,nrow=m+1 )
  rownames(beta1)<-c(paste("beta1",0:m,sep=""))
  colnames(beta1)<-c("Est","SE","Stat","Pvalue")

  beta1[,1]=pars.est[(4):(4+m)]
  beta1[,2]=SE[(4):(4+m)]
  beta1[,3]=statist[(4):(4+m)]^2
  beta1[,4]=p.value[(4):(4+m)]

  ## estimates of beta20 and beta21
  beta2<-matrix(NA,ncol=4,nrow=m+1 )
  rownames(beta2)<-c(paste("beta2",0:m,sep=""))
  colnames(beta2)<-c("Est","SE","Stat","Pvalue")

  beta2[,1]=pars.est[(4+m+1):(4+m+1+m)]
  beta2[,2]=SE[(4+m+1):(4+m+1+m)]
  beta2[,3]=statist[(4+m+1):(4+m+1+m)]^2
  beta2[,4]=p.value[(4+m+1):(4+m+1+m)]

  df2.est <- pars.est[4+m+1+m+1];

  par2.est=NA
  maf.est=NA
  h2.est=NA

  if(cop=="Student"){
    par2.est <- pars.est[4+m+1+m+2];
  }

  if(lik=="retros" || lik=="pros"){
    lglgc.maf.est <- pars.est[4+m+1+m+2]
    maf.est <- exp( lglgc.maf.est )/(1+exp( lglgc.maf.est ))
    h2.est= ( 2*beta2[2,1]^2*maf.est*(1-maf.est) )/(1+2*beta2[2,1]^2*maf.est*(1-maf.est))

    if(cop=="Student")
      par2.est <- pars.est[4+m+1+m+2+1];
  }

  tau.est=copPar2KT(alpha=alpha.est,cop=cop)

  list(beta1=beta1, beta2=beta2, sigmay1=sigma.y1.est, sigmay2=sigma.y2.est, dep=c(alpha.est,tau=tau.est), df2=df2.est,
       h2=h2.est, DF=par2.est, maf=maf.est, AIC=2*length(pars.init)-2*maxloglik)
}




# Only for simulation
# this function performe a single secondary phenotype analysis for a single SNP
# gwas_cop_et_snps<-function(y1,y2,marker,covar,
#                            yub, ylb, lik, cop){
#
#     outfit<-try( cop.fit.et.spl(y1=y1,y2=y2,marker=marker,covar=covar,
#                               yub=yub, ylb=ylb,  cop=cop, lik=lik), silent = TRUE)
#     if(class(outfit)=="try-error"){
#       res<-c(alpha=NA,tau=NA, b11=NA, b21=NA, se1=NA, se2=NA,
#              Pvalue1=NA, Pvalue2=NA, mlog10Pvalue1=NA,mlog10Pvalue2=NA, AIC=NA,maf=NA,h2=NA)
#     }else{
#     res<-c(outfit$dep, b11=outfit$beta1[2,1], b21=outfit$beta2[2,1], se1=outfit$beta1[2,2], se2=outfit$beta2[2,2],
#            Pvalue1=outfit$beta1[2,4], Pvalue2=outfit$beta2[2,4], mlog10Pvalue1=-log10(outfit$beta1[2,4]),
#            mlog10Pvalue2=-log10(outfit$beta2[2,4]), AIC=outfit$AIC,maf=outfit$maf,h2=outfit$h2) }
#   names(res)=c("alpha","tau","b11","b21","se1","se2","pval1","pval2","mlog10p1","mlog10p2","AIC","maf","h2")
#   res
# }


gwas_cop_et_snps<-function(y1,y2,marker,covar,
                           yub, ylb, lik, cop, p.lu){

  outfit<-cop.fit.et.spl(y1=y1,y2=y2,marker=marker,covar=covar,
                         yub=yub, ylb=ylb,  cop=cop, lik=lik, p.lu = p.lu)

  res<-list(dep=outfit$dep, b11=outfit$beta1[2,1], b21=outfit$beta2[2,1], se1=outfit$beta1[2,2], se2=outfit$beta2[2,2],
            Pvalue1=outfit$beta1[2,4], Pvalue2=outfit$beta2[2,4], sigmay1=outfit$sigmay1, sigmay2=outfit$sigmay2, AIC=outfit$AIC, maf=outfit$maf,df2=outfit$df2,
            b01=outfit$beta1[1,1], se01=outfit$beta1[1,2], b02=outfit$beta2[1,1], se02=outfit$beta2[1,2])
  #names(res)=c("alpha","tau","b11","b21","se1","se2","pval1","pval2","sigmay1","sigmay2","AIC","maf","df2","b01","se01","b02","se02")
  res
}

## this function performe a matrix of several secondary phenotype analysis for a single SNP
## only for real data analysis

gwas_cop_et_geno <- function(pri.trait, pheno, geno, covar,
                             yub, ylb, lik, cop) {
  ### Tounkara Fode
  ### Created: 31 august 2016
  ### Last modified: december, 7,  2016
  ###
  ### For each of the secondary traits in pheno, perform MTA on one or more
  ### allelic models, using the genotype data provided in geno

  stopifnot( nrow(pheno) == ncol(geno) )

  ## Get phenotype names
  if(is.null(colnames(pheno))) { traits <- paste("trait",1:ncol(pheno),sep="") } else { traits <- colnames(pheno) }
  num.traits <- length(traits)

  ## Get marker names
  if(is.null(rownames(geno))) { markers <- paste("marker",1:nrow(geno),sep="") } else { markers <- rownames(geno) }
  num.markers <- length(markers)

  ## get sample names
  pheno.samples <- rownames(pheno)
  geno.samples <- colnames(geno)
  if(is.null(pheno.samples) | is.null(pheno.samples)) {
    if(nrow(pheno) == ncol(geno)) {
      samples <- paste("sample",1:nrow(pheno),sep="")
    } else {
      stop("Different sample numbers in pheno and geno!\n")
    }
  } else {
    if(!all(is.element(pheno.samples, geno.samples)) ) {
      stop("Different samples in pheno and geno!\n")
    } else {
      if(!all(pheno.samples == geno.samples) ) {   #samples in different order between geno & pheno
        pheno <- pheno[match(geno.samples, pheno.samples),]
      }
      samples <- geno.samples
    }
  }
  rm(geno.samples, pheno.samples)

  res <- list()

  ## Perform GWAS on each secondary trait
  for(i in 1:num.traits)
  {
    cat(traits[i],"\n")
    cur.phval <- as.vector(pheno[,i])
    inds <- which(!is.na(cur.phval))
    if (length(inds)<3) { next }  #skip trait if fewer than 2 datapoints
    cur.phval <- cur.phval[inds]


    res.mat <- matrix(NA, ncol=13, nrow=num.markers,
                      dimnames=list(markers,c("alpha", "tau","b11","b21", "se1","se2",
                                              "pval1", "pval2", "-log10pval1", "-log10pval2","AIC","maf","h2")) )

    for(j in 1:num.markers)
    {
      cur.geval <- as.numeric(geno[j,])
      if(nlevels(as.factor(cur.geval))<=1) {next}

      res.mat[j,]<-gwas_cop_et_snps(y1=pri.trait, y2=cur.phval, marker=cur.geval,
                                    covar=covar, yub=yub, ylb=ylb, lik=lik, cop=cop)
    }


    res[[traits[i]]] <- res.mat

  }
  res

}


#----------------
# last update 19/03/2019 : change compared to the version "SPA_copST_cc_samples-KO-5.R" is
# output of the gwas_cop_cc_snps() is now a list instead of a vector
#----------------
# last update 11/02/2019 : change compared to the version "SPA_copST_cc_samples-KO-4.R" is
# 1) the output which now provides beta01 and beat02; reporting bias for b02 was required by Reviewer 1
# 2) the covariates can be now adjusted for
# 3) the prospectrive likelihood is now fitted with fixed prevalance at ist true value
#----------------
#----------------
# last update and solved issues: May 1st
# L-BFGS-G method needs bounderies for alpha the cop. par. which is now fixed with lower.cop.par and upper.cop.par in optim
# the intercept in the y1 marginal model causes problems of convergence in optim(), so its domain also is constrained to the Interval (-5,5)
# initial value calculation of DF parameter for Student causes problems/bugs also, it is fixed now in optim see lines 275-276
#----------------
# update april 27th
# KO: now the prev is modeled, but constrained optimisation is used in optim "L-BFGS-G" method
# KO: In this version (lik="pros"), prevalances are fixed at their true values 0.1 at 0.06112 and all the remaining parameters are estimated

### Copyright 20016 - 2017 Fode Tounkara
### Created: December 2016
### Last modified:


## glm for the primary phenotype (case-control and multiple-trait studies)
## eta : the linear predictor eta= xiT*beta
#--------- changed by KO ----------#
meanOfY1.KO<-function(eta, link=c("probit","logit","cloglog") )
{

  if(link=="probit")
    res<-pnorm(eta)
  else if(link=="logit")
    res<- 1/(1+exp(-eta)) # KO: changed here also
  else if(link=="cloglog")
    res<-1-exp(-exp(eta))
  res
}


### This function implement the log-likelihood

cop.negloglik.cc.spl<-function(pars,y1,y2,marker,covar, link, cop, lik, prev){

  ### *** inputs ***
  ###
  ###  pars : vector of model parameters
  ###  y1,y2 : vector of primary (y1) and  secondary (y2) traits values
  ###  marker : vector of genotype value  : {0,1,2,NA}
  ###  covar: matrix, or data.frame, of covariables values: one covariable per column
  ### yub, ylb: the upper (yub) and the lower (ylb) primary trait thresholds
  ### prev : the prevelence for seletive samples
  ### lik : vector of three ascertainment corrections : naive, retrospective and prospective
  ###
  ### *** output ***
  ###
  ### the negative log-likelihood


  ## if covariable is available, we combine it with matrix of SNPs
  covar.mat<-as.matrix(cbind(marker))
  m<-ncol(covar.mat)

  if(length(covar)>0){
    covar.mat<-as.matrix(cbind(marker,covar))
    m<-ncol(covar.mat)  }


  # copula dependance parameter
  alpha=pars[1];

  if ( abs(alpha)> 1 &  (cop=="Gaussian" | cop=="Student") )
    return(Inf)

  if ( ( alpha == 0 | abs(alpha)> 100 ) & (cop == "Frank") )
    return(Inf)

  if( (alpha < 0.0000001 |  alpha>99 ) &   (cop=="Clayton" | cop=="Gumbel"  | cop=="Joe") )
    return(Inf)

  # the variance of Y2
  log.sigma.y2 <- pars[2];   sigma.y2 <- exp(log.sigma.y2);

  # regression parameters for Y1
  beta.10<-pars[3]; beta.1gc<-pars[(3+1):(3+m)];

  # regression parameters for Y2
  beta.20<-pars[3+m+1];  beta.2gc<-pars[(3+m+2):(3+m+1+m)];
  df2<- pars[3+m+1+m+1];

  # degree of freedom parameter for student copula
  par2=0
  if(cop=="Student")
    par2<- pars[3+m+1+m+2];

  if((lik=="retros")||(lik=="pros")){
    lglgc.maf <- pars[3+m+1+m+2]
    maf <- exp( lglgc.maf )/(1+exp( lglgc.maf ))
    par2=0
    if(cop=="Student")
      par2<- pars[3+m+1+m+2+1];
  }


  if( (par2!=0) & (par2 <= 2 ) )
    return(Inf)

  if(df2 <=2 ) return(Inf)

  ## GLM for Y1
  yc<-0
  eta1<-beta.10+covar.mat%*%beta.1gc
  mu1<-meanOfY1.KO(eta=as.vector((yc - eta1)),link=link)  # function difined above

  sigma22<-sqrt( sigma.y2^2*(df2-2)/df2 )
  if (sigma22 <= 0) sigma22 = 1e-4

  ## student model for Y2

  mu2<-beta.20+covar.mat%*%beta.2gc
  F.y2.given.gc<-pst(y2, mu2, sigma=sigma22, nu=df2);
  f.y2.given.gc<-dst(y2, mu2, sigma=sigma22, nu=df2);

  ## Joint probability : P(y1,y2|marker)
  #------ KO: change here also
  div.C.y2 <- Hfunc(u1=mu1, u2=F.y2.given.gc, alpha=alpha, cop=cop,par2=par2)
  joint_density<-f.y2.given.gc*( (y1==0)* div.C.y2 + (y1==1)*(1-div.C.y2) )

  ## sampling ascertainment corrections
  Likelihood <- rep(NA, length(y1))

  if(lik=="naive")
    Likelihood<-joint_density
  else if(lik=="retros") {

    # Pr(Y1=1)
    g<-0:2
    Pg<-c((1-maf)^2,2*maf*(1-maf),maf^2)
    mu0g <- meanOfY1.KO(eta=(yc-(beta.10+g*beta.1gc[1])),link=link)
    mu1g<- 1 - mu0g # KO:change here also
    proby1.eq1.given.g<-mu1g;
    pi.un<-sum( proby1.eq1.given.g*Pg )

    prob.y1 <- (y1==0)*(1-pi.un) + (y1==1)*pi.un

    ## Genotype frequency
    prob.marker<-(marker==0)*(1-maf)^2+ (marker==1)*2*maf*(1-maf) + (marker==2)*maf^2

    Likelihood <- joint_density*prob.marker/prob.y1

  } else if(lik=="pros") {

    # 11/02/2019: set the prev to its true value, and omit the lines 152-->157
    #prev <- 0.1
    #g<-0:2
    #Pg<-c((1-maf)^2,2*maf*(1-maf),maf^2)
    #mu0g <- meanOfY1.KO(eta=(yc-(beta.10+g*beta.1gc[1])),link=link)
    #mu1g<- 1 - mu0g # KO:change here also
    #proby1.eq1.given.g<-mu1g;
    #prev<-sum( proby1.eq1.given.g*Pg[g+1] )

    ncase<-length(y1[y1==1])
    ncon<-length(y1[y1==0])
    ncc<-ncase+ncon

    ## sampling probabilities (for prospective likelihood only)
    prob.spl.given.y1.eq0<-ncon/(ncc*(1-prev))
    prob.spl.given.y1.eq1<-ncase/(ncc*prev)

    prob.spl.given.y1<-(y1==0)*prob.spl.given.y1.eq0 + (y1==1)*prob.spl.given.y1.eq1

    ## conditional probabilities given marker (for prospective likelihood only)
    ## P(Y1=0 | g)
    proby1.eq0.given.marker<-mu1

    ## P(Y1=1| g)
    proby1.eq1.given.marker<-1-proby1.eq0.given.marker

    num<-prob.spl.given.y1*joint_density
    den<-prob.spl.given.y1.eq0*proby1.eq0.given.marker + prob.spl.given.y1.eq1*proby1.eq1.given.marker

    Likelihood<-num/den     # prospective likelihood
  }

  #print(-sum(log(Likelihood)))
  return( -sum(log(Likelihood)) )
}


### This function estimates model parmeters and perform the secondary phenotype test

cop.fit.cc.spl<-function(y1, y2, marker, covar, link=c("probit","logit","cloglog"),
                         cop=c("Gaussian","Student","Clayton","Gumbel","Frank","Joe"),
                         lik=c("naive","retros","pros"), prev){
  ### *** inputs ***
  ###  y1,y2 : vector of primary (y1) and  secondary (y2) traits values
  ###  marker : vector of genotype value  : {0,1,2,NA}
  ###  covar: matrix, or data.frame, of covariables values: one covariable per column
  ### yub, ylb: the upper (yub) and the lower (ylb) primary trait thresholds
  ### pct.lower, pct.upper: the upper (yub) and the lower (ylb) percentages for seletive samples
  ### lik : vector of three ascertainment corrections : naive, retrospective and prospective
  ###
  ### *** output ***
  ### the list of the estimates of regression parameters for markers/y1, markers/y2 and covariable/y2,
  ###, the estimates of the variance of the primary and seconadary trait distribution
  ### the estimates of correlation between y1 and y2, and the estimates of the estimates of AIV values

  if(length(cop)>1)
    cop="Gaussian"

  if(length(lik)>1)
    lik="pros"

  if(length(link)>1)
    link="probit"

  covar.mat<-as.matrix(cbind(marker))
  m<-ncol(covar.mat)

  ## Get marker names
  colnames(covar.mat) <- paste("SNP",1:m,sep="")

  if(!is.null(covar)){
    covar.mat<-as.matrix(cbind(marker,covar))
    m<-ncol(covar.mat)
    # colnames(covar.mat) <- paste("SNP",1:(m-ncol(covar)),sep="") ####---> KO (11 Feb 2019) moved this line since it showed bugs when covariates are present
  }

  inds=which(  is.na(y1) | is.na(y2) )

  if(length(inds)>0) {
    warning(length(inds), " samples are missing primary or secondary or SNP information; they are ignored.\n")
    y1 <- y1[-inds]
    y2 <- y2[-inds]
    covar.mat<-covar.mat[-inds,]
  }

  alpha.init=cor( cbind(y1,y2) )[1,2]
  ###### Initial values for :

  ## correltaion parameter
  ## DF parameter for Student copula
  if(cop=="Student"){
    t.cop <- tCopula(dim=2)
    m1 <- pobs(as.matrix(cbind(y1,y2)))
    tau.init=cor(m1, method="kendall")[1,2]
    start.vals = c(tau.init, 3)
    names(start.vals) = c("rho","df")
    fit <- fitCopula(t.cop,m1,
                     method="mpl",
                     start=start.vals, optim.method="L-BFGS-B",
                     lower=c(-0.985, 2.5),
                     upper=c(0.985, 100))
    alpha.init=coef(fit)[1]
    par2.init= coef(fit)[2]

  }

  ## the variance for secondary phenotype
  log.sigma.y2.init=0;

  #----------------------------------------------#
  ## regression parameters for the Y1
  glmcoef<-glm(y1~covar.mat+1,family=binomial(link=link))$coefficients #---> KO (11 Feb 2019): change code in this line to adujust for other covar
  beta1.init<-glmcoef

  #----------------------------------------------#
  ## regression parameters for the Y2
  beta2.init<-lm(y2~covar.mat+1)$coefficients

  ## df for the Y2
  nloglik=function(pars1){
    df2=pars1[1]
    beta.20<-pars1[2];  beta.21.vec<-pars1[(2+1):(2+m)];
    mu2<-beta.20+covar.mat%*%beta.21.vec
    sigma22<-sqrt( (df2-2)/df2 )
    -sum( log( dst(y2, mu2, sigma=sigma22, nu=df2) ) )
  }

  pary2.init=try( optim(c(10,beta2.init), nloglik, # KO added constrained bonderies here
                        method = "L-BFGS-B",lower=c(3, rep(-Inf, m+1)), upper = c(Inf, rep(Inf, m+1)))$par, silent=TRUE)

  if(class(pary2.init)=="try-error"){
    fit.t = fitdistr(y2, densfun="t")
    df2.init = coef(fit.t)["df"]
    pary2.init=c(coef(fit.t)["df"],beta2.init)
  }

  df2.init =pary2.init[1]
  beta2.init=pary2.init[-1]

  #--- cop param domain for L-BFGS-G method
  if (cop=="Gaussian" | cop=="Student") {
    lower.par.cop = -.98
    upper.par.cop = .98
  }

  if (cop=="Clayton" | cop=="Gumbel"  | cop=="Joe") {
    lower.par.cop = .001
    upper.par.cop = 10
  }

  if ( ( alpha.init > 0 ) & (cop == "Frank") ) {
    lower.par.cop = .001
    upper.par.cop = 10
  }

  if ( ( alpha.init < 0 ) & (cop == "Frank") ) {
    lower.par.cop = -.001
    upper.par.cop = -10
  }

  #----------------------------------------#
  if((lik=="retros") ||(lik=="pros")){
    maf.init <- mean(covar.mat[,1])/2 # KO: added this initial value for MAF
    lglgc.maf.init <- log(maf.init/(1-maf.init))

    if(cop=="Student"){
      pars.init<- c(alpha.init,log.sigma.y2.init, beta1.init, beta2.init, df2.init, lglgc.maf.init, par2.init);
      outoptim<- try(optim(pars.init, cop.negloglik.cc.spl,
                           y1=y1, y2=y2, marker=marker, covar=covar, link=link,  cop=cop, lik=lik, prev =prev,
                           method = "L-BFGS-B", lower=c(lower.par.cop,-Inf,rep(-Inf,2*(m+1)),3,-Inf,3), upper=c(upper.par.cop,Inf,rep(Inf,2*(m+1)),Inf,Inf,Inf),
                           hessian=T),TRUE);
    }
    if(cop!="Student"){
      pars.init<- c(alpha.init,log.sigma.y2.init, beta1.init, beta2.init, df2.init, lglgc.maf.init);
      outoptim<- try(optim(pars.init, cop.negloglik.cc.spl,
                           y1=y1, y2=y2, marker=marker, covar=covar, link=link, cop=cop, lik=lik, prev = prev,
                           method = "L-BFGS-B", lower=c(lower.par.cop,-Inf,-3,rep(-Inf,2*m +1),3,-Inf), upper=c(upper.par.cop,Inf,3,rep(Inf,2*m +1),Inf,Inf), # domaine of intercept of y1 has limited to (-3,3)
                           hessian=T),TRUE);
    }
  }

  #print(outoptim)
  maxloglik <- try(- outoptim$value,TRUE)
  pars.est=try(outoptim$par,TRUE)
  # if(cop=="Student"){ mvar<-ginv(outoptim$hes) solve(outoptim$hes[1:(2*m+4),1:(2*m+4)])} # August 31-2017: change lines 342-325, change solve with g
  # if(cop!="Student"){ mvar<-solve(outoptim$hes)}
  mvar<-ginv(outoptim$hes)
  SE<-sqrt(diag(mvar))[1:(2*m+4)]
  statist=pars.est[1:(2*m+4)]/SE
  p.value <- 1-pchisq(statist^2,df=1)

  alpha.est<- pars.est[1];

  tau.est=copPar2KT(alpha=alpha.est,cop=cop)

  log.sigma.y2.est <- pars.est[2]; sigma.y2.est <- exp(log.sigma.y2.est);
  y1c.est <- pars.est[3]
  ## estimates of beta10 and beta11
  beta1<-matrix(NA,ncol=4,nrow=m+1 )
  rownames(beta1)<-c(paste("beta1",0:m,sep=""))
  colnames(beta1)<-c("Est","SE","Stat","Pvalue")

  beta1[,1]=pars.est[(3):(3+m)]
  beta1[,2]=SE[(3):(3+m)]
  beta1[,3]=statist[(3):(3+m)]^2
  beta1[,4]=p.value[(3):(3+m)]

  ## estimates of beta20 and beta21
  beta2<-matrix(NA,ncol=4,nrow=m+1 )
  rownames(beta2)<-c(paste("beta2",0:m,sep=""))
  colnames(beta2)<-c("Est","SE","Stat","Pvalue")

  beta2[,1]=pars.est[(3+m+1):(3+m+1+m)]
  beta2[,2]=SE[(3+m+1):(3+m+1+m)]
  beta2[,3]=statist[(3+m+1):(3+m+1+m)]^2
  beta2[,4]=p.value[(3+m+1):(3+m+1+m)]

  df2.est <- pars.est[3+m+1+m+1];

  par2.est=NA
  maf.est=NA
  h2.est=NA

  if(cop=="Student"){
    par2.est <- pars.est[3+m+1+m+2];
  }

  if(lik=="retros" || lik=="pros"){
    lglgc.maf.est <- pars.est[3+m+1+m+2]
    maf.est <- exp( lglgc.maf.est )/(1+exp( lglgc.maf.est ))
    h2.est= ( 2*beta2[2,1]^2*maf.est*(1-maf.est) )/(1+2*beta2[2,1]^2*maf.est*(1-maf.est))

    if(cop=="Student"){
      par2.est <- pars.est[3+m+1+m+2+1];
    }
  }

  list(beta1=beta1, beta2=beta2, sigmay2=sigma.y2.est, dep=c(alpha.est,tau=tau.est), h2=0,
       DF=par2.est, maf=.3, log.l = 2*maxloglik, AIC=2*length(pars.init)-2*maxloglik, df2 = df2.est)

}


#---------------------------------------------------------------------------------#
#                                                                                 #
#---------------------------------------------------------------------------------#
gwas_cop_cc_snps<-function(y1,y2,marker,covar, link, cop, lik, prev){


  outfit<-cop.fit.cc.spl(y1=y1,y2=y2,marker=marker,covar=covar,
                         link=link,  cop=cop, lik=lik, prev = prev)

  res<-list(dep = outfit$dep, b11=outfit$beta1[2,1], b12=outfit$beta2[2,1], se1=outfit$beta1[2,2], se2=outfit$beta2[2,2],
            Pvalue1=outfit$beta1[2,4], Pvalue2=outfit$beta2[2,4], sigmay1=1, sigma2=outfit$sigmay2, AIC=outfit$AIC,maf=outfit$maf, df2=outfit$df2,
            b01=outfit$beta1[1,1], se01=outfit$beta1[1,2], b02=outfit$beta2[1,1], se02=outfit$beta2[1,2])
  #  names(res)=c("alpha","tau","b11","b21","se1","se2","pval1","pval2","sigmay1","sigmay2","AIC","maf", "df2","b01","se01","b02","se02")
  res
}

#---------------------------------------------------------------------------------#
#                                                                                 #
#---------------------------------------------------------------------------------#
gwas_cop_cc_geno <- function(pri.trait, pheno, geno, covar, link,  lik, cop) {

  stopifnot( nrow(pheno) == ncol(geno) )

  ## Get phenotype names
  if(is.null(colnames(pheno))) { traits <- paste("trait",1:ncol(pheno),sep="") } else { traits <- colnames(pheno) }
  num.traits <- length(traits)

  ## Get marker names
  if(is.null(rownames(geno))) { markers <- paste("marker",1:nrow(geno),sep="") } else { markers <- rownames(geno) }
  num.markers <- length(markers)

  ## get sample names
  pheno.samples <- rownames(pheno)
  geno.samples <- colnames(geno)
  if(is.null(pheno.samples) | is.null(pheno.samples)) {
    if(nrow(pheno) == ncol(geno)) {
      samples <- paste("sample",1:nrow(pheno),sep="")
    } else {
      stop("Different sample numbers in pheno and geno!\n")
    }
  } else {
    if(!all(is.element(pheno.samples, geno.samples)) ) {
      stop("Different samples in pheno and geno!\n")
    } else {
      if(!all(pheno.samples == geno.samples) ) {   #samples in different order between geno & pheno
        pheno <- pheno[match(geno.samples, pheno.samples),]
      }
      samples <- geno.samples
    }
  }
  rm(geno.samples, pheno.samples)

  res <- list()

  ## Perform GWAS on each secondary trait
  for(i in 1:num.traits)
  {
    cat(traits[i],"\n")
    cur.phval <- as.vector(pheno[,i])
    inds <- which(!is.na(cur.phval))
    if (length(inds)<3) { next }  #skip trait if fewer than 2 datapoints
    cur.phval <- cur.phval[inds]


    res.mat <- matrix(NA, ncol=13, nrow=num.markers,
                      dimnames=list(markers,c("alpha", "tau","b11","b21", "se1","se2",
                                              "pval1", "pval2", "-log10pval1", "-log10pval2","AIC","maf","h2")) )

    for(j in 1:num.markers)
    {
      cur.geval <- as.numeric(geno[j,])
      if(nlevels(as.factor(cur.geval))<=1) {next}

      cur.geval <- as.numeric(geno[j,])
      res.mat[j,]<-gwas_cop_cc_snps(y1=pri.trait, y2=cur.phval, marker=cur.geval,
                                    covar=covar, link=link, lik=lik, cop=cop)
    }


    res[[traits[i]]] <- res.mat

  }
  res

}


#----------------
# last update 19/03/2019 : change compared to the version "SPA_copST_mt_samples-KO-5.R" is
# output of the gwas_cop_mt_snps() is now a list instead of a vector
#----------------
# last update 11/02/2019 : change compared to the version "SPA_LCcopST_MT_samples-KO-4.R" is
# 1) the output now provides also beta01 and beat02; reporting bias for b02 was required by Reviewer 1
# 2) the covariates can be now adjusted for using prospective likelihood
#----------------
#----------------
# last update: May 1st
# KO: In this version (lik="pros"), prevalances are fixed at their true values 0.1 at 0.06112 and all the remaining parameters are estimated
# last update and solved issues: May 1st
# L-BFGS-G method needs bounderies for alpha the cop. par. which is now fixed with lower.cop.par and upper.cop.par in optim
# the intercept in the y1 marginal model causes problems of convergence in optim(), so its domain also is constrained to the Interval (-5,5)
# initial value calculation of DF parameter for Student causes problems/bugs also, it is fixed now in optim see lines 275-276
#----------------

### Copyright 20016 - 2017 Fode Tounkara
### Created: December 2016
### Last modified:
###


## glm for the primary phenotype (case-control and multiple-trait studies)
## eta : the linear predictor eta= xiT*beta
#--------- changed by KO ----------#
meanOfY1.KO<-function(eta, link=c("probit","logit","cloglog") )
{

  if(link=="probit")
    res<-pnorm(eta)
  else if(link=="logit")
    res<- 1/(1+exp(-eta)) # KO: changed here also
  else if(link=="cloglog")
    res<-1-exp(-exp(eta))
  res
}


### This function implement the log-likelihood

cop.negloglik.mt.spl<-function(pars,y1,y2,marker,covar, y2ub, link, cop, lik, prev, prev2){

  ### *** inputs ***
  ###
  ###  pars : vector of model parameters
  ###  y1,y2 : vector of primary (y1) and  secondary (y2) traits values
  ###  marker : vector of genotype value  : {0,1,2,NA}
  ###  covar: matrix, or data.frame, of covariables values: one covariable per column
  ### yub, ylb: the upper (yub) and the lower (ylb) primary trait thresholds
  ### prev : the prevelence for seletive samples
  ### lik : vector of three ascertainment corrections : naive, retrospective and prospective
  ###
  ### *** output ***
  ###
  ### the negative log-likelihood


  ## if covariable is available, we combine it with matrix of SNPs
  covar.mat<-as.matrix(cbind(marker))
  m<-ncol(covar.mat)

  if(length(covar)>0){
    covar.mat<-as.matrix(cbind(marker,covar))
    m<-ncol(covar.mat)
  }

  # copula dependance parameter
  alpha=pars[1];

  if ( abs(alpha)> 1 &  (cop=="Gaussian" | cop=="Student") )
    return(Inf)

  if ( ( alpha == 0 | abs(alpha)> 100 ) & (cop == "Frank") )
    return(Inf)

  if( (alpha < 0.0000001 |  alpha>99 ) &   (cop=="Clayton" | cop=="Gumbel"  | cop=="Joe") )
    return(Inf)

  # the variance of Y2
  log.sigma.y2 <- pars[2];   sigma.y2 <- exp(log.sigma.y2);

  # regression parameters for Y1
  beta.10<-pars[3]; beta.1gc<-pars[(3+1):(3+m)];

  # regression parameters for Y2
  beta.20<-pars[3+m+1];  beta.2gc<-pars[(3+m+2):(3+m+1+m)];
  df2<- pars[3+m+1+m+1];

  # degree of freedom parameter for student copula
  par2=0
  if(cop=="Student")
    par2<- pars[3+m+1+m+2];

  if(lik=="retros"){
    lglgc.maf <- pars[3+m+1+m+2]
    maf <- exp( lglgc.maf )/(1+exp( lglgc.maf ))
    par2=0
    if(cop=="Student")
      par2<- pars[3+m+1+m+2+1];
  }

  if( (par2!=0) & (par2 <= 2 ) )
    return(Inf)

  if(df2 <=2 ) return(Inf)

  ## GLM for Y1
  yc<-0
  eta1<-beta.10+covar.mat%*%beta.1gc
  mu1<-meanOfY1.KO(eta=as.vector((yc - eta1)),link=link)  # function difined above

  sigma22<-sqrt( sigma.y2^2*(df2-2)/df2 )
  if (sigma22 <= 0) sigma22 = 1e-4

  ## student model for Y2

  mu2<-beta.20+covar.mat%*%beta.2gc
  F.y2.given.gc<-pst(y2, mu2, sigma=sigma22, nu=df2);
  f.y2.given.gc<-dst(y2, mu2, sigma=sigma22, nu=df2);

  ## Joint probability : P(y1,y2|marker)
  #------ KO: change here also
  div.C.y2 <- Hfunc(u1=mu1, u2=F.y2.given.gc, alpha=alpha, cop=cop, par2=par2)
  joint_density<-f.y2.given.gc*( (y1==0)* div.C.y2 + (y1==1)*(1-div.C.y2) )

  ## sampling ascertainment corrections
  Likelihood <- rep(NA, length(y1))

  if(lik=="naive")
    Likelihood<-joint_density
  else if(lik=="retros") {

    # Pr(Y1=1)
    #g<-0:2
    #Pg<-c((1-maf)^2,2*maf*(1-maf),maf^2)
    #proby1.eq1.given.g<-meanOfY1(eta=beta.10+g*beta.1gc[1],link=link)
    #mu1g<- proby1.eq1.given.g
    #pi.un<-1-sum( proby1.eq1.given.g*Pg )
     
    g<-0:2
    Pg<-c((1-maf)^2,2*maf*(1-maf),maf^2)
    mu0g <- meanOfY1.KO(eta=(yc-(beta.10+g*beta.1gc[1])),link=link)
    mu1g<- 1 -  mu0g # KO:change here also
    proby1.eq1.given.g<-mu1g;
    pi.un<-sum( proby1.eq1.given.g*Pg )
    
    
    
    ## P(Y1=1, Y2>y2ub) (prev2)
    mu2g<-beta.20+g*beta.2gc[1]
    F.y2ub.given.g <- pst(y2ub, mu2g, sigma=sigma22, nu=df2);
    proby1.eq0.y2.lower.given.g<-copCDF(u1=mu0g, u2=F.y2ub.given.g, alpha=alpha,cop=cop,par2=par2)
    proby1.eq1.y2.upper.given.g<- (1-F.y2ub.given.g) - mu0g + proby1.eq0.y2.lower.given.g
    pi.af<- sum( proby1.eq1.y2.upper.given.g*Pg )
    
    prob.y1.y2 <- (y1==0)*(1-pi.un) + (y1==1)*pi.af

    ## Genotype frequency
    prob.marker<-(marker==0)*(1-maf)^2+ (marker==1)*2*maf*(1-maf) + (marker==2)*maf^2

    Likelihood<-joint_density*prob.marker/prob.y1.y2

  } else if(lik=="pros") {

    #g<-0:2
    #Pg<-c((1-maf)^2,2*maf*(1-maf),maf^2)
    #mu0g <- meanOfY1.KO(eta=(yc-(beta.10+g*beta.1gc[1])),link=link)
    #mu1g<- 1 - mu0g # KO:change here also
    #proby1.eq1.given.g<-mu1g;
    #pi.un<-sum( proby1.eq1.given.g*Pg[g+1] )
    #pi.un<-0.1

    ## P(Y1=1, Y2>y2ub) (prev2)
    #mu2g<-beta.20+g*beta.2gc[1]
    #F.y2ub.given.g<-pnorm( (y2ub-mu2g)/sqrt(sigma.y2) )
    #proby1.eq0.y2.lower.given.g<-copCDF(u1=mu0g, u2=F.y2ub.given.g, alpha=alpha,cop=cop,par2=par2)
    #proby1.eq1.y2.upper.given.g<- (1-F.y2ub.given.g) - mu0g + proby1.eq0.y2.lower.given.g
    #pi.af<- sum( proby1.eq1.y2.upper.given.g*Pg[g+1] )
    #pi.af<- 0.06112

    naff<-length(y1[y1==1])
    nunaff<-length(y1[y1==0])
    nmt<-nunaff+naff

    ## sampling probabilities (for prospective likelihood only)
    prob.spl.given.y1.eq0<-nunaff/(nmt*(1-prev))
    prob.spl.given.y1.eq1<-naff/(nmt*prev2)

    prob.spl.given.y1.y2<-(y1==0)*prob.spl.given.y1.eq0 + (y1==1)*prob.spl.given.y1.eq1

    ## conditional probabilities given marker (for prospective likelihood only)
    ## P(Y1=0 | g)
    proby1.eq0.given.marker<-mu1

    ## P(Y1=1, Y2>y2ub | g)
    F.y2ub.given.marker <- pst(y2ub, mu2, sigma=sigma22, nu=df2);
    proby1.eq0.y2.lower.given.marker <- copCDF(u1=mu1, u2=F.y2ub.given.marker, alpha=alpha,cop=cop, par2=par2)
    proby1.eq1.y2.upper.given.marker <- (1 - F.y2ub.given.marker) - mu1 + proby1.eq0.y2.lower.given.marker

    num<-prob.spl.given.y1.y2*joint_density
    den<-prob.spl.given.y1.eq0*proby1.eq0.given.marker + prob.spl.given.y1.eq1*proby1.eq1.y2.upper.given.marker

    Likelihood<-num/den     # prospective likelihood

  }  else stop("Unrecognized likelihood")

  #print(-sum(log(Likelihood)))
  return( -sum(log(Likelihood)) )
}


### This function estimates model parmeters and perform the secondary phenotype test

cop.fit.mt.spl<-function(y1, y2, marker, covar,y2ub, link=c("probit","logit","cloglog"),
                         cop=c("Gaussian","Student","Clayton","Gumbel","Frank","Joe"),
                         lik=c("naive","retros","pros"), prev, prev2){
  ### *** inputs ***
  ###  y1,y2 : vector of primary (y1) and  secondary (y2) traits values
  ###  marker : vector of genotype value  : {0,1,2,NA}
  ###  covar: matrix, or data.frame, of covariables values: one covariable per column
  ### yub, ylb: the upper (yub) and the lower (ylb) primary trait thresholds
  ### pct.lower, pct.upper: the upper (yub) and the lower (ylb) percentages for seletive samples
  ### lik : vector of three ascertainment corrections : naive, retrospective and prospective
  ###
  ### *** output ***
  ### the list of the estimates of regression parameters for markers/y1, markers/y2 and covariable/y2,
  ###, the estimates of the variance of the primary and seconadary trait distribution
  ### the estimates of correlation between y1 and y2, and the estimates of the estimates of AIV values

  if(length(cop)>1)
    cop="Gaussian"

  if(length(lik)>1)
    lik="pros"

  if(length(link)>1)
    link="probit"

  covar.mat<-as.matrix(cbind(marker))
  m<-ncol(covar.mat)

  ## Get marker names
  colnames(covar.mat) <- paste("SNP",1:m,sep="")

  if(!is.null(covar)){
    covar.mat<-as.matrix(cbind(marker,covar))
    m<-ncol(covar.mat)
    # colnames(covar.mat) <- paste("SNP",1:(m-ncol(covar)),sep="") ### KO 12/02/2019: omit this line
  }

  inds=which(  is.na(y1) | is.na(y2) )

  if(length(inds)>0) {
    warning(length(inds), " samples are missing primary or secondary or SNP information; they are ignored.\n")
    y1 <- y1[-inds]
    y2 <- y2[-inds]
    covar.mat<-covar.mat[-inds,]
  }

  alpha.init=cor( cbind(y1,y2) )[1,2]
  ###### Initial values for :

  ## correltaion parameter
  ## DF parameter for Student copula
  if(cop=="Student"){
    t.cop <- tCopula(dim=2)
    m1 <- pobs(as.matrix(cbind(y1,y2)))
    tau.init=cor(m1, method="kendall")[1,2]
    start.vals = c(tau.init, 3)
    names(start.vals) = c("rho","df")
    fit <- fitCopula(t.cop,m1,
                     method="mpl",
                     start=start.vals, optim.method="L-BFGS-B",
                     lower=c(-0.98, 2.5),
                     upper=c(0.98, 100))
    alpha.init=coef(fit)[1]
    par2.init= coef(fit)[2]

  }

  ## the variance for secondary phenotype
  log.sigma.y2.init=0;

  #----------------------------------------------#
  ## regression parameters for the Y1
  glmcoef<-glm(y1~covar.mat+1,family=binomial(link=link))$coefficients
  beta1.init<-glmcoef

  #----------------------------------------------#
  ## regression parameters for the Y2
  beta2.init<-lm(y2~covar.mat+1)$coefficients

  ## df for the Y2
  nloglik=function(pars1){
    df2=pars1[1]
    beta.20<-pars1[2];  beta.21.vec<-pars1[(2+1):(2+m)];
    mu2<-beta.20+covar.mat%*%beta.21.vec
    sigma22<-sqrt( (df2-2)/df2 )
    -sum( log( dst(y2, mu2, sigma=sigma22, nu=df2) ) )
  }

  pary2.init=try( optim(c(10, beta2.init), nloglik, # KO added constrained bonderies here
                        method = "L-BFGS-B",lower=c(3,rep(-Inf, (m+1))), upper = c(Inf, rep(Inf, (m+1))) )$par, silent=TRUE)

  if(class(pary2.init)=="try-error"){
    fit.t = fitdistr(y2, densfun="t")
    df2.init = coef(fit.t)["df"]
    pary2.init=c(coef(fit.t)["df"],beta2.init)
  }

  df2.init =pary2.init[1]
  beta2.init=pary2.init[-1]

  #--- cop param domain for L-BFGS-G method
  if (cop=="Gaussian" | cop=="Student") {
    lower.par.cop = -.98
    upper.par.cop = .98
  }

  if (cop=="Clayton" | cop=="Gumbel"  | cop=="Joe") {
    lower.par.cop = .001
    upper.par.cop = 10
  }

  if ( ( alpha.init > 0 ) & (cop == "Frank") ) {
    lower.par.cop = .001
    upper.par.cop = 10
  }

  if ( ( alpha.init < 0 ) & (cop == "Frank") ) {
    lower.par.cop = -.001
    upper.par.cop = -10
  }

  #----------------------------------------#
  if(lik=="retros"){
    maf.init <- mean(covar.mat[,1])/2 # KO: added this initial value for MAF
    lglgc.maf.init <- log(maf.init/(1-maf.init))

    if(cop=="Student"){
      pars.init<- c(alpha.init,log.sigma.y2.init, beta1.init, beta2.init, df2.init, lglgc.maf.init, par2.init);
      outoptim<- try(optim(pars.init, cop.negloglik.mt.spl,
                           y1=y1, y2=y2, marker=marker, covar=covar, y2ub=y2ub, link=link,  cop=cop, lik=lik, prev = prev, prev2 = prev2,
                           method = "L-BFGS-B", lower=c(lower.par.cop,-Inf,rep(-Inf, 2*(m+1)),3,-Inf,3), upper=c(upper.par.cop,Inf,rep(Inf, 2*(m+1)),Inf,Inf,Inf),
                           hessian=T),TRUE);
    }
    if(cop!="Student"){
      pars.init<- c(alpha.init,log.sigma.y2.init, beta1.init, beta2.init, df2.init, lglgc.maf.init);
      #print(pars.init)
      outoptim<- try(optim(pars.init, cop.negloglik.mt.spl,
                           y1=y1, y2=y2, marker=marker, covar=covar, y2ub=y2ub, link=link, cop=cop, lik=lik, prev = prev, prev2 = prev2,
                           method = "L-BFGS-B", lower=c(lower.par.cop,-Inf,-3,rep(-Inf, (2*m+1)),3,-Inf), upper=c(upper.par.cop,Inf,3,rep(Inf, (2*m+1)),Inf,Inf),
                           hessian=T),TRUE);
    }
  }

  if(lik=="pros"){
    if(cop=="Student"){
      pars.init<- c(alpha.init,log.sigma.y2.init, beta1.init, beta2.init,df2.init, par2.init);
      outoptim<- try(optim(pars.init, cop.negloglik.mt.spl,
                           y1=y1, y2=y2, marker=marker, covar=covar,y2ub=y2ub , link=link,  cop=cop, lik=lik, prev = prev, prev2 = prev2,
                           method = "L-BFGS-B", lower=c(lower.par.cop,-Inf,rep(-Inf,2*(m+1)),3,3), upper=c(upper.par.cop,Inf,rep(Inf, 2*(m+1)),Inf,Inf),
                           hessian=T),TRUE);
    }

    if(cop!="Student"){
      pars.init<- c(alpha.init,log.sigma.y2.init, beta1.init, beta2.init,df2.init);
      outoptim<- try(optim(pars.init, cop.negloglik.mt.spl,
                           y1=y1, y2=y2, marker=marker, covar=covar, y2ub=y2ub, link=link,  cop=cop, lik=lik, prev = prev, prev2 = prev2,
                           method = "L-BFGS-B", lower=c(lower.par.cop,-Inf,-3,rep(-Inf, (2*m+1)),3), upper=c(upper.par.cop,Inf,3,rep(Inf, (2*m+1)),Inf), # !!! need upper-lower for each copula
                           hessian=T),TRUE);
    }
  }

  #print("#-------------------------------------#")
  #print(outoptim)
  maxloglik <- try(- outoptim$value,TRUE)
  pars.est=try(outoptim$par,TRUE)
  if(cop=="Student"){ mvar<-solve(outoptim$hes[1:(2*m+4),1:(2*m+4)])}
  if(cop!="Student"){ mvar<-solve(outoptim$hes)}
  SE<-sqrt(diag(mvar))[1:(2*m+4)]
  statist=pars.est[1:(2*m+4)]/SE
  p.value <- 1-pchisq(statist^2,df=1)

  alpha.est<- pars.est[1];

  tau.est=copPar2KT(alpha=alpha.est,cop=cop)

  log.sigma.y2.est <- pars.est[2]; sigma.y2.est <- exp(log.sigma.y2.est);
  y1c.est <- pars.est[3]
  ## estimates of beta10 and beta11
  beta1<-matrix(NA,ncol=4,nrow=m+1 )
  rownames(beta1)<-c(paste("beta1",0:m,sep=""))
  colnames(beta1)<-c("Est","SE","Stat","Pvalue")

  beta1[,1]=pars.est[(3):(3+m)]
  beta1[,2]=SE[(3):(3+m)]
  beta1[,3]=statist[(3):(3+m)]^2
  beta1[,4]=p.value[(3):(3+m)]

  ## estimates of beta20 and beta21
  beta2<-matrix(NA,ncol=4,nrow=m+1 )
  rownames(beta2)<-c(paste("beta2",0:m,sep=""))
  colnames(beta2)<-c("Est","SE","Stat","Pvalue")

  beta2[,1]=pars.est[(3+m+1):(3+m+1+m)]
  beta2[,2]=SE[(3+m+1):(3+m+1+m)]
  beta2[,3]=statist[(3+m+1):(3+m+1+m)]^2
  beta2[,4]=p.value[(3+m+1):(3+m+1+m)]

  df2.est <- pars.est[3+m+1+m+1];

  par2.est=NA
  maf.est=NA
  h2.est=NA

  if(cop=="Student"){
    par2.est <- pars.est[3+m+1+m+2];
  }

  if(lik=="retros"){
    lglgc.maf.est <- pars.est[3+m+1+m+2]
    maf.est <- exp( lglgc.maf.est )/(1+exp( lglgc.maf.est ))
    h2.est= ( 2*beta2[2,1]^2*maf.est*(1-maf.est) )/(1+2*beta2[2,1]^2*maf.est*(1-maf.est))

    if(cop=="Student"){
      par2.est <- pars.est[3+m+1+m+2+1];
    }
  }

  list(beta1=beta1, beta2=beta2, sigmay2=sigma.y2.est, dep=c(alpha.est,tau=tau.est), h2=0,
       DF=par2.est, maf=.3, log.l = 2*maxloglik, AIC=2*length(pars.init)-2*maxloglik, df2 = df2.est)

}


#---------------------------------------------------------------------------------#
#                                                                                 #
#---------------------------------------------------------------------------------#
gwas_cop_mt_snps<-function(y1,y2,marker,covar, y2ub, link, cop, lik, prev, prev2){


  outfit<-cop.fit.mt.spl(y1=y1,y2=y2,marker=marker,covar=covar, y2ub=y2ub,
                         link=link,  cop=cop, lik=lik, prev = prev, prev2 = prev2)

  res<-list(dep = outfit$dep, b11=outfit$beta1[2,1], b21=outfit$beta2[2,1], se1=outfit$beta1[2,2], se2=outfit$beta2[2,2],
            Pvalue1=outfit$beta1[2,4], Pvalue2=outfit$beta2[2,4], sigmay1=1, sigma2=outfit$sigmay2, AIC=outfit$AIC,maf=outfit$maf, df2=outfit$df2,
            b01=outfit$beta1[1,1], se01=outfit$beta1[1,2], b02=outfit$beta2[1,1], se02=outfit$beta2[1,2])
  # names(res)=c("alpha","tau","b11","b21","se1","se2","pval1","pval2","sigmay1", "sigmay2","AIC","maf", "df2","b01","se01","b02","se02")
  res
}

#---------------------------------------------------------------------------------#
#                                                                                 #
#---------------------------------------------------------------------------------#
gwas_cop_mt_geno <- function(pri.trait, pheno, geno, covar, link,  lik, cop) {

  stopifnot( nrow(pheno) == ncol(geno) )

  ## Get phenotype names
  if(is.null(colnames(pheno))) { traits <- paste("trait",1:ncol(pheno),sep="") } else { traits <- colnames(pheno) }
  num.traits <- length(traits)

  ## Get marker names
  if(is.null(rownames(geno))) { markers <- paste("marker",1:nrow(geno),sep="") } else { markers <- rownames(geno) }
  num.markers <- length(markers)

  ## get sample names
  pheno.samples <- rownames(pheno)
  geno.samples <- colnames(geno)
  if(is.null(pheno.samples) | is.null(pheno.samples)) {
    if(nrow(pheno) == ncol(geno)) {
      samples <- paste("sample",1:nrow(pheno),sep="")
    } else {
      stop("Different sample numbers in pheno and geno!\n")
    }
  } else {
    if(!all(is.element(pheno.samples, geno.samples)) ) {
      stop("Different samples in pheno and geno!\n")
    } else {
      if(!all(pheno.samples == geno.samples) ) {   #samples in different order between geno & pheno
        pheno <- pheno[match(geno.samples, pheno.samples),]
      }
      samples <- geno.samples
    }
  }
  rm(geno.samples, pheno.samples)

  res <- list()

  ## Perform GWAS on each secondary trait
  for(i in 1:num.traits)
  {
    cat(traits[i],"\n")
    cur.phval <- as.vector(pheno[,i])
    inds <- which(!is.na(cur.phval))
    if (length(inds)<3) { next }  #skip trait if fewer than 2 datapoints
    cur.phval <- cur.phval[inds]


    res.mat <- matrix(NA, ncol=13, nrow=num.markers,
                      dimnames=list(markers,c("alpha", "tau","b11","b21", "se1","se2",
                                              "pval1", "pval2", "-log10pval1", "-log10pval2","AIC","maf","h2")) )

    for(j in 1:num.markers)
    {
      cur.geval <- as.numeric(geno[j,])
      if(nlevels(as.factor(cur.geval))<=1) {next}

      cur.geval <- as.numeric(geno[j,])
      res.mat[j,]<-gwas_cop_cc_snps(y1=pri.trait, y2=cur.phval, marker=cur.geval,
                                    covar=covar, link=link, lik=lik, cop=cop)
    }


    res[[traits[i]]] <- res.mat

  }
  res

}






