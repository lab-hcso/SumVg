require(locfdr)
require(sfsmisc)

# zall = vector of z-statistics 
# method = one of "jack", "paraboot", "fdrboot1", "fdrboot2"; please refer to the paper for further details 
# d = number of observations to delete in delete-d-jackknife 
# repl= no. of bootstrap or  jackknife iterations
# out="unconditional" or "conditional" for type of SumVg estimator
# SE =  standard error of the logistic regression coefficients (a vector)
# K = prevalence of disease
# MAF = list of minor allele frequencies for the corresponding SNPs (does not need to match to the effect allele)


SumVg.binary <- function (zall, method="paraboot", d=1, repl=200, out="unconditional", SE, K, MAF) 
                          
{

#  adapted from https://github.com/alexgillett/scale_transformation/blob/master/scale_transformation_code_HUMAN_HEREDITY_PAPER.R 
#   please refer to https://doi.org/10.1101/385740 for details
### Scale transformation functions
### Equation 3 in paper is equivalent to equation 4 for a binary RF
### Eqns 4, 5 and 6 function:
logit2liab.f <- function(OR, K, prob){  #prob is MAF 
  
  functionb0.e4 <- function(p){
    r <- rep(NA, length(p))
    r[1] <- ((1-prob)*(1/(1 + exp(-p[1])))) + (prob*(1/(1 + exp(-p[1])*exp(-log(OR))))) - K
    r
  }
  root.out <- uniroot(functionb0.e4, c(-100,100))
  b0.e4 <- root.out$root
  pDlogit.RF1.e4 <- 1/(1 + (exp(-b0.e4)*exp(-log(OR))))
  pDlogit.RF0.e4 <- 1/(1 + (exp(-b0.e4)))
  b0.e6 <- log(K/(1-K))
  pDlogit.RF1.e6 <- 1/(1 + (exp(-b0.e6)*exp(-log(OR))))
  pDlogit.RF0.e6 <- 1/(1 + (exp(-b0.e6)))
  varRF <- 2*prob*(1-prob)  #amended by me; original is    prob*(1-prob)
  a <- qnorm(pDlogit.RF1.e4) - qnorm(pDlogit.RF0.e4)
  tau.e4 <- a/sqrt(1 + (a^2)*varRF)
  tau.e5 <- a
  
  tau.e6 <- qnorm(plogis(log(K/(1-K)) + log(OR))) - qnorm(K)
  
  
  out <- c(b0.e4, b0.e6, pDlogit.RF0.e4, pDlogit.RF1.e4, pDlogit.RF0.e6, pDlogit.RF1.e6, tau.e4, tau.e5, tau.e6)
  out
}

 
ztoVg.func <-function(z, SE, K, MAF) {   

   beta = z*SE 
  OR = exp(beta)
  res = mapply(logit2liab.f,  OR,  K, prob=MAF)
  tau = res[7,]  #we use equation 4 result from https://doi.org/10.1101/385740 as default
  tau_standard = sqrt(2*MAF*(1-MAF)) * tau    # SD of genotype X * tau to produce standardized coefficient
  Vg = sum(tau_standard^2) 
  return(Vg)
}



#******************************************************************************
#  Modified Selection bias correction by Efron's empirical Bayes method                *
#******************************************************************************
bias.corr.kernel<- function (zz,xout,bw="nrd0")

{
     density.obj=density(zz,bw=bw)
     fz.func = splinefun(density.obj$x, density.obj$y)
    Psi.z = log(  fz.func(zz) /dnorm(zz)  )
    truez = D1ss(x= zz, y=Psi.z ,xout=xout)
return(truez)
}



#*************************************************************************

#source("D:\\MRes\\Alz gene\\exponential distribution of effect sizes\\Vg to Z function.R")
#source("D:\\MRes\\Alz gene\\exponential distribution of effect sizes\\Z-value to Vg.R") ##uses lnOR/varlnOR to cal z
zall.corr.ker =  bias.corr.kernel(zall,xout=zall)
#________________________________________________________________________________________________________________
SumVg.cont <- function (zall, SE, K, MAF) 
{
zall.corr.ker =  bias.corr.kernel(zall, xout=zall)
Vg.corr.ker = mapply( ztoVg.func, z=zall.corr.ker, SE=SE, K=K, MAF=MAF)
sum.kernel = sum(Vg.corr.ker)
return ( sum.kernel)
}

#****************************************************************************************************************************************************
SumVg.cont.givenH1 <- function (zall, SE, K, MAF) 

{   
zall.corr.ker =  bias.corr.kernel(zall,xout=zall)
fdr = locfdr(zall,bre=500,df=10,nulltype=0,plot=0)$fdr
largefdr.ind = which (fdr>=0.95) 
truez.given.H1 = zall.corr.ker[-largefdr.ind]/(1-fdr[-largefdr.ind])   ##note that we don't consider any z-values with fdr>0.95   
Vg.corr = mapply( ztoVg.func, z=truez.given.H1, SE=SE, K=K, MAF=MAF )*(1-fdr[-largefdr.ind]) 
sum.kernel.givenH1 = sum(Vg.corr)

return ( sum.kernel.givenH1)
}



#__________________________________________________________________________________________
if (out=="unconditional") 
{
Est.SumVg = SumVg.cont(zall, SE=SE, K=K, MAF=MAF)  
#******************************************************
#delete-d-jackknife                                   *
#******************************************************
if (method=="jack") 
{

totalNoofZ= length(zall)
u =NULL

jack.del.d <- function(x, theta, d, repl, SE, K, MAF) {
  n <- length(x)

    for (i in 1:repl) {
        indices = sample (1:totalNoofZ, size= d)
        u[i] <- theta(x[-indices], SE=SE[-indices], K=K, MAF=MAF[-indices])
        }
    thetahat <- theta(x, SE=SE, K=K, MAF=MAF)
    jack.bias <- (n - d) * (mean(u) - thetahat)    ##need to look up
    jack.se <- sqrt(        ((n - d)/d/repl) * sum((u - mean(u))^2)            )   ##formula from Shao 1988 consistency of jackknife variance estimators
    return(list(jack.se = jack.se, jack.bias = jack.bias, jack.values = u
        ))
   }
   
   t1=proc.time()
   jackobj =   jack.del.d(x=zall, SumVg.cont, d=d, repl=repl, SE=SE, K=K, MAF=MAF)
   SE.SumVg = jackobj$jack.se
   resamp_val = jackobj$jack.values
}
   
   #
#***********************
#****************************
#  "parametric" bootstrap ver 2 
# each sample New Z statistic = N(corrected z stastic,1)  '
#*********************************************************** 
if (method=="paraboot") 
{

SumVg.cont.paraboot= NULL
SumVg.cont.givenH1.paraboot = NULL

n <- length(zall)
boot.para<- function (listOfZ, repl) {
   
for (i in 1:repl) {
        x = rnorm (n, mean= listOfZ, sd=1)
        SumVg.cont.paraboot[i]  = SumVg.cont (x, SE=SE, K=K, MAF=MAF) 
                  }
        
        return ( list(boot.para.se=sd(SumVg.cont.paraboot),
        bootvalues1=SumVg.cont.paraboot   )  )   
                                      }
                                      
boot.para.obj.corrZ = boot.para(listOfZ= zall.corr.ker , repl=repl)
SE.SumVg = boot.para.obj.corrZ$boot.para.se
resamp_val = boot.para.obj.corrZ$bootvalues1
}


#___________________________________
 #****************************
#  Amended weighted "parametric" bootstrap 
# each sample New Z statistic = N(observed Z stastic,1) with probability Pr(H1) 
#***********************************************************
if (method=="fdrboot1") 
{

SumVg.cont.paraboot= NULL
SumVg.cont.givenH1.paraboot = NULL

fdrobj = locfdr(zall,bre=120,df=10,nulltype=2,plot=0)
fdr = fdrobj$fdr

n <- length(zall)
z.sim = length(n) 

repl=100
boot.para<- function (listOfZ, repl) {
   
for (i in 1:repl) {
        z.sim = length(n) 
              for (j in 1:n)  {
                      dummy = runif(1, min=0,max=1); 
                      if (dummy>fdr[j]) z.sim[j]= rnorm(1,mean=zall[j],sd=1 )  
                      else  ( z.sim[j]=rnorm(1, mean=0,sd=1)  )  
                              }
                              
        
        SumVg.cont.paraboot[i]  = SumVg.cont (z.sim, SE=SE, K=K, MAF=MAF)  } 
      
      return ( list(boot.para.se=sd(SumVg.cont.paraboot), bootvalues1=SumVg.cont.paraboot) )       }
        

boot.para.obj.wt1 = boot.para(listOfZ=zall, repl=repl)
SE.SumVg= boot.para.obj.wt1$boot.para.se
resamp_val = boot.para.obj.wt1$bootvalues1 

}

 #****************************
#  Amended weighted "parametric" bootstrap method 2 
# each sample New Z statistic = N(corrected Z stastic,1) with probability Pr(H1)'
#***********************************************************
if (method=="fdrboot2")
{ 
SumVg.cont.paraboot= NULL
SumVg.cont.givenH1.paraboot = NULL

fdrobj = locfdr(zall,bre=120,df=10,nulltype=2,plot=0)
fdr = fdrobj$fdr

zall.corr.ker =  bias.corr.kernel(zall,xout=zall)


n <- length(zall)
z.sim = length(n) 

repl=100
boot.para<- function (listOfZ, repl) {
   
for (i in 1:repl) {
        z.sim = length(n) 
              for (j in 1:n)  {
                      dummy = runif(1, min=0,max=1); 
                      if (dummy>fdr[j]) z.sim[j]= rnorm(1,mean= zall.corr.ker[j],sd=1 )  
                      else  ( z.sim[j]=rnorm(1, mean=0,sd=1)  )  
                              }
                              
        
        SumVg.cont.paraboot[i]  = SumVg.cont (z.sim, SE=SE, K=K, MAF=MAF)  } 
      
      return ( list(boot.para.se=sd(SumVg.cont.paraboot), bootvalues1=SumVg.cont.paraboot) )       }
        

boot.para.obj.wt2 = boot.para(listOfZ=zall, repl=repl)
SE.SumVg= boot.para.obj.wt2$boot.para.se
resamp_val = boot.para.obj.wt2$bootvalues1

}

} #end of the part '  if (out="unconditional") '
#********************************************************

























if (out=="conditional")
{
Est.SumVg = SumVg.cont.givenH1(zall=zall, SE=SE, K=K, MAF=MAF)  
#******************************************************
#delete-d-jackknife                                   *
#******************************************************
if (method=="jack") 
{

totalNoofZ= length(zall)
u =NULL

jack.del.d <- function(x, theta, d, repl, ...) {
  n <- length(x)

    for (i in 1:repl) {
        indices = sample (1:totalNoofZ, size= d)
        u[i] <- theta(x[-indices], ...)
        }
    thetahat <- theta(x, ...)
    jack.bias <- (n - d) * (mean(u) - thetahat)    ##need to look up
    jack.se <- sqrt(        ((n - d)/d/repl) * sum((u - mean(u))^2)            )   ##formula from Shao 1988 consistency of jackknife variance estimators
    return(list(jack.se = jack.se, jack.bias = jack.bias, jack.values = u
        ))
   }
   
   t1=proc.time()
   jackobj =   jack.del.d(x=zall, SumVg.cont.givenH1, d=d, repl=repl, SE=SE, K=K, MAF=MAF)
   SE.SumVg = jackobj$jack.se
   resamp_val = jackobj$jack.values
}
   
   #
#***********************
#****************************
#  "parametric" bootstrap ver 2 
# each sample New Z statistic = N(corrected z stastic,1)  '
#*********************************************************** 
if (method=="paraboot") 
{

SumVg.cont.paraboot= NULL
SumVg.cont.givenH1.paraboot = NULL
zall.corr.ker =  bias.corr.kernel(zall,xout=zall)
n <- length(zall)
boot.para<- function (listOfZ, repl) {
   
for (i in 1:repl) {
        x = rnorm (n, mean= listOfZ, sd=1)
        SumVg.cont.paraboot[i]  = SumVg.cont.givenH1 (x, SE=SE, K=K, MAF=MAF) 
                  }
        
        return ( list(boot.para.se=sd(SumVg.cont.paraboot),
        bootvalues1=SumVg.cont.paraboot   )  )   
                                      }
boot.para.obj.corrZ = boot.para(listOfZ= zall.corr.ker , repl=repl)
SE.SumVg = boot.para.obj.corrZ$boot.para.se
resamp_val = boot.para.obj.corrZ$bootvalues1
}


#___________________________________
 #****************************
#  Amended weighted "parametric" bootstrap 
# each sample New Z statistic = N(observed Z stastic,1) with probability Pr(H1) 
#***********************************************************
if (method=="fdrboot1") 
{

SumVg.cont.paraboot= NULL
SumVg.cont.givenH1.paraboot = NULL

fdrobj = locfdr(zall,bre=120,df=10,nulltype=2,plot=0)
fdr = fdrobj$fdr

n <- length(zall)
z.sim = length(n) 

repl=100
boot.para<- function (listOfZ, repl) {
   
for (i in 1:repl) {
        z.sim = length(n) 
              for (j in 1:n)  {
                      dummy = runif(1, min=0,max=1); 
                      if (dummy>fdr[j]) z.sim[j]= rnorm(1,mean=zall[j],sd=1 )  
                      else  ( z.sim[j]=rnorm(1, mean=0,sd=1)  )  
                              }
                              
        
        SumVg.cont.paraboot[i]  = SumVg.cont.givenH1 (z.sim, SE=SE, K=K, MAF=MAF)  } 
      
      return ( list(boot.para.se=sd(SumVg.cont.paraboot), bootvalues1=SumVg.cont.paraboot) )       }
        

boot.para.obj.wt1 = boot.para(listOfZ=zall, repl=repl)
SE.SumVg= boot.para.obj.wt1$boot.para.se
resamp_val = boot.para.obj.wt1$bootvalues1

}

 #****************************
#  Amended weighted "parametric" bootstrap method 2 
# each sample New Z statistic = N(corrected Z stastic,1) with probability Pr(H1)'
#***********************************************************
if (method=="fdrboot2")
{ 
SumVg.cont.paraboot= NULL
SumVg.cont.givenH1.paraboot = NULL

fdrobj = locfdr(zall,bre=120,df=10,nulltype=2,plot=0)
fdr = fdrobj$fdr

zall.corr.ker =  bias.corr.kernel(zall,xout=zall)


n <- length(zall)
z.sim = length(n) 

repl=100
boot.para<- function (listOfZ, repl) {
   
for (i in 1:repl) {
        z.sim = length(n) 
              for (j in 1:n)  {
                      dummy = runif(1, min=0,max=1); 
                      if (dummy>fdr[j]) z.sim[j]= rnorm(1,mean= zall.corr.ker[j],sd=1 )  
                      else  ( z.sim[j]=rnorm(1, mean=0,sd=1)  )  
                              }
                              
        
        SumVg.cont.paraboot[i]  = SumVg.cont.givenH1 (z.sim, SE=SE, K=K, MAF=MAF)  } 
      
      return ( list(boot.para.se=sd(SumVg.cont.paraboot), bootvalues1=SumVg.cont.paraboot) )       }
        

boot.para.obj.wt2 = boot.para(listOfZ=zall, repl=repl)
SE.SumVg= boot.para.obj.wt2$boot.para.se
resamp_val = boot.para.obj.wt2$bootvalues1
}

} #end of the part '  if (out="conditional") '
#************************************************************************************************************************************************


return( list(Est.SumVg= Est.SumVg, SE.SumVg= SE.SumVg) ) 


} ##end of SumVg.binary function