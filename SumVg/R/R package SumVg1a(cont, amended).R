require(locfdr)
require(sfsmisc)

 
SumVg <- function (zall, totalN, method="paraboot", d=1, repl=50,out="unconditional" ) 

{


#***************************************************************************
#                 function to cal. Z from a given Vg			         *
#***************************************************************************

VgtoZ.func <- function(varexp,samp=totalN){
sqrt( varexp/(1-varexp)*(samp-2) )
}


#***************************************************************************
#                 function to cal. Vg from a given z			                 *
#***************************************************************************

ztoVg.func <-function(Z,samp=totalN){
Z^2 /(samp-2 + Z^2)
}

#******************************************************************************
#  Modified Selection bias correction by Efron's empirical Bayes method       *
#******************************************************************************
bias.corr.kernel<- function (zz,xout,bw="nrd0")

{
   density.obj=density(zz,bw=bw)
   fz.func = splinefun(density.obj$x, density.obj$y)
   Psi.z = log(  fz.func(zz) /dnorm(zz)  )
   truez = D1ss(x= zz, y=Psi.z ,xout=xout)
return(truez)
}

#________________________________________________________________________________________________________________
SumVg.cont <- function (zall) 
{

zall.corr.ker =  bias.corr.kernel(zall,xout=zall)
Vg.corr.ker = ztoVg.func(zall.corr.ker )
sum.kernel = sum(Vg.corr.ker)
return ( sum.kernel)
}

#****************************************************************************************************************************************************
SumVg.cont.givenH1 <- function (zall)

{
#Vgall = ztoVg.func(zall)
zall.corr.ker =  bias.corr.kernel(zall,xout=zall)
fdrobj2 = locfdr(zall,bre=120,df=7,nulltype=2,plot=0)
fdr2 = fdrobj2$fdr
largefdr.ind = which (fdr2>=0.95)
truez.given.H1 = zall.corr.ker[-largefdr.ind]/(1-fdr2[-largefdr.ind])   ##note that we don't consider any z-values with fdr>0.95
Vg.corr = ztoVg.func(truez.given.H1)*(1-fdr2[-largefdr.ind ])

sum.kernel.givenH1 = sum(Vg.corr)

return ( sum.kernel.givenH1)
}

#__________________________________________________________________________________________
if (out=="unconditional") 
{
Est.SumVg = SumVg.cont(zall=zall)  
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
   jackobj =   jack.del.d(x=zall, SumVg.cont, d=d, repl=repl)
   SE.SumVg = jackobj$jack.se
  
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
        SumVg.cont.paraboot[i]  = SumVg.cont (x) 
                  }
        
        return ( list(boot.para.se=sd(SumVg.cont.paraboot),
        bootvalues1=SumVg.cont.paraboot   )  )   
                                      }
boot.para.obj.corrZ = boot.para(listOfZ= zall.corr.ker , repl=repl)
SE.SumVg = boot.para.obj.corrZ$boot.para.se

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
                              
        
        SumVg.cont.paraboot[i]  = SumVg.cont (z.sim)  } 
      
      return ( list(boot.para.se=sd(SumVg.cont.paraboot), bootvalues1=SumVg.cont.paraboot) )       }
        

boot.para.obj.wt1 = boot.para(listOfZ=zall, repl=repl)
SE.SumVg= boot.para.obj.wt1$boot.para.se


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
                              
        
        SumVg.cont.paraboot[i]  = SumVg.cont (z.sim)  } 
      
      return ( list(boot.para.se=sd(SumVg.cont.paraboot), bootvalues1=SumVg.cont.paraboot) )       }
        

boot.para.obj.wt2 = boot.para(listOfZ=zall, repl=repl)
SE.SumVg= boot.para.obj.wt2$boot.para.se

}

} #end of the part '  if (out="unconditional") '
#********************************************************































if (out=="conditional")
{
Est.SumVg = SumVg.cont.givenH1(zall=zall)  
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
   jackobj =   jack.del.d(x=zall, SumVg.cont.givenH1, d=d, repl=repl)
   SE.SumVg = jackobj$jack.se
  
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
        SumVg.cont.paraboot[i]  = SumVg.cont.givenH1 (x) 
                  }
        
        return ( list(boot.para.se=sd(SumVg.cont.paraboot),
        bootvalues1=SumVg.cont.paraboot   )  )   
                                      }
boot.para.obj.corrZ = boot.para(listOfZ= zall.corr.ker , repl=repl)
SE.SumVg = boot.para.obj.corrZ$boot.para.se

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
                              
        
        SumVg.cont.paraboot[i]  = SumVg.cont.givenH1 (z.sim)  } 
      
      return ( list(boot.para.se=sd(SumVg.cont.paraboot), bootvalues1=SumVg.cont.paraboot) )       }
        

boot.para.obj.wt1 = boot.para(listOfZ=zall, repl=repl)
SE.SumVg= boot.para.obj.wt1$boot.para.se


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
                              
        
        SumVg.cont.paraboot[i]  = SumVg.cont.givenH1 (z.sim)  } 
      
      return ( list(boot.para.se=sd(SumVg.cont.paraboot), bootvalues1=SumVg.cont.paraboot) )       }
        

boot.para.obj.wt2 = boot.para(listOfZ=zall, repl=repl)
SE.SumVg= boot.para.obj.wt2$boot.para.se

}

} #end of the part '  if (out="conditional") '
#************************************************************************************************************************************************


return( list(Est.SumVg= Est.SumVg, SE.SumVg= SE.SumVg) ) 


} ##ned of Sum.Vg function