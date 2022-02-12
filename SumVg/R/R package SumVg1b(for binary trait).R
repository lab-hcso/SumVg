require(locfdr)
require(sfsmisc)

 
SumVg.binary <- function (zall, method="paraboot", d=1, repl=50,out="unconditional",caseNo, ctrlNo, K ) 

{

muaa=0 ##fixed
preval = K ##to be changed for another disease

sum.kernel = numeric(repl)
sum.kernel.givenH1 = numeric(repl)
sum.uncorr = numeric(repl)
true.sum = numeric(repl)


func.RR <- function(RR1,Vg,PA=0.5){
RR2=RR1^2
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2

faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa 
T = qnorm(1-faa) #muaa is set to 0 and residual var set to 1
muAa = T-qnorm(1-fAa)
muAA = T-qnorm(1-fAA)
##overall mean
mean.all= PAa*muAa+ PAA*muAA
expVg= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
expVg2 = expVg/(1+expVg)
return( (expVg2-Vg)^2 ) 
}

####function to cal. the resulting RR from a given Vg##########
resRR.func <-function(varexp){
optimize(func.RR,c(1,1000),Vg=varexp)$minimum}


##function to cal. power ###############
power.func <- function(RR1,RR2=RR1^2, PA=0.5, K, case.size=caseNo,ctrl.size=ctrlNo,alpha=5e-05){
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2
muaa=0 #fixed
faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa 

Paa.case = faa*Paa/K
PAa.case = fAa*PAa/K
PAA.case = fAA*PAA/K

Paa.ctrl = (1-faa)*Paa/(1-K)
PAa.ctrl = (1-fAa)*PAa/(1-K)
PAA.ctrl = (1-fAA)*PAA/(1-K)

PA.case = PAa.case/2 + PAA.case
Pa.case = 1- PA.case
PA.ctrl = PAa.ctrl/2 + PAA.ctrl 
Pa.ctrl = 1- PA.ctrl

ORall = PA.case*Pa.ctrl/Pa.case/PA.ctrl
										
VarlnOR = 1/(2*case.size)*(1/PA.case + 1/Pa.case) +  1/(2*ctrl.size)*(1/PA.ctrl +1/Pa.ctrl)
Z = log(ORall)/sqrt(VarlnOR)

critR = qnorm(1-alpha/2)
critL = qnorm(alpha/2)

power = 1-pnorm(critR,mean=Z,sd=1)+pnorm(critL,mean=Z,sd=1)
res.power.func <- list()
res.power.func$Z = Z
res.power.func$power = power
return(res.power.func) 
}

muaa=0 #fixed

##function to be optimized###########################
func.RR <- function(RR1,Vg,PA=0.5){
RR2=RR1^2
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2

faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa 
T = qnorm(1-faa) #muaa is set to 0 and residual var set to 1
muAa = T-qnorm(1-fAa)
muAA = T-qnorm(1-fAA)
##overall mean
mean.all= PAa*muAa+ PAA*muAA
expVg= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
expVg2 = expVg/(1+expVg)
return( (expVg2-Vg)^2 ) 
}

#***************************************************************************
#                 function to cal. Z from a given Vg			         *
#***************************************************************************

VgtoZ.func <- function(varexp,PA=0.5, K=0.001, case.size=caseNo,ctrl.size=ctrlNo){

RR1 = optimize(func.RR,c(1,100000),Vg=varexp)$minimum
#mapply(resRR.func,c(0.01,0.02)) 

RR2=RR1^2
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2
muaa=0 #fixed
faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa 

Paa.case = faa*Paa/K
PAa.case = fAa*PAa/K
PAA.case = fAA*PAA/K

Paa.ctrl = (1-faa)*Paa/(1-K)
PAa.ctrl = (1-fAa)*PAa/(1-K)
PAA.ctrl = (1-fAA)*PAA/(1-K)

PA.case = PAa.case/2 + PAA.case
Pa.case = 1- PA.case
PA.ctrl = PAa.ctrl/2 + PAA.ctrl 
Pa.ctrl = 1- PA.ctrl


Acase = 2*case.size*PA.case
acase = 2*case.size*Pa.case
Actrl = 2*ctrl.size*PA.ctrl
actrl = 2*ctrl.size*Pa.ctrl
chisqmat = matrix(c(Acase,acase,Actrl,actrl),nrow=2, byrow=F)
Zsq = chisq.test(chisqmat)$statistic
Z = sqrt(Zsq)

#Z = log(ORall)/sqrt(VarlnOR) ##old method; will cause problem when varexp> ~0.2
return(Z) 
}

##must specify caseNo and ctrlNo at the start of your function
##function to be optimized
power.optim <- function(RR1,obsZ,RR2=RR1^2,case.size=caseNo,ctrl.size=ctrlNo,PA=0.5, K=preval){
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2
muaa=0 #fixed
faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa

Paa.case = faa*Paa/K
PAa.case = fAa*PAa/K
PAA.case = fAA*PAA/K

Paa.ctrl = (1-faa)*Paa/(1-K)
PAa.ctrl = (1-fAa)*PAa/(1-K)
PAA.ctrl = (1-fAA)*PAA/(1-K)

PA.case = PAa.case/2 + PAA.case
Pa.case = 1- PA.case
PA.ctrl = PAa.ctrl/2 + PAA.ctrl
Pa.ctrl = 1- PA.ctrl

ORall = PA.case*Pa.ctrl/Pa.case/PA.ctrl
							
VarlnOR = 1/(2*case.size)*(1/PA.case + 1/Pa.case) +  1/(2*ctrl.size)*(1/PA.ctrl +1/Pa.ctrl)
Z = log(ORall)/sqrt(VarlnOR)



return( (Z-obsZ)^2 )
}

##########################################################################################
ztoVg.func <-function(obsZ,PA=0.5,K=preval){
#RR1=optimize(power.optim,c(0,5),obsZ=obsZ,case.size=caseNo,ctrl.size=ctrlNo,PA=0.5,K=preval)$minimum
RR1=nlminb(1.13,power.optim,obsZ=obsZ,case.size=caseNo,ctrl.size=ctrlNo,PA=0.5,K=preval)$par
RR2=RR1^2

Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2
muaa=0 #fixed

faa= K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa

T = qnorm(1-faa) #muaa is set to 0 and residual var set to 1
muAa = T-qnorm(1-fAa)
muAA = T-qnorm(1-fAA)

##overall mean
mean.all= PAa*muAa+ PAA*muAA
Vg= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
return( Vg/(1+Vg) )
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
SumVg.cont <- function (zall) 
{
zall.corr.ker =  bias.corr.kernel(zall,xout=zall)
Vg.corr.ker = mapply( ztoVg.func, obsZ= zall.corr.ker  )
sum.kernel = sum(Vg.corr.ker)
return ( sum.kernel)
}

#****************************************************************************************************************************************************
SumVg.cont.givenH1 <- function (zall)

{   
zall.corr.ker =  bias.corr.kernel(zall,xout=zall)
fdr = locfdr(zall,bre=500,df=10,nulltype=0,plot=0)$fdr
largefdr.ind = which (fdr>=0.95) 
truez.given.H1 = zall.corr.ker[-largefdr.ind]/(1-fdr[-largefdr.ind])   ##note that we don't consider any z-values with fdr>0.95   
Vg.corr = mapply( ztoVg.func, obsZ= truez.given.H1 )*(1-fdr[-largefdr.ind]) 
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


} ##end of SumVg.binary function