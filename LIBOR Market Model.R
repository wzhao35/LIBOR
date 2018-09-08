setwd("~/Desktop")
#LMM
data<-read.csv('eta.csv',header=T)

#number of dates selected, <=60
maxnum<-60

#observation dates
obserdate<-data$Maturity
obserdate<-obserdate[2:61] 

#tau, time between obserdate
tau<-data$Year.fraction..Tau
tau<-tau[2:61]

#F(0), all forward rates at time 0
fwdrate_0<-data$forward.rate
fwdrate_0<-fwdrate_0[1:60]

#instant vol, from Excel
eta<-data$eta
eta<-eta[2:61]  

#instant correlation
rho_beta<-0.5
rho_alpha<-0.5

rho<-function(ti,tj)
{
  return (rho_alpha+(1-rho_alpha)*exp(-rho_beta*(abs(ti-tj))))
}

#covariance matrix
covrc<-matrix(0,nrow = maxnum,ncol = maxnum)
for(i in 1:maxnum)
{
  for(j in 1:maxnum)
  {
    covrc[i,j]<-eta[i]*eta[j]*rho(obserdate[i],obserdate[j])
  }
}
#cholensky
chol_cov<-t(chol(covrc))

#swaption strikes
swptn_strike<-c(1.04,1.12,1.21,1.29,1.37,1.45,1.53,1.59,1.65,1.70)
swptn_strike<-swptn_strike/100



#simulation starts here, LONG JUMP, Lecture11 Page28
numsim<-50
#store all the security prices
price_euroswaption<-rep(0,10)
price_bermudan<-0

for(lalala in 1:numsim)
{

#fwd rates matrix
fwd_matrix<-matrix(0,nrow = maxnum,ncol = maxnum+1)
fwd_matrix[,1]<-fwdrate_0

#F(t) for the longest F
for(i in 2:(maxnum+1))
{
  phi<-cbind(rnorm(maxnum,0,1))
  fwd_matrix[maxnum,i]<-fwd_matrix[maxnum,i-1]*exp(-0.5*covrc[maxnum,maxnum]*tau[i-1]+chol_cov[maxnum,]%*%phi)
}
#lecture11, page29
for(i in (maxnum-1):1)
{
  for(j in 2:(i+1))
  {
    drift<-0 #for drift, k->i, j->k, i->j
    for(k in (i+1):maxnum)
    {
      predcrect<-fwd_matrix[k,(j-1)]*tau[k]/(1+fwd_matrix[k,(j-1)]*tau[k])+fwd_matrix[k,j]*tau[k]/(1+fwd_matrix[k,j]*tau[k])
      drift<-drift-0.5*covrc[k,i]*predcrect
    }
    
    phi<-cbind(rnorm(i,0,1))
    fwd_matrix[i,j]<-fwd_matrix[i,(j-1)]*exp(drift-0.5*covrc[i,i]*tau[j-1]+chol_cov[i,1:i]%*%phi)
  }
}

#write.csv(fwd_matrix,'fwd_rates.csv')

#transfer forward rates to ZCB price under current simulated fwd_matrix
zcbp<-function(i,j)  #notice here t->i,j=1, T0->i,j=2, T(M-1)->only j =M+1
{
  multiplication<-1
  for (k in 1:j)
  {
    multiplication<-multiplication*(1+fwd_matrix[k,i]*tau[k])
  }
  return(1/multiplication)
}

#European swaption prices
euroswaption<-function(opt_maturity,swap_maturity,strike)
{
  swaprate<-1-zcbp(opt_maturity,swap_maturity)
  for(k in (opt_maturity+1):swap_maturity)
  {
    swaprate<-swaprate-strike*tau[k]*zcbp(opt_maturity,k)
  }
  swaprate<-pnorm(swaprate,0,1)
  return (zcbp(1,swap_maturity)/zcbp(opt_maturity,swap_maturity)*max(swaprate,0))
}

for (i in 1:10)
{
  price_euroswaption[i]<-(lalala-1)/lalala*price_euroswaption[i]+euroswaption(4*i+1,4*i+5,swptn_strike[i])/lalala
}

#Bermudan swaption price
bermudan<-function(opt_maturity,swap_maturity,strike)
{
  #bermudan swaption price for this path
  path_bermudan<-0
  
  #at option maturity, value of the swaption
  intermediate<-0
  for (j in opt_maturity:swap_maturity)
  {
    intermediate<-intermediate+tau[j]*(fwd_matrix[j,opt_maturity]-strike)*zcbp(opt_maturity,j)
  }
  
  #envolve backwards in time
  for (i in (opt_maturity-1):1)
  {
    #early exercise value at time i
    early_exercise<-0
    for(j in (i+1):swap_maturity)
    {
      early_exercise<-early_exercise+tau[j]*(fwd_matrix[j,i]-strike)*zcbp(i,j)
    }
    
    #regression value at time i    T^T
    regression_value<-0
    numreg<-10
    dep<-rep(0,numreg)
    indep1<-rep(0,numreg)
    
    for(k in 1:numreg)
    {
      fk_fwd<-matrix(runif(swap_maturity*swap_maturity,0,0.09), swap_maturity)
      indep1[k]<-mean(fk_fwd[,i])
      
      fk_zcbp<-function(a,b)  #transfer fwd to zcbp under each senario
      {
        multiplication<-1
        for (c in 1:b)
        {
          multiplication<-multiplication*(1+fk_fwd[c,a]*tau[c])
        }
        return(1/multiplication)
      }
      
      dep[k]<-0
      for (j in opt_maturity:swap_maturity)
      {
        dep[k]<-dep[k]+tau[j]*(fk_fwd[j,opt_maturity]-strike)*fk_zcbp(opt_maturity,j)
      }
      dep[k]<-dep[k]*fk_zcbp(i,opt_maturity)
    }
    
    indep2<-indep1*indep1
    #only regress on those ITM senarios
    itm<-(indep1>strike)
    dependent<-dep[itm]
    independent1<-indep1[itm]
    independent2<-indep2[itm]
    linear_regress<-lm(dependent~independent1+independent2)
    reg_coeff<-summary(linear_regress)$coefficients
    regression_value<-reg_coeff[1,1]+reg_coeff[2,1]*mean(fwd_matrix[,i])+reg_coeff[3,1]*mean(fwd_matrix[,i])
    
    #compare early exercise to regression value to decide the path value at time i 
    #continuation value if EE<regression value, EE if else
    if(regression_value<=early_exercise)
    {
      path_bermudan<-early_exercise
    }
    else{
      path_bermudan<-zcbp(i,opt_maturity)*intermediate
    }
    
  }#now reach time 1
 
  return (abs(path_bermudan))
  
}#end of bermudan function

#bermudan swaption price update  ??????????????? need change inputs to real data
price_bermudan<-(lalala-1)/lalala*price_bermudan+bermudan(2,4,0.005)
}#simulation ends here

help(lm)

bermudan(2,4,0.05)
