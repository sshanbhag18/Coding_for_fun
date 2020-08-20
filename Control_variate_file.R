library(Sim.DiffProc)
library(LSMonteCarlo)
geometric_brownianmotion=GBM(100,1000,100,t0=0,T=1,theta = 0.01,sigma = 0.1)
View(geometric_brownianmotion)
final_prices=geometric_brownianmotion[101,]
final_prices
strike_price=105
payoffinitial=final_prices-strike_price
payoff=ifelse(payoffinitial<0,0,payoffinitial)
payoff
beta=(cov(final_prices,payoff,method = "pearson")/var(final_prices))
beta  
discounted_initialprice=strike_price*exp(0.05*0.5)
discounted_initialprice
a=(final_prices-discounted_initialprice)*beta
num=payoff-a
num
control_variate=mean(num)
control_variate
call_price=(mean(payoff)*exp(-0.01*0.1))
call_price
call_price_BS=EuPutBS(40,0.2,40,0.05,0,0.5)
call_price_BS
######variance reduction due to control variates######
rho=cor(payoff,final_prices,method = "pearson")
rho
p=var(payoff)
b=var(final_prices)
c=sd(payoff)
d=sd(final_prices)
variance_controlvariate=(p+(beta^2*b)-2*beta*rho*c*d)/10
variance_controlvariate  
var(num)
########European put option###
geometric_brownianmotion2=GBM(100,10,40,t0=0,T=0.5,theta = 0.05,sigma = 0.2)
initial_price2=40
View(geometric_brownianmotion2)
final_prices2=geometric_brownianmotion[101,]
final_prices2
strike_price2=40
payoffinitial2=final_prices2-strike_price2
payoffinitial2
payoffinitial
payoff2=ifelse(payoffinitial2<0,0,payoffinitial2)
beta2=(cov(final_prices2,payoff2,method = "pearson")/var(final_prices2))
beta2  
discounted_initialprice2=strike_price2*exp(0.05*0.5)
discounted_initialprice2
p2=(final_prices2-discounted_initialprice2)*beta2
num2=payoff2-p2
num2
control_variate2=mean(num2)
control_variate2
put_price=(mean(payoff2)*exp(-0.05*0.5))
put_price
######variance reduction due to control variates######
rho2=cor(payoff2,final_prices2,method = "pearson")
rho2
a2=var(payoff2)
b2=var(final_prices2)
c2=sd(payoff2)
d2=sd(final_prices2)
variance_controlvariate2=(a2+((beta2)^2*b2)-2*beta2*rho2*c2*d2)/10
variance_controlvariate2  
var(num2)
#########American option########
library(LSMonteCarlo)########Package for monte carlo simulations
library(Sim.DiffProc)#####Package to generate GBM
AmerPut=function (Spot,sigma,n,m,Strike,r,dr,mT){
  GBM=matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    GBM[i,]=Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  print(GBM)
  X=ifelse(GBM<Strike,GBM,NA)
  Payoff=matrix(pmax(0,Strike-GBM),nrow=n,ncol=m)
  Price=X[,-m]
  Price_square=Price*Price
  Dis_payoffs=Payoff*exp(-1*r*(mT/m))
  Dis_payoffs_2=cbind((matrix(NA, nrow=n, ncol=m-1)), Dis_payoffs[,m])
  Cont_val=matrix(NA, nrow=n, ncol=m-1)
  for (i in m-1:1) {
    reg1=lm(Dis_payoffs_2[,i+1]~Price[,i]+Price_square[,i])
    Cont_val[,i]=(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Dis_payoffs [,i])+((matrix(reg1$coefficients)[3,1])*Dis_payoffs_2[,i])
    Cont_val[,i]=(ifelse(is.na(Cont_val[,i]),0,Cont_val[,i]))
    Dis_payoffs_2[,i]<-ifelse(Payoff[,i]>Cont_val[,i], Dis_payoffs[,i], Dis_payoffs_2[,i+1]*exp(-1*r*(mT/m)))
  }
  Cont_val=ifelse(is.na(Cont_val),0,Cont_val)######To calculate the continuation value##
  Cont_val2=cbind(Cont_val,(matrix(0,nrow = n,ncol = 1)))
  Payoff_initial=ifelse(Cont_val2>Payoff,0,Payoff)
  Payoff_final=firstValueRow(Payoff_initial)#####GBM value########
  dPayoff_final=matrix(NA,nrow=n,ncol=m)
  for (i in 1:m)
  {
    dPayoff_final[,i]=Payoff_final[,i]*exp(-1*mT/m*r*i)
    
  }
  Put_price<-mean(rowSums(dPayoff_final))
  res<- list(price=(Put_price), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res)<-"AmerPut"######American Put######
  return(res)
}
AmerPut(Spot=100, sigma=0.1, n=1000, m=100, Strike=105, r=0.01, dr=0.0, mT=1)
########American Call############
AmerCall=function (Spot,sigma,n,m,Strike,r,dr,mT){
  GBM=matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    GBM[i,]=Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X=ifelse(GBM<Strike,GBM,NA)
  Payoff=matrix(pmax(0,GBM-Strike),nrow=n,ncol=m)#####Payoff
  Price=X[,-m]
  Price_square=Price*Price
  Dis_payoffs=Payoff*exp(-1*r*(mT/m))
  Dis_payoffs_2=cbind((matrix(NA, nrow=n, ncol=m-1)), Dis_payoffs[,m])
  Cont_val=matrix(NA, nrow=n, ncol=m-1)
  for (i in m-1:1) {
    reg1=lm(Dis_payoffs_2[,i+1]~Price[,i]+Price_square[,i])
    Cont_val[,i]=(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Dis_payoffs [,i])+((matrix(reg1$coefficients)[3,1])*Dis_payoffs_2[,i])
    Cont_val[,i]=(ifelse(is.na(Cont_val[,i]),0,Cont_val[,i]))
    Dis_payoffs_2[,i]<-ifelse(Payoff[,i]>Cont_val[,i], Dis_payoffs[,i], Dis_payoffs_2[,i+1]*exp(-1*r*(mT/m)))
  }
  Cont_val=ifelse(is.na(Cont_val),0,Cont_val)
  Cont_val2=cbind(Cont_val,(matrix(0,nrow = n,ncol = 1)))#####continuation value###
  Payoff_initial=ifelse(Cont_val2>Payoff,0,Payoff)
  Payoff_final=firstValueRow(Payoff_initial)
  dPayoff_final=matrix(NA,nrow=n,ncol=m)
  for (i in 1:m)
  {
    dPayoff_final[,i]=Payoff_final[,i]*exp(-1*mT/m*r*i)
    
  }
  Call_price<-mean(rowSums(dPayoff_final))
  res<- list(price=(Call_price), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res)<-"AmerCall"
  return(res)
}
Amer_call=AmerCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=1)

Amer_call$price #########American Call price###

#####Control variates for american valuePut value###
AmerPutLSM_CV=function (Spot,sigma,n,m,Strike,r,dr,mT){
  GBM=matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    GBM[i,]=Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X=ifelse(GBM<Strike,GBM,NA)
  Payoff=matrix(pmax(0,Strike-GBM),nrow=n,ncol=m)
  Price=X[,-m]
  Price_square=Price*Price
  Dis_payoffs=Payoff*exp(-1*r*(mT/m))
  Dis_payoffs_2=cbind((matrix(NA, nrow=n, ncol=m-1)), Dis_payoffs[,m])
  Cont_val=matrix(NA, nrow=n, ncol=m-1)
  for (i in m-1:1) {
    reg1=lm(Dis_payoffs_2[,i+1]~Price[,i]+Price_square[,i])
    Cont_val[,i]=(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Dis_payoffs [,i])+((matrix(reg1$coefficients)[3,1])*Dis_payoffs_2[,i])
    Cont_val[,i]=(ifelse(is.na(Cont_val[,i]),0,Cont_val[,i]))
    Dis_payoffs_2[,i]<-ifelse(Payoff[,i]>Cont_val[,i], Dis_payoffs[,i], Dis_payoffs_2[,i+1]*exp(-1*r*(mT/m)))
  }
  Cont_val=ifelse(is.na(Cont_val),0,Cont_val)######To calculate the continuation value##
  Cont_val2=cbind(Cont_val,(matrix(0,nrow = n,ncol = 1)))
  Payoff_initial=ifelse(Cont_val2>Payoff,0,Payoff)
  Payoff_final=firstValueRow(Payoff_initial)#####GBM value########
  dPayoff_final=matrix(NA,nrow=n,ncol=m)
  for (i in 1:m)
  {
    dPayoff_final[,i]=Payoff_final[,i]*exp(-1*mT/m*r*i)
    
  }
  Put_price=mean(rowSums(dPayoff_final))
  Dis_payoffs=Payoff[,m]*exp(-1*r*mT)
  Eursimul=mean(Dis_payoffs)
  EuBS=EuPutBS(Spot,sigma,Strike,r,dr,mT)
  Put_price_CV=Put_price-(Eursimul-EuBS)
  res= list(price=(Put_price_CV), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res)="AmerPut_CV"######American Put######
  return(res)
}
AmerPutLSM_CV(Spot=1, sigma=0.2, n=1000, m=100, Strike=1.1, r=0.06, dr=0.0, mT=1)$price
########Price for Amer_call_CV##########
AmerCallLSM_CV=function (Spot,sigma,n,m,Strike,r,dr,mT){
  GBM=matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    GBM[i,]=Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X=ifelse(GBM<Strike,GBM,NA)
  Payoff=matrix(pmax(0,GBM-Strike),nrow=n,ncol=m)
  Price=X[,-m]
  Price_square=Price*Price
  Dis_payoffs=Payoff*exp(-1*r*(mT/m))
  Dis_payoffs_2=cbind((matrix(NA, nrow=n, ncol=m-1)), Dis_payoffs[,m])
  Cont_val=matrix(NA, nrow=n, ncol=m-1)
  for (i in m-1:1) {
    reg1=lm(Dis_payoffs_2[,i+1]~Price[,i]+Price_square[,i])
    Cont_val[,i]=(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Dis_payoffs [,i])+((matrix(reg1$coefficients)[3,1])*Dis_payoffs_2[,i])
    Cont_val[,i]=(ifelse(is.na(Cont_val[,i]),0,Cont_val[,i]))
    Dis_payoffs_2[,i]<-ifelse(Payoff[,i]>Cont_val[,i], Dis_payoffs[,i], Dis_payoffs_2[,i+1]*exp(-1*r*(mT/m)))
  }
  Cont_val=ifelse(is.na(Cont_val),0,Cont_val)######To calculate the continuation value##
  Cont_val2=cbind(Cont_val,(matrix(0,nrow = n,ncol = 1)))
  Payoff_initial=ifelse(Cont_val2>Payoff,0,Payoff)
  Payoff_final=firstValueRow(Payoff_initial)#####GBM value########
  dPayoff_final=matrix(NA,nrow=n,ncol=m)
  for (i in 1:m)
  {
    dPayoff_final[,i]=Payoff_final[,i]*exp(-1*mT/m*r*i)
    
  }
  Put_price=mean(rowSums(dPayoff_final))
  Dis_payoffs=Payoff[,m]*exp(-1*r*mT)
  Eursimul=mean(Dis_payoffs)
  EuBS=EuPutBS(Spot,sigma,Strike,r,dr,mT)
  Call_price_CV=Amer_call$price-(Eursimul-EuBS)
  res= list(price=(Call_price_CV), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res)="AmerCall_CV"######American_Call_CV######
  return(res)
}
AmerCallLSM_CV(Spot=1, sigma=0.2, n=1000, m=100, Strike=1.1, r=0.06, dr=0.0, mT=1)$price
######Bermuda_PutPrice########### 
BermPut=function (Spot,sigma,n,m,Strike,r,dr,mT){
  GBM=matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    GBM[i,]=Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X=ifelse(GBM<Strike,GBM,NA)
  Payoff=matrix(pmax(0,Strike-GBM),nrow=n,ncol=m)
  Price=X[,-m]
  Price_square=Price*Price
  Dis_payoffs=Payoff*exp(-1*r*(mT/m))
  Dis_payoffs_2=cbind((matrix(NA, nrow=n, ncol=m-1)), Dis_payoffs[,m])
  Cont_val=matrix(NA, nrow=n, ncol=m-1)
  for (i in m-1:1) {
    reg1=lm(Dis_payoffs_2[,i+1]~Price[,i]+Price_square[,i])
    Cont_val[,i]=(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Dis_payoffs [,i])+((matrix(reg1$coefficients)[3,1])*Dis_payoffs_2[,i])
    Cont_val[,i]=(ifelse(is.na(Cont_val[,i]),0,Cont_val[,i]))
    Dis_payoffs_2[,i]<-ifelse(Payoff[,i]>Cont_val[,i], Dis_payoffs[,i], Dis_payoffs_2[,i+1]*exp(-1*r*(mT/m)))
  }
  Cont_val=ifelse(is.na(Cont_val),0,Cont_val)
  Cont_val2=cbind(Cont_val,(matrix(0,nrow = n,ncol = 1)))
  Payoff_initial=ifelse(Cont_val2>Payoff,0,Payoff)
  Payoff_final=firstValueRow(Payoff_initial)
  dPayoff_final=matrix(NA,nrow=n,ncol=m)
  for (i in 1:m)
  {
    dPayoff_final[,i]=Payoff_final[,i]*exp(-1*mT/m*r*i)
  }
  Put_price<-mean(rowSums(dPayoff_final))
  res<- list(price=(Put_price), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res)<-"BermPut"
  return(res)
}
a=BermPut(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.1)$price*(exp(-0.06*0.1))
b=BermPut(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.2)$price*(exp(-0.06*0.1))
c=BermPut(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0,mT=0.3)$price*(exp(-0.06*0.1))
d=BermPut(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.4)$price*(exp(-0.06*0.1))
e=BermPut(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.5)$price*(exp(-0.06*0.1))
f=BermPut(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.6)$price*(exp(-0.06*0.1))
g=BermPut(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.7)$price*(exp(-0.06*0.1))
h=BermPut(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.8)$price*(exp(-0.06*0.1))
k=BermPut(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.9)$price*(exp(-0.06*0.1))
put_price_bermuda=mean(cbind(a,b,c,d,e,f,g,h,k))
put_price_bermuda
#########BermudaCall_price#######
BermCall=function (Spot,sigma,n,m,Strike,r,dr,mT){
  GBM=matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    GBM[i,]=Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X=ifelse(GBM<Strike,GBM,NA)
  Payoff=matrix(pmax(0,GBM-Strike),nrow=n,ncol=m)
  Price=X[,-m]
  Price_square=Price*Price
  Dis_payoffs=Payoff*exp(-1*r*(mT/m))
  Dis_payoffs_2=cbind((matrix(NA, nrow=n, ncol=m-1)), Dis_payoffs[,m])
  Cont_val=matrix(NA, nrow=n, ncol=m-1)
  for (i in m-1:1) {
    reg1=lm(Dis_payoffs_2[,i+1]~Price[,i]+Price_square[,i])
    Cont_val[,i]=(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Dis_payoffs [,i])+((matrix(reg1$coefficients)[3,1])*Dis_payoffs_2[,i])
    Cont_val[,i]=(ifelse(is.na(Cont_val[,i]),0,Cont_val[,i]))
    Dis_payoffs_2[,i]<-ifelse(Payoff[,i]>Cont_val[,i], Dis_payoffs[,i], Dis_payoffs_2[,i+1]*exp(-1*r*(mT/m)))
  }
  Cont_val=ifelse(is.na(Cont_val),0,Cont_val)
  Cont_val2=cbind(Cont_val,(matrix(0,nrow = n,ncol = 1)))
  Payoff_initial=ifelse(Cont_val2>Payoff,0,Payoff)
  Payoff_final=firstValueRow(Payoff_initial)
  dPayoff_final=matrix(NA,nrow=n,ncol=m)
  for (i in 1:m)
  {
    dPayoff_final[,i]=Payoff_final[,i]*exp(-1*mT/m*r*i)
  }
  Call_price<-mean(rowSums(dPayoff_final))
  res<- list(price=(Call_price), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res)<-"BermCall"
  return(res)
}
a=BermCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.1)$price*(exp(-0.06*0.1))
b=BermCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.2)$price*(exp(-0.06*0.1))
c=BermCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0,mT=0.3)$price*(exp(-0.06*0.1))
d=BermCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.4)$price*(exp(-0.06*0.1))
e=BermCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.5)$price*(exp(-0.06*0.1))
f=BermCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.6)$price*(exp(-0.06*0.1))
g=BermCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.7)$price*(exp(-0.06*0.1))
h=BermCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.8)$price*(exp(-0.06*0.1))
k=BermCall(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.9)$price*(exp(-0.06*0.1))
call_price_bermuda=mean(cbind(a,b,c,d,e,f,g,h,k))
call_price_bermuda
#######Control_variate_for_Bermuda####
BermPutLSM_CV=function (Spot,sigma,n,m,Strike,r,dr,mT){
  GBM=matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    GBM[i,]=Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X=ifelse(GBM<Strike,GBM,NA)
  Payoff=matrix(pmax(0,Strike-GBM),nrow=n,ncol=m)
  Price=X[,-m]
  Price_square=Price*Price
  Dis_payoffs=Payoff*exp(-1*r*(mT/m))
  Dis_payoffs_2=cbind((matrix(NA, nrow=n, ncol=m-1)), Dis_payoffs[,m])
  Cont_val=matrix(NA, nrow=n, ncol=m-1)
  for (i in m-1:1) {
    reg1=lm(Dis_payoffs_2[,i+1]~Price[,i]+Price_square[,i])
    Cont_val[,i]=(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Dis_payoffs [,i])+((matrix(reg1$coefficients)[3,1])*Dis_payoffs_2[,i])
    Cont_val[,i]=(ifelse(is.na(Cont_val[,i]),0,Cont_val[,i]))
    Dis_payoffs_2[,i]<-ifelse(Payoff[,i]>Cont_val[,i], Dis_payoffs[,i], Dis_payoffs_2[,i+1]*exp(-1*r*(mT/m)))
  }
  Cont_val=ifelse(is.na(Cont_val),0,Cont_val)######To calculate the continuation value##
  Cont_val2=cbind(Cont_val,(matrix(0,nrow = n,ncol = 1)))
  Payoff_initial=ifelse(Cont_val2>Payoff,0,Payoff)
  Payoff_final=firstValueRow(Payoff_initial)#####GBM value########
  dPayoff_final=matrix(NA,nrow=n,ncol=m)
  for (i in 1:m)
  {
    dPayoff_final[,i]=Payoff_final[,i]*exp(-1*mT/m*r*i)
    
  }
  Put_price=mean(rowSums(dPayoff_final))
  Dis_payoffs=Payoff[,m]*exp(-1*r*mT)
  Eursimul=mean(Dis_payoffs)
  EuBS=EuPutBS(Spot,sigma,Strike,r,dr,mT)
  Put_price_CV=Put_price-(Eursimul-EuBS)
  res= list(price=(Put_price_CV), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res)="BermPutLSM_CV"######BermudaCV######
  return(res)
}
a=BermPutLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.1)$price*(exp(-0.06*0.1))
b=BermPutLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.2)$price*(exp(-0.06*0.1))
c=BermPutLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0,mT=0.3)$price*(exp(-0.06*0.1))
d=BermPutLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.4)$price*(exp(-0.06*0.1))
e=BermPutLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.5)$price*(exp(-0.06*0.1))
f=BermPutLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.6)$price*(exp(-0.06*0.1))
g=BermPutLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.7)$price*(exp(-0.06*0.1))
h=BermPutLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.8)$price*(exp(-0.06*0.1))
k=BermPutLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.9)$price*(exp(-0.06*0.1))
put_price_bermuda_cv=mean(cbind(a,b,c,d,e,f,g,h,k))
put_price_bermuda_cv
#######Control_variate for Bermuda_call###
BermCallLSM_CV=function (Spot,sigma,n,m,Strike,r,dr,mT){
  GBM=matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    GBM[i,]=Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X=ifelse(GBM<Strike,GBM,NA)
  Payoff=matrix(pmax(0,GBM-Strike),nrow=n,ncol=m)
  Price=X[,-m]
  Price_square=Price*Price
  Dis_payoffs=Payoff*exp(-1*r*(mT/m))
  Dis_payoffs_2=cbind((matrix(NA, nrow=n, ncol=m-1)), Dis_payoffs[,m])
  Cont_val=matrix(NA, nrow=n, ncol=m-1)
  for (i in m-1:1) {
    reg1=lm(Dis_payoffs_2[,i+1]~Price[,i]+Price_square[,i])
    Cont_val[,i]=(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Dis_payoffs [,i])+((matrix(reg1$coefficients)[3,1])*Dis_payoffs_2[,i])
    Cont_val[,i]=(ifelse(is.na(Cont_val[,i]),0,Cont_val[,i]))
    Dis_payoffs_2[,i]<-ifelse(Payoff[,i]>Cont_val[,i], Dis_payoffs[,i], Dis_payoffs_2[,i+1]*exp(-1*r*(mT/m)))
  }
  Cont_val=ifelse(is.na(Cont_val),0,Cont_val)######To calculate the continuation value##
  Cont_val2=cbind(Cont_val,(matrix(0,nrow = n,ncol = 1)))
  Payoff_initial=ifelse(Cont_val2>Payoff,0,Payoff)
  Payoff_final=firstValueRow(Payoff_initial)#####GBM value########
  dPayoff_final=matrix(NA,nrow=n,ncol=m)
  for (i in 1:m)
  {
    dPayoff_final[,i]=Payoff_final[,i]*exp(-1*mT/m*r*i)
    
  }
  Put_price=mean(rowSums(dPayoff_final))
  Dis_payoffs=Payoff[,m]*exp(-1*r*mT)
  Eursimul=mean(Dis_payoffs)
  EuBS=EuPutBS(Spot,sigma,Strike,r,dr,mT)
  Put_price_CV=Put_price-(Eursimul-EuBS)
  res= list(price=(Put_price_CV), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res)="BermPutLSM_CV"######BermudaCV######
  return(res)
}
a=BermCallLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.1)$price*(exp(-0.06*0.1))
b=BermCallLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.2)$price*(exp(-0.06*0.1))
c=BermCallLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0,mT=0.3)$price*(exp(-0.06*0.1))
d=BermCallLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.4)$price*(exp(-0.06*0.1))
e=BermCallLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.5)$price*(exp(-0.06*0.1))
f=BermCallLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.6)$price*(exp(-0.06*0.1))
g=BermCallLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.7)$price*(exp(-0.06*0.1))
h=BermCallLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.8)$price*(exp(-0.06*0.1))
k=BermCallLSM_CV(Spot=40, sigma=0.2, n=1000, m=100, Strike=40, r=0.06, dr=0.0, mT=0.9)$price*(exp(-0.06*0.1))
call_price_bermuda_cv=mean(cbind(a,b,c,d,e,f,g,h,k))
call_price_bermuda_cv

