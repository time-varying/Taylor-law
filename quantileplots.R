library(ggplot2)
library(dplyr)
library(tidyr)
library(pracma)
library(Rlab)


#The following function compute the normalized statistics (see Figure 1 of the paper).

ic2 <- function(t,n,lbd){
 
   #First case: square of Gaussian.
  param=c(log(2),2)
 mat <- (matrix(rnorm(n*t,0,1), nrow=n, ncol=t)^2)*matrix(lbd,n,t,byrow=T)

  #Second case: Poisson variables
     # param=c(0,1)   
 #mat <- matrix(NA, nrow=n, ncol=t)
  #mat <- t(apply(mat,1,function(x) rpois(t,lbd)))
  
  # sample mean/variance
  moy <- apply(mat, 2, mean)
  moy2 <- apply(mat^2, 2, mean)
  sigma <- apply(mat, 2, var)
  
    ### 1- Dt ###
  dt_hat <- matrix(c(1,(1/t)*sum(log(moy)),(1/t)*sum(log(moy)),(1/t)*sum(log(moy)^2)),2,2,byrow=T)
  design=inv(dt_hat)
  estimators=inv(dt_hat)%*%c(mean(log(moy2-moy^2)),mean(log(moy)*log(moy2-moy^2)))
 theta1=estimators[1]
 theta2=estimators[2]
   # -----  mat var-covar
  var_covar <- list()
 for(i in 1:length(sigma)){
    var_covar[[i]] <- matrix(c(sigma[i],cov(mat[,i],mat[,i]^2),cov(mat[,i],mat[,i]^2),var(mat[,i]^2)),2,2,byrow=T)
 }
  
 
  ### 2- Ct ###
  # ----- mat jacobienne
  phi <- function(x){
    coef=c(theta1,theta2)
    c(log(x[2]-x[1]^2)-coef[1]-coef[2]*log(x[1]),log(x[1])*log(x[2]-x[1]^2)-coef[1]*log(x[1])-coef[2]*(log(x[1]))^2)
  }
  
    j <- list()
    for(i in 1:length(moy)){
   j[[i]] <- jacobian(phi,c(moy[i],moy2[i]))
    }
  phi2 <- function(x){
    coef=c(theta1,theta2)
    log(x[1])*log(x[2]-x[1]^2)-coef[1]*log(x[1])-coef[2]*(log(x[1]))^2
  }
  Jphi1 <- function(x,y){
    coef=c(theta1,theta2)
    a=-2/(y-x^2)-4*(x^2)/(y-x^2)^2+coef[2]/x^2
    b=2*x/(y-x^2)^2
    c=-1/(y-x^2)^2
  matrix(c(a,b,b,c),2,2)
  }
  library(pracma)
 
  j <- list()
  jjj<-list()
  
  for(i in 1:length(moy)){
    j[[i]] <- jacobian(phi,c(moy[i],moy2[i]))
    jjj[[i]]<-hessian(phi2,c(moy[i],moy2[i]))
  
    }
  # -----  mat var-covar

# ----- variance
  mult_mat <- list()
  CORBIAS1=list()
  CORBIAS2=list()
  for (i in 1: length(j)){
    mult_mat[[i]] <- j[[i]]%*%var_covar[[i]]%*%t(j[[i]])
    CORBIAS1[[i]]<-Jphi1(moy[i],moy2[i])%*%var_covar[[i]]
    CORBIAS2[[i]]<-jjj[[i]]%*%var_covar[[i]]
    }
  
  cible <- Reduce("+", mult_mat)/t # somme des matrices
  
  CT=inv(sqrtm(cible)$B)
  CORB1<-Reduce("+",CORBIAS1)/(2*n*t)
  CORB2<-Reduce("+",CORBIAS2)/(2*n*t)
  CORB1=CORB1[1,1]+CORB1[2,2]
  CORB2=CORB2[1,1]+CORB2[2,2]
 
StatCOR=sqrt(n*t)*CT%*%dt_hat%*%matrix(c(theta1-param[1],theta2-param[2]),2,1)
 
  H=design%*%cible%*%design
 # m=sqrt(n*t)*design%*%matrix(c(CORB1,CORB2),2,1)
# m=matrix(c(0,0),2,1)
 # ren1=(sqrt(n*t)*(theta1-param[1])-m[1])/sqrt(H[1,1])
 # ren2=(sqrt(n*t)*(theta2-param[2])-m[2])/sqrt(H[2,2])
 # Finalp=StatCOR
  #Below the version with correction
Finalp=StatCOR-sqrt(n*t)*CT%*%matrix(c(CORB1,CORB2),2,1)

    return(Finalp)
  }

#quantile/quantile plot to check Gaussianity.

n =1000
t=n

#Choice of the time-varying mean.

lbd =exp(cos(1:t))+1+sin(1:t)
#lbd=3+cos(1:t)


repet <-100

testing<-replicate(repet,ic2(t,n,lbd))  

par(mar=c(5.1,5.1,4.1,2.1))
qqnorm(testing[1,1,],cex.main=2,cex.lab=2,cex.axis=2,tck=-0.03,lwd=3)
abline(0,1,lwd=3,col='red')
qqnorm(testing[2,1,],cex.main=2,cex.lab=2,cex.axis=2,tck=-0.03,lwd=3)
abline(0,1,lwd=3,col='red')
