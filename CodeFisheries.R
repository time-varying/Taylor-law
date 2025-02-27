
library(ggplot2)
library(dplyr)
library(tidyr)
library(pracma)
library(Rlab)
#_______________________________________________________________________________
#                         I) DATA NS-IBTS ----
#_______________________________________________________________________________

# Standard Species : pas comme dans l'article

#setwd("")
data_ns <- read.csv("D:/The Best/recherche stats/TaylorLaw/CPUE per length per subarea.csv", head=T)

attach(data_ns)

data_ns$Quarter <- as.factor(Quarter)
data_ns$Area <- as.factor(Area)
data_ns$SubArea <- as.factor(SubArea)
data_ns$AphiaID <- as.factor(AphiaID)
data_ns$Species <- as.factor(Species)

data_ns <- data_ns[,2:9]



yy <- as.factor(Year)
unique(yy) # 39 ans
y <- summary(yy)
mean(y) # 16137.82 obs par années en moy
median(y) # 16053 

summary(Species) # unique(AphiaID)
unique(Species) # 9 espèces
unique(SubArea) #195 subarea



year <- unique(data_ns$Year)

#_____________________________
# 1.1) Groupe par espèces ----
#_____________________________

# Liste du tableau de données initial par espèces
data_esp <- split(data_ns, f = data_ns$Species)  

class(data_ns$Species) 

# Liste des séries temporelles par espèces
library(reshape2)
library(stringr)

serie_esp <- lapply(data_esp, function(x){ dcast(x,x$SubArea ~ x$Year,value.var = "CPUE_number_per_hour",fun.aggregate = sum)})
serie_esp <- lapply(serie_esp, function(x){ row.names(x)<-x[,1]; x})
serie_esp <-  lapply(serie_esp,function(x) x[,-1])


# Plot une série temporelle
cpue_serie_temp <- t(serie_esp[[1]][1,])
serie_temp <- cbind(year,cpue_serie_temp)
serie_temp <- as.data.frame(serie_temp)

ggplot(serie_temp, aes(year, cpue_serie_temp)) +  
  geom_point() +
  geom_line() +
  xlab("Years") +
  ylab("CPUE") +
  scale_x_discrete(limits=year)+
  theme(axis.text.x = element_text(angle=90))



# ----------------  Spatial TL ---------------- 

#Tableaux des moy et variance pour chaque espèces / années en fonction des subArea
mt <- sapply(serie_esp,function(x) apply(x,2,mean))
mt2 <- sapply(serie_esp, function(x) apply(x^2,2,mean))
st <- sapply(serie_esp,function(x) apply(x,2,var))
st=sqrt(st)



#checking significance of the three first autocorrelations.
i=1
m=matrix(mt[,i],39,195)
m=t(m)
sigma=matrix(st[,i],39,195)
sigma=t(sigma)
A=as.matrix(serie_esp[[i]])
A=A-m
A=A/sigma
CR=matrix(0,195,3)
for (i in 1:195)
{CF=acf(A[i,])
CF=CF[[1]]
  for (j in 1:3){
  r=CF[j+1]
CR[i,j]=(tanh(atanh(r)-1.96/sqrt(39))<0)&(0<tanh(atanh(r)+1.96/sqrt(39)))  
  }
}
apply(CR,2,mean)



#Log/log plots for the variance/mean of a given species (here i=9).
i=9
reg=lm(log(st[,i])~log(mt[,i]))
plot(log(mt[,i]),log(st[,i]),xlab="Log(Mean)", ylab="Log(Variance)")
abline(reg, col="red")








#Consider species 1, normalize the abundance and compute the confidence bands with or 
#without correction for the slope parameter.
i=1
A=as.matrix(serie_esp[[i]])
A=A/mean(A)
IC(39,195,A)


#Function IC for computing the confidence intervals of the slope.


IC <- function(t,n,donnees){
  
  
  moy <- apply(donnees,2,mean)
  moy2 <- apply(donnees^2,2,mean)
  sigma <- moy2-moy^2
  mat=donnees
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
  
  H=design%*%cible%*%design
  m=sqrt(n*t)*design%*%matrix(c(CORB1,CORB2),2,1)
  imoins=(sqrt(H[2,2])*1.96+m[2])/sqrt(n*t)
  iplus=(sqrt(H[2,2])*1.96-m[2])/sqrt(n*t)
  Finalp=matrix(0,2,3)
  Finalp[1,]=c(theta2,theta2-imoins,theta2+iplus)
  imoins=(sqrt(H[2,2])*1.96)/sqrt(n*t)
  iplus=(sqrt(H[2,2])*1.96)/sqrt(n*t)
  Finalp[2,]=c(theta2,theta2-imoins,theta2+iplus)
  
  return(Finalp)
}




