###open source data:
library(dlnm)
library(data.table)
library(splines)
library(mgcv)
library(forecast)

data.cv<-chicagoNMMAPS
summary(chicagoNMMAPS)


par(mfrow=c(4,1))
plot(data.cv$time, data.cv$resp,type="l")
##remove extreme value?? is it caused by air pollutent??

plot(data.cv$time, data.cv$pm10,type="l")
plot(data.cv$time, data.cv$o3,type="l")
plot(data.cv$time, data.cv$temp,type="l")

par(mfrow=c(1,2))
acf(data.cv$o3)
acf(data.cv$pm10)


##demo

cb1.pm <- crossbasis(data.cv$pm10[1:1000], lag=5, argvar=list(fun="lin"),
                       arglag=list(fun="poly",degree=4))
cb1.temp <- crossbasis(data.cv$temp, lag=3, argvar=list(df=5),
                         arglag=list(fun="strata",breaks=1))
model1 <- glm(death ~ cb1.pm + ns(temp, 4) + ns(time, 20) + as.numeric(dow),
              family=quasipoisson(), chicagoNMMAPS[1:1000,])

sss<-summary(model1)
sss

pred1.pm <- crosspred(cb1.pm, model1, at=10, bylag=1, cumul=TRUE)

plot(pred1.pm, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2), main="Association with a 10-unit increase in PM10")
plot(pred1.pm, "slices", var=10, col=2, cumul=TRUE, ylab="Cumulative RR",
       main="Cumulative association with a 10-unit increase in PM10")

##Extract single day effet and CIs?
dim(pred1.pm$matRRfit)
pred1.pm$matRRfit-1
pred1.pm$matRRlow
pred1.pm$matRRhigh


par(mfrow=c(1,2))
o3<-data.cv$o3
mean(abs(acf(o3)$acf[2:6]))
##missing value in pm10
which(is.na(data.cv$pm10))
data.cv$pm10<-na.interp(data.cv$pm10)
mean(abs(acf(pm10)$acf[2:6]))


##the ac level can range from 0.1 to 0.6 or higher...

###====================================================================

auto.arima(o3-smooth.spline(o3,df=50)$y)

tr<- 8*sin(time/58-1)+30
X<-8*arima.sim(list(ar=0.31, ma=0.23),1000)+tr

plot(X,type="l",ylim=c(0,70))
mean(acf(X)$acf[2:6])




###============================PREDICTOR=====================================

### now we are going to randomly generate time series for 3000 times.
simsize=1000
ar.coef<-seq(0.2,0.8,by=0.1)
simdat<-NULL

tr<-8*sin(time[1:1000]/58)+40

simdat<-NULL
simac<-NULL

for (i in 1:length(ar.coef)){
  temp<-NULL
  tempac<-NULL
  for (j in 1:simsize){
    temp<-rbind(temp,8*arima.sim(list(ar=ar.coef[i],ma=0.1),1000)+tr)
    tempac<-c(tempac, mean(abs(acf(temp[j,])$acf[2:6])))
  }
  simdat<-rbind(simdat, temp)
  simac<-c(simac,tempac)
}

data<-as.data.frame(simdat)
data$ac<-simac
dim(data)
summary(data$ac)
hist(data$ac)


data<-data[order(data$ac),]


###==========================COVARAITES========================================
time<-seq(1,1000,by=1)
s<-ns(time,df=20)
coef_time<-as.vector(c(-0.106890934, -0.109039657,-0.226018013,-0.066980289,-0.103221086,
                       0.043865213,-0.004445378, 0.030004012,-0.166819518,-0.108588736,
                       -0.259338882,0.107178911, -0.118470791,-0.187391080, 0.145316375,
                       -0.355723728,0.090670258, -0.284979146, -0.220898118,-0.023576640))
sm_time<-as.matrix(s)%*%coef_time

DOW <- c(rep(c(1,2,3,4,5,6,7),150))[1:1000]

temp<-3*rnorm(1000,0,1)+8*sin(time/60)+70
s<-ns(temp,df=4)
coef_temp<-as.vector(c(0.12404438, 0.03116323, 0.30611978, 0.13411584))

sm_temp<-as.matrix(s)%*%coef_temp





###=========================3 DAY LAG EFFECT+===================================
###================================MODELS=======================================
aic.glm<-NULL;coef.glm<-NULL; CI.glm<-NULL; sd.glm<-NULL
aic.gam<-NULL; coef.gam<-NULL; CI.gam<-NULL; sd.gam<-NULL
aic.dlnm<-NULL; coef.dlnm<-NULL; CI.dlnm<-NULL; sd.dlnm<-NULL

N<-length(data$ac)
N


for (i in 1:N){
  print(i)
  X<-as.numeric(data[i,1:1000])
  X_1<-c(X[1:999],mean(X[1:999]))
  X_2<-c(X[1:998], mean(X[1:998]), mean(X[1:998]))
  mu<-exp(0.005*X+0.001*X_1+0.002*X_2+sm_time+sm_temp-0.01*DOW+2.68)
  Y<-rpois(1000,mu)


  ##glm
  fit.glm<-glm(Y~X+ns(time,df=20)+ns(temp,df=4)+DOW,family=poisson(link="log"))
  sum.1<-summary(fit.glm)$coefficients
  coef.glm<-c(coef.glm, sum.1[2,1])
  sd.glm<-c(sd.glm, sum.1[2,2])
  CI.glm<-rbind(CI.glm, c(sum.1[2,1]-1.96*sum.1[2,2], sum.1[2,1]+1.96*sum.1[2,2]))
  aic.glm<-c(aic.glm, AIC(fit.glm))

  ##gam
  s=mgcv:::s
  fit.gam<-gam( Y ~ X+s(time, k=20)+s(temp, k=4)+DOW, family=poisson(link = "log"))
  sum.2<-summary(fit.gam)
  coef.gam<-c(coef.gam, sum.2$p.coeff[2])
  sd.gam<-c(sd.gam, sum.2$se[2])
  CI.gam<-rbind(CI.gam, c(sum.2$p.coeff[2]-1.96*sum.2$se[2], sum.2$p.coeff[2]+1.96*sum.2$se[2]))
  aic.gam<-c(aic.gam, AIC(fit.gam))
  
  ##dlnm
  basis.X<-crossbasis(X,lag=10,argvar=list(fun="lin"), arglag=list(fun="poly",degree=3))
  fit.dlnm <-gam(Y ~ s(time,k=20)+s(temp,k=4)+DOW+basis.X, family = poisson(link="log"))
  sum.3<-summary(fit.dlnm)
  
  pred.dlnm<-crosspred(basis.X, fit.dlnm, at=1, bylag=1, cumul=TRUE)
  coef.dlnm<-rbind(coef.dlnm, pred.dlnm$matfit[1:3])
  sd.dlnm<-rbind(sd.dlnm,pred.dlnm$matse[1:3])
  CI.dlnm<-rbind(CI.dlnm, log(c(pred.dlnm$matRRlow[1],pred.dlnm$matRRhigh[1],
                                pred.dlnm$matRRlow[2],pred.dlnm$matRRhigh[2], 
                                pred.dlnm$matRRlow[3],pred.dlnm$matRRhigh[3])))
  aic.dlnm<-c(aic.dlnm, AIC(fit.dlnm))
  
}




par(mfrow=c(1,3))


#glm
length(coef.dlnm[,1])
length(data$ac)


plot(x=data$ac,y=smooth.spline(coef.glm[1:length(data$ac)],df=700)$y,type="l",ylim=c(-0.001,0.015),ylab="beta",xlab="auto-correlation index",main="GLM Performance")
lines(x=data$ac,y=smooth.spline(CI.glm[1:length(data$ac),1],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.glm[1:length(data$ac),2],df=700)$y,lwd=0.5,lty=2)
abline(a=0.005,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta0", "95% CI for beta0", "true beta0", "no effect beta0=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)


#gam
plot(x=data$ac,y=smooth.spline(coef.gam[1:length(data$ac)],df=700)$y,type="l",ylim=c(-0.001,0.015),ylab="beta",xlab="auto-correlation index",main="GAM Performance")
lines(x=data$ac,y=smooth.spline(CI.gam[1:length(data$ac),1],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.gam[1:length(data$ac),2],df=700)$y,lwd=0.5,lty=2)
abline(a=0.005,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta0", "95% CI for beta0", "true beta0", "no effect beta0=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)

par(mfrow=c(1,3))
#dlnm
plot(x=data$ac,y=smooth.spline(coef.dlnm[1:length(data$ac),1],df=700)$y,type="l",ylim=c(-0.001,0.01),ylab="beta for lag0",xlab="auto-correlation index",main="DLNM Performance")
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),1],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),2],df=700)$y,lwd=0.5,lty=2)
abline(a=0.005,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta1", "95% CI for beta0", "true beta0", "no effect beta0=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)


plot(x=data$ac,y=smooth.spline(coef.dlnm[1:length(data$ac),2],df=700)$y,type="l",ylim=c(-0.001,0.01),ylab="beta for lag1",xlab="auto-correlation index",main="DLNM Performance")
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),3],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),4],df=700)$y,lwd=0.5,lty=2)
abline(a=0.001,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta1", "95% CI for beta1", "true beta1", "no effect beta1=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)

plot(x=data$ac,y=smooth.spline(coef.dlnm[1:length(data$ac),3],df=700)$y,type="l",ylim=c(-0.001,0.01),ylab="beta for lag2",xlab="auto-correlation index",main="DLNM Performance")
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),5],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),6],df=700)$y,lwd=0.5,lty=2)
abline(a=0.002,b=0,lty=1,col="red")
abline(0,0,lty=1,col="blue")
legend("topright", c("Estimated beta2", "95% CI for beta2", "true beta2", "no effect beta2=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)

par(mfrow=c(1,2))
##variance compare

plot(x=data$ac,y=smooth.spline(sd.glm[1:length(data$ac)],df=700)$y,type="l",ylim=c(0,0.001),ylab="sd(beta)", xlab="auto-correlation index",main="Standard Error of Beta's")
lines(x=data$ac,y=smooth.spline(sd.gam[1:length(data$ac)],df=700)$y,col="red",lty=2)
lines(x=data$ac,y=smooth.spline(sd.dlnm[1:length(data$ac),1],df=700)$y,col="green",lty=1)
lines(x=data$ac, y=smooth.spline(sd.dlnm[1:length(data$ac),2],df=700)$y, col="cyan",lty=2)
lines(x=data$ac, y=smooth.spline(sd.dlnm[1:length(data$ac),3],df=700)$y, col="dark green",lty=2)

legend("topright",c("SE(beta) of GLM", "SE(beta) of GAM", "SE(beta0) of DLNM", "SE(beta1) of DLNM", "SE(beta2) of DLNM"), 
       col=c("black","red","green","cyan","dark green"), lty=c(1,2,1,2,2),adj=0,cex=0.7)

##AIC compare
plot(x=data$ac,y=smooth.spline(aic.glm[1:length(data$ac)],df=700)$y,type="l",ylim=c(5700,6000),ylab="AIC",xlab="auto-correlation index", main="AIC of Each model")
lines(x=data$ac,y=smooth.spline(aic.gam[1:length(data$ac)],df=700)$y,col="red",lty=2)
lines(x=data$ac,y=smooth.spline(aic.dlnm[1:length(data$ac)],df=700)$y,col="green",lty=2)
legend("topright",c("AIC of GLM", "AIC of GAM", "AIC of DLNM"),lty=c(1,2,2), col=c("black","red","green"),cex=0.7)

##Summary Statistics
result<-as.data.frame(cbind(data$ac, coef.glm, CI.glm, sd.glm, coef.gam, CI.gam, sd.gam, coef.dlnm, CI.dlnm, sd.dlnm))
colnames(result)<-c("ac", "beta.glm", "LCB.glm", "UCB.glm", "sd.glm", 'beta.gam', 'LCB.gam','UCB.gam', "sd.gam", 
                    "beta0.dlnm","beta1.dlnm","beta2.dlnm","LCB0.dlnm","UCB0.dlnm","LCB1.dlnm","UCB1.dlnm","LCB2.dlnm","UCB2.dlnm","sd0.dlnm","sd1.dlnm","sd2.dlnm")
rownames(result)<-seq(1,7000,by=1)
attach(result)

head(result)


difference<-matrix(c(mean(beta.glm[result$ac<0.4]-0.005),
                     mean(beta.gam[result$ac<0.4]-0.005),
                     mean(beta0.dlnm[result$ac<0.4]-0.005),
                     mean(beta1.dlnm[result$ac<0.4]-0.001),
                     mean(beta2.dlnm[result$ac<0.4]-0.002),
                     mean(beta.glm[result$ac>0.4 & result$ac<0.6]-0.005),
                     mean(beta.gam[result$ac>0.4 & result$ac<0.6]-0.005),
                     mean(beta0.dlnm[result$ac>0.4 & result$ac<0.6]-0.005),
                     mean(beta1.dlnm[result$ac>0.4 & result$ac<0.6]-0.001),
                     mean(beta2.dlnm[result$ac>0.4 & result$ac<0.6]-0.002),
                     mean(beta.glm[result$ac>0.6]-0.005),
                     mean(beta.gam[result$ac>0.6]-0.005),
                     mean(beta0.dlnm[result$ac>0.6]-0.005),
                     mean(beta1.dlnm[result$ac>0.6]-0.001),
                     mean(beta2.dlnm[result$ac>0.6]-0.002)),nrow=5,ncol=3)


colnames(difference)<-c("low","mid","high")
rownames(difference)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")



coverage<-matrix(c(mean(LCB.glm[result$ac<0.4]<0.005 & UCB.glm[result$ac<0.4]>0.005),
                 mean(LCB.gam[result$ac<0.4]<0.005 & UCB.gam[result$ac<0.4]>0.005),
                 mean(LCB0.dlnm[result$ac<0.4]<0.005 & UCB0.dlnm[result$ac<0.4]>0.005),
                 mean(LCB1.dlnm[result$ac<0.4]<0.001 & UCB1.dlnm[result$ac<0.4]>0.001),
                 mean(LCB2.dlnm[result$ac<0.4]<0.002 & UCB2.dlnm[result$ac<0.4]>0.002),
                 mean(LCB.glm[result$ac>0.4 & result$ac<0.6]<0.005 & UCB.glm[result$ac>0.4 & result$ac<0.6]>0.005),
                 mean(LCB.gam[result$ac>0.4 & result$ac<0.6]<0.005 & UCB.gam[result$ac>0.4 & result$ac<0.6]>0.005),
                 mean(LCB0.dlnm[result$ac>0.4 & result$ac<0.6]<0.005 & UCB0.dlnm[result$ac>0.4 & result$ac<0.6]>0.005),
                 mean(LCB1.dlnm[result$ac>0.4 & result$ac<0.6]<0.001 & UCB1.dlnm[result$ac>0.4 & result$ac<0.6]>0.001),
                 mean(LCB2.dlnm[result$ac>0.4 & result$ac<0.6]<0.002 & UCB2.dlnm[result$ac>0.4 & result$ac<0.6]>0.002),
                 mean(LCB.glm[result$ac>0.6]<0.005 & UCB.glm[result$ac>0.6]>0.005),
                 mean(LCB.gam[result$ac>0.6]<0.005 & UCB.gam[result$ac>0.6]>0.005),
                 mean(LCB0.dlnm[result$ac>0.6]<0.005 & UCB0.dlnm[result$ac>0.6]>0.005),
                 mean(LCB1.dlnm[result$ac>0.6]<0.001 & UCB1.dlnm[result$ac>0.6]>0.001),
                 mean(LCB2.dlnm[result$ac>0.6]<0.002 & UCB2.dlnm[result$ac>0.6]>0.002)),nrow=5,ncol=3,byrow=FALSE)

colnames(coverage)<-c("low","mid","high")
rownames(coverage)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")


temp<-matrix(c(mean(LCB.glm[result$ac<0.4]<0 & UCB.glm[result$ac<0.4]>0),
                       mean(LCB.gam[result$ac<0.4]<0 & UCB.gam[result$ac<0.4]>0),
                       mean(LCB0.dlnm[result$ac<0.4]<0 & UCB0.dlnm[result$ac<0.4]>0),
                       mean(LCB1.dlnm[result$ac<0.4]<0 & UCB1.dlnm[result$ac<0.4]>0),
                       mean(LCB2.dlnm[result$ac<0.4]<0 & UCB2.dlnm[result$ac<0.4]>0),
                       mean(LCB.glm[result$ac>0.4 & result$ac<0.6]<0 & UCB.glm[result$ac>0.4 & result$ac<0.6]>0),
                       mean(LCB.gam[result$ac>0.4 & result$ac<0.6]<0 & UCB.gam[result$ac>0.4 & result$ac<0.6]>0),
                       mean(LCB0.dlnm[result$ac>0.4 & result$ac<0.6]<0 & UCB0.dlnm[result$ac>0.4 & result$ac<0.6]>0),
                       mean(LCB1.dlnm[result$ac>0.4 & result$ac<0.6]<0 & UCB1.dlnm[result$ac>0.4 & result$ac<0.6]>0),
                       mean(LCB2.dlnm[result$ac>0.4 & result$ac<0.6]<0 & UCB2.dlnm[result$ac>0.4 & result$ac<0.6]>0),
                       mean(LCB.glm[result$ac>0.6]<0 & UCB.glm[result$ac>0.6]>0),
                       mean(LCB.gam[result$ac>0.6]<0 & UCB.gam[result$ac>0.6]>0),
                       mean(LCB0.dlnm[result$ac>0.6]<0 & UCB0.dlnm[result$ac>0.6]>0),
                       mean(LCB1.dlnm[result$ac>0.6]<0 & UCB1.dlnm[result$ac>0.6]>0),
                       mean(LCB2.dlnm[result$ac>0.6]<0 & UCB2.dlnm[result$ac>0.6]>0)),nrow=5,ncol=3,byrow=FALSE)

significance<-1-temp

colnames(significance)<-c("low","mid","high")
rownames(significance)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")

beta<-matrix(c(mean(beta.glm[result$ac<0.4]),mean(beta.gam[result$ac<0.4]),
             mean(beta0.dlnm[result$ac<0.4]),mean(beta1.dlnm[result$ac<0.4]),
             mean(beta2.dlnm[result$ac<0.4]),
             mean(beta.glm[result$ac>0.4 & result$ac<0.6]),mean(beta.gam[result$ac>0.4 & result$ac<0.6]),
             mean(beta0.dlnm[result$ac>0.4 & result$ac<0.6]),mean(beta1.dlnm[result$ac>0.4 & result$ac<0.6]),
             mean(beta2.dlnm[result$ac>0.4 & result$ac<0.6]),
             mean(beta.glm[result$ac>0.6]),mean(beta.gam[result$ac>0.6]),
             mean(beta0.dlnm[result$ac>0.6]),mean(beta1.dlnm[result$ac>0.6]),
             mean(beta2.dlnm[result$ac>0.6]),0.005,0.005,0.005,0.001,0.002),nrow=5,ncol=4)
colnames(beta)<-c("low","mid","high","TRUE beta")
rownames(beta)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")

SD<-matrix(c(mean(sd.glm[result$ac<0.4]),mean(sd.gam[result$ac<0.4]),
             mean(sd0.dlnm[result$ac<0.4]),mean(sd0.dlnm[result$ac<0.4]),
             mean(sd2.dlnm[result$ac<0.4]),
             mean(sd.glm[result$ac>0.4 & result$ac<0.6]),mean(sd.gam[result$ac>0.4 & result$ac<0.6]),
             mean(sd0.dlnm[result$ac>0.4 & result$ac<0.6]),mean(sd0.dlnm[result$ac>0.4 & result$ac<0.6]),
             mean(sd2.dlnm[result$ac>0.4 & result$ac<0.6]),
             mean(sd.glm[result$ac>0.6]),mean(sd.gam[result$ac>0.6]),
             mean(sd0.dlnm[result$ac>0.6]),mean(sd0.dlnm[result$ac>0.6]),
             mean(sd2.dlnm[result$ac>0.6])),nrow=5,ncol=3)

colnames(SD)<-c("low","mid","high")
rownames(SD)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")


beta
difference
SD
coverage
significance



###What Wendy Asked For
###Pick a random trial, look into all the coefficients:

X<-as.numeric(data[2501,1:1000])
X_1<-c(X[1:999],mean(X[1:999]))
X_2<-c(X[1:998], mean(X[1:998]), mean(X[1:998]))
mu<-exp(0.005*X+0.001*X_1+0.002*X_2+sm_time+sm_temp-0.01*DOW+2.68)
Y<-rpois(1000,mu)

fit.glm<-glm(Y~X+ns(time,df=20)+ns(temp,df=4)+DOW,family=poisson(link="log"))
summary(fit.glm)


s=mgcv:::s
fit.gam<-gam( Y ~ X+s(time, k=20)+s(temp, k=4)+DOW, family=poisson(link = "log"))
summary(fit.gam)

basis.X<-crossbasis(X,lag=10,argvar=list(fun="lin"), arglag=list(fun="poly",degree=3))
fit.dlnm <-gam(Y ~ s(time,k=20)+s(temp,k=4)+DOW+basis.X, family = poisson(link="log"))
summary(fit.dlnm)









###============================SINGLE DAY EFFECT====================================
###================================MODELS=======================================

aic.glm<-NULL;coef.glm<-NULL; CI.glm<-NULL; sd.glm<-NULL
aic.gam<-NULL; coef.gam<-NULL; CI.gam<-NULL; sd.gam<-NULL
aic.dlnm<-NULL; coef.dlnm<-NULL; CI.dlnm<-NULL; sd.dlnm<-NULL

N<-length(data$ac)


for (i in 1:N){
  print(i)
  X<-as.numeric(data[i,1:1000])
  X_1<-c(X[1:999],mean(X[1:999]))
  X_2<-c(X[1:998], mean(X[1:998]), mean(X[1:998]))
  mu<-exp(0.005*X+sm_time+sm_temp-0.01*DOW+2.68)
  Y<-rpois(1000,mu)
  
  ##glm
  fit.glm<-glm(Y~X+ns(time,df=20)+ns(temp,df=4)+DOW,family=poisson(link="log"))
  sum.1<-summary(fit.glm)$coefficients
  coef.glm<-c(coef.glm, sum.1[2,1])
  sd.glm<-c(sd.glm, sum.1[2,2])
  CI.glm<-rbind(CI.glm, c(sum.1[2,1]-1.96*sum.1[2,2], sum.1[2,1]+1.96*sum.1[2,2]))
  aic.glm<-c(aic.glm, AIC(fit.glm))
  
  ##gam
  fit.gam<-gam(Y~X+s(time,k=20)+s(temp,k=4)+DOW, family=poisson(link = "log"))
  sum.2<-summary(fit.gam)
  coef.gam<-c(coef.gam, sum.2$p.coeff[2])
  sd.gam<-c(sd.gam, sum.2$se[2])
  CI.gam<-rbind(CI.gam, c(sum.2$p.coeff[2]-1.96*sum.2$se[2], sum.2$p.coeff[2]+1.96*sum.2$se[2]))
  aic.gam<-c(aic.gam, AIC(fit.gam))
  
  ##dlnm
  basis.X<-crossbasis(X,lag=5,argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  fit.dlnm <-gam(Y ~ s(time,k=20)+s(temp,k=4)+DOW+basis.X, family = poisson(link="log"))
  sum.3<-summary(fit.dlnm)
  
  pred.dlnm<-crosspred(basis.X, fit.dlnm, at=1, bylag=1, cumul=TRUE)
  coef.dlnm<-rbind(coef.dlnm, pred.dlnm$matfit[1:3])
  sd.dlnm<-rbind(sd.dlnm,pred.dlnm$matse[1:3])
  CI.dlnm<-rbind(CI.dlnm, log(c(pred.dlnm$matRRlow[1],pred.dlnm$matRRhigh[1],
                                pred.dlnm$matRRlow[2],pred.dlnm$matRRhigh[2], 
                                pred.dlnm$matRRlow[3],pred.dlnm$matRRhigh[3])))
  aic.dlnm<-c(aic.dlnm, AIC(fit.dlnm))
}



par(mfrow=c(1,3))


#glm
length(coef.dlnm[,1])
length(data$ac)


plot(x=data$ac,y=smooth.spline(coef.glm[1:length(data$ac)],df=700)$y,type="l",ylim=c(-0.001,0.015),ylab="beta",xlab="auto-correlation index",main="GLM Performance")
lines(x=data$ac,y=smooth.spline(CI.glm[1:length(data$ac),1],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.glm[1:length(data$ac),2],df=700)$y,lwd=0.5,lty=2)
abline(a=0.005,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta0", "95% CI for beta0", "true beta0", "no effect beta0=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)


#gam
plot(x=data$ac,y=smooth.spline(coef.gam[1:length(data$ac)],df=700)$y,type="l",ylim=c(-0.001,0.015),ylab="beta",xlab="auto-correlation index",main="GAM Performance")
lines(x=data$ac,y=smooth.spline(CI.gam[1:length(data$ac),1],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.gam[1:length(data$ac),2],df=700)$y,lwd=0.5,lty=2)
abline(a=0.005,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta0", "95% CI for beta0", "true beta0", "no effect beta0=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)

par(mfrow=c(1,3))
#dlnm
plot(x=data$ac,y=smooth.spline(coef.dlnm[1:length(data$ac),1],df=700)$y,type="l",ylim=c(-0.0025,0.01),ylab="beta for lag0",xlab="auto-correlation index",main="DLNM Performance")
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),1],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),2],df=700)$y,lwd=0.5,lty=2)
abline(a=0.005,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta1", "95% CI for beta0", "true beta0", "no effect beta0=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)


plot(x=data$ac,y=smooth.spline(coef.dlnm[1:length(data$ac),2],df=700)$y,type="l",ylim=c(-0.0025,0.01),ylab="beta for lag1",xlab="auto-correlation index",main="DLNM Performance")
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),3],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),4],df=700)$y,lwd=0.5,lty=2)
abline(a=0.001,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta1", "95% CI for beta1", "true beta1", "no effect beta1=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)

plot(x=data$ac,y=smooth.spline(coef.dlnm[1:length(data$ac),3],df=700)$y,type="l",ylim=c(-0.0025,0.01),ylab="beta for lag2",xlab="auto-correlation index",main="DLNM Performance")
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),5],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),6],df=700)$y,lwd=0.5,lty=2)
abline(a=0.002,b=0,lty=1,col="red")
abline(0,0,lty=1,col="blue")
legend("topright", c("Estimated beta2", "95% CI for beta2", "true beta2", "no effect beta2=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)

par(mfrow=c(1,2))
##variance compare

plot(x=data$ac,y=smooth.spline(sd.glm[1:length(data$ac)],df=700)$y,type="l",ylim=c(0.0004,0.0013),ylab="sd(beta)", xlab="auto-correlation index",main="Standard Error of Beta's")
lines(x=data$ac,y=smooth.spline(sd.gam[1:length(data$ac)],df=700)$y,col="red",lty=2)
lines(x=data$ac,y=smooth.spline(sd.dlnm[1:length(data$ac),1],df=700)$y,col="green",lty=1)


legend("topright",c("SE(beta) of GLM", "SE(beta) of GAM", "SE(beta0) of DLNM", "SE(beta1) of DLNM", "SE(beta2) of DLNM"), 
       col=c("black","red","green","cyan","dark green"), lty=c(1,2,1,2,2),adj=0,cex=0.7)

##AIC compare
plot(x=data$ac,y=smooth.spline(aic.glm[1:length(data$ac)],df=700)$y,type="l",ylim=c(5600,5900),ylab="AIC",xlab="auto-correlation index", main="AIC of Each model")
lines(x=data$ac,y=smooth.spline(aic.gam[1:length(data$ac)],df=700)$y,col="red",lty=2)
lines(x=data$ac,y=smooth.spline(aic.dlnm[1:length(data$ac)],df=700)$y,col="green",lty=2)
legend("topright",c("AIC of GLM", "AIC of GAM", "AIC of DLNM"),lty=c(1,2,2), col=c("black","red","green"),cex=0.7)

##Summary Statistics
result<-as.data.frame(cbind(data$ac, coef.glm, CI.glm, coef.gam, CI.gam, coef.dlnm, CI.dlnm))
colnames(result)<-c("ac", "beta.glm", "LCB.glm", "UCB.glm", 'beta.gam', 'LCB.gam','UCB.gam',"beta0.dlnm","beta1.dlnm","beta2.dlnm","LCB0.dlnm","UCB0.dlnm","LCB1.dlnm","UCB1.dlnm","LCB2.dlnm","UCB2.dlnm")
rownames(result)<-seq(1,7000,by=1)
attach(result)

head(result)


difference<-matrix(c(mean(beta.glm[result$ac<0.4]-0.005),
                     mean(beta.gam[result$ac<0.4]-0.005),
                     mean(beta0.dlnm[result$ac<0.4]-0.005),
                     mean(beta1.dlnm[result$ac<0.4]-0.001),
                     mean(beta2.dlnm[result$ac<0.4]-0.002),
                     mean(beta.glm[result$ac>0.4 & result$ac<0.6]-0.005),
                     mean(beta.gam[result$ac>0.4 & result$ac<0.6]-0.005),
                     mean(beta0.dlnm[result$ac>0.4 & result$ac<0.6]-0.005),
                     mean(beta1.dlnm[result$ac>0.4 & result$ac<0.6]-0.001),
                     mean(beta2.dlnm[result$ac>0.4 & result$ac<0.6]-0.002),
                     mean(beta.glm[result$ac>0.6]-0.005),
                     mean(beta.gam[result$ac>0.6]-0.005),
                     mean(beta0.dlnm[result$ac>0.6]-0.005),
                     mean(beta1.dlnm[result$ac>0.6]-0.001),
                     mean(beta2.dlnm[result$ac>0.6]-0.002)),nrow=5,ncol=3)


colnames(difference)<-c("low","mid","high")
rownames(difference)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")



coverage<-matrix(c(mean(LCB.glm[result$ac<0.4]<0.005 & UCB.glm[result$ac<0.4]>0.005),
                   mean(LCB.gam[result$ac<0.4]<0.005 & UCB.gam[result$ac<0.4]>0.005),
                   mean(LCB0.dlnm[result$ac<0.4]<0.005 & UCB0.dlnm[result$ac<0.4]>0.005),
                   mean(LCB1.dlnm[result$ac<0.4]<0.001 & UCB1.dlnm[result$ac<0.4]>0.001),
                   mean(LCB2.dlnm[result$ac<0.4]<0.002 & UCB2.dlnm[result$ac<0.4]>0.002),
                   mean(LCB.glm[result$ac>0.4 & result$ac<0.6]<0.005 & UCB.glm[result$ac>0.4 & result$ac<0.6]>0.005),
                   mean(LCB.gam[result$ac>0.4 & result$ac<0.6]<0.005 & UCB.gam[result$ac>0.4 & result$ac<0.6]>0.005),
                   mean(LCB0.dlnm[result$ac>0.4 & result$ac<0.6]<0.005 & UCB0.dlnm[result$ac>0.4 & result$ac<0.6]>0.005),
                   mean(LCB1.dlnm[result$ac>0.4 & result$ac<0.6]<0.001 & UCB1.dlnm[result$ac>0.4 & result$ac<0.6]>0.001),
                   mean(LCB2.dlnm[result$ac>0.4 & result$ac<0.6]<0.002 & UCB2.dlnm[result$ac>0.4 & result$ac<0.6]>0.002),
                   mean(LCB.glm[result$ac>0.6]<0.005 & UCB.glm[result$ac>0.6]>0.005),
                   mean(LCB.gam[result$ac>0.6]<0.005 & UCB.gam[result$ac>0.6]>0.005),
                   mean(LCB0.dlnm[result$ac>0.6]<0.005 & UCB0.dlnm[result$ac>0.6]>0.005),
                   mean(LCB1.dlnm[result$ac>0.6]<0.001 & UCB1.dlnm[result$ac>0.6]>0.001),
                   mean(LCB2.dlnm[result$ac>0.6]<0.002 & UCB2.dlnm[result$ac>0.6]>0.002)),nrow=5,ncol=3,byrow=FALSE)

colnames(coverage)<-c("low","mid","high")
rownames(coverage)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")


temp<-matrix(c(mean(LCB.glm[result$ac<0.4]<0 & UCB.glm[result$ac<0.4]>0),
               mean(LCB.gam[result$ac<0.4]<0 & UCB.gam[result$ac<0.4]>0),
               mean(LCB0.dlnm[result$ac<0.4]<0 & UCB0.dlnm[result$ac<0.4]>0),
               mean(LCB1.dlnm[result$ac<0.4]<0 & UCB1.dlnm[result$ac<0.4]>0),
               mean(LCB2.dlnm[result$ac<0.4]<0 & UCB2.dlnm[result$ac<0.4]>0),
               mean(LCB.glm[result$ac>0.4 & result$ac<0.6]<0 & UCB.glm[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB.gam[result$ac>0.4 & result$ac<0.6]<0 & UCB.gam[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB0.dlnm[result$ac>0.4 & result$ac<0.6]<0 & UCB0.dlnm[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB1.dlnm[result$ac>0.4 & result$ac<0.6]<0 & UCB1.dlnm[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB2.dlnm[result$ac>0.4 & result$ac<0.6]<0 & UCB2.dlnm[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB.glm[result$ac>0.6]<0 & UCB.glm[result$ac>0.6]>0),
               mean(LCB.gam[result$ac>0.6]<0 & UCB.gam[result$ac>0.6]>0),
               mean(LCB0.dlnm[result$ac>0.6]<0 & UCB0.dlnm[result$ac>0.6]>0),
               mean(LCB1.dlnm[result$ac>0.6]<0 & UCB1.dlnm[result$ac>0.6]>0),
               mean(LCB2.dlnm[result$ac>0.6]<0 & UCB2.dlnm[result$ac>0.6]>0)),nrow=5,ncol=3,byrow=FALSE)

significance<-1-temp

colnames(significance)<-c("low","mid","high")
rownames(significance)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")


beta<-matrix(c(mean(beta.glm[result$ac<0.4]),mean(beta.gam[result$ac<0.4]),
               mean(beta0.dlnm[result$ac<0.4]),mean(beta1.dlnm[result$ac<0.4]),
               mean(beta2.dlnm[result$ac<0.4]),
               mean(beta.glm[result$ac>0.4 & result$ac<0.6]),mean(beta.gam[result$ac>0.4 & result$ac<0.6]),
               mean(beta0.dlnm[result$ac>0.4 & result$ac<0.6]),mean(beta1.dlnm[result$ac>0.4 & result$ac<0.6]),
               mean(beta2.dlnm[result$ac>0.4 & result$ac<0.6]),
               mean(beta.glm[result$ac>0.6]),mean(beta.gam[result$ac>0.6]),
               mean(beta0.dlnm[result$ac>0.6]),mean(beta1.dlnm[result$ac>0.6]),
               mean(beta2.dlnm[result$ac>0.6]),0.005,0.005,0.005,0.001,0.002),nrow=5,ncol=4)
colnames(beta)<-c("low","mid","high","TRUE beta")
rownames(beta)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")

SD<-matrix(c(mean(sd.glm[result$ac<0.4]),mean(sd.gam[result$ac<0.4]),
             mean(sd0.dlnm[result$ac<0.4]),mean(sd0.dlnm[result$ac<0.4]),
             mean(sd2.dlnm[result$ac<0.4]),
             mean(sd.glm[result$ac>0.4 & result$ac<0.6]),mean(sd.gam[result$ac>0.4 & result$ac<0.6]),
             mean(sd0.dlnm[result$ac>0.4 & result$ac<0.6]),mean(sd0.dlnm[result$ac>0.4 & result$ac<0.6]),
             mean(sd2.dlnm[result$ac>0.4 & result$ac<0.6]),
             mean(sd.glm[result$ac>0.6]),mean(sd.gam[result$ac>0.6]),
             mean(sd0.dlnm[result$ac>0.6]),mean(sd0.dlnm[result$ac>0.6]),
             mean(sd2.dlnm[result$ac>0.6])),nrow=5,ncol=3)

colnames(SD)<-c("low","mid","high")
rownames(SD)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")


beta
SD
difference
coverage
significance









###===============================NO EFFECT=====================================
###================================MODELS=======================================
aic.glm<-NULL;coef.glm<-NULL; CI.glm<-NULL; sd.glm<-NULL
aic.gam<-NULL; coef.gam<-NULL; CI.gam<-NULL; sd.gam<-NULL
aic.dlnm<-NULL; coef.dlnm<-NULL; CI.dlnm<-NULL; sd.dlnm<-NULL

N<-length(data$ac)



for (i in 1:N){
  print(i)
  X<-as.numeric(data[i,1:1000])
  X_1<-c(X[1:999],mean(X[1:999]))
  X_2<-c(X[1:998], mean(X[1:998]), mean(X[1:998]))
  mu<-exp(sm_time+sm_temp-0.01*DOW+2.68)
  Y<-rpois(1000,mu)
  
  ##glm
  fit.glm<-glm(Y~X+ns(time,df=20)+ns(temp,df=4)+DOW,family=poisson(link="log"))
  sum.1<-summary(fit.glm)$coefficients
  coef.glm<-c(coef.glm, sum.1[2,1])
  sd.glm<-c(sd.glm, sum.1[2,2])
  CI.glm<-rbind(CI.glm, c(sum.1[2,1]-1.96*sum.1[2,2], sum.1[2,1]+1.96*sum.1[2,2]))
  aic.glm<-c(aic.glm, AIC(fit.glm))
  
  ##gam
  fit.gam<-gam(Y~X+s(time,k=20)+s(temp,k=4)+DOW, family=poisson(link = "log"))
  sum.2<-summary(fit.gam)
  coef.gam<-c(coef.gam, sum.2$p.coeff[2])
  sd.gam<-c(sd.gam, sum.2$se[2])
  CI.gam<-rbind(CI.gam, c(sum.2$p.coeff[2]-1.96*sum.2$se[2], sum.2$p.coeff[2]+1.96*sum.2$se[2]))
  aic.gam<-c(aic.gam, AIC(fit.gam))
  
  ##dlnm
  basis.X<-crossbasis(X,lag=5,argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  fit.dlnm <-gam(Y ~ s(time,k=20)+s(temp,k=4)+DOW+basis.X, family = poisson(link="log"))
  sum.3<-summary(fit.dlnm)
  
  pred.dlnm<-crosspred(basis.X, fit.dlnm, at=1, bylag=1, cumul=TRUE)
  coef.dlnm<-rbind(coef.dlnm, pred.dlnm$matfit[1:3])
  sd.dlnm<-rbind(sd.dlnm,pred.dlnm$matse[1:3])
  CI.dlnm<-rbind(CI.dlnm, log(c(pred.dlnm$matRRlow[1],pred.dlnm$matRRhigh[1],
                                pred.dlnm$matRRlow[2],pred.dlnm$matRRhigh[2], 
                                pred.dlnm$matRRlow[3],pred.dlnm$matRRhigh[3])))
  aic.dlnm<-c(aic.dlnm, AIC(fit.dlnm))
}


par(mfrow=c(1,3))


#glm
length(coef.dlnm[,1])
length(data$ac)


plot(x=data$ac,y=smooth.spline(coef.glm[1:length(data$ac)],df=700)$y,type="l",ylim=c(-0.001,0.015),ylab="beta",xlab="auto-correlation index",main="GLM Performance")
lines(x=data$ac,y=smooth.spline(CI.glm[1:length(data$ac),1],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.glm[1:length(data$ac),2],df=700)$y,lwd=0.5,lty=2)
abline(a=0.005,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta0", "95% CI for beta0", "true beta0", "no effect beta0=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)


#gam
plot(x=data$ac,y=smooth.spline(coef.gam[1:length(data$ac)],df=700)$y,type="l",ylim=c(-0.001,0.015),ylab="beta",xlab="auto-correlation index",main="GAM Performance")
lines(x=data$ac,y=smooth.spline(CI.gam[1:length(data$ac),1],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.gam[1:length(data$ac),2],df=700)$y,lwd=0.5,lty=2)
abline(a=0.005,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta0", "95% CI for beta0", "true beta0", "no effect beta0=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)

par(mfrow=c(1,3))
#dlnm
plot(x=data$ac,y=smooth.spline(coef.dlnm[1:length(data$ac),1],df=700)$y,type="l",ylim=c(-0.001,0.01),ylab="beta for lag0",xlab="auto-correlation index",main="DLNM Performance")
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),1],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),2],df=700)$y,lwd=0.5,lty=2)
abline(a=0.005,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta1", "95% CI for beta0", "true beta0", "no effect beta0=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)


plot(x=data$ac,y=smooth.spline(coef.dlnm[1:length(data$ac),2],df=700)$y,type="l",ylim=c(-0.001,0.01),ylab="beta for lag1",xlab="auto-correlation index",main="DLNM Performance")
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),3],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),4],df=700)$y,lwd=0.5,lty=2)
abline(a=0.001,b=0,lty=1,col="red")
abline(0,0,col="blue")
legend("topright", c("Estimated beta1", "95% CI for beta1", "true beta1", "no effect beta1=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)

plot(x=data$ac,y=smooth.spline(coef.dlnm[1:length(data$ac),3],df=700)$y,type="l",ylim=c(-0.001,0.01),ylab="beta for lag2",xlab="auto-correlation index",main="DLNM Performance")
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),5],df=700)$y,lwd=0.5,lty=2)
lines(x=data$ac,y=smooth.spline(CI.dlnm[1:length(data$ac),6],df=700)$y,lwd=0.5,lty=2)
abline(a=0.002,b=0,lty=1,col="red")
abline(0,0,lty=1,col="blue")
legend("topright", c("Estimated beta2", "95% CI for beta2", "true beta2", "no effect beta2=0"),
       col=c("black","black","red","blue"),lty=c(1,2,1,1),cex=0.8)

par(mfrow=c(1,2))
##variance compare

plot(x=data$ac,y=smooth.spline(sd.glm[1:length(data$ac)],df=700)$y,type="l",ylim=c(0,0.001),ylab="sd(beta)", xlab="auto-correlation index",main="Standard Error of Beta's")
lines(x=data$ac,y=smooth.spline(sd.gam[1:length(data$ac)],df=700)$y,col="red",lty=2)
lines(x=data$ac,y=smooth.spline(sd.dlnm[1:length(data$ac),1],df=700)$y,col="green",lty=1)
lines(x=data$ac, y=smooth.spline(sd.dlnm[1:length(data$ac),2],df=700)$y, col="cyan",lty=2)
lines(x=data$ac, y=smooth.spline(sd.dlnm[1:length(data$ac),3],df=700)$y, col="dark green",lty=2)

legend("topright",c("SE(beta) of GLM", "SE(beta) of GAM", "SE(beta0) of DLNM", "SE(beta1) of DLNM", "SE(beta2) of DLNM"), 
       col=c("black","red","green","cyan","dark green"), lty=c(1,2,1,2,2),adj=0,cex=0.7)

##AIC compare
plot(x=data$ac,y=smooth.spline(aic.glm[1:length(data$ac)],df=700)$y,type="l",ylim=c(5700,6000),ylab="AIC",xlab="auto-correlation index", main="AIC of Each model")
lines(x=data$ac,y=smooth.spline(aic.gam[1:length(data$ac)],df=700)$y,col="red",lty=2)
lines(x=data$ac,y=smooth.spline(aic.dlnm[1:length(data$ac)],df=700)$y,col="green",lty=2)
legend("topright",c("AIC of GLM", "AIC of GAM", "AIC of DLNM"),lty=c(1,2,2), col=c("black","red","green"),cex=0.7)

##Summary Statistics
result<-as.data.frame(cbind(data$ac, coef.glm, CI.glm, coef.gam, CI.gam, coef.dlnm, CI.dlnm))
colnames(result)<-c("ac", "beta.glm", "LCB.glm", "UCB.glm", 'beta.gam', 'LCB.gam','UCB.gam',"beta0.dlnm","beta1.dlnm","beta2.dlnm","LCB0.dlnm","UCB0.dlnm","LCB1.dlnm","UCB1.dlnm","LCB2.dlnm","UCB2.dlnm")
rownames(result)<-seq(1,7000,by=1)
attach(result)

head(result)


difference<-matrix(c(mean(beta.glm[result$ac<0.4]-0.005),
                     mean(beta.gam[result$ac<0.4]-0.005),
                     mean(beta0.dlnm[result$ac<0.4]-0.005),
                     mean(beta1.dlnm[result$ac<0.4]-0.001),
                     mean(beta2.dlnm[result$ac<0.4]-0.002),
                     mean(beta.glm[result$ac>0.4 & result$ac<0.6]-0.005),
                     mean(beta.gam[result$ac>0.4 & result$ac<0.6]-0.005),
                     mean(beta0.dlnm[result$ac>0.4 & result$ac<0.6]-0.005),
                     mean(beta1.dlnm[result$ac>0.4 & result$ac<0.6]-0.001),
                     mean(beta2.dlnm[result$ac>0.4 & result$ac<0.6]-0.002),
                     mean(beta.glm[result$ac>0.6]-0.005),
                     mean(beta.gam[result$ac>0.6]-0.005),
                     mean(beta0.dlnm[result$ac>0.6]-0.005),
                     mean(beta1.dlnm[result$ac>0.6]-0.001),
                     mean(beta2.dlnm[result$ac>0.6]-0.002)),nrow=5,ncol=3)


colnames(difference)<-c("low","mid","high")
rownames(difference)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")



coverage<-matrix(c(mean(LCB.glm[result$ac<0.4]<0.005 & UCB.glm[result$ac<0.4]>0.005),
                   mean(LCB.gam[result$ac<0.4]<0.005 & UCB.gam[result$ac<0.4]>0.005),
                   mean(LCB0.dlnm[result$ac<0.4]<0.005 & UCB0.dlnm[result$ac<0.4]>0.005),
                   mean(LCB1.dlnm[result$ac<0.4]<0.001 & UCB1.dlnm[result$ac<0.4]>0.001),
                   mean(LCB2.dlnm[result$ac<0.4]<0.002 & UCB2.dlnm[result$ac<0.4]>0.002),
                   mean(LCB.glm[result$ac>0.4 & result$ac<0.6]<0.005 & UCB.glm[result$ac>0.4 & result$ac<0.6]>0.005),
                   mean(LCB.gam[result$ac>0.4 & result$ac<0.6]<0.005 & UCB.gam[result$ac>0.4 & result$ac<0.6]>0.005),
                   mean(LCB0.dlnm[result$ac>0.4 & result$ac<0.6]<0.005 & UCB0.dlnm[result$ac>0.4 & result$ac<0.6]>0.005),
                   mean(LCB1.dlnm[result$ac>0.4 & result$ac<0.6]<0.001 & UCB1.dlnm[result$ac>0.4 & result$ac<0.6]>0.001),
                   mean(LCB2.dlnm[result$ac>0.4 & result$ac<0.6]<0.002 & UCB2.dlnm[result$ac>0.4 & result$ac<0.6]>0.002),
                   mean(LCB.glm[result$ac>0.6]<0.005 & UCB.glm[result$ac>0.6]>0.005),
                   mean(LCB.gam[result$ac>0.6]<0.005 & UCB.gam[result$ac>0.6]>0.005),
                   mean(LCB0.dlnm[result$ac>0.6]<0.005 & UCB0.dlnm[result$ac>0.6]>0.005),
                   mean(LCB1.dlnm[result$ac>0.6]<0.001 & UCB1.dlnm[result$ac>0.6]>0.001),
                   mean(LCB2.dlnm[result$ac>0.6]<0.002 & UCB2.dlnm[result$ac>0.6]>0.002)),nrow=5,ncol=3,byrow=FALSE)

colnames(coverage)<-c("low","mid","high")
rownames(coverage)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")


temp<-matrix(c(mean(LCB.glm[result$ac<0.4]<0 & UCB.glm[result$ac<0.4]>0),
               mean(LCB.gam[result$ac<0.4]<0 & UCB.gam[result$ac<0.4]>0),
               mean(LCB0.dlnm[result$ac<0.4]<0 & UCB0.dlnm[result$ac<0.4]>0),
               mean(LCB1.dlnm[result$ac<0.4]<0 & UCB1.dlnm[result$ac<0.4]>0),
               mean(LCB2.dlnm[result$ac<0.4]<0 & UCB2.dlnm[result$ac<0.4]>0),
               mean(LCB.glm[result$ac>0.4 & result$ac<0.6]<0 & UCB.glm[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB.gam[result$ac>0.4 & result$ac<0.6]<0 & UCB.gam[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB0.dlnm[result$ac>0.4 & result$ac<0.6]<0 & UCB0.dlnm[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB1.dlnm[result$ac>0.4 & result$ac<0.6]<0 & UCB1.dlnm[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB2.dlnm[result$ac>0.4 & result$ac<0.6]<0 & UCB2.dlnm[result$ac>0.4 & result$ac<0.6]>0),
               mean(LCB.glm[result$ac>0.6]<0 & UCB.glm[result$ac>0.6]>0),
               mean(LCB.gam[result$ac>0.6]<0 & UCB.gam[result$ac>0.6]>0),
               mean(LCB0.dlnm[result$ac>0.6]<0 & UCB0.dlnm[result$ac>0.6]>0),
               mean(LCB1.dlnm[result$ac>0.6]<0 & UCB1.dlnm[result$ac>0.6]>0),
               mean(LCB2.dlnm[result$ac>0.6]<0 & UCB2.dlnm[result$ac>0.6]>0)),nrow=5,ncol=3,byrow=FALSE)

significance<-1-temp

colnames(significance)<-c("low","mid","high")
rownames(significance)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")


beta<-matrix(c(mean(beta.glm[result$ac<0.4]),mean(beta.gam[result$ac<0.4]),
               mean(beta0.dlnm[result$ac<0.4]),mean(beta1.dlnm[result$ac<0.4]),
               mean(beta2.dlnm[result$ac<0.4]),
               mean(beta.glm[result$ac>0.4 & result$ac<0.6]),mean(beta.gam[result$ac>0.4 & result$ac<0.6]),
               mean(beta0.dlnm[result$ac>0.4 & result$ac<0.6]),mean(beta1.dlnm[result$ac>0.4 & result$ac<0.6]),
               mean(beta2.dlnm[result$ac>0.4 & result$ac<0.6]),
               mean(beta.glm[result$ac>0.6]),mean(beta.gam[result$ac>0.6]),
               mean(beta0.dlnm[result$ac>0.6]),mean(beta1.dlnm[result$ac>0.6]),
               mean(beta2.dlnm[result$ac>0.6]),0.005,0.005,0.005,0.001,0.002),nrow=5,ncol=4)
colnames(beta)<-c("low","mid","high","TRUE beta")
rownames(beta)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")

SD<-matrix(c(mean(sd.glm[result$ac<0.4]),mean(sd.gam[result$ac<0.4]),
             mean(sd0.dlnm[result$ac<0.4]),mean(sd0.dlnm[result$ac<0.4]),
             mean(sd2.dlnm[result$ac<0.4]),
             mean(sd.glm[result$ac>0.4 & result$ac<0.6]),mean(sd.gam[result$ac>0.4 & result$ac<0.6]),
             mean(sd0.dlnm[result$ac>0.4 & result$ac<0.6]),mean(sd0.dlnm[result$ac>0.4 & result$ac<0.6]),
             mean(sd2.dlnm[result$ac>0.4 & result$ac<0.6]),
             mean(sd.glm[result$ac>0.6]),mean(sd.gam[result$ac>0.6]),
             mean(sd0.dlnm[result$ac>0.6]),mean(sd0.dlnm[result$ac>0.6]),
             mean(sd2.dlnm[result$ac>0.6])),nrow=5,ncol=3)

colnames(SD)<-c("low","mid","high")
rownames(SD)<-c("glm","gam","dlnm_0","dlnm_1","dlnm_2")


beta
SD
difference
coverage
significance




####================================================================================
####==============================Cross Validation in Chicago Data==================


data.cv<-chicagoNMMAPS[1:4000,]
summary(chicagoNMMAPS)
data.cv$pm10<-na.interp(data.cv$pm10)


par(mfrow=c(4,1))
plot(data.cv$time, data.cv$resp,type="l",xlab="",ylab="death by respiratory disease ",main="Time Seires Plot")
##remove extreme value?? is it caused by air pollutent??
plot(data.cv$time, data.cv$o3,type="l",xlab="",ylab="O3 concentration")
plot(data.cv$time, data.cv$temp,type="l",xlab="",ylab="temperature")
plot(data.cv$time, data.cv$rhum,type="l",xlab="time",ylab="humidity")

par(mfrow=c(1,2))
hist(data.cv$resp,main="Histogram of Respiratory Death Data", xlab="count")
acf(data.cv$o3, main="Series O3")
mean(abs(acf(data.cv$o3)$acf[2:6]))


sd(data.cv$resp)^2
mean(data.cv$resp)

data.cv$o3<-data.cv$o3/10
###models
fit1<-glm(resp~o3+ns(time, df=7*11)+ns(rhum, df=6)+ns(temp, df=7)+dow, family = poisson(),data=data.cv)
fit2<-gam(resp~o3+s(time,k=7*11)+s(rhum,k=6)+s(temp, k=7)+dow, family= poisson(), data=data.cv)


cb1.o3 <- crossbasis(data.cv$o3, lag=15, argvar=list(fun="lin"), arglag=list(fun="poly", degree=4))
fit3 <- gam(resp~cb1.o3+s(time,k=7*11)+s(rhum,k=6)+s(temp, k=7)+dow, family= poisson(), data=data.cv)
pred3 <- crosspred(cb1.o3, fit3, at=1, bylag=1, cumul=TRUE)

pred3$matfit
pred3$matse
log(c(pred3$matRRlow[1],pred3$matRRhigh[1]))


cb2.o3 <- crossbasis(data.cv$o3, lag=10, argvar=list(fun="lin"), arglag=list(fun="poly", degree=4))
fit4 <- gam(resp~cb2.o3+s(time,k=7*11)+s(rhum,k=6)+s(temp, k=7)+dow, family= poisson(), data=data.cv)
pred4 <- crosspred(cb2.o3, fit4, at=1, bylag=1, cumul=TRUE)


pred4$matfit
pred4$matse
log(c(pred4$matRRlow[1],pred4$matRRhigh[1]))



plot(pred3, "slices", var=1, col=3, ylab="RR", ci.arg=list(density=15, lwd=2),main="Lag-response curve for 
                                                                                    10-unit increase in o3")
plot(pred3, "slices", var=1, col=2, ylab="Cumulative RR", ci.arg=list(density=15, lwd=2),main="Lag-response curve of 
                                                                              incremental cumulative effects")

plot(pred4, "slices", var=1, col=3, ylab="RR", ci.arg=list(density=15, lwd=2),main="Lag-response curve for 
                                                                                    10-unit increase in o3")
plot(pred4, "slices", var=1, col=2, ylab="Cumulative RR", ci.arg=list(density=15, lwd=2),main="Lag-response curve of 
     incremental cumulative effects")





summary(fit1)
confint(fit1,"o3",level=0.95)

summary(fit2)

c(-0.006231-1.96*0.009080, -0.006231+1.96*0.009080)






###PM10
fit1<-glm(resp~pm10+ns(time, df=7*11)+ns(rhum, df=6)+ns(temp, df=7)+dow, family = poisson(),data=data.cv)
fit2<-gam(resp~pm10+s(time,k=7*11)+s(rhum,k=6)+s(temp, k=7)+dow, family= poisson(), data=data.cv)

summary(fit1)
confint(fit1,"pm10",level=0.95)

s<-summary(fit2)
s$se
c(summary(fit2)$p.coef[2]-1.96*summary(fit2)$se[2], summary(fit2)$p.coef[2]+1.96*summary(fit2)$se[2])

cb1.pm10 <- crossbasis(pm10, lag=15, argvar=list(fun="lin"), arglag=list(fun="poly", degree=4))
fit3 <- gam(resp~cb1.pm10+s(time,k=7*11)+s(rhum,k=6)+s(temp, k=7)+dow, family= poisson(), data=data.cv)
pred3 <- crosspred(cb1.pm10, fit3, at=1, bylag=1, cumul=TRUE)

pred3$matfit
pred3$matse
log(c(pred3$matRRlow[1],pred3$matRRhigh[1]))


cb2.pm10 <- crossbasis(data.cv$pm10, lag=10, argvar=list(fun="lin"), arglag=list(fun="poly", degree=4))
fit4 <- gam(resp~cb2.pm10+s(time,k=7*11)+s(rhum,k=6)+s(temp, k=7)+dow, family= poisson(), data=data.cv)
pred4 <- crosspred(cb2.pm10, fit4, at=1, bylag=1, cumul=TRUE)


pred4$matfit
pred4$matse
log(c(pred4$matRRlow[1],pred4$matRRhigh[1]))




