f1 = function (theta,x) {
  theta/(1+theta^2)*exp(-(x-theta)^2)/2 }
f2 = function (theta,x) {
  1/(1+theta^2)*exp(-(x-theta)^2)/2 }

TH<- seq(-8, 8, length=1000)
Y1 <-f1(theta=TH,x=1)
Y2 <-f2(theta=TH,x=1)

par(mfrow = c(1, 1))
plot(range(TH),range(c(Y1,Y2)), type="n",ylab="Integrand value",xlab="Theta",
     main="Integrand functions of theta for x=1")
lines(TH,Y1,col="red",lwd=2.5)
lines(TH,Y2,col="blue",lwd=2.5,lty=2)
legend("topright",c("f1: Numerator integrand","f2: Denominator integrand"), 
       lty=c(1,2), col=c("red","blue"),lwd=c(2.5,2.5))


TH<- seq(-12, 12, length=1000)
par(mfrow = c(3, 3))
for (i in 0:7) {
  Y1 <-f1(theta=TH,x=2*i-6)
  Y2 <-f2(theta=TH,x=2*i-6)
  plot(range(TH),range(c(Y1,Y2)), type="n",ylab="Integrand",xlab="Theta",
       main=paste("Integrands of theta for x=",toString(2*i-6)))
  lines(TH,Y1,col="red",lwd=2.5)
  lines(TH,Y2,col="blue",lwd=2.5,lty=2)
  abline(v=2*i-6,col="black")
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,c("f1","f2"), lty=c(1,2), col=c("red","blue"),lwd=c(2.5,2.5),horiz=TRUE)


install.packages("plot3D")
library(plot3D)

l<-500
TH<-seq(-4, 4, length=l)
X<-seq(-2,2, length=l)
z <- outer(TH, X, f1)
par(mfrow = c(1,1))
persp3D(x=TH,y=X,z=z,xlab="Theta",ylab="X",zlab="",clab = "Numerator integrand",
        phi=30,theta=-30,ticktype="detailed",axes=TRUE,expand=0.6)
par(mfrow = c(1,1))
X<-seq(-1,1, length=l)
z <- outer(TH, X, f1)
persp3D(x=TH,y=X,z=z,xlab="Theta",ylab="X",zlab="",clab = "Numerator integrand",
        phi=30,theta=-30,ticktype="detailed",axes=TRUE,expand=0.6)


boxmuller=function(n)  {
  u1<-runif(n)
  u2<-runif(n)
  x<-c()
  for (i in 1:n)  {x<-c(x,sqrt(-2*log(u1[i]))*cos(2*pi*u2[i]))}
  return(x)}

delta<-function(x,n){
  theta1<-x+boxmuller(n)
  theta2<-x+boxmuller(n)
  f1_numerator<-mean(sqrt(2*pi)*theta1/(1+theta1^2))
  f2_denominator<-mean(sqrt(2*pi)*1/(1+theta2^2))
  return(f1_numerator/f2_denominator)}

X<-seq(-10,10, length=100)
Y<-sapply(X,FUN=function(x){delta(x,1000)})
par(mfrow = c(1, 1))
plot(X,Y,type="n",main="Monte Carlo estimation of delta(x) for 1000 simulations",
     xlab="x",ylab="delta(x)")
lines(X,Y,col="grey",lwd=5)

f1_n<-function(x,n){
  theta<-x+boxmuller(n)
  return(mean(sqrt(2*pi)*theta/(1+theta^2)))}
f2_n<-function(x,n){
  theta<-x+boxmuller(n)
  return(mean(sqrt(2*pi)*1/(1+theta^2)))}
sigma1_n<-function(x,n){
  theta<-x+boxmuller(n)
  return(var(sqrt(2*pi)*theta/(1+theta^2)))}
sigma2_n<-function(x,n){
  theta<-x+boxmuller(n)
  return(var(sqrt(2*pi)*1/(1+theta^2)))}
omega_n<-function(x,n)  {
  return(sigma1_n(x,n)/(f2_n(x,n)^2)+sigma2_n(x,n)*f1_n(x,n)^2/(f2_n(x,n)^4))}
nlim<-function(x,n,digits){
  return(floor(16*(qchisq(.95, df=1)^2)*omega_n(x,n)*10^(2*digits)))}

N<-seq(0,10000, length=100)
Y<-sapply(N,FUN=function(n){omega_n(x=1,n)/n})
par(mfrow = c(1, 1))
plot(N,Y,type="n",main="Asymptotic variance estimation for x=1",
     xlab="Number of simulations",ylab="omega_n(x=1,n)")
lines(N,Y,col="red",lwd=3)

X<-seq(-10,10,length=10)
mean(sapply(X,FUN=function(x){nlim(x,n=10^5,digits=3)}))
