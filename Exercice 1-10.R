#Creates the Student density function
student = function (x,p,sigma,theta) {
  1/sigma*(1+((x-theta)^2/(p*sigma^2)))^(-(p+1)/2)  }

#Particular case where p=1 and sigma=1 (cauchy)
case_cauchy=function(x,theta)  {
  student(x,p = 1,sigma = 1,theta)  }

par(mfrow = c(1, 1))
X <- seq(-5, 5, length=1000)
Y <-case_cauchy(x=X,theta=1)
plot(X,Y, lwd=2, col="blue",type="n",ylab="Density",main="Cauchy density function for theta=1")
lines(X,Y,col="black")

#Vector of observed variables
observed<-c(0,5,9)

#Likelihood of observations
likelihood=function(theta, obs){
  l<-1
  for (i in obs)  {
    l<-l*case_cauchy(x=i,theta)
  }
  return(l)
}
par(mfrow = c(1, 1))
X <- seq(-5, 15, length=10000)
Y <-likelihood(X,observed)
plot(X,Y, lwd=2, col="blue",type="n",xlab="Theta", 
     ylab="Likelihood",main="Likelihood of theta for the 3 given observations")
lines(X,Y,col="blue")
abline(v=observed)

#What happens with a fourth observation which varies?

X <- seq(-5, 15, length=10000)
par(mfrow = c(4, 3))
for (i in 1:12) {
  Y <-likelihood(X,c(observed,i))
  plot(X,Y, lwd=2, col="blue",type="n",xlab="Theta", ylab="Likelihood",
       main=paste("Likelihood of theta with X4=",toString(i)))
  lines(X,Y,col="blue")
  abline(v=i,col="red",)
}