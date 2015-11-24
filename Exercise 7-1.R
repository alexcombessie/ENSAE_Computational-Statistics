# Exercise 7.1

#Calculate the mean of a Gamma(4.3,6.2) random variable using

#(a) Accept-Reject with a Gamma(4,7) candidate.

X <- seq(-2,100,length=100000)
Gamma_var1<-dgamma(X, shape = 4.3, scale = 6.2)
Gamma_var2<-dgamma(X, shape = 4, scale = 7)
par(mfrow = c(1,1))
plot(range(X), c(0,max(Gamma_var1)*1.2),type="n",xlab="",ylab="Density")
lines(X,Gamma_var1,col="blue",lty=1,lwd=1.5)
lines(X,Gamma_var2,col="red",lty=2,lwd=1.5)
legend("topright",c("Density of a Gamma(4.3,6.2)","Density of a Gamma(4,7)"), 
       lty=c(1,2), col=c("blue","red"),lwd=c(2,2))
Ratio<-Gamma_var1/Gamma_var2
M<-max(Ratio,na.rm=TRUE)


AcceptReject=function(nb_simulations, max_ratio=M){
  n<-ceiling(nb_simulations*max_ratio)
  candidate<-rgamma(n, shape = 4, scale = 7)
  uniform<-runif(n,max=max_ratio)
  result<-candidate[uniform<=dgamma(candidate, shape = 4.3, scale = 6.2)
                                /dgamma(candidate, shape = 4, scale = 7)]
  if(length(result)>=nb_simulations)  {
    return(result[1:nb_simulations])
  }
  else {
    return(c(result,AcceptReject(nb_simulations-length(result))))
  }
}
par(mfrow = c(1, 1))
hist(AcceptReject(1000000), prob=T, col='grey')
mean(AcceptReject(1000000))

real_mean_value<-4.3*6.2
real_variance_value<-4.3*6.2^2
N<-seq(0,100000, length=1000)
Y_acceptreject<-sapply(N,FUN=function(n){mean(AcceptReject(n))})
Y2_acceptreject<-sapply(N,FUN=function(n){var(AcceptReject(n))})
par(mfrow = c(2, 1))
plot(N,Y_acceptreject,type="n",main="Convergence of the AcceptReject algorithm",
     xlab="Number of simulations",ylab="Value of the estimator")
lines(N,Y_acceptreject,col="red",lwd=2)
abline(h=real_mean_value,col="black")
plot(N,Y2_acceptreject,type="n",main="Variance of the AcceptReject algorithm",
     xlab="Number of simulations",ylab="Variance of the estimator")
lines(N,Y2_acceptreject,col="blue",lwd=2)
abline(h=real_variance_value,col="black")

#(b) Metropolis–Hastings with a Gamma(4,7) candidate.

Candidate_shape<-4
Candidate_scale<-7
target_density=function(x){
  dgamma(x, shape = 4.3, scale = 6.2)}
candidate_density=function(x){
  dgamma(x, shape = Candidate_shape, scale = Candidate_scale)}
rho_probability=function(x,y){
  min(1,(target_density(y)/target_density(x))*(candidate_density(x)/candidate_density(y)))}

metropolis_hastings_1step=function(x) {
  y<-rgamma(1,shape=Candidate_shape, scale=Candidate_scale)
  rho<-rho_probability(x,y)
  bin<-rbinom(n=1, size=1, prob=rho)
  if (bin==1)  {return(y)}
  else {return(x)}}

metropolis_hastings_full=function(nb_simulations, starting_point=4){
  x<-rep(NA,nb_simulations)
  x[1]<-starting_point
  for(i in 2:nb_simulations){
    x[i]<-metropolis_hastings_1step(x[i-1])}
  return(x)}

hist(metropolis_hastings_full(100000), prob=T, col='grey')

N<-seq(2,10000, length=1000)
Y_mh1<-sapply(N,FUN=function(n){mean(metropolis_hastings_full(n))})
Y2_mh1<-sapply(N,FUN=function(n){var(metropolis_hastings_full(n))})
par(mfrow = c(2, 1))
plot(N,Y_mh1,type="n",main="Convergence of the Metropolis-Hastings algorithm - Candidate 1",
     xlab="Number of simulations",ylab="Value of the estimator")
lines(N,Y_mh1,col="red",lwd=2)
abline(h=real_mean_value,col="black")
plot(N,Y2_mh1,type="n",main="Variance of the Metropolis-Hastings algorithm - Candidate 1",
     xlab="Number of simulations",ylab="Variance of the estimator")
lines(N,Y2_mh1,col="blue",lwd=2)
abline(h=real_variance_value,col="black")

#(c) Metropolis–Hastings with a Gamma(5,6) candidate.

Candidate_shape<-5
Candidate_scale<-6

target_density=function(x){
  dgamma(x, shape = 4.3, scale = 6.2)}
candidate_density=function(x){
  dgamma(x, shape = Candidate_shape, scale = Candidate_scale)}
rho_probability=function(x,y){
  min(1,(target_density(y)/target_density(x))*(candidate_density(x)/candidate_density(y)))}

metropolis_hastings_1step=function(x) {
  y<-rgamma(1,shape=Candidate_shape, scale=Candidate_scale)
  rho<-rho_probability(x,y)
  bin<-rbinom(n=1, size=1, prob=rho)
  if (bin==1)  {return(y)}
  else {return(x)}
}

metropolis_hastings_full=function(nb_simulations, starting_point=4){
  x<-rep(NA,nb_simulations)
  x[1]<-starting_point
  for(i in 2:nb_simulations){
    x[i]<-metropolis_hastings_1step(x[i-1])}
  return(x)
}

hist(metropolis_hastings_full(100000), prob=T, col='grey')

N<-seq(2,10000, length=1000)
Y_mh2<-sapply(N,FUN=function(n){mean(metropolis_hastings_full(n))})
Y2_mh2<-sapply(N,FUN=function(n){var(metropolis_hastings_full(n))})
par(mfrow = c(2, 1))
plot(N,Y_mh2,type="n",main="Convergence of the Metropolis-Hastings algorithm - Candidate 2",
     xlab="Number of simulations",ylab="Value of the estimator")
lines(N,Y_mh2,col="red",lwd=2)
abline(h=real_mean_value,col="black")
plot(N,Y2_mh2,type="n",main="Variance of the Metropolis-Hastings algorithm - Candidate 2",
     xlab="Number of simulations",ylab="Variance of the estimator")
lines(N,Y2_mh2,col="blue",lwd=2)
abline(h=real_variance_value,col="black")

#In each case monitor the convergence

par(mfrow = c(1, 1))
plot(N,Y_mh2,type="n",main="Comparative convergence of algorithms",
     xlab="Number of simulations",ylab="Value of the estimator",
     ylim=c(25,28))
lines(N,Y_acceptreject,col="red",lwd=1)
lines(N,Y_mh1,col="green",lwd=1)
lines(N,Y_mh2,col="blue",lwd=1)
legend("topright",c("Accept-Reject with Candidate Gamma(4,7)",
                    "Metropolis-Hastings with Candidate Gamma(4,7)",
                    "Metropolis-Hastings with Candidate Gamma(5,6)"),
                    col=c("red","green","blue"),lwd=c(1,1,1))
abline(h=real_mean_value,col="black")

precision=function(algorithm,num_decimals)  {
  n<-1000
  while(abs(mean(algorithm(n))-real_mean_value)>10^(-num_decimals)){
    n<-n+100}
  return(n)}
set.seed(1)
precision(AcceptReject,3)
precision(metropolis_hastings_full,3)


