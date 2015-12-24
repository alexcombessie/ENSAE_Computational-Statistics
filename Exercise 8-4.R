#Exercise 8.4

TargetDensity<-function(x,d)  exp(-x^d)

set.seed(42)

SliceSampler_1step<-function(x_t, d)    {
  u<-runif(1, min=0, max=TargetDensity(x_t,d))
  x_t_plus1<-runif(1, min=0, max=(-log(u))^(1/d))
  return(x_t_plus1)}

SliceSampler_full<-function(nb_simulations, starting_point, d)    {
  X<-rep(NA,nb_simulations)
  X[1]<-starting_point
  for(t in 2:nb_simulations) X[t]<-SliceSampler_1step(X[t-1], d)
  return(X)}

d_list<-c(0.15,0.25,0.4)
par(mfrow = c(3,1))
for (d in d_list) {
  S<-SliceSampler_full(1000000, starting_point = 0, d)
  hist(S, prob=T, col='grey', breaks="Scott", xlim=c(0,0.1*max(S)),
       main=paste("Slice Sampler histograms for d=",toString(d)),xlab="")
  x_range<-seq(0,0.1*max(S),length=1000)
  normalizing_constant<-1/integrate(function(x) TargetDensity(x,d), 0, Inf)$value
  lines(x_range,TargetDensity(x_range,d)*normalizing_constant, add=TRUE, col="red", lwd=2)}

N<-seq(2,10000, length=1000)
d_list<-c(0.1,0.25,0.4)
par(mfrow = c(3,1))
for (d in d_list) {
  Y<-sapply(N,FUN=function(n){mean(SliceSampler_full(n,d=d,0))})
  plot(N,Y,type="n",main=paste("Slice Sampler convergence for d=",toString(d)),
       xlab="Number of simulations",ylab="Mean")
  lines(N,Y,col="red",lwd=1)}