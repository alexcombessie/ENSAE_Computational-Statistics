TargetDensity<-function(x)  exp(-x^2/2)
set.seed(42)
SliceSampler_1step<-function(x_t)    {
  omega<-runif(1, min=0, max=TargetDensity(x_t))
  x_t_plus1<-runif(1, min=-sqrt(-2*log(omega)), max=sqrt(-2*log(omega)))
  return(x_t_plus1)}
SliceSampler_full<-function(nb_simulations, starting_point=0)    {
  X<-rep(NA,nb_simulations)
  X[1]<-starting_point
  for(t in 2:nb_simulations) X[t]<-SliceSampler_1step(X[t-1])
  return(X)}

quantilevalue_list<-c(0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9999)
cdfvalue_list<-qnorm(quantilevalue_list)
N<-seq(2,1000, length=500)

par(mfrow = c(3,3), xpd=FALSE)
for (i in 1:length(cdfvalue_list)){
  Y1<-sapply(N, FUN=function(n){mean(SliceSampler_full(n,0)<=cdfvalue_list[i])})
  Y2<-sapply(N, FUN=function(n){mean(rnorm(n)<=cdfvalue_list[i])})
  plot(N, Y1, type="n", xlab="Number of simulations", ylab="Quantile", ylim=c(0,1),
       main=paste("Convergence for the quantile ",toString(quantilevalue_list[i])))
  lines(N, Y1, col="red", lty=1)
  lines(N, Y2, col="blue", lty=2)
  abline(h=quantilevalue_list[i], col="black")
  legend("bottomright", c("Slice Sampler", "R iid Sampler"),
         col=c("red","blue"), lty=c(1,2), bty="n",cex=0.7)}

set.seed(42)
precision<-function(algorithm,cdfvalue,quantile,num_decimals=4,num_simulations=10)  {
  n<-rep(1000,num_simulations)
  for(i in 1:num_simulations){
    while(abs(mean(algorithm(n[i])<=cdfvalue)-quantile)>10^(-num_decimals)){
      n[i]<-n[i]+500}}
  return(mean(n))}
M<-length(cdfvalue_list)
precision4rnorm<-rep(NA,M)
precision4slice<-rep(NA,M)
for (i in 1:M){
  precision4rnorm[i]<-precision(rnorm,cdfvalue_list[i],quantilevalue_list[i])
  precision4slice[i]<-precision(SliceSampler_full,cdfvalue_list[i],quantilevalue_list[i])}

precisionmatrix<-matrix(data=c(precision4slice,mean(precision4slice),
                               precision4rnorm,mean(precision4rnorm)),
                        byrow=TRUE, nrow=2,ncol=M+1,
                        dimnames=list(c("Slice Sampler","R iid Sampler"),
                                      c(quantilevalue_list,"Mean")))
par(mfrow = c(1,1))
barplot(precisionmatrix, main="Number of simulations needed to reach 4-decimal precision",
        xlab="Number of simulations", col=c("red","blue"),
        legend = row.names(precisionmatrix), beside=TRUE)

