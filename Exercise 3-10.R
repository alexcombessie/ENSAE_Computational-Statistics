boxmuller=function(n)  {
  u1<-runif(n)
  u2<-runif(n)
  x<-c()
  for (i in 1:n)  {
    x<-c(x,sqrt(-2*log(u1[i]))*cos(2*pi*u2[i]))}
  return(x)}

p1=function(n)  {
  normal<-boxmuller(n)
  ind<-normal>=1&normal<=2
  return(1/n*sum(ind))}

h=function(x) {
  return(1/sqrt(2*pi)*exp(-x^2/2))}

p2=function(n)  {
  u<-runif(n,min=1,max=2)
  h_u<-sapply(u,FUN=h)
  return(1/n*sum(h_u))}

h_inv=function(x) {
  return(sqrt(2*pi)*exp(x^2/2))}

p3=function(n)  {
  u<-runif(n,min=1,max=2)
  h_inv_u<-sapply(u,FUN=h_inv)
  return(n/sum(h_inv_u))}


n<-100000
p1_n<-c()
p2_n<-c()
p3_n<-c()
n_range<-seq.int(1,n,n/10)
for (i in n_range)  {
  p1_n <-c(p1_n,p1(i))
  p2_n <-c(p2_n,p2(i))
  p3_n <-c(p3_n,p3(i))}
par(mfrow = c(1, 1))
plot(range(n_range),c(0.1,0.18),type="n",xlab="Number of simulations", ylab="Value of the estimator",
     main="Convergence of the estimators")
lines(n_range,p1_n, col="red",lwd=2)
lines(n_range,p2_n, col="blue", lty=2,lwd=2)
lines(n_range,p3_n, col="green",lty=3,lwd=2)
abline(h=pnorm(2)-pnorm(1),col="black")
legend("topright",c("p1","p2","p3"), 
       lty=c(1,2,3), col=c("red","blue","green"),lwd=c(2,2,2))
text(x =85000 ,y=pnorm(2)-pnorm(1),pos=1,"Real value")


step_precision=function(p,nb_decimals)  {
  target<-pnorm(2)-pnorm(1)
  i<-1000
  while(abs(p(i)-target)>1/2*10^-nb_decimals)
    i<-i+100
  return(i)}
p1_step_precision4<-mean(replicate(100, step_precision(p1,4)))
p2_step_precision4<-mean(replicate(100, step_precision(p2,4)))


