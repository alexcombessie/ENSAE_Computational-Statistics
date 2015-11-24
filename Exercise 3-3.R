p_n=function(n)  {
  normal<-rnorm(n)
  ind<-normal>2.5
  return(1/n*sum(ind))}
p_n(10^6)
1-pnorm(2.5)
16*(1-pnorm(2.5))*pnorm(2.5)*10^6

q_n=function(n)  {
  gamma<-rgamma(n,1,1)
  ind<-gamma>5.3
  return(1/n*sum(ind))}
q_n(10^6)

floor(2*qchisq(0.995,df=1)*10^3-qchisq(0.995,df=1))