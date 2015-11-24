standard_boxmuller=function(n)  {
  u1<-runif(n)
  u2<-runif(n)
  x1<-c()
  x2<-c()
  for (i in 1:n){
    x1<-c(x1,sqrt(-2*log(u1[i]))*cos(2*pi*u2[i]))
    x2<-c(x2,sqrt(-2*log(u1[i]))*sin(2*pi*u2[i]))
  }
  return(x1)
}

alternative_boxmuller=function(n)  {
  x1<-c()
  x2<-c()
  counter<-1
  while(counter<=n){
    u1<-runif(1,min=-1,max=1)
    u2<-runif(1,min=-1,max=1)
    s<-u1^2+u2^2
    if (s<=1) {
      counter<-counter+1
      x1<-c(x1,sqrt(-2*log(s)/s)*u1)
      x2<-c(x2,sqrt(-2*log(s)/s)*u2)
    }
  }
  return(x1)
}

n<-20000
time_standard<-c()
time_alternative<-c()
for (i in seq.int(1,n,n/10))  {
  time_standard<-c(time_standard,system.time(standard_boxmuller(i))[3])
  time_alternative<-c(time_alternative,system.time(alternative_boxmuller(i))[3])
}
plot(c(0,n),range(time_alternative),xlab="Number of generations",ylab="Execution time")
lines(seq.int(1,n,n/10),time_standard, col="red")
lines(seq.int(1,n,n/10),time_alternative, col="blue")
legend("topleft",c("Standard algorithm","Alternative algorithm"), 
       lty=c(1,1), col=c("red","blue"))
