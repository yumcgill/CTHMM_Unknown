library(expm)
library(abind)
library(parallel)
options(mc.cores = 20)
set.seed(2332)
state<-4


ini<-c(0.35,0.25,0.2,0.2)

Q<-matrix(c(0,2,1,0,1,0,0.75,0.05,0.15,0.55,0,0.35,0.0,0.25,0.4,0),byrow=T,nrow=state)
r<-rowSums(Q)
tranQ<-Q/r
diag(Q)<--r
be<-matrix(c(-1.28,-0.88,0.70,-0.55,1.15,0.68,-1.05,1.36,-1.12,0.99,1.73,-1.20),byrow=F,ncol=state)



simdata<-NULL
N.sim<-1000
T.obs<-15
sigma<-1
for(m in 1:N.sim){
  x<-NULL
  y<-NULL
  n.obs<-round(runif(1,20,60))
  z1<-rbinom(n.obs,1,0.6)
  z2<-rnorm(n.obs,-1,1)
  x[1]<-sample(1:state,1,prob=ini)
  y[1]<-rnorm(1,cbind(1,z1[1],z2[1])%*%be[,x[1]],sigma)
  delt<-NULL
  i = 2
  t<-0
  while(max(t)<=T.obs){
    
    #tran<-NULL
    x[i]<-sample(1:state,1,prob=tranQ[x[i-1],])
    delt[i-1]<-rexp(1,rate=r[x[i-1]])
    #tran<-expm(Q*delt[i-1])
    #x[i]<-sample(1:state,1,prob=tran[x[i-1],])
    i<-i+1
    t<-c(t,sum(delt))
  }
  t[length(t)]<-T.obs
  tim.g<-c(0,sort(runif(n.obs-1,0,T.obs)))
  
  del<-0
  xh<-rep(0,n.obs)
  xh[1]<-x[1]
  for (j in 2:n.obs){
    ind<-max(which(t<=tim.g[j]))
    xh[j]<-x[ind]
    y[j]<-rnorm(1,cbind(1,z1[j],z2[j])%*%be[,xh[j]],sigma)
    del<-c(del,tim.g[j]-tim.g[j-1])
  }
  
  
  id<-rep(m,n.obs)
  sub<-cbind(id,y,z1,z2,tim.g,del,xh)
  simdata<-as.data.frame(rbind(simdata,sub))
}
simdata$int<-1
