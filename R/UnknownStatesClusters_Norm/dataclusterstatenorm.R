set.seed(22221)
#state<-3
ini1<-c(0.5,0.4,0.1)
Q1<-matrix(c(-2.5,2,0.5,0.5,-1.5,1,0.1,0.9,-1),byrow=T,nrow=3)
r1<-c(2.5,1.5,1)
tranQ1<-matrix(c(0,2/2.5,0.5/2.5,0.5/1.5,0,1/1.5,0.1,0.9,0),byrow=T,nrow=3)
be1<-c(-3,0,2)

ini2<-c(0.6,0.4)
Q2<-matrix(c(-1.2,1.2,0.25,-0.25),byrow=T,nrow=2)
r2<-c(1.2,0.25)
tranQ2<-matrix(c(0,1,1,0),byrow=T,nrow=2)
be2<-c(-3.7,3.9)


ini3<-c(0.45,0.45,0.1)
Q3<-matrix(c(-0.5,0.49,0.01,0.25,-0.3,0.05,0.01,0.1,-0.11),byrow=T,nrow=3)
r3<-c(0.5,0.3,0.11)
tranQ3<-matrix(c(0,0.49/0.5,0.01/0.5,0.25/0.3,0,0.05/0.3,0.01/0.11,0.1/0.11,0),byrow=T,nrow=3)
be3<-c(-3.8,1,4)

ini4<-c(0.35,0.25,0.2,0.2)
Q4<-matrix(c(0,2,1,0,1,0,0.75,0.05,0.15,0.55,0,0.35,0.0,0.25,0.4,0),byrow=T,nrow=4)
r4<-rowSums(Q4)
tranQ4<-Q4/r4
diag(Q4)<--r4
be4<-c(-4,-1.2,1.7,3.8)

simdata<-NULL

T.obs<-15
N.sim1<-400
N.sim2<-500
N.sim3<-450
N.sim4<-550
sigma<-1

state1<-nrow(Q1)
state2<-nrow(Q2)
state3<-nrow(Q3)
state4<-nrow(Q4)
for(m in 1:N.sim1){
  x<-NULL
  y<-NULL
  n.obs<-round(runif(1,20,60))
 
  z<-rnorm(n.obs,1,1)
  x[1]<-sample(1:state1,1,prob=ini1)
  
  y[1]<-rnorm(1,be1[x[1]],sigma)
  
  delt<-NULL
  i = 2
  t<-0
  while(max(t)<=T.obs){
    
    #tran<-NULL
    x[i]<-sample(1:state1,1,prob=tranQ1[x[i-1],])
    delt[i-1]<-rexp(1,rate=r1[x[i-1]])
    
    i<-i+1
    t<-c(t,sum(delt))
  }
  
  
  tim.g<-c(0,sort(runif(n.obs-1,0,T.obs)))
  
  del<-NULL
  del[1]<-0
  xh<-rep(0,n.obs)
  xh[1]<-x[1]
  #th<-0
  #xh<-x[1]
  for (j in 2:n.obs){
    ind<-max(which(t<=tim.g[j]))
    xh[j]<-x[ind]
    #th<-c(th,t[ind])
    y[j]<-rnorm(1,be1[xh[j]],sigma)
    #y[j]<-rgamma(1,exp(cbind(1,z[j])%*%be1[,xh[j]]),sigma)
    del[j]<-tim.g[j]-tim.g[j-1]
  }
  
  id<-rep(m,n.obs)
  sub<-cbind(id,y,z,tim.g,del,xh)
  simdata<-as.data.frame(rbind(simdata,sub))
  
} 


for(m in 1:N.sim2){
  x<-NULL
  y<-NULL
  n.obs<-round(runif(1,20,60))
  
  z<-rnorm(n.obs,1,1)
  x[1]<-sample(1:state2,1,prob=ini2)
  
  y[1]<-rnorm(1,be2[x[1]],sigma)
  
  delt<-NULL
  i = 2
  t<-0
  while(max(t)<=T.obs){
    
    
    x[i]<-sample(1:state2,1,prob=tranQ2[x[i-1],])
    delt[i-1]<-rexp(1,rate=r2[x[i-1]])
   
    i<-i+1
    t<-c(t,sum(delt))
  }
  
  
  tim.g<-c(0,sort(runif(n.obs-1,0,T.obs)))
  
  del<-NULL
  del[1]<-0
  xh<-rep(0,n.obs)
  xh[1]<-x[1]
  
  for (j in 2:n.obs){
    ind<-max(which(t<=tim.g[j]))
    xh[j]<-x[ind]
    
    y[j]<-rnorm(1,be2[xh[j]],sigma)
   
    del[j]<-tim.g[j]-tim.g[j-1]
  }
  xh<-xh+state1
  id<-rep(m+N.sim1,n.obs)
  sub<-cbind(id,y,z,tim.g,del,xh)
  simdata<-as.data.frame(rbind(simdata,sub))
  
} 

for(m in 1:N.sim3){
  x<-NULL
  y<-NULL
  n.obs<-round(runif(1,20,60))
  
  z<-rnorm(n.obs,1,1)
  x[1]<-sample(1:state3,1,prob=ini3)
  
  y[1]<-rnorm(1,be3[x[1]],sigma)
  
  delt<-NULL
  i = 2
  t<-0
  while(max(t)<=T.obs){
    
    
    x[i]<-sample(1:state3,1,prob=tranQ3[x[i-1],])
    delt[i-1]<-rexp(1,rate=r3[x[i-1]])
    
    i<-i+1
    t<-c(t,sum(delt))
  }
  
  
  tim.g<-c(0,sort(runif(n.obs-1,0,T.obs)))
  
  del<-NULL
  del[1]<-0
  xh<-rep(0,n.obs)
  xh[1]<-x[1]
 
  for (j in 2:n.obs){
    ind<-max(which(t<=tim.g[j]))
    xh[j]<-x[ind]
    
    y[j]<-rnorm(1,be3[xh[j]],sigma)
   
    del[j]<-tim.g[j]-tim.g[j-1]
  }
  xh<-xh+state2+state1
  id<-rep(m+N.sim1+N.sim2,n.obs)
  sub<-cbind(id,y,z,tim.g,del,xh)
  simdata<-as.data.frame(rbind(simdata,sub))
  
} 


for(m in 1:N.sim4){
  x<-NULL
  y<-NULL
  n.obs<-round(runif(1,20,60))
  
  z<-rnorm(n.obs,1,1)
  x[1]<-sample(1:state4,1,prob=ini4)
  
  y[1]<-rnorm(1,be4[x[1]],sigma)
 
  delt<-NULL
  i = 2
  t<-0
  while(max(t)<=T.obs){
    
    
    x[i]<-sample(1:state4,1,prob=tranQ4[x[i-1],])
    delt[i-1]<-rexp(1,rate=r4[x[i-1]])
    
    i<-i+1
    t<-c(t,sum(delt))
  }
  
  
  tim.g<-c(0,sort(runif(n.obs-1,0,T.obs)))
  
  del<-NULL
  del[1]<-0
  xh<-rep(0,n.obs)
  xh[1]<-x[1]
 
  for (j in 2:n.obs){
    ind<-max(which(t<=tim.g[j]))
    xh[j]<-x[ind]
    
    y[j]<-rnorm(1,be4[xh[j]],sigma)
    
    del[j]<-tim.g[j]-tim.g[j-1]
  }
  xh<-xh+state3+state2+state1
  id<-rep(m+N.sim1+N.sim2+N.sim3,n.obs)
  sub<-cbind(id,y,z,tim.g,del,xh)
  simdata<-as.data.frame(rbind(simdata,sub))
  
} 


simdata$int<-1
