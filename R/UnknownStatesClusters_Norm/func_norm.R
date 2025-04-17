library(MCMCpack)
library(ECctmc)
library(mvtnorm)
library(expm)
library(abind)
library(castor)
library(parallel)


options(mc.cores = 22)

###compute the jump matrix and holding time
jump <- function(n.stat, observed_jumps,observed_time){
  N <- array(0, dim=c(n.stat,n.stat))
  R <- array(0, dim=n.stat)
  for(i in 1:(length(observed_jumps)-1)){
    i_ind <- observed_jumps[i]
    j_ind <- observed_jumps[i+1]
    R[i_ind]<-R[i_ind]+(observed_time[i+1]-observed_time[i])
    N[i_ind,j_ind] <- N[i_ind,j_ind] + 1
  }
  return(list(R=R,N=N))
}

markov.sim<-function(subdat,a,Qe){
  n.stat<-ncol(Qe)
  time.R<-array(0,c(n.stat))
  time.N<-array(0,c(n.stat,n.stat))
  tim.g<-cumsum(subdat$del)
  #for( j in 1:(dim(subdat)[1]-1)) {
  #pp<-NULL
  pp<-lapply(1:(dim(subdat)[1]-1),function(j) sample_path_mr(a=a[[j]][1],b=a[[j]][2],t0=tim.g[j],t1=tim.g[j+1],Q=Qe))
  #pp<-sample_path_mr(a=a[[j]][1],b=a[[j]][2],t0=tim.g[j],t1=tim.g[j+1],Q=Qe)
  #full.path<-rbind(full.path,pp)
  N.R<-lapply(1:length(pp),function(t) {jump(n.stat,pp[[t]][,2],pp[[t]][,1])})
  #N.R<-jump(n.stat,pp[,2],pp[,1])
  time.R<-Reduce("+",lapply(N.R, function(t){t$R}))
  time.N<-Reduce("+",lapply(N.R, function(t){t$N}))
  #}
  return(list(time.N=time.N,time.R=time.R))
}



forward.backward <- function(Qe, pi, beta, subdat) {
  n.stat <- length(pi)
  n.cov <- ncol(beta)
  N <- nrow(subdat)
  
  beta.mult <- t(as.matrix(subdat[c("int")]) %*% beta)
  def <- mapply(dnorm, subdat$y, as.list(data.frame(beta.mult)),sd=sig)
  expm.Q <- lapply(subdat$del, function(x) expm(Qe*x))
  
  alp <- array(dim = c(N, n.stat))
  alp[1,] <- pi * def[,1]
  
  for (j in 2:N) {
    sec <- c(alp[j-1,] %*% expm.Q[[j]])
    alp[j,] <- def[,j] * sec
  }
  
  gam <- array(dim = c(N, n.stat))
  gam[N,] <- rep(1, n.stat)
  left <- array(dim = c(n.stat, N))
  
  for (m in (N-1):1) {
    left[,m+1] <- def[,m+1] * gam[m+1,]
    val <- expm.Q[[m+1]] %*% left[,m+1]
    gam[m,] <- c(val)
  }
  
  # b <- array(dim=c(N-1, n.stat, n.stat))
  b <- array(0, dim=c(N, n.stat, n.stat))
  dem <- sum(alp[1,]*gam[1,])
  for (m in (N-1):1) {
    right <- t(expm.Q[[m+1]] * alp[m,])
    #dem <- sum(left[,m+1] %*% right)
    for (k in 1:n.stat)
      b[m,k,] <- left[,m+1] * right[,k] / dem
  }
  
  # a <- t(apply(b, 1, rowSums))
  a <- t(apply(b[-N,,], 1, rowSums))
  a <- rbind(a, alp[N,] / sum(alp[N,]))
  list(a=a, b=b)
}

marlik.obs1 <- function(Qe, pi, beta, subdat) {
  n.stat <- length(pi)
  n.cov <- ncol(beta)
  N <- nrow(subdat)
  
  beta.mult <- t(as.matrix(subdat[c("int")]) %*% beta)
  def <- mapply(dnorm, subdat$y, as.list(data.frame(beta.mult)),sig)
  if(n.stat>1){
    expm.Q <- lapply(subdat$del, function(x) expm(Qe*x))
    
    alp <- array(dim = c(N, n.stat))
    alp[1,] <- pi * def[,1]
    
    for (j in 2:N) {
      sec <- c(alp[j-1,] %*% expm.Q[[j]])
      alp[j,] <- def[,j] * sec
    }
    #marlik<-log(apply(alp,1,sum))
    #marloglik.obs<-c(marlik[1],diff(marlik))
    marloglik.obs<-log(sum(alp[N,]))
  } else {
    marloglik.obs<-sum(log(def)) 
  }
  
  return(marloglik.obs)
}


### Set up initial values and hyper parameters

N.iter<-10000
burn<-0
n.cov<-1
data<-simdata
Q.s<-list()
Q.s[[1]]<-as.matrix(0)
label.c<-array(0,c(burn+N.iter,length(unique(data$id))))
Qem<-list()
Qem[[1]]<-Q.s



glm(y~1,data=data)
beta<-list()
beta[[1]]<-list((as.matrix(0.3469)))
ini.p<-list()
ini.p[[1]]<-list(c(1))

label.c[1,]<-rep(1,length(unique(data$id)))
fir.obs1<-c(1,cumsum(as.vector(table(data$id)))+1)
fir.obs1<-fir.obs1[1:(length(fir.obs1)-1)]
las.obs1<-c(fir.obs1[-1]-1,dim(data)[1])


prior.n0<-2
prior.r0<-4

prior.n<-1
prior.r<-2


a.pl<-0.1
prior.ini<-1
prior.mean<-0


split.datatt<-function(data){
  data <- data[order(data$id, data$tim.g),]
  split.data <- split(data, data$id)
  split.data <- 
    lapply(split.data, function(x) data.frame(x, obs.num = 1:nrow(x)))
  return(split.data)
}

split.data<-split.datatt(data)



sig <-1

prior.loglik.q0<-function(q){
  sum(dgamma(q,shape=prior.n0,rate=prior.r0,log=TRUE))
}
prior.loglik.q<-function(q){
  sum(dgamma(q,shape=prior.n,rate=prior.r,log=TRUE))
}
prior.loglik.ini<-function(ini.dis){
  log(ddirichlet(ini.dis, rep(prior.ini,length(ini.dis))))
}
dpois0<-function(xv,lam,log=F){ #Poisson pmf with zero removed.
  fv<-dpois(xv,lam,log=log)-log(1-exp(-lam))
  if(log){return(fv)
  }else{return(exp(fv))
  }
}

dpois1<-function(xv,lam,log=F){ #Poisson pmf with zero removed.
  fv<-dpois(xv,lam,log=log)-log(1-exp(-lam))
  if(log){return(fv)
  }else{return(exp(fv))
  }
}
lamd<-2
lamd.cl<-3.5
