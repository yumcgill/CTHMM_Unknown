library(MCMCpack)
library(ECctmc)
library(mvtnorm)
library(expm)

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

###generate the number of jumps and holding time for a given trajectory
markov.sim<-function(subdat,a,Qe){
  n.stat<-ncol(Qe)
  time.R<-array(0,c(n.stat))
  time.N<-array(0,c(n.stat,n.stat))
  tim.g<-cumsum(subdat$del)
  pp<-lapply(1:(dim(subdat)[1]-1),function(j) sample_path_mr(a=a[[j]][1],b=a[[j]][2],t0=tim.g[j],t1=tim.g[j+1],Q=Qe))
  N.R<-lapply(1:length(pp),function(t) {jump(n.stat,pp[[t]][,2],pp[[t]][,1])})
  time.R<-Reduce("+",lapply(N.R, function(t){t$R}))
  time.N<-Reduce("+",lapply(N.R, function(t){t$N}))
  return(list(time.N=time.N,time.R=time.R))
}


###forward and backward algorithm
forward.backward <- function(Qe, pi, beta, subdat) {
  n.stat <- length(pi)
  n.cov <- ncol(beta)
  N <- nrow(subdat)
  
  beta.mult <- t(as.matrix(subdat[c("int","z1","z2")]) %*% beta)
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
  
  b <- array(0, dim=c(N, n.stat, n.stat))
  dem <- sum(alp[1,]*gam[1,])
  for (m in (N-1):1) {
    right <- t(expm.Q[[m+1]] * alp[m,])
    for (k in 1:n.stat)
      b[m,k,] <- left[,m+1] * right[,k] / dem
  }
  a <- t(apply(b[-N,,], 1, rowSums))
  a <- rbind(a, alp[N,] / sum(alp[N,]))
  list(a=a, b=b)
}


### calculate the marginal likelihood
marlik.obs <- function(Qe, pi, beta, subdat) {
  n.stat <- length(pi)
  n.cov <- ncol(beta)
  N <- nrow(subdat)
  
  beta.mult <- t(as.matrix(subdat[c("int","z1","z2")]) %*% beta)
  def <- mapply(dnorm, subdat$y, as.list(data.frame(beta.mult)),sig)
  if(n.stat>1){
    expm.Q <- lapply(subdat$del, function(x) expm(Qe*x))
    
    alp <- array(dim = c(N, n.stat))
    alp[1,] <- pi * def[,1]
    
    for (j in 2:N) {
      sec <- c(alp[j-1,] %*% expm.Q[[j]])
      alp[j,] <- def[,j] * sec
    }
    marloglik.obs<-log(sum(alp[N,]))
  } else {
    marloglik.obs<-sum(log(def)) 
  }
  
  return(marloglik.obs)
}


#####setting initial values and hyper parameters
N.iter<-20000
burn<-0
n.cov<-3
data<-simdata
beta<-list()
beta[[1]]<-as.matrix(be[,3])
ini.p<-list()
ini.p[[1]]<-1

Qem<-list()
Qem[[1]]<-as.matrix(0)

fir.obs<-c(1,cumsum(as.vector(table(data$id)))+1)
fir.obs<-fir.obs[1:(length(fir.obs)-1)]
las.obs<-c(fir.obs[-1]-1,dim(data)[1])

prior.n<-1
prior.r<-2

a.pl<-0.1
prior.ini<-1
prior.mean<-0


split.datatt<-function(data){
  split.data <- split(data, data$id)
  split.data <- 
    lapply(split.data, function(x) data.frame(x, obs.num = 1:nrow(x)))
  return(split.data)
}

split.data<-split.datatt(data)

sig <-1
