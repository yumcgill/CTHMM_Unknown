# loglik prior for Q
prior.loglik.q<-function(q){
  sum(dgamma(q,shape=prior.n,rate=prior.r,log=TRUE))
}
#loglik prior for pi
prior.loglik.ini<-function(ini.dis){
  log(ddirichlet(ini.dis, rep(prior.ini,length(ini.dis))))
}
dpois0<-function(xv,lam,log=F){ #Poisson pmf with zero removed.
  fv<-dpois(xv,lam,log=log)-log(1-exp(-lam))
  if(log){return(fv)
  }else{return(exp(fv))
  }
}
### prior for number of states
lamd<-3.5

ss<-rep(1,nrow(data))
library(castor)
## sigma for the proposal
p.be.sq<-0.25

for (i in 2:(burn+N.iter)){
  Q.dem<-ncol(Qem[[i-1]])
  
  ###Update K
  if(runif(1)<0.5 || Q.dem==1){
    label<-sample(1:(Q.dem),1)
    
    ###proposal for Q
    if (Q.dem>1){
      Qe.pp<-array(0,c((Q.dem+1),(Q.dem+1)))
      Qe.pp[1:Q.dem,1:Q.dem]<-Qem[[i-1]]
      Qe.pp[Q.dem+1,1:Q.dem]<-Qem[[i-1]][label,]
      Qe.pp[1:Q.dem,Q.dem+1]<-Qem[[i-1]][,label]
      w.Q<-rbeta(Q.dem+1,2,2)
      w.Q[c(label,Q.dem+1)]<-0
      Qe.pp[,label]<-Qe.pp[,label]*(w.Q)
      Qe.pp[,Q.dem+1]<-Qe.pp[,Q.dem+1]*(1-w.Q)
      Qe.pp[label,Q.dem+1]<-rgamma(1,shape=prior.n,rate=prior.r)
      Qe.pp[Q.dem+1,label]<-rgamma(1,shape=prior.n,rate=prior.r)
  
      diag(Qe.pp)<-0
      diag(Qe.pp) <- -rowSums(Qe.pp)

      
      prior.o<-prior.loglik.q(c(Qem[[i-1]][Qem[[i-1]]>0]))+prior.loglik.ini(ini.p[[i-1]])+dpois0(Q.dem,lamd,log=T)
      prop.w<-sum(log(prod(Qem[[i-1]][-label,label])))-sum(dbeta(w.Q[w.Q>0],2,2,log=TRUE))
      
    } else{
      
      Qe.pp<-array(rgamma((Q.dem+1)*(Q.dem+1),shape=prior.n,rate=prior.r),c((Q.dem+1),(Q.dem+1)))
      diag(Qe.pp)<-0
      diag(Qe.pp) <- -rowSums(Qe.pp)
      prior.o<-dpois0(Q.dem,lamd,log=T)
      
      prop.w<-0
      
    }
    q.val1<-c(Qe.pp[label,Q.dem+1],Qe.pp[Q.dem+1,label])
   
   ### proposal for pi
    phi<-rbeta(1,2,2)
    ini.pp<-c(ini.p[[i-1]],0)
    ini.pp[c(label,length(ini.pp))]<-ini.p[[i-1]][label]*c(phi,1-phi)
    ###proposal for beta
    bet.pp<-cbind(beta[[i-1]],0)
    bet.pp[,Q.dem+1]<-bet.pp[,label]

    u<-rnorm(1,0,p.be.sq)
    bet.pp[1,Q.dem+1]<-beta[[i-1]][1,label]+u
  
  ##calculate the marginal likelihood
    q.lik.o<-sum(unlist(mclapply(split.data, marlik.obs, 
                                 Qe=Qem[[i-1]], pi=ini.p[[i-1]], beta=beta[[i-1]])))
    
    q.lik.s1<-sum(unlist(mclapply(split.data, marlik.obs, 
                                  Qe=Qe.pp, pi=ini.pp, beta=bet.pp)))
    prior.s1<-prior.loglik.q(c(Qe.pp[Qe.pp>0]))+prior.loglik.ini(ini.pp)+dpois0(Q.dem+1,lamd,log=T)
    
    
    like.r.s<-(q.lik.s1-q.lik.o)+(prior.s1-prior.o)
    
    p.split<-log(Q.dem+1)+log(ini.p[[i-1]][label])-dbeta(phi,2,2,log=TRUE)-
      sum(dgamma(q.val1,shape=prior.n,rate=prior.r,log=TRUE))+prop.w-
      sum(dnorm(bet.pp[1,Q.dem+1],mean=beta[[i-1]][1,label],sd=p.be.sq,log=TRUE))
   
    
    jmp <- min(like.r.s+p.split, 0)
    
    if(log(runif(1)) <= jmp){
      Qe.ini<-Qe.pp
      beee<-bet.pp
      ini.sam<-ini.pp
      
    } else {
      Qe.ini<-Qem[[i-1]]
      beee<-beta[[i-1]]
      ini.sam<-ini.p[[i-1]]
    }
    
  } else {
    
    if (Q.dem==2){
      
      labelm<-sample(c(1:2),2)
      
      m.int<-mean(beta[[i-1]][1,])
      
      bet.pp<-beta[[i-1]]
      bet.pp[1,labelm[1]]<-m.int
      bet.pp<-bet.pp[,-labelm[2]]
      
      q.val1<--diag(Qem[[i-1]])
      Qe.pp<-as.matrix(0)
      
      
      orgw<-1
      phi<-ini.p[[i-1]][labelm[1]]/sum(ini.p[[i-1]][labelm])
      ini.pp<-1
      q.lik.m<-sum(unlist(mclapply(split.data, marlik.obs, 
                                   Qe=Qe.pp, pi=ini.pp, beta=bet.pp)))
      
      prior.m<-dpois0(Q.dem-1,lamd,log=T)
      
      #prop.q2m<-0
      prop.wm<-0
      
      
    } else {
      sam<-sample(1:(Q.dem),1)
      cut.mean<-tapply(data$y,ss,mean)
      sam2<-order(abs(cut.mean-cut.mean[sam]))[2]
     
      labelm<-sample(c(sam,sam2),2)
      
      ini.pp<-ini.p[[i-1]]
      ini.pp[labelm[1]]<-ini.p[[i-1]][labelm[1]]+ini.p[[i-1]][labelm[2]]
      
      
      phi<-ini.p[[i-1]][labelm[1]]/ini.pp[labelm[1]]
      orgw<-ini.pp[labelm[1]]
      pi.stat<-get_stationary_distribution(Qem[[i-1]])[labelm]
      ini.pp<-ini.pp[-labelm[2]]
      Qe.pp<-Qem[[i-1]]
      q.val1<-c(Qe.pp[labelm[1],labelm[2]],Qe.pp[labelm[2],labelm[1]])
     
      Qe.pp[,labelm[1]]<- Qe.pp[,labelm[1]]+ Qe.pp[,labelm[2]]
      Qe.pp[labelm[1],]<-pi.stat[1]/sum(pi.stat)*Qe.pp[labelm[1],]+pi.stat[2]/sum(pi.stat)*Qe.pp[labelm[2],]
      W.qm<-Qem[[i-1]][,labelm[1]]/Qe.pp[,labelm[1]]
      W.qm[c(labelm[1],labelm[2])]<-0
      
      
      prop.wm<-sum(log(prod(Qe.pp[-c(labelm[1],labelm[2]),labelm[1]])))-sum(dbeta(W.qm[W.qm>0],2,2,log=TRUE))
      
      Qe.pp<-Qe.pp[-labelm[2],-labelm[2]]
      diag(Qe.pp) <- 0
      diag(Qe.pp) <- -rowSums(Qe.pp)
      
      
      bet.pp<-beta[[i-1]]
      bet.pp[,labelm[1]]<-beta[[i-1]][,labelm[1]]*(pi.stat[1]/sum(pi.stat))+beta[[i-1]][,labelm[2]]*(pi.stat[2]/sum(pi.stat))
      
      bet.pp<-bet.pp[,-labelm[2]]
      
      q.lik.m<-sum(unlist(mclapply(split.data, marlik.obs, 
                                   Qe=Qe.pp, pi=ini.pp, beta=bet.pp)))
      
      prior.m<-prior.loglik.ini(ini.pp)+prior.loglik.q(c(Qe.pp[Qe.pp>0]))+dpois0(Q.dem-1,lamd,log=T)
      
    }
    
    
    q.lik.m.s1<-sum(unlist(mclapply(split.data, marlik.obs, 
                                    Qe=Qem[[i-1]], pi=ini.p[[i-1]], beta=beta[[i-1]])))
    
    prior.m1<-prior.loglik.q(c(Qem[[i-1]][Qem[[i-1]]>0]))+prior.loglik.ini(ini.p[[i-1]])+dpois0(Q.dem,lamd,log=T)
    
    
    like.r.m<-(q.lik.m-q.lik.m.s1)+(prior.m-prior.m1)
    
    p.merge<--log(Q.dem)-log(orgw)+dbeta(phi,2,2,log=TRUE)+
      sum(dgamma(q.val1,shape=prior.n,rate=prior.r,log=TRUE))-prop.wm
   
    jmp <- min(like.r.m+p.merge, 0)
    if(log(runif(1)) <= jmp){
      Qe.ini<-Qe.pp
      beee<-bet.pp
      ini.sam<-ini.pp
      
    } else {
      Qe.ini<-Qem[[i-1]]
      beee<-beta[[i-1]]
      ini.sam<-ini.p[[i-1]]
    } 
    
  }
  
  
  
  
  
  ###update inital 
  a.and.b <- 
    mclapply(split.data, forward.backward, 
             Qe=Qe.ini, pi=ini.sam, beta=beee)
  
  
  split.a <- lapply(a.and.b, `[[`, "a")
  aa <- do.call(rbind, split.a)
  nstat<-length(ini.sam)
  ss<-NULL
  ss<-apply(apply(aa,1,function(x) rmultinom(1,1,prob=x)),2,which.max)
  ss<-factor(ss,levels=c(1:nstat))
  
  ### update beta
  
  B<-3
  bet.tt<-array(NA,c(n.cov,nstat))
  for (k in 1:nstat){
    bet<-array(NA,c(n.cov))
    bet.cov<-array(NA,c(n.cov,n.cov))
    gen.fit<-glm(y~z1+z2,data=data[ss==k,],family = gaussian)
    bet<-gen.fit$coefficients
    bet.cov<-vcov(gen.fit)
    ###  Metropolis-Hastings
    mh.bet<-array(NA,c(B+1,n.cov))
    mh.bet[1,]<-bet
    for (mh in 2:(B+1)){
      u<-rmvnorm(1,mean=mh.bet[mh-1,],sigma=bet.cov)
      rr<-sum(dnorm(data$y[ss==k],cbind(1,data$z1[ss==k],data$z2[ss==k])%*%t(u),sigma,log=TRUE))-
        sum(dnorm(data$y[ss==k],cbind(1,data$z1[ss==k],data$z2[ss==k])%*%mh.bet[mh-1,],sigma,log=TRUE))
      rej <- min(exp(rr), 1)
      if(runif(1) <= rej) {
        mh.bet[mh,] <- t(u)
      } else {
        mh.bet[mh,] <- mh.bet[mh-1,]
      }
    }
    bet.tt[,k]<-mh.bet[B+1,]
  }
  
  beta[[i]]<-bet.tt
  
  ###update inital 
  
  
  ini.p[[i]]<-rdirichlet(1,1+prior.ini+table(factor(ss[fir.obs], levels = 1:nstat)))
  
  ### update Q
  if (nstat>1){
    
    split.b <- mclapply(a.and.b, `[[`, "b") 
   
    bb <- do.call(function(...) abind(..., along = 1), split.b)
    doubless<-apply(bb,1,function(x) 
    {if (sum(x)!=0){return(which(
      matrix(rmultinom(1,1,prob=x),byrow=F,ncol=ncol(Qe.ini))==1,arr.ind=T))} else{
        return(0)}
    })
    doubless <- split(doubless, data$id)
    
    time.RQ.Q<-mclapply(1:length(split.data),function(tt){markov.sim(split.data[[tt]],a=doubless[[tt]],Qe=Qe.ini)})
    N.n.Q<-lapply(time.RQ.Q, function(x) x$time.N)
    R.n.Q<-lapply(time.RQ.Q, function(x) x$time.R)
    N.n1<-Reduce(`+`, N.n.Q)
    R.n1<-Reduce(`+`, R.n.Q)
    Q.n<-array(0,c(nstat,nstat))
    for (mmm in 1:nstat){
      Q.n[mmm,-mmm]<-rgamma(nstat-1,shape=N.n1[mmm,-mmm]+prior.n+1,rate=R.n1[mmm]+prior.r)
    }
    diag(Q.n) <- -rowSums(Q.n)
    Qem[[i]]<-Q.n 
    
  } else {  Qem[[i]]<-as.matrix(0)   }
  
  plot(unlist(mclapply(ini.p,length)),type="s") 
} 

