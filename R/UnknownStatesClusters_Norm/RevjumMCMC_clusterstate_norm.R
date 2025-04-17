library(plyr)
library(castor)
p.be.sq<-0.5
p.be.sq1<-0.3

varpi<-list()
varpi[[1]]<-c(1)
marlik.obs <- function(vap,Qe, pi, beta, subdat) {
  n.stat <- length(pi)
  n.cov <- ncol(beta)
  N <- nrow(subdat)
  
  beta.mult <- t(as.matrix(subdat[c("int")]) %*% beta)
  def <- mapply(dnorm, subdat$y, as.list(data.frame(beta.mult)),sig)
  expm.Q <- lapply(subdat$del, function(x) expm(Qe*x))
  
  alp <- array(dim = c(N, n.stat))
  alp[1,] <- pi * def[,1]
  
  for (j in 2:N) {
    sec <- c(alp[j-1,] %*% expm.Q[[j]])
    alp[j,] <- def[,j] * sec
  }
 
  marloglik.obs<-(vap*sum(alp[N,]))
  
  return(marloglik.obs)
}

for (i in 2:(burn+N.iter)){
  num.cluster<-length(unique(label.c[i-1,]))
  ###update number of states in each cluster
  Q.stat.list<-list()
  ini.stat.list<-list()
  beta.stat.list<-list()
  for (cl in 1:num.cluster){
    kset.s<-which(label.c[i-1,]==cl)
    Qe.st<-Qem[[i-1]][[cl]]
    ini.st<-ini.p[[i-1]][[cl]]
    beta.st<-beta[[i-1]][[cl]]
    
    Q.dem<-ncol(Qe.st)
    if(runif(1)<0.5 || Q.dem==1){
      label<-sample(1:(Q.dem),1)
      
      if (Q.dem>1){
        Qe.pp<-array(0,c((Q.dem+1),(Q.dem+1)))
        Qe.pp[1:Q.dem,1:Q.dem]<-Qe.st
        Qe.pp[Q.dem+1,1:Q.dem]<-Qe.st[label,]
        Qe.pp[1:Q.dem,Q.dem+1]<-Qe.st[,label]
        w.Q<-rbeta(Q.dem+1,2,2)
        w.Q[c(label,Q.dem+1)]<-0
        Qe.pp[,label]<-Qe.pp[,label]*(w.Q)
        Qe.pp[,Q.dem+1]<-Qe.pp[,Q.dem+1]*(1-w.Q)
        Qe.pp[label,Q.dem+1]<-rgamma(1,shape=prior.n,rate=prior.r)
        Qe.pp[Q.dem+1,label]<-rgamma(1,shape=prior.n,rate=prior.r)
        
        diag(Qe.pp)<-0
        diag(Qe.pp) <- -rowSums(Qe.pp)
        
        q.val2<-c(Qe.pp[Q.dem+1,-c(label,Q.dem+1)])
        
        prior.o<-prior.loglik.q0(c(Qe.st[Qe.st>0]))+prior.loglik.ini(ini.st)
        
        
        prop.w<-log(prod(Qe.st[-label,label]))-dbeta(w.Q[w.Q>0],2,2,log=TRUE)
      } else{
        
        Qe.pp<-array(rgamma((Q.dem+1)*(Q.dem+1),shape=prior.n,rate=prior.r),c((Q.dem+1),(Q.dem+1)))
        diag(Qe.pp)<-0
        diag(Qe.pp) <- -rowSums(Qe.pp)
        prior.o<-dpois0(Q.dem,lamd,log=T)
       
        prop.w<-0
        
      }
      q.val1<-c(Qe.pp[label,Q.dem+1],Qe.pp[Q.dem+1,label])
      phi<-rbeta(1,2,2)
      ini.pp<-c(ini.st,0)
      ini.pp[c(label,length(ini.pp))]<-ini.st[label]*c(phi,1-phi)
      
      bet.pp<-c(beta.st,0)
      bet.pp[Q.dem+1]<-bet.pp[label]
      
      #u<-rnorm(1,0,p.be.sq)
      bet.pp[Q.dem+1]<-bet.pp[Q.dem+1]+rnorm(1,0,p.be.sq)
      
      q.lik.o<-sum(unlist(mclapply(split.data[kset.s], marlik.obs1, 
                                   Qe=Qe.st, pi=ini.st, beta=beta.st)))
      
      q.lik.s1<-sum(unlist(mclapply(split.data[kset.s], marlik.obs1, 
                                    Qe=Qe.pp, pi=ini.pp, beta=bet.pp)))
      prior.s1<-prior.loglik.q0(c(Qe.pp[Qe.pp>0]))+prior.loglik.ini(ini.pp)+dpois0(Q.dem+1,lamd,log=T)
      
      
      like.r.s<-(q.lik.s1-q.lik.o)+(prior.s1-prior.o)
      
      p.split<-log(Q.dem+1)+log(ini.st[label])-dbeta(phi,2,2,log=TRUE)-
        sum(dgamma(q.val1,shape=prior.n,rate=prior.r,log=TRUE))+prop.w-
        sum(dnorm(bet.pp[Q.dem+1],mean=beta.st[label],sd=p.be.sq,log=TRUE))
     
      
      jmp <- min(like.r.s+p.split, 0)
      
      if(log(runif(1)) <= jmp){
        bet.pp<-bet.pp[order(bet.pp)]
        Qe.pp<-Qe.pp[order(bet.pp),order(bet.pp)]
        ini.pp<-ini.pp[order(bet.pp)]
        
        a.and.b.s <- 
          mclapply(split.data[kset.s],forward.backward,
                   Qe=Qe.pp, pi=ini.pp, 
                   beta=bet.pp)
        split.as <- mclapply(a.and.b.s, `[[`, "a")
        aas <- do.call(rbind, split.as)
        split.bs <- mclapply(a.and.b.s, `[[`, "b") 
        ss1<-NULL
        ss1<-apply(apply(aas,1,function(x) rmultinom(1,1,prob=x)),2,which.max)
        ss1<-factor(ss1,levels=c(1:(Q.dem+1)))
        ###update inital 
        fir.obs<-c(1,cumsum(as.vector(table(data$id[data$id%in%kset.s])))+1)
        fir.obs<-fir.obs[1:(length(fir.obs)-1)]
        
        
        ini.st<-rdirichlet(1,1+prior.ini+table(factor(ss1[fir.obs], levels = c(1:(Q.dem+1)))))
        
        ### update beta
        obs.mean<-tapply(data$y[data$id%in%kset.s],ss1,mean)
        obs.count<-table(ss1)
        post.mean<-obs.mean/(1+1/obs.count)+(1/obs.count)*prior.mean/(1+1/obs.count)
        post.mean[obs.count==0]<-prior.mean
        beta.st<-rnorm(length(obs.count),post.mean,sig*sqrt((1/(1+obs.count))))
        
        ###update Q
        
        doubless<-mclapply(split.bs,function(y) 
        {apply(y,1, function (x) {if (sum(x)!=0){return(which(
          matrix(rmultinom(1,1,prob=x),byrow=F,ncol=(Q.dem+1))==1,arr.ind=T))} else{
            return(0)} })
        })
        
        time.RQ.Q<-mclapply(1:length(split.data[kset.s]),function(tt){markov.sim(split.data[[kset.s[tt]]],a=doubless[[tt]],Qe=Qe.pp)})
        N.n.Q<-mclapply(time.RQ.Q, function(x) x$time.N)
        R.n.Q<-mclapply(time.RQ.Q, function(x) x$time.R)
        N.n1<-Reduce(`+`, N.n.Q)
        R.n1<-Reduce(`+`, R.n.Q)
        Q.n<-array(0,c(Q.dem+1,Q.dem+1))
        for (mmm in 1:(Q.dem+1)){
          Q.n[mmm,-mmm]<-rgamma(Q.dem,shape=N.n1[mmm,-mmm]+prior.n+1,rate=R.n1[mmm]+prior.r)
        }
        diag(Q.n) <- -rowSums(Q.n)
        Qe.st<-Q.n 
        
        
      } 
      
    } else {
      
      if (Q.dem==2){
        
        labelm<-sample(c(1:2),2)
        
        m.int<-mean(beta.st[1])
        
        bet.pp<-beta.st
        bet.pp[labelm[1]]<-m.int
        bet.pp<-bet.pp[-labelm[2]]
        
        q.val1<--diag(Qe.st)
        Qe.pp<-as.matrix(0)
        
        
        orgw<-1
        phi<-ini.st[labelm[1]]/sum(ini.st[labelm])
        ini.pp<-1
        q.lik.m<-sum(unlist(mclapply(split.data, marlik.obs1, 
                                     Qe=Qe.pp, pi=ini.pp, beta=bet.pp)))
        
        prior.m<-dpois0(Q.dem-1,lamd,log=T)
        
        #prop.q2m<-0
        prop.wm<-0
        
        
      } else {
        
        diff.m<-which.min(diff(beta.st))
        order.label<-order(beta.st)
        pair.mm<-rbind(order.label,c(order.label[-1],order.label[1]))[,-length(order.label)]
        labelm<-sample((pair.mm[,diff.m]),2)
        
        
        
        ini.pp<-ini.st
        ini.pp[labelm[1]]<-ini.st[labelm[1]]+ini.st[labelm[2]]
        
        
        phi<-ini.st[labelm[1]]/ini.pp[labelm[1]]
        orgw<-ini.pp[labelm[1]]
        pi.stat<-get_stationary_distribution(Qe.st)[labelm]
        ini.pp<-ini.pp[-labelm[2]]
        Qe.pp<-Qe.st
        q.val1<-c(Qe.pp[labelm[1],labelm[2]],Qe.pp[labelm[2],labelm[1]])
        
        Qe.pp[,labelm[1]]<- Qe.pp[,labelm[1]]+ Qe.pp[,labelm[2]]
        Qe.pp[labelm[1],]<-pi.stat[1]/sum(pi.stat)*Qe.pp[labelm[1],]+pi.stat[2]/sum(pi.stat)*Qe.pp[labelm[2],]
        W.qm<-Qe.st[,labelm[1]]/Qe.pp[,labelm[1]]
        W.qm[c(labelm[1],labelm[2])]<-0
        
        prop.wm<-sum(log(prod(Qe.pp[-c(labelm[1],labelm[2]),labelm[1]])))-sum(dbeta(W.qm[W.qm>0],2,2,log=TRUE))
        
        Qe.pp<-Qe.pp[-labelm[2],-labelm[2]]
        diag(Qe.pp) <- 0
        diag(Qe.pp) <- -rowSums(Qe.pp)
        
        
        bet.pp<-beta.st
        bet.pp[labelm[1]]<-beta.st[labelm[1]]*(pi.stat[1]/sum(pi.stat))+beta.st[labelm[2]]*pi.stat[2]/sum(pi.stat)
        
        bet.pp<-bet.pp[-labelm[2]]
        
        q.lik.m<-sum(unlist(mclapply(split.data[kset.s], marlik.obs1, 
                                     Qe=Qe.pp, pi=ini.pp, beta=bet.pp)))
        
        prior.m<-prior.loglik.ini(ini.pp)+prior.loglik.q0(c(Qe.pp[Qe.pp>0]))+dpois0(Q.dem-1,lamd,log=T)
        
      }
      
      
      
      q.lik.m.s1<-sum(unlist(mclapply(split.data[kset.s], marlik.obs1, 
                                      Qe=Qe.st, pi=ini.st, beta=beta.st)))
      
      prior.m1<-prior.loglik.q0(c(Qe.st[Qe.st>0]))+prior.loglik.ini(ini.st)+dpois0(Q.dem,lamd,log=T)
      
      
      like.r.m<-(q.lik.m-q.lik.m.s1)+(prior.m-prior.m1)
      
      p.merge<--log(Q.dem)-log(orgw)+dbeta(phi,2,2,log=TRUE)+
        sum(dgamma(q.val1,shape=prior.n,rate=prior.r,log=TRUE))-prop.wm
      
      
      jmp <- min(like.r.m+p.merge, 0)
      if(log(runif(1)) <= jmp){
        bet.pp<-bet.pp[order(bet.pp)]
        Qe.pp<-Qe.pp[order(bet.pp),order(bet.pp)]
        ini.pp<-ini.pp[order(bet.pp)]
        
        if (Q.dem==2){
          ini.st<-1
          Qe.pp<-as.matrix(0)
          bet.pp<-rnorm(1,mean(data$y[data$id%in%kset.s]),sd=sig)
          
        } else{
          
          a.and.b.s <- 
            mclapply(split.data[kset.s],forward.backward,
                     Qe=Qe.pp, pi=ini.pp, 
                     beta=bet.pp)
          split.as <- mclapply(a.and.b.s, `[[`, "a")
          aas <- do.call(rbind, split.as)
          split.bs <- mclapply(a.and.b.s, `[[`, "b") 
          ss1<-NULL
          ss1<-apply(apply(aas,1,function(x) rmultinom(1,1,prob=x)),2,which.max)
          ss1<-factor(ss1,levels=c(1:(Q.dem-1)))
          ###update inital 
          fir.obs<-c(1,cumsum(as.vector(table(data$id[data$id%in%kset.s])))+1)
          fir.obs<-fir.obs[1:(length(fir.obs)-1)]
          
          
          ini.st<-rdirichlet(1,1+prior.ini+table(factor(ss1[fir.obs], levels = c(1:(Q.dem-1)))))
          
          ### update beta
          obs.mean<-tapply(data$y[data$id%in%kset.s],ss1,mean)
          obs.count<-table(ss1)
          post.mean<-obs.mean/(1+1/obs.count)+(1/obs.count)*prior.mean/(1+1/obs.count)
          post.mean[obs.count==0]<-prior.mean
          beta.st<-rnorm(length(obs.count),post.mean,sig*sqrt((1/(1+obs.count))))
          
          ###update Q
          
          doubless<-mclapply(split.bs,function(y) 
          {apply(y,1, function (x) {if (sum(x)!=0){return(which(
            matrix(rmultinom(1,1,prob=x),byrow=F,ncol=(Q.dem-1))==1,arr.ind=T))} else{
              return(0)} })
          })
          
          time.RQ.Q<-mclapply(1:length(split.data[kset.s]),function(tt){markov.sim(split.data[[kset.s[tt]]],a=doubless[[tt]],Qe=Qe.pp)})
          N.n.Q<-mclapply(time.RQ.Q, function(x) x$time.N)
          R.n.Q<-mclapply(time.RQ.Q, function(x) x$time.R)
          N.n1<-Reduce(`+`, N.n.Q)
          R.n1<-Reduce(`+`, R.n.Q)
          Q.n<-array(0,c(Q.dem-1,Q.dem-1))
          for (mmm in 1:(Q.dem-1)){
            Q.n[mmm,-mmm]<-rgamma(Q.dem-2,shape=N.n1[mmm,-mmm]+prior.n+1,rate=R.n1[mmm]+prior.r)
          }
          diag(Q.n) <- -rowSums(Q.n)
          Qe.st<-Q.n 
        }
      } 
      
    }
    
    
    Q.stat.list[[cl]]<-Qe.st
    beta.stat.list[[cl]]<-beta.st
    ini.stat.list[[cl]]<-ini.st
  } 
  
  ###update number of clusters
  
  vap.vec<-varpi[[i-1]]
  n.stat<-unlist(lapply(ini.stat.list,length))
  if(runif(1)<0.5 ||  num.cluster==1 || !any(duplicated(n.stat))){
    
    label<-sample(sort(unique(label.c[i-1,]))[table(label.c[i-1,])!=1],1)
    kset.s<-which(label.c[i-1,]==label)
    
    cl.stat<-n.stat[label]
    
    Qe.ss<-c(Q.stat.list,Q.stat.list[label])
    
    
    phi<-rbeta(1,2,2)
    vap.ss<-c(vap.vec,0)
    vap.ss[label]<-vap.vec[label]*phi
    vap.ss[num.cluster+1]<-vap.vec[label]*(1-phi)
    
    ini.ss<-c(ini.stat.list,ini.stat.list[label])
    
    bet.pp<- beta.stat.list[[label]]
    u<-rnorm(cl.stat,0,p.be.sq1)
    bet.pp<-bet.pp+u
    bet.ss<-c(beta.stat.list,list(bet.pp))
    
    if(num.cluster>1){
      thre<-sort(unlist(mclapply(beta.stat.list, function (x) sum(abs(x-beta.stat.list[[label]])))))[2]} else{
        thre<-Inf
      }
    if (sum(abs(u))>thre) {
      jmp<--Inf
    } else {
      
      prior.o<-dpois0(num.cluster,lamd.cl,log=T)
      q.lik.o<-sum(unlist(mclapply(1:length(split.data), function (tt) { lik<-NULL
      for (jjj in 1:num.cluster) {lik<-c(lik,marlik.obs(vap=vap.vec[jjj],Qe=Q.stat.list[[jjj]], pi=ini.stat.list[[jjj]], 
                                                        beta=beta.stat.list[[jjj]],subdat=split.data[[tt]]))}
      return(log(sum((lik))))})))
      
      
      prior.s1<-prior.loglik.q(c(Q.stat.list[[label]][Q.stat.list[[label]]>0]))+prior.loglik.ini(ini.stat.list[[label]])+dpois0(num.cluster+1,lamd.cl,log=T)
      q.lik.s1<-sum(unlist(mclapply(1:length(split.data), function (tt) { lik<-NULL
      for (jjj in 1:(num.cluster+1)) {lik<-c(lik,marlik.obs(vap=vap.ss[jjj],Qe=Qe.ss[[jjj]], pi=ini.ss[[jjj]], 
                                                            beta=bet.ss[[jjj]],subdat=split.data[[tt]]))}
      return(log(sum((lik))))})))
      
      like.r.s<-(q.lik.s1-q.lik.o)+(prior.s1-prior.o)
      
      p.split<-log(num.cluster+1)+log(vap.vec[label])-sum(dbeta(phi,2,2,log=TRUE))-sum(dnorm(u,mean=0,sd=p.be.sq1,log=TRUE))
      
      jmp <- min(like.r.s+p.split, 0)
      
    }
    
    if(log(runif(1)) <= jmp){
      Q.stat.list[[num.cluster+1]]<-Q.stat.list[[label]]
      beta.stat.list[[num.cluster+1]]<-bet.pp
      ini.stat.list[[num.cluster+1]]<-ini.stat.list[[label]]
      vap.vec<-vap.ss
      num.cluster<-num.cluster+1
    }
    
  } else {
    dups<-n.stat[duplicated(n.stat)]
    if(length(dups)==1){
      m.stat<-dups
    } else{
      m.stat<-sample(dups,1)
    }
    merge.cl<-which(n.stat==m.stat)
    sam<-sample(merge.cl,1)
    cut.mean<-do.call(rbind,beta.stat.list[merge.cl])
    cl.sam<-order(apply(cut.mean,1, function (x) sum(abs(x-cut.mean[which(merge.cl==sam),]))))[2]
    sam2<-merge.cl[cl.sam]
    
    labelm<-sort(c(sam,sam2))
    #labelm<-label.c[i-1,ran.sm]
    phi<-vap.vec[labelm[1]]/sum(vap.vec[labelm])
    orgw<-sum(vap.vec[labelm])
    pi.weight<-get_stationary_distribution(Q.stat.list[[labelm[1]]])/
      (get_stationary_distribution(Q.stat.list[[labelm[1]]])+get_stationary_distribution(Q.stat.list[[labelm[2]]]))
    
    ini.mm<-pi.weight*ini.stat.list[[labelm[1]]]+ (1-pi.weight)*ini.stat.list[[labelm[2]]]
    ini.mm<-ini.mm/sum(ini.mm)
    ini.m<-ini.stat.list
    ini.m[[labelm[1]]]<-ini.mm
    ini.m[[labelm[2]]]<-ini.mm
    
    Q.merge<- t(pi.weight*t(Q.stat.list[[labelm[1]]])+ (1-pi.weight)*t(Q.stat.list[[labelm[2]]]))
    Q.m<-Q.stat.list
    Q.m[[labelm[1]]]<-Q.merge
    Q.m[[labelm[2]]]<-Q.merge
    
    be.merge<-(pi.weight*(beta.stat.list[[labelm[1]]])+(1-pi.weight)*(beta.stat.list[[labelm[2]]]))
    be.m<-beta.stat.list
    be.m[[labelm[1]]]<-be.merge
    be.m[[labelm[2]]]<-be.merge
    
    
    q.lik.m<-sum(unlist(mclapply(1:length(split.data), function (tt) { lik<-NULL
    for (jjj in 1:num.cluster) {lik<-c(lik,marlik.obs(vap=vap.vec[jjj],Qe=Q.m[[jjj]], pi=ini.m[[jjj]], 
                                                      beta=be.m[[jjj]],subdat=split.data[[tt]]))}
    return(log(sum((lik))))})))
    prior.m<-prior.loglik.ini(ini.mm)+prior.loglik.q(c(Q.merge[Q.merge>0]))+dpois0(num.cluster-1,lamd.cl,log=T)
    
    
    q.lik.m.o<-sum(unlist(mclapply(1:length(split.data), function (tt) { lik<-NULL
    for (jjj in 1:num.cluster) {lik<-c(lik,marlik.obs(vap=vap.vec[jjj],Qe=Q.stat.list[[jjj]], pi=ini.stat.list[[jjj]], 
                                                      beta=beta.stat.list[[jjj]],subdat=split.data[[tt]]))}
    return(log(sum((lik))))})))
    
    Qe.org<-as.matrix(bdiag(Q.stat.list[[labelm[1]]],Q.stat.list[[labelm[2]]]))
    prior.m.o<-prior.loglik.q(c(Qe.org[Qe.org>0]))+prior.loglik.ini(ini.stat.list[[labelm[1]]])+prior.loglik.ini(ini.stat.list[[labelm[2]]])+
      dpois0(num.cluster,lamd.cl,log=T)
    
    like.r.m<-(q.lik.m-q.lik.m.o)+(prior.m-prior.m.o)
    
    p.merge<--log(num.cluster-1)-sum(log(orgw))+sum(dbeta(phi,2,2,log=TRUE))
    
    jmp <- min(like.r.m+p.merge, 0)
    
    if(log(runif(1)) <= jmp){
      
      Q.stat.list[[labelm[1]]]<-Q.merge
      beta.stat.list[[labelm[1]]]<-be.merge
      ini.stat.list[[labelm[1]]]<-ini.mm
      vap.vec[labelm[1]]<-vap.vec[labelm[1]]+vap.vec[labelm[2]]
      
      Q.stat.list<-Q.stat.list[-labelm[2]]
      beta.stat.list<- beta.stat.list[-labelm[2]]
      ini.stat.list<-ini.stat.list[-labelm[2]]
      num.cluster<-num.cluster-1
    }
    
  }
  ###update cluster membership
  n.stat<-unlist(mclapply(ini.stat.list,length))
  margin.lik<-
    mclapply(1:length(split.data), function (tt) { lik<-NULL
    for (jjj in 1:num.cluster) {lik<-c(lik,marlik.obs1(Qe=Q.stat.list[[jjj]], pi=ini.stat.list[[jjj]], 
                                                       beta=beta.stat.list[[jjj]],subdat=split.data[[tt]]))}
    return(lik)})
  
  post.prob<-do.call(rbind,mclapply(margin.lik, function(x) exp(x-max(x))))
  post.prob[is.nan(post.prob)]<-0
  label.c[i,]<-apply(post.prob,1,function(t){
    sample(1:num.cluster,1,prob = t)})
  
  emp.cl<-which(table(factor(label.c[i,],levels=1:num.cluster))==0)
  if( length(emp.cl)>0){
    num.cluster<-num.cluster-length(emp.cl)
    label.c[i,]<- mapvalues(label.c[i,], from = names(table(label.c[i,])), to = c(1:(num.cluster)))
  }
  
  ###
  varpi[[i]]<-rdirichlet(1,table(label.c[i,])+1)
  
  
  
  #####update parameters associated with each cluster
  Q.n.list<-list()
  ini.list<-list()
  beta.list<-list()
  
  a.and.b <- 
    mclapply(1:length(split.data),function(tt){ forward.backward(
      Qe=Q.stat.list[[label.c[i,tt]]], pi=ini.stat.list[[label.c[i,tt]]], 
      beta=beta.stat.list[[label.c[i,tt]]],subdat=split.data[[tt]])})
  
  for (jj in 1:num.cluster){
    mm11<-which(label.c[i,]==jj)
    split.a <- mclapply(a.and.b[mm11], `[[`, "a")
    aa <- do.call(rbind, split.a)
    split.b <- mclapply(a.and.b[mm11], `[[`, "b") 
    ss<-NULL
    ss<-apply(apply(aa,1,function(x) rmultinom(1,1,prob=x)),2,which.max)
    ss<-factor(ss,levels=c(1:n.stat[jj]))
    ###update inital 
    fir.obs<-c(1,cumsum(as.vector(table(data$id[data$id%in%mm11])))+1)
    fir.obs<-fir.obs[1:(length(fir.obs)-1)]
    
    
    ini.list[[jj]]<-rdirichlet(1,1+prior.ini+table(factor(ss[fir.obs], levels = c(1:n.stat[jj]))))
    
    ### update beta
    obs.mean<-tapply(data$y[data$id%in%mm11],ss,mean)
    obs.count<-table(ss)
    post.mean<-obs.mean/(1+1/obs.count)+(1/obs.count)*prior.mean/(1+1/obs.count)
    post.mean[obs.count==0]<-prior.mean
    beta.list[[jj]]<-rnorm(length(obs.count),post.mean,sig*sqrt((1/(1+obs.count))))
    
    ###update Q
    
    doubless<-mclapply(split.b,function(y) 
    {apply(y,1, function (x) {if (sum(x)!=0){return(which(
      matrix(rmultinom(1,1,prob=x),byrow=F,ncol=n.stat[jj])==1,arr.ind=T))} else{
        return(0)} })
    })
    
    time.RQ.Q<-mclapply(1:length(split.data[mm11]),function(tt){markov.sim(split.data[[mm11[tt]]],a=doubless[[tt]],Qe=Q.stat.list[[jj]])})
    N.n.Q<-mclapply(time.RQ.Q, function(x) x$time.N)
    N.n.Q<-mclapply(N.n.Q, function(x) {
      x[is.na(x)]<-0; return(x)})
    R.n.Q<-mclapply(time.RQ.Q, function(x) x$time.R)
    R.n.Q<-mclapply(R.n.Q, function(x) {
      x[is.na(x)]<-0; return(x)})
    N.n1<-Reduce(`+`, N.n.Q)
    R.n1<-Reduce(`+`, R.n.Q)
    Q.n<-array(0,c(n.stat[jj],n.stat[jj]))
    for (mmm in 1:(n.stat[jj])){
      Q.n[mmm,-mmm]<-rgamma((n.stat[jj])-1,shape=N.n1[mmm,-mmm]+prior.n+1,rate=R.n1[mmm]+prior.r)
    }
    diag(Q.n) <- -rowSums(Q.n)
    Q.n.list[[jj]]<-Q.n 
    
  }
  ini.p[[i]]<-ini.list
  Qem[[i]]<-Q.n.list
  beta[[i]]<-beta.list
  
  print(list(num.cluster=table(label.c[i,]),n.stat=n.stat))
  plot(unlist(mclapply(ini.p,length)),type="s") 
}
