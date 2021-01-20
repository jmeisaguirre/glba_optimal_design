library(foreach)
library(doParallel)
registerDoParallel(cores=20)

load('rec_work_re_y_13.RData')

#### score function
calc.q=function(x){
  y=sum((x-mean(x))^2)/length(x)
  return(y)
}

###### reduce to 1000 designs
No.designs=1000
Transects.2017.l=Transects.2017.l[1:No.designs]


###
### get indexes for imputed data
###

set.seed(23) # get the same indexes for parallel
imps=sample(1:dim(ImputedData)[2],100)[1:10] # submitted 10 jobs to get all 1:100


###
### containers
###

N.tot=list()
n.iter.t=dim(ImputedData)[2]
N.est=matrix(,n.iter.t,No.designs)
p.save=matrix(,n.iter.t,No.designs)
theta.save=matrix(,n.iter.t,No.designs)
kappa.save=matrix(,n.iter.t,No.designs)
gamma.save=matrix(,n.iter.t,No.designs)
scores=matrix(,length(imps),No.designs)

###
### parallelized across imputed data
###

scores=foreach(j=1:length(imps),.packages=c("raster","sp")) %dopar% { ## loop over imputed data
  ####
  set.seed(233+j)
  ####
  for(l in (1):(No.designs)){
    
    ###
    ### starting values
    ###
    
    t.ind=which(Transects.2017.l[[l]][which(BoundaryInf[]==1)]==1)
    lambda=lambda.y[,1]
    p=p.y[1]
    kappa=kappa.y[1]
    theta=theta.y[1]
    gamma=gamma.y[1]
    N=ImputedN[,1]
    p.star=p
    
    Y.imp=Transects.2017.l[[l]][which(BoundaryInf[]==1)]*ImputedData[,imps[j]]
    
    ###
    ### filter MCMC sample for each possible design
    ###
    
    for(k in 1:n.iter.t){
      
      ##
      ## sample (with replacement) which iteration to use as proposal
      ##
      
      kk=sample(1:n.iter.t,1,TRUE)
      
      ##
      ## recursive update
      ##
      
      N.star=ImputedN[,kk]*Transects.2017.l[[l]][which(BoundaryInf[]==1)]
      p.star=p.y[kk]
      mh1.tmp1=dbinom(Y.imp[t.ind],N.star[t.ind],p.star,log=TRUE)
      mh2.tmp1=dbinom(Y.imp[t.ind],N[t.ind],p,log=TRUE)
      mh1.tmp=mh1.tmp1[mh2.tmp1!=-Inf&mh1.tmp1!=-Inf]  # discard sites where Y>N
      mh2.tmp=mh2.tmp1[mh2.tmp1!=-Inf&mh1.tmp1!=-Inf]
      mh1=sum(mh1.tmp)
      mh2=sum(mh2.tmp)
      mh=exp(mh1-mh2)
      if(mh>runif(1)){
        N=N.star
        p=p.star
        lambda=lambda.y[,kk]
        kappa=kappa.y[kk]
        theta=theta.y[kk]
        gamma=gamma.y[kk]
      }
      
      if(k==1){N.tot[[l]]=lambda/n.iter.t}
      if(k>1){N.tot[[l]]=N.tot[[l]]+lambda/n.iter.t} # running mean for plotting
      ## total abudance in 2013 estimated with design l
      N.est[k,l]=sum(lambda,na.rm=TRUE)
      p.save[k,l]=p
      theta.save[k,l]=theta
      kappa.save[k,l]=kappa
      gamma.save[k,l]=gamma
      
    }
    if(l%%10==0){print(l)}
  }
  
  ###
  ### get design scores
  ###
  
  apply(N.est,2,calc.q)
  
}

save(scores,file='scores1.RData')
