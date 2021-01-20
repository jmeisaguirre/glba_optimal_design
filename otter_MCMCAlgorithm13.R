### MCMC Algorithm to estimate ecological diffusion
### of sea otters in Glacier Bay, Alaska written by
### Perry Williams 01Feb2016


run.mcmc=function(n.iter,seed,name,
                    alpha.start,beta.start,theta.start,kappa.start,
                    q.p,r.p,q.r,r.r,mu.kappa,sigma2.kappa,
                    mu.theta,sigma2.theta,mu.beta,sigma2.beta,
                    q.alpha,r.alpha,
                    beta.tune,alpha.tune,theta.tune,kappa.tune,
                    C,CISU,N,Y,X,W,cell,ISU,BoundaryNA,
                    keep,dt,timesteps,us.fact,d,
                    ho){

    set.seed(seed)

    ## Libraries
    library(raster)
    library(fields)
    library(RColorBrewer)
    library(rasterVis)
    library(gridExtra)
    library(fBasics)
    library(coda)
    library(truncnorm)
    library(inline)
    library(RcppArmadillo)
    source("HomoginizationFunctions1.R")

    ## Dimensions
    years=length(1993:2017)
    q=dim(cell)[1]*dim(cell)[2]
    Bound=rep(BoundaryNA[],years)
    burn=0

    ## Remove observations outside the park
    ## For inference
    BoundaryInf=BoundaryNA
    BoundaryInf[120:170,0:80]=NA
    BoundaryInf[144:170,0:140]=NA
    BoundaryInf[134:170,116:140]=NA
    BoundaryInf.vec=rep(BoundaryInf[],length(keep))

    ## Hold out sample
    ho.ind=matrix(0,nr=length(Y),nc=1)
    ## ho.ind[ho,]=1  # uncomment to hold out data

    ## Create indicator vector for ISU cells
    ind=1:length(N)
    ISU.ind=matrix(0,nr=length(Y),nc=1)
    ISU.ind[ISU,]=1
    ## ISU.ind[ho,]=0

    ## Identify which cells were: 1) observed to be occupied
    ##                            2) not in ISU sample
    ##                            3) not in hold-out sample
    N.occ.tmp=ind[Y==1&ISU.ind==0& ho.ind==0]
    N.occ=N.occ.tmp[!is.na(N.occ.tmp)]

    ## Identify which cells were: 1) observed to be unoccupied
    ##                            2) not in ISU sample
    ##                            3) not in hold-out sample
    N.abb.tmp=ind[Y==0&ISU.ind==0&ho.ind==0]
    N.abb=N.abb.tmp[!is.na(N.abb.tmp)]

    ## Identify which sites were surveyed and not in hold-out
    SurvSites=ind[!is.na(Y)&ho.ind==0]

    # ## C++ sampler for projection
    # code <- '
    # arma::mat Hmat = Rcpp::as<arma::mat>(H);
    # arma::mat c0mat = Rcpp::as<arma::mat>(c0);
    # int steps = Rcpp::as<int>(timesteps);
    # int n = steps;
    # int k = Hmat.n_rows;
    # arma::mat call(k,n);
    # call.col(0)=c0mat;
    # for(int i = 1; i < n; ++i){
    #   call.col(i)=Hmat*call.col(i-1);
    # }
    # return Rcpp::wrap(call);
    # '
    ## C++ sampler for projection (sparse matrix ops)
    code <- '
    arma::mat Hmat = Rcpp::as<arma::mat>(H);
    arma::mat c0mat = Rcpp::as<arma::mat>(c0);
    int steps = Rcpp::as<int>(timesteps);
    int n = steps;
    int k = Hmat.n_rows;
    arma::sp_mat temp(Hmat);
    arma::mat call(k,n);
    call.col(0)=c0mat;
    for(int i = 1; i < n; ++i){
      call.col(i)=temp*call.col(i-1);
    }
    return Rcpp::wrap(call);
  '
    calcc = cxxfunction(signature(H="numeric",c0="numeric",
                                  timesteps="numeric"),
                        body=code, plugin="RcppArmadillo")

    rtbinom=function(N,size,prob){
        qbinom(runif(N,dbinom(0,size,prob),1),size,prob)
    }

    ## Broad-scale grid cells
    us.cells=aggregate(cell,fact=us.fact,fun=mean)
    us.cells[]=1:length(us.cells[])
    us.res=res(us.cells)

    ## Fine-scale grid cells
    ds.cells=cell

    ## First-order neighborhood matrix
    NN=neighborhood(us.cells,Boundary.us)

    ## Distance matrix for initial state in km
    D=rdist(data.frame(SpatialPoints(cell)),
              matrix(d,1,2))/1000

    ## Containers
    p.save=matrix(,n.iter,1)
    theta.save=matrix(,n.iter,1)
    kappa.save=matrix(,n.iter,1)
    alpha.save=matrix(,n.iter,length(alpha.start))
    beta.save=matrix(,n.iter,length(beta.start))
    tune=c(beta.tune,alpha.tune,theta.tune,kappa.tune)
    accept=matrix(0,n.iter,length(tune))
    lam.tot=matrix(,n.iter,years)
    b0ind=b1ind=b2ind=b3ind=b4ind=0
    H.ind.save=matrix(0,nr=n.iter,length(beta.start))
    N.holdout=matrix(,nr=length(Y),nc=1)
    Y.ppd.save=matrix(,sum(!is.na(BoundaryInf[])),n.iter)
    N.save=matrix(,sum(!is.na(BoundaryInf[])),n.iter)
    lam.save=matrix(,sum(!is.na(BoundaryInf[])),n.iter)

    ## Starting Values
    alpha=alpha.start
    beta=beta.start
    theta=theta.start
    kappa=kappa.start
    delta=exp(X%*%beta)
    gamma=W%*%alpha
    ds.cells$delta=delta
    ds.cells$gamma=gamma
    delta.bar=aggregate(ds.cells$delta,fact=us.fact,
                        fun=function(x,na.rm){
                            (1/mean(1/x,na.rm=TRUE))
                        })
    gamma.bar=delta.bar*aggregate(ds.cells$gamma/ds.cells$delta,
                                  fact=us.fact,fun=function(x,na.rm){
                                      (mean(x,na.rm=TRUE))
                                  })
    H=propagator.plainZF(NN,delta.bar[],gamma.bar[],
                       dx=us.res[1],dy=us.res[2],dt=dt)
    ds.cells$lambda0=(exp(-D^2/kappa^2)/sum(exp(-D^2/kappa^2))*
                      theta)
    us.cells$c0=extract(ds.cells$delta*ds.cells$lambda0,
                        SpatialPoints(us.cells))
    c.all=brick(nrows=dim(us.cells$c0)[1],ncols=dim(us.cells$c0)[2],
                xmn=extent(cell)[1],xmx=extent(cell)[2],
                ymn=extent(cell)[3],ymx=extent(cell)[4])
    c.all=setValues(c.all,calcc(H,vec(us.cells$c0[]),timesteps)[,keep])
    lambda.all=disaggregate(c.all,us.fact)/ds.cells$delta
    lambda=vec(lambda.all[])*Bound
    k=1

    ## Begin MCMC loop
    for (k in 1:n.iter){
        if (k%%100==0)
            cat(k,"")

        ## Sample p
        p=rbeta(1,sum(CISU,na.rm =TRUE)+
                  q.p,sum(N[ISU]-CISU,na.rm=TRUE)+r.p)

        ## Sample N.tilde
        ## sites observed occupied (does not sample ISU units)
        N[N.occ]=rpois(length(N[N.occ]),lambda[N.occ]*(1-p))+C[N.occ]
        ## sites not observed occupied (does not sample ISU units)
        N[N.abb]=rpois(length(N[N.abb]),lambda[N.abb]*(1-p))

        ## Sample beta0
        beta0.star=rnorm(1,mean=beta[1],sd=tune[1])
        beta.star=c(beta0.star,beta[2:5])
        delta.star=exp(X%*%beta.star)
        ds.cells$delta.star=delta.star
        delta.bar.star=aggregate(ds.cells$delta.star,fact=us.fact,
                                 fun=function(x,na.rm){(1/mean(1/x,na.rm=TRUE))})
        gamma.bar.star=delta.bar.star*aggregate(ds.cells$gamma/ds.cells$delta.star,
                                                fact=us.fact,fun=function(x,na.rm){mean(x,na.rm=TRUE)})
        H.star=propagator.plainZF(NN,delta.bar.star[],gamma.bar.star[],
                                    dx=us.res[1],dy=us.res[2],dt=dt)
        if(min(range(H.star,na.rm=TRUE))>=0){
            us.cells$c0.star =extract(ds.cells$delta.star*ds.cells$lambda0,
                                       SpatialPoints(us.cells))
            c.all.star=setValues(c.all,calcc(H.star,vec(us.cells$c0.star[]),timesteps)[,keep])
            lambda.all.star=disaggregate(c.all.star,us.fact)/ds.cells$delta.star
            lambda.star=vec(lambda.all.star[])*Bound
            mh1=sum(dpois(N[SurvSites],lambda.star[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta.star,mu.beta,10,log=TRUE))
            mh2=sum(dpois(N[SurvSites],lambda[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta,mu.beta,10,log=TRUE))
            mh=exp(mh1-mh2)
            if(mh>runif(1)){
                beta=beta.star
                delta=delta.star
                ds.cells$delta=ds.cells$delta.star
                delta.bar=delta.bar.star
                gamma.bar=gamma.bar.star
                H=H.star
                us.cells$c0=us.cells$c0.star
                c.all=c.all.star
                lambda.all=lambda.all.star
                lambda=lambda.star
                accept[k,1]=1
            }
        }else{b0ind=b0ind+1}

        ## Sample beta1
        beta1.star=rnorm(1,beta[2],tune[2])
        beta.star=c(beta[1],beta1.star,beta[3:5])
        delta.star=exp(X%*%beta.star)
        ds.cells$delta.star=delta.star
        delta.bar.star=aggregate(ds.cells$delta.star,fact=us.fact,
                                 fun=function(x,na.rm){(1/mean(1/x,na.rm=TRUE))})
        gamma.bar.star=delta.bar.star*aggregate(ds.cells$gamma/ds.cells$delta.star,
                                                fact=us.fact,fun=function(x,na.rm){mean(x,na.rm=TRUE)})
        H.star=propagator.plainZF(NN,delta.bar.star[],gamma.bar.star[],
                                  dx=us.res[1],dy=us.res[2],dt=dt)
        if(min(range(H.star,na.rm=TRUE))>=0){
            us.cells$c0.star=extract(ds.cells$delta.star*ds.cells$lambda0,
                                       SpatialPoints(us.cells))
            c.all.star=setValues(c.all,calcc(H.star,vec(us.cells$c0.star[]),timesteps)[,keep])
            lambda.all.star=disaggregate(c.all.star,us.fact)/ds.cells$delta.star
            lambda.star=vec(lambda.all.star[])*Bound
            mh1=sum(dpois(N[SurvSites],lambda.star[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta.star,mu.beta,sigma2.beta^0.5,log=TRUE))
            mh2=sum(dpois(N[SurvSites],lambda[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta,mu.beta,sigma2.beta^0.5,log=TRUE))
            mh=exp(mh1-mh2)
            if(mh>runif(1)){
                beta=beta.star
                delta=delta.star
                ds.cells$delta=ds.cells$delta.star
                delta.bar=delta.bar.star
                gamma.bar=gamma.bar.star
                H=H.star
                us.cells$c0=us.cells$c0.star
                c.all=c.all.star
                lambda.all=lambda.all.star
                lambda=lambda.star
                accept[k,2]=1
            }
        }else{b1ind=b1ind+1}

        ## Sample beta2
        beta2.star=rnorm(1,beta[3],tune[3])
        beta.star=c(beta[1:2],beta2.star,beta[4:5])
        delta.star=exp(X%*%beta.star)
        ds.cells$delta.star=delta.star
        delta.bar.star=aggregate(ds.cells$delta.star,fact=us.fact,
                                 fun=function(x,na.rm){(1/mean(1/x, na.rm = TRUE))})
        gamma.bar.star=delta.bar.star*aggregate(ds.cells$gamma/ds.cells$delta.star,
                                                fact=us.fact,fun=function(x,na.rm){mean(x,na.rm=TRUE)})
        H.star = propagator.plainZF(NN, delta.bar.star[],gamma.bar.star[],
                                    dx=us.res[1],dy=us.res[2],dt=dt)
        if(min(range(H.star,na.rm=TRUE))>=0){
            us.cells$c0.star = extract(ds.cells$delta.star * ds.cells$lambda0,
                                       SpatialPoints(us.cells))
            c.all.star=setValues(c.all,calcc(H.star, vec(us.cells$c0.star[]),timesteps)[,keep])
            lambda.all.star=disaggregate(c.all.star, us.fact)/ds.cells$delta.star
            lambda.star=vec(lambda.all.star[])*Bound
            mh1=sum(dpois(N[SurvSites],lambda.star[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta.star, mu.beta,sigma2.beta^0.5,log=TRUE))
            mh2=sum(dpois(N[SurvSites],lambda[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta, mu.beta,sigma2.beta^0.5,log=TRUE))
            mh=exp(mh1-mh2)
            if(mh>runif(1)){
                beta=beta.star
                delta=delta.star
                ds.cells$delta=ds.cells$delta.star
                delta.bar=delta.bar.star
                gamma.bar=gamma.bar.star
                H=H.star
                us.cells$c0=us.cells$c0.star
                c.all=c.all.star
                lambda.all=lambda.all.star
                lambda=lambda.star
                accept[k,3]=1
            }
        }else{b2ind=b2ind+1}

        ## Sample beta3
        beta3.star=rnorm(1,beta[4],tune[4])
        beta.star=c(beta[1:3],beta3.star,beta[5])
        delta.star=exp(X%*%beta.star)
        ds.cells$delta.star=delta.star
        delta.bar.star=aggregate(ds.cells$delta.star,fact=us.fact,
                                 fun=function(x,na.rm){(1/mean(1/x, na.rm = TRUE))})
        gamma.bar.star=delta.bar.star*aggregate(ds.cells$gamma/ds.cells$delta.star,
                                                fact=us.fact,fun=function(x,na.rm){mean(x,na.rm=TRUE)})
        H.star = propagator.plainZF(NN, delta.bar.star[],gamma.bar.star[],
                                    dx=us.res[1],dy=us.res[2],dt=dt)
        if(min(range(H.star,na.rm=TRUE))>=0){
            us.cells$c0.star = extract(ds.cells$delta.star * ds.cells$lambda0,
                                       SpatialPoints(us.cells))
            c.all.star=setValues(c.all,calcc(H.star, vec(us.cells$c0.star[]),timesteps)[,keep])
            lambda.all.star=disaggregate(c.all.star, us.fact)/ds.cells$delta.star
            lambda.star=vec(lambda.all.star[])*Bound
            mh1=sum(dpois(N[SurvSites],lambda.star[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta.star, mu.beta,sigma2.beta^0.5,log=TRUE))
            mh2=sum(dpois(N[SurvSites],lambda[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta, mu.beta,sigma2.beta^0.5,log=TRUE))
            mh=exp(mh1-mh2)
            if(mh>runif(1)){
                beta=beta.star
                delta=delta.star
                ds.cells$delta=ds.cells$delta.star
                delta.bar=delta.bar.star
                gamma.bar=gamma.bar.star
                H=H.star
                us.cells$c0=us.cells$c0.star
                c.all=c.all.star
                lambda.all=lambda.all.star
                lambda=lambda.star
                accept[k,4]=1
            }
        }else{b3ind=b3ind+1}

        ## Sample beta4
        beta4.star=rnorm(1,beta[5],tune[5])
        beta.star=c(beta[1:4],beta4.star)
        delta.star=exp(X%*%beta.star)
        ds.cells$delta.star=delta.star
        delta.bar.star=aggregate(ds.cells$delta.star,fact=us.fact,
                                 fun=function(x,na.rm){(1/mean(1/x, na.rm = TRUE))})
        gamma.bar.star=delta.bar.star*aggregate(ds.cells$gamma/ds.cells$delta.star,
                                                fact=us.fact,fun=function(x,na.rm){mean(x,na.rm=TRUE)})
        H.star = propagator.plainZF(NN, delta.bar.star[],gamma.bar.star[],
                                    dx=us.res[1],dy=us.res[2],dt=dt)
        if(min(range(H.star,na.rm=TRUE))>=0){
            us.cells$c0.star = extract(ds.cells$delta.star * ds.cells$lambda0,
                                       SpatialPoints(us.cells))
            c.all.star=setValues(c.all,calcc(H.star, vec(us.cells$c0.star[]),timesteps)[,keep])
            lambda.all.star=disaggregate(c.all.star, us.fact)/ds.cells$delta.star
            lambda.star=vec(lambda.all.star[])*Bound
            mh1=sum(dpois(N[SurvSites],lambda.star[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta.star, mu.beta,sigma2.beta^0.5,log=TRUE))
            mh2=sum(dpois(N[SurvSites],lambda[SurvSites],log=TRUE),na.rm=TRUE)+
                sum(dnorm(beta, mu.beta,sigma2.beta^0.5,log=TRUE))
            mh=exp(mh1-mh2)
            if(mh>runif(1)){
                beta=beta.star
                delta=delta.star
                ds.cells$delta=ds.cells$delta.star
                delta.bar=delta.bar.star
                gamma.bar=gamma.bar.star
                H=H.star
                us.cells$c0=us.cells$c0.star
                c.all=c.all.star
                lambda.all=lambda.all.star
                lambda=lambda.star
                accept[k,5]=1
            }
        }else{b4ind=b4ind+1}

        ## Sample alpha
        alpha.star=rnorm(1,mean=alpha,sd=tune[6])
        gamma.star=W%*%alpha.star
        ds.cells$gamma.star=gamma.star
        gamma.bar.star=delta.bar*aggregate(ds.cells$gamma.star/ds.cells$delta,
                                                fact=us.fact,fun=function(x,na.rm){mean(x,na.rm=TRUE)})
        H.star = propagator.plainZF(NN, delta.bar[],gamma.bar.star[],
                                  dx=us.res[1],dy=us.res[2],dt=dt)
        c.all.star=setValues(c.all,calcc(H.star,vec(us.cells$c0[]),timesteps)[,keep])
        lambda.all.star=disaggregate(c.all.star,us.fact)/ds.cells$delta
        lambda.star=vec(lambda.all.star[])*Bound
        mh1=sum(dpois(N[SurvSites],lambda.star[SurvSites],log=TRUE),na.rm=TRUE)+
            sum(dbeta(alpha.star,q.alpha,r.alpha,log=TRUE))
        mh2=sum(dpois(N[SurvSites],lambda[SurvSites],log=TRUE),na.rm=TRUE)+
            sum(dbeta(alpha,q.alpha,r.alpha,log=TRUE))
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            alpha=alpha.star
            gamma=gamma.star
            ds.cells$gamma=ds.cells$gamma.star
            gamma.bar=gamma.bar.star
            H=H.star
            c.all=c.all.star
            lambda.all=lambda.all.star
            lambda=lambda.star
            accept[k,6]=1
        }

        ## Sample theta
        theta.star=rnorm(1, theta, tune[7])
        if(theta.star>0){
            ds.cells$lambda0.star = exp(-D^2/kappa^2)/sum(exp(-D^2/kappa^2))*theta.star
            us.cells$c0.star = extract(ds.cells$delta * ds.cells$lambda0.star,
                                       SpatialPoints(us.cells))
            c.all.star=setValues(c.all,calcc(H, vec(us.cells$c0.star[]),timesteps)[,keep])
            lambda.all.star=disaggregate(c.all.star, us.fact)/ds.cells$delta
            lambda.star=vec(lambda.all.star[])*Bound
            mh1=sum(dpois(N[SurvSites],lambda.star[SurvSites],log=TRUE),na.rm=TRUE)+
                dnorm(theta.star,mu.theta,sigma2.theta^0.5,log=TRUE)
            mh2=sum(dpois(N[SurvSites],lambda[SurvSites],log=TRUE),na.rm=TRUE)+
                dnorm(theta,mu.theta,sigma2.theta^0.5,log=TRUE)
            mh=exp(mh1-mh2)
            if(mh>runif(1)){
                theta=theta.star
                ds.cells$lambda0=ds.cells$lambda0.star
                us.cells$c0=us.cells$c0.star
                c.all=c.all.star
                lambda.all=lambda.all.star
                lambda=lambda.star
                accept[k,7]=1
            }
        }

        ## Sample kappa
        kappa.star = rnorm(1,kappa, tune[8])
        if(kappa.star>1){  ## 1 instead of 0 for computational stability
            ds.cells$lambda0.star = exp(-D^2/kappa.star^2)/sum(exp(-D^2/kappa.star^2))*theta
            us.cells$c0.star = extract(ds.cells$delta*ds.cells$lambda0.star,
                                       SpatialPoints(us.cells))
            c.all.star=setValues(c.all,calcc(H, vec(us.cells$c0.star[]),timesteps)[,keep])
            lambda.all.star=disaggregate(c.all.star, us.fact)/ds.cells$delta
            lambda.star=vec(lambda.all.star[])*Bound
            mh1=sum(dpois(N[SurvSites],lambda.star[SurvSites],log=TRUE),na.rm=TRUE)+
                dnorm(kappa.star,mu.kappa,sigma2.kappa^0.5,log=TRUE)
            mh2=sum(dpois(N[SurvSites],lambda[SurvSites],log=TRUE),na.rm=TRUE)+
                dnorm(kappa,mu.kappa,sigma2.kappa^0.5,log=TRUE)
            mh=exp(mh1-mh2)
            if(mh>runif(1)){
                kappa=kappa.star
                ds.cells$lambda0=ds.cells$lambda0.star
                us.cells$c0=us.cells$c0.star
                c.all=c.all.star
                lambda.all=lambda.all.star
                lambda=lambda.star
                accept[k,8]=1
            }
        }

        #future.lam=lambda[571201:595000] # 2017
        future.lam=lambda[476001:499800] # 2013
        N.sim=rpois(q,future.lam)
        Y.ppd=rbinom(q,N.sim*BoundaryInf[],p)


        ## Save Samples
        H.ind.save[k,]=c(b0ind,b1ind,b2ind,b3ind,b4ind)
        p.save[k,]=p
        beta.save[k,]=beta
        alpha.save[k,]=alpha
        theta.save[k,]=theta
        kappa.save[k,]=kappa
        lam.tot[k,]=colSums(lambda.all[]*BoundaryInf[],na.rm=TRUE)
        Y.ppd.save[,k]=Y.ppd[!is.na(Y.ppd)]
        lam.save[,k]=(future.lam*BoundaryInf[])[!is.na(Y.ppd)]
        N.save[,k]=(N.sim*BoundaryInf[])[!is.na(Y.ppd)]

        if (k%%1000==0){
            out=list(accept,
                     p.save,
                     beta.save,
                     alpha.save,
                     theta.save,
                     kappa.save,
                     lam.tot,
                     H.ind.save,
                     Y.ppd.save,
                     lam.save,
                     N.save
                     )
            save(out,file=name)
        }
        cat(" ")
    }
}
