rm(list=ls())

###
### Script used for optimal design
###

## Settings to change for multiple chains
name="MCMCOutput_recursive13_1.RData"
seed=1

library(fBasics)
source("otter_MCMCAlgorithm13.R")
load("RealData.RData")

###  Settings
n.iter=60000
burn=0
dt=1/100
timesteps=1/dt*length(1993:2017)
keep=seq(1,timesteps,1/dt)
us.fact=10

## up-scaled boundary layer surrounded by zeros
Boundary.us=aggregate(Boundary,fact=us.fact,na.rm=TRUE,fun=max)
Boundary.us[1,]=0
Boundary.us[,1]=0
Boundary.us[dim(Boundary.us)[1],]=0
Boundary.us[,dim(Boundary.us)[2]]=0
q=dim(cell)[1]*dim(cell)[2]
N=vec(N[])

## Covariates
load("DepthCov.RData")
load("DistCov.RData")
load("SlopeCov.RData")
load("ShoreCov.RData")
X=cbind(1,DepthCov[],DistCov[],SlopeCov[]*DepthCov[],ShoreCov[])
W=matrix(1,nr=q,nc=1)

## Remove large groups
N[N>65]=NA
C[C>65]=NA

## Tuning parameters
beta.tune=c(0.01667718,
            0.01984703,
            0.01796850,
            0.01312200,
            0.01403076)
alpha.tune=0.00166617
theta.tune=54.45
kappa.tune=3.1185

## Starting Values
p.start=0.75
beta.start=c(17.8060824,
             -1.3800455,
             0.4898585,
             -0.3665777,
             0.6549915)
alpha.start=0.2259298
theta.start=2592.491
kappa.start=6.446541

### Priors
q.p=1
r.p=1
q.r=1
r.r=1
mu.kappa=5
sigma2.kappa=1
mu.theta=500
sigma2.theta=250^2
mu.beta=0
sigma2.beta=1^2  ## optimal value based on cross-validation
q.alpha=1
r.alpha=1

## Cross validation hold out samples
ind=1:length(Y)
set.seed(2017)
ho.ind=round(runif(length(Y),.5,8.5))
ho=ind[ho.ind==seed]

## Run mcmc
sys.time=Sys.time()
run.mcmc(n.iter,seed,name,
         alpha.start,beta.start,theta.start,kappa.start,
         q.p,r.p,q.r,r.r,mu.kappa,sigma2.kappa,
         mu.theta,sigma2.theta,mu.beta,sigma2.beta,
         q.alpha,r.alpha,
         beta.tune,alpha.tune,theta.tune,kappa.tune,
         C,CISU,N,Y,X,W,cell,ISU,BoundaryNA,
         keep,dt,timesteps,us.fact,d,
         ho)
Sys.time()-sys.time
warnings()
