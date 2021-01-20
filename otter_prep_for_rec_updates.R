
###
### gather MCMC output
###

load("MCMCOutput_recursive13_1.RData")
chain1=out
load("MCMCOutput_recursive13_2.RData")
chain2=out
load("MCMCOutput_recursive13_3.RData")
chain3=out
out=NULL
gc()


###
### extract forecast for 2013
###

burn=10001
n.iter=60000
lambda=cbind(chain1[[10]][,burn:n.iter],chain2[[10]][,burn:n.iter],chain3[[10]][,burn:n.iter])
ImputedData=cbind(chain1[[9]][,burn:n.iter],chain2[[9]][,burn:n.iter],chain3[[9]][,burn:n.iter])
ImputedN=cbind(chain1[[11]][,burn:n.iter],chain2[[11]][,burn:n.iter],chain3[[11]][,burn:n.iter])
p=c(chain1[[2]][burn:n.iter],chain2[[2]][burn:n.iter],chain3[[2]][burn:n.iter])
theta=c(chain1[[5]][burn:n.iter],chain2[[5]][burn:n.iter],chain3[[5]][burn:n.iter])
kappa=c(chain1[[6]][burn:n.iter],chain2[[6]][burn:n.iter],chain3[[6]][burn:n.iter])
gamma=c(chain1[[4]][burn:n.iter],chain2[[4]][burn:n.iter],chain3[[4]][burn:n.iter])

gc()
dim(ImputedData)


###
### generate surveys
###

##
## Set up for optimal design in 2013
##

load('RealData.RData')

BoundaryInf=BoundaryNA
BoundaryInf[120:170,0:80]=NA
BoundaryInf[144:170,0:140]=NA
BoundaryInf[134:170,116:140]=NA

row.ind=matrix(,170,140)
for(i in 1:170){
  row.ind[i,]=ifelse(!is.na(BoundaryInf[i,]),i,NA)
}
designs=sort(row.ind[!is.na(row.ind)])

## proposed optimal transects
No.designs=3000
Transects.2017.l=list()  ## actually 2013
NumOptTransects=5
set.seed(2017)
for(l in 1:No.designs){
  design.tmp.opt=sort(sample(
    unique(designs),NumOptTransects))
  t.m=matrix(,170,140)
  t.m[design.tmp.opt,]=1
  Transects.2017.l[[l]]=c(t(t.m))
}

n.thin=3

lambda.y=lambda[,seq(1,dim(lambda)[2],n.thin)]

p.y=p[seq(1,length(p),n.thin)]

theta.y=theta[seq(1,length(p),n.thin)]
kappa.y=kappa[seq(1,length(p),n.thin)]
gamma.y=gamma[seq(1,length(p),n.thin)]

ImputedData=ImputedData[,seq(1,dim(lambda)[2],n.thin)]
ImputedN=ImputedN[,seq(1,dim(lambda)[2],n.thin)]

imp=NULL
chains=NULL
N.ppd=NULL
lambda=NULL
chain1=NULL
chain2=NULL
chain3=NULL
p=NULL
gc()
save.image('rec_work_re_y_13.RData')


