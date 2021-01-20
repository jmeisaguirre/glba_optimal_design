
load('RealData.RData')
load('scores.m.RData')

BoundaryInf=BoundaryNA
BoundaryInf[120:170,0:80]=NA
BoundaryInf[144:170,0:140]=NA
BoundaryInf[134:170,116:140]=NA
row.ind=matrix(,170,140)
for(i in 1:170){
  row.ind[i,]=ifelse(!is.na(BoundaryInf[i,]),i,NA)
}
designs=sort(row.ind[!is.na(row.ind)])

### proposed designs
No.designs=1000
Transects.2017.l=list()
NumOptTransects=20
set.seed(2017)
for(l in 1:No.designs){
  design.tmp.opt=sort(sample(
    unique(designs),NumOptTransects))
  t.m=matrix(,170,140)
  t.m[design.tmp.opt,]=1
  Transects.2017.l[[l]]=c(t(t.m))
}

sd=scores.m
hist(sd,breaks=No.designs)
plot(sd,MCMC.error)
summary(lm(sd~MCMC.error))
which(sd==min(sd,na.rm=TRUE))

## 561
foo=matrix(Transects.2017.l[[561]],
           170,140,byrow=TRUE)
OptimalTransects=which(foo[,1]==1)  ## optimal transects
