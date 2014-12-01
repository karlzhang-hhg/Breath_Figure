setwd("C:/My Life Style/Courses of Statistics/STA 561 Probability Machine Learning/Project/R-code")

source("Evol.R")


h=5
w=5
nstep=1
lambda=100
alpha=0.1
rcr=0.5
dr=0.01
dt=0.1
nuclr=0.0001


t=0
tseq=t
coalseq=0
ncoal=0
coalall=NULL
depart.list=list()
bfobj=NULL
bfallobj=NULL
coalall=NULL

Evol(h,w,nstep,lambda,alpha,rcr,dr,dt,nuclr)
