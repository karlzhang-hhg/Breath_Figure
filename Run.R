source("Evol.R")

#setwd("/home/vis/kz29/ZKG folder/STA561/Breath_Figure")

h=1
w=1
nstep=200L
lambda=20
alpha=0.08
rcr=0.03
dr=1.05
dt=0.01
nuclr=0.0001


t=0
tseq=t
coalseq=NULL
ncoal=0L
coalall=NULL
depart.list=list()
bfobj=NULL
bfallobj=NULL
coalall=NULL
nucl.num=0

Evol(h,w,nstep,lambda,alpha,rcr,dr,dt,nuclr)
