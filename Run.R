source("Evol.R")

#setwd("/home/vis/kz29/ZKG folder/STA561/Breath_Figure")

h=1
w=1
nstep=4000L
lambda=30
alpha=0.15
rcr=0.03
dr=0.001
dt=0.01
nuclr=0.0001
#the upper ratio for jumping happen
#in order for jumping happening, the merged drop size should smaller
#than up_ratio times of smaller droplet
up_ratio=3

phypar=list(height=h,width=w,num_step=nstep,nucleation_density=lambda,
            condensing_rate=alpha,critical_radius=rcr,repulsive_length=dr,
            time_step=dt,nucleate_radius=nuclr,upper_ratio=up_ratio)

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

saveim=0

Evol(h,w,nstep,lambda,alpha,rcr,dr,dt,nuclr,saveim)

old.wd=getwd()
setwd(paste0(old.wd,"/Dat/Run8"))
save(phypar,file="phypar")
save(bfallobj,file="bfallobj")
save(bfobj,file="bfobj")
save(tseq,file="tseq")
save(depart.list,file="depart.list")
save(coalall,file="coalall")
save(coalseq,file="coalseq")
save(nucl.num,file="nucl.num")
setwd(old.wd)