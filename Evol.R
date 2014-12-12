source("Nucl.R")
source("Min_dist.R")
source("Grow_Chk_coal_depa.R")
source("Plot_bf.R")

######################################################################
#Evolve in between two calculated time
#nstep: total step of evolution
Evol<-function(h,w,nstep,lambda,alpha,rcr,dr,dt,nuclr,saveim)
{
  ######################################################################
  #Nucleation (routine)
  search.pool.min=NULL
  search.pool.min=Nucl(t,lambda,h,w,dr,nuclr,search.pool.min)
  Plot_bf(0,t,h,w,dt)
  #ncoal=ncoal+Neg_dist(t,ncoal,h,w,search.pool)
  ######################################################################
  #Find initial min_dist
  min_dist0=matrix(c(-1,-1,0,0,max(h,w)),nrow=1)
  min_dist=min_dist0
  Min_dist_info<-Min_dist(t,min_dist,search.pool.min)
  min_dist=Min_dist_info[[1]]
  #search.pool.min=Min_dist_info[2]
  ######################################################################
  bfallobj[names(bfobj)]<<-bfobj
  ######################################################################
  for (i in 1:nstep)
  {
    ######################################################################
    #Grow (condense) and check if there is coalescence happen between [t,t+dtp]
    coalinfo<-Grow_Chk_coal_depa(alpha,ncoal,t,dt,h,w,min_dist,rcr,search.pool.min)
    print(1)
    #print(c("coalinfo:",coalinfo))
    coalseq<<-c(coalseq,coalinfo[[1]])
    ncoal<<-ncoal+coalinfo[[1]]
    dtp=coalinfo[[2]]
    min_dist=coalinfo[[3]]
    depart.list<<-c(depart.list,coalinfo[[4]])
    print(c("jump:",coalinfo[[4]]))
    search.pool.min=coalinfo[[5]]
    ######################################################################
    #Nucleation (routine)
    t<<-t+dtp
    tseq<<-c(tseq,t)
    if (t>=nucl.num*dt)
    {
      search.pool.min=Nucl(t,lambda,h,w,dr,nuclr,search.pool.min)
      nucl.num<<-nucl.num+1
    }
    ######################################################################
    #Find min_dist after nucleation for preparation of next round Grow_Chk_coal
    #only need to go through those droplets in preivous min_dist and newly nucleated
    #but the minimum distance may be higher, so that should start with min_dist0
    min_dist0=matrix(c(-1,-1,0,0,max(h,w)),nrow=1)
    min_dist=min_dist0
    Min_dist_info<-Min_dist(t,min_dist,search.pool.min)
    min_dist=Min_dist_info[[1]]
    print(2)
    search.pool.min=Min_dist_info[[2]]
    ######################################################################
    #store the updtated information into bfallobj
    bfallobj[names(bfobj)]<<-bfobj
    if (saveim==1)
    {
      Plot_bf(i,t,h,w,dt)
    }
    print(c("i&nucl.num&t:",i,nucl.num,t))
  }
}