######################################################################
#Physical parameters:
#lambda: nucleation rate or average nucleates per area
#alpha: condensation rate or the constant factor of radius 
#       as function of time (assume constant incoming volume rate)
#h: height of edge of simulation
#w: width of edge of simulation 
#rcr: critical radius for droplets depart
#dr: repulsive distance for nucleation
#nuclr: nucleation radius (nuclr<dr)
#theta: contact angle
#dt: default time step for evolving

######################################################################
#Relation between physical parameters:
#The alpha and c(lambda,nuclr,dt) should have some constriction, 
#   because the growth rate by nucleation or condensation should
#   be governed by some physical laws
#The rcr should much smaller than c(h,w)
#The nuclr should be a very small number
#The dr should smaller than rcr but larger than nuclr

######################################################################
# Functions:
# Crea_drop<-function(ind,t,x,y,inir)
# Dist<-function(p1,p2,h,w)
# Nucl<-function(t,lambda, h, w, dr, nuclr,search.pool)
# Update_coal<-function(ncoal,t,small,big,edge.circ,ds,h,w)
# Grow_Chk_coal_depa<-function(alpha,ncoal,t,dt,h,w,min_dist,rcr)
# Neg_dist<-function(t,ncoal,h,w,depart.list,search.pool)
# Min_dist<-function(t,min_dist,search.pool.min)
# Plot_bf<-function(i,t,h,w)
# Evol<-function(h,w,nstep,lambda,alpha,rcr,dr,dt,nuclr)

######################################################################
# Structure of node, bfobj, and bfallobj:
# node1=list(birt.t=double(),#birth time
#            deat.t=double(),#death time
#            r=double(),#current radius
#            posi=as.vector(double()),#center position
#            log.tr=as.data.frame(as.integer()),#log of time and radius
#            abso.drop=as.vector(),#indices of droplets absorbed to generate this droplet
#            coal=list(coal1=list(ind=as.integer(),#index of coalescence
#                                 coal.t=as.double(),#time of coalescence
#                                 delt.s=as.double(),#drifting distance of coalescence
#                                 invo.drop=as.vector(),#indices of droplets involved in this coalescence
#            ))
# )
# 
# bfallobj<-list(node1=list(birt.t=double(),#birth time
#                           deat.t=double(),#death time
#                           r=double(),#current radius
#                           posi=as.vector(double()),#center position
#                           log.tr=as.data.frame(),#log of time and radius
#                           abso.drop=as.vector(),#indices of droplets absorbed to generate this droplet
#                           coal=list(coal1=list(ind=as.integer(),#index of coalescence
#                                                coal.t=as.double(),#time of coalescence
#                                                delt.s=as.double(),#drifting distance of coalescence
#                                                invo.drop=as.vector(),#indices of droplets involved in this coalescence
#                                                ))
#                           )
#               )
# 
# bfobj<-list(node1=list(   birt.t=double(),#birth time
#                           deat.t=double(),#death time
#                           r=double(),#current radius
#                           posi=as.vector(double()),#center position
#                           log.tr=as.data.frame(),#log of time and radius
#                           abso.drop=as.vector(),#indices of droplets absorbed to generate this droplet
#                           coal=list(coal1=list(ind=as.integer(),#index of coalescence
#                                                coal.t=as.double(),#time of coalescence
#                                                delt.s=as.double(),#drifting distance of coalescence
#                                                invo.drop=as.vector(),#indices of droplets involved in this coalescence
#                           ))
# )
# )


library(spatstat)
library(plotrix)
#project2set
#rHardcore
#summary.ppp
#nndist

######################################################################
#create a element of bfobj or bfallobj as a list
Crea_drop<-function(ind,t,x,y,inir)
{
  stopifnot(is.integer(ind))
  node=list( ind=ind, #index of droplet
             birt.t=t,#birth time
             deat.t=-1,#death time
             jump.t=-1,#jump time
             r=inir,#current radius
             posi=c(x,y),#center position
             log.tr=as.data.frame("t"=t,"x"=x,"y"=y,"inir"=inir),#log of time, position, and radius
             abso.drop=as.vector(as.integer()),#indices of droplets absorbed to generate this droplet
             coal=list(coal1=list(ind=as.integer(),#index of coalescence
                                  coal.t=as.double(),#time of coalescence
                                  delt.s=as.double(),#drifting distance of coalescence
                                  invo.drop=as.vector(as.double())#indices of droplets involved in this coalescence
             ))
  )
  print(node)
  return(node)           
}


######################################################################
#Check distance of two points in a rectangle of w*l
#The distance should be periodically defined: meaning
#if the point near the opposite edges of the window 
#should also be checked in terms of center distance
#By using different Dist function, we can obmit the edge effect. The distribution would be different.
Dist<-function(p1,p2,h,w)
{
  dist0<-sqrt(sum((p1-p2)^2))
  distw<-sqrt(sum((c(w,0)-(p1-p2))^2))
  disth<-sqrt(sum((c(0,h)-(p1-p2))^2))
  if (dist0<distw & dist0<disth) return(c(0,0,dist0))
  else
  {
    if (distw<disth) return(c(1,0,distw))
    else return(c(0,1,disth))
  }
}

######################################################################
#Nucleation function: Generate nucleation and 
#check if nucleated droplets fall in repulsive 
#ranges of previous droplets
#the search.pool as input is used for fast calculating min_dist
#?what does stationary distribution mean in the function?
Nucl<-function(t,lambda, h, w, dr, nuclr,search.pool)
{
  #bfobj is the list of points currently existing
  #ndropall is the number of droplets in bfallobj
  #ndrop is the number of droplets in bfobj
  ndrop<-length(bfobj)
  ndropall<-length(bfallobj)
  #this can also change to Poisson process
  nucl<-rMaternII(lambda,dr,win=owin(c(0,w),c(0,h)),stationary = T)
  plot(nucl)
  if (length(bfobj)==0)
  {
    ndrop=ndrop+1L
    ndropall=ndropall+1L
    search.pool=c(search.pool,ndropall)
    #The name of Crea_drop wouldn't be passed, so that first asign element; then change name
    bfobj<<-c(bfobj,Crea_drop(ndropall,t,nucl$x[1],nucl$y[1],nuclr))
    names(bfobj)[ndrop]<<-as.character(ndropall)
  }
  for (i in 2:nucl$n)
  {
    ndrop<-length(bfobj)
    flag.thin=0 #thin or not
    for (j in 1:ndrop)
    {
      if (Dist(c(nucl$x[i],nucl$y[i]),bfobj[j]$posi,h,w)[3]<=nuclr+dr+bfobj[j]$r)
      {
        flag.thin=1
        break
      }
    }
    if (flag.thin==0)
    {
      ndrop=ndrop+1L
      ndropall=ndropall+1L
      search.pool=c(search.pool,ndropall)
      #The name of Crea_drop wouldn't be passed, so that first asign element; then change name
      bfobj[ndrop]<<-Crea_drop(ndropall,t,nucl$x[i],nucl$y[i],nuclr)
      names(bfobj)[ndrop]<<-as.character(ndropall)
    }
  }
  return(search.pool)
}

######################################################################
#Update information after one coalescence between two droplets
#Delete the smaller droplet
#Update information for the big droplet in bfobj
#Updata information in bfallobj
Update_coal<-function(ncoal,t,small,big,edge.circ,ds,h,w)
{
  stopifnot(is.integer(ncoal))
  stopifnot(is.integer(small))
  stopifnot(is.integer(big))
  #convert the index number to character
  small=as.character(small)
  big=as.character(big)
  #update information for the small droplet after a coalescence
  bfobj[small]$deat.t<<-t
  bfobj[small]$log.tr<<-rbind(bfobj[small]$log.tr,c(t,bfobj[small]$r))
  coalall[ncoal]<<-list(ind=ncoal,t=t,invo.drop=rbind(invo.drop,c(i,j)),delt.s=ds)#invo.drop initially is NULL
  names(coalall)[ncoal]<<-as.character(ncoal)
  bfobj[small]$coal[length(bfobj[small]$coal)+1L]<<-coalall[ncoal]
  #update information for the big droplet after a coalescence
  bfobj[big]$r<<-((bfobj[small]$r+bfobj[big]$r)^3)^(1/3)
  #Consider the circlic effect
  if (edge.circ[1]==0 & edge.circ[2]==0)
  {
    #mass weighted average center position as the new position
    bfobj[big]$posi<<-bfobj[big]$posi+(bfobj[small]$posi-bfobj[big]$posi)*bfobj[small]$r^3/(bfobj[small]$r^3+bfobj[big]$r^3)
  }
  else
  {
    if (edge.circ[1]==1)
    {
      bfobj[big]$posi<<-(bfobj[big]$posi+(bfobj[small]$posi-bfobj[big]$posi-c(w,0)*sign(bfobj[small]$posi[1]-bfobj[big]$posi[1]))*bfobj[small]$r^3/(bfobj[small]$r^3+bfobj[big]$r^3))%%w
    }
    else #edge.circ[2]==1
    {
      bfobj[big]$posi<<-(bfobj[big]$posi+(bfobj[small]$posi-bfobj[big]$posi-c(0,h)*sign(bfobj[small]$posi[2]-bfobj[big]$posi[2]))*bfobj[small]$r^3/(bfobj[small]$r^3+bfobj[big]$r^3))%%h
    }
  }
  bfobj[big]$log.tr<<-rbind(bfobj[big]$log.tr,t=t,x=bfobj[big]$posi[1],y=bfobj[big]$posi[2],bfobj[big]$r)
  bfobj[big]$coal[length(bfobj[big]$coal)+1L]<<-coalall[ncoal]
#   bfallobj[names(bfobj)]<<-bfobj
#   bfobj<<-bfobj[!(names(bfobj) %in% as.character(delete.list))]
}

######################################################################
#Check if there is any coalescence should have happened
#in time interval [t,t+dt]: if does, return a new dt, the
#time elapsed until the first coalescence in this time
#interval, since the last calculated time t.
#First grow; then check coalescence; then check departure. 
#Update the bfobj (growth & coalescence).
Grow_Chk_coal_depa<-function(alpha,ncoal,t,dt,h,w,min_dist,rcr)
{
  ncoal.old=ncoal
  ######################################################################
  #for potential coalescence happening in time interval [t,t+dt]
  if (min_dist[1,5]>2*alpha*dt) #there is no coalescence in time interval [t,t+dt]
  {
    for (i in 1:length(bfobj))
    {
      #Note that condensation couldn't result in departure
      bfobj[i]$r<<-bfobj[i]$r+alpha*dt
      bfobj[i]$log.rt<<-rbind(bfobj[i]$log.rt,c(t+dt,bfobj[i]$r))
    }
    for (i in 1:nrow(min_dist))
    {
      min_dist[i,5]=min_dist[i,5]-2*alpha*dt
    }
    search.pool=as.vector(min_dist[,1:2])
    return(list(num_coal=ncoal-ncoal.old,dt,min_dist,search.pool))
  }
  else #there is some coalescence in time interval [t,t+dt]; update dt to dtp; update droplet radius to time t+dtp
  {
    search.pool.new=NULL #only droplets coalesced need to be checked (put into search.pool.new)
    delete.list=NULL #those absorbed droplets need to be deleted from bfobj
    depart.list=NULL #if, after coalesecence, radius is larger than rcr, the droplet departs from the substrate
    #Because every droplet can only appear once in min_dist, so that I don't need worry about
    #whether previous coalescence would result in the case of unphysical later coalescence
    #where droplet-pair distance is larger than zero
    for (i in 1:nrow(min_dist))
    {      
      ncoal=ncoal+1L
      dtp=min_dist[1,5]/2/alpha
      bfobj[as.character(min_dist[i,1])]$r<<-bfobj[as.character(min_dist[i,1])]$r+alpha*dtp
      bfobj[as.character(min_dist[i,2])]$r<<-bfobj[as.character(min_dist[i,2])]$r+alpha*dtp
      #The larger drop before coalescence stay. The smaller one is regarded as absorbed
      if (bfobj[as.character(min_dist[i,1])]$r<bfobj[as.character(min_dist[i,2])]$r)
      {
        small=min_dist[i,1]
        big=min_dist[i,2]         
      }
      else
      {
        small=min_dist[i,2]
        big=min_dist[i,1]
      }
      #update information for both droplets
      Update_coal(ncoal,t=t+dtp,small,big,min_dist[i,3:4],ds = 0,h,w)
      #Check if the radius after coalescence is larger than rcr
      #If does, the drop jump away
      #Otherwise, it stays
      if (bfobj[as.character(big)]$r>rcr)
      {
        depart.list=c(depart.list,bfobj[as.character(big)]$ind) #included those departed
        search.pool.new=c(search.pool.new,big)
        delete.list=c(delete.list,c(small)) #included those absorbed
      }
      else
      {
        search.pool.new=c(search.pool.new,big)
        delete.list=c(delete.list,small)
      }
    }
    bfallobj[names(bfobj)]<<-bfobj
    bfobj<<-bfobj[!(names(bfobj) %in% as.character(delete.list))]
    ######################################################################
    #Check if there is still some coalescence because of previous drop-pair coalescence
    #need to happen using Neg_dist
    search.pool=search.pool.new
    Neg_dist_info=Neg_dist(t+dtp,ncoal,h,w,search.pool)
    ncoal=ncoal+Neg_dist_info[1]
    depart.list=Neg_dist_info[2]
    ######################################################################
    #Find min_dist after coalescence (a new breath figure)
    #Once coalescence happen, min_dist should be calculated over all existing droplets
    min_dist0=matrix(c(-1,-1,0,0,max(h,w)),nrow=1)
    min_dist=min_dist0
    Min_dist_info<-Min_dist(t,min_dist,as.integer(names(bfobj)))
    min_dist=Min_dist_info[1]
    #search.pool.min is a search pool for finding minimum gap in droplet pattern
    search.pool.min=Min_dist_info[2]
    ######################################################################
    return(list(num_coal=ncoal-ncoal.old,dtp,min_dist,depart.list,search.pool.min))
  }
}

######################################################################
#Find neg_dist at time t. neg_dist indicate coalescence
#Update bfobj if there is some coalescence
#Stop when no coalescence happens at time t
Neg_dist<-function(t,ncoal,h,w,depart.list,search.pool)
{
  ncoal.old=ncoal
  repeat{
    ndrop=length(bfobj)
    #only those coalesced droplets need to be checked in the next loop (search.pool.new)
    search.pool.new=NULL
    neg_dist=NULL
    #check if there is negative distance
    for (i in 1:ndrop)
    {
      for (j in search.pool[as.character(search.pool) %in% names(bfobj)(i+1L:ndrop)])
      {
        #if the two droplets for distance calculation are the same, move to the next iteration
        if (bfobj[i]$ind==bfobj[as.character(j)]$ind) {next}
        cal_dist=Dist(bfobj[i]$posi,bfobj[as.character(j)]$posi,h,w)
        gap=cal_dist[3]-bfobj[i]$r-bfobj[as.character(j)]$r
        if(gap<=0)
        {
          neg_dist=rbind(neg_dist,c(bfobj[i]$ind,bfobj[as.character(j)]$ind,cal_dist))
        }
      }
    }
    
    if (is.null(neg_dist)) #a valid breath figure; nothing needs to change
    {
      break
    }
    else #some coalescence need to happen
    {
      ######################################################################
      #Sort the neg_dist so that droplets can only appear once in search.pool.neg or neg_dist
      search.pool.neg=as.vector(neg_dist[1,1:2])
      delete.pair=NULL
      for (i in 2:nrow(neg_dist))
      {
        if (sum(neg_dist[i,1:2] %in% search.pool.neg))
        {
          delete.pair=c(delete.pair,i)
        }
        else
        {
          search.pool.neg=c(search.pool.neg,as.vector(neg_dist[i,1:2]))
        }
      }
      if (!is.null(delete.pair))
      {
        neg_dist=neg_dist[-delete.pair,]
      }
      ######################################################################
      delete.list=NULL
      for (i in 1:nrow(neg_dist))
      {
        ncoal=ncoal+1L
        if (bfobj[as.character(neg_dist[i,1])]$r<bfobj[as.character(neg_dist[i,2])]$r)
        {
          small=neg_dist[i,1]
          big=neg_dist[i,2]         
        }
        else
        {
          small=neg_dist[i,2]
          big=neg_dist[i,1]
        }
        #update information for both droplets
        Update_coal(ncoal,t=t+dtp,small,big,neg_dist[i,3:4],ds = 0,h,w)
        #Check if the radius after coalescence is larger than rcr
        #If does, the drop jump away
        #Otherwise, it stays
        if (bfobj[as.character(big)]$r>rcr)
        {
          #Only include the droplet jump at the last time
          if (!((small %in% depart.list)|(big %in% depart.list)))
          {
            depart.list=c(depart.list,big) #included those departed
          }
          else
          {
            if (small %in% depart.list)
            {
              stopifnot(big %in% depart.list)
              depart.list=depart.list[depart.list[]!=small]
            }
            else
            {
              if(!(big %in% depart.list))
              {
                depart.list=c(depart.list,big)
              }
            }
          }
        }
        search.pool.new=c(search.pool.new,big)
        delete.list=c(delete.list,c(small)) #included those absorbed
      }
      search.pool=search.pool.new
      stopifnot(!is.null(delete.list))
      bfallobj[names(bfobj)]<<-bfobj
      bfobj<<-bfobj[!(names(bfobj) %in% as.character(delete.list))]
    }
  }
  return(list(num_coal=ncoal-ncoal.old,depart.list))
}

######################################################################
#Find min_dist after nucleation
Min_dist<-function(t,min_dist,search.pool.min)
{
  ndrop=length(bfobj)
  #find the minimum non-negative distance among the droplet pattern
  for (i in 1:ndrop)
  {
    for (j in search.pool.min[as.character(search.pool.min) %in% names(bfobj)(i+1L:ndrop)])
    {
      if (bfobj[i]$ind==bfobj[as.character(j)]$ind) {next}
      cal_dist=Dist(bfobj[i]$posi,bfobj[as.character(j)]$posi,h,w)
      gap=cal_dist[3]-bfobj[i]$r-bfobj[as.character(j)]$r
      stopifnot(gap>=0)
      if (gap<min_dist[1,5])
      {
        min_dist=matrix(c(bfobj[i]$ind,bfobj[as.character(j)]$ind,cal.dist[1:2],gap),nrow=1)
      }
      else
      {
        if (gap==min_dist[1,5])
        {
          min_dist=rbind(min_dist,c(bfobj[i]$ind,bfobj[as.character(j)]$ind,cal.dist[1:2],gap))
        }
      }
    }
  }
  #Sort the min_dist so that droplets can only appear once in search.pool.min or min_dist
  search.pool.min=as.vector(min_dist[1,1:2])
  delete.pair=NULL
  for (i in 2:nrow(min_dist))
  {
    if (sum(min_dist[i,1:2] %in% search.pool.min))
    {
      delete.pair=c(delete.pair,i)
    }
    else
    {
      search.pool.min=c(search.pool.min,as.vector(min_dist[i,1:2]))
    }
  }
  if (!is.null(delete.pair))
  {
    min_dist=min_dist[-delete.pair,]
  }
  search.pool.min=as.vector(min_dist[,1:2])
  return(list(min_dist,search.pool.min))
}

######################################################################
#Plot the frame of current time t
Plot_bf<-function(i,t,h,w)
{
  setwd("C:/My Life Style/Courses of Statistics/STA 561 Probability Machine Learning/Project/R-code/image")
  pdf(paste0(i,"-",t,".pdf"))
  plot.new()
  plot.window(xlim=c(0,w),ylim=c(0,h),xaxs="i",yaxs="i")
  box()
  title(mian=paste0(i,"-",t))
  drop.info=NULL
  for (i in 1:length(bfobj))
  {
    drop.info.i=c(bfobj[i]$ind,bfobj[i]$posi,bfobj[i]$r)
    #if droplet exceed left and right boundaries
    if (bfobj[i]$posi[1]-bfobj[i]$r<0)
    {
      drop.info.i=rbind(drop.info.i,c(bfobj[i]$ind,bfobj[i]$posi+c(w,0),bfobj[i]$r))
    }
    else
    {
      if (bfobj[i]$posi[1]+bfobj[i]$r>w)
      {
        drop.info.i=rbind(drop.info.i,c(bfobj[i]$ind,bfobj[i]$posi+c(-w,0),bfobj[i]$r))
      }
    }
    
    if (bfobj[i]$posi[2]-bfobj[i]$r<0)
    {
      drop.info.i=rbind(drop.info.i,c(bfobj[i]$ind,bfobj[i]$posi+c(0,h),bfobj[i]$r))
    }
    else
    {
      if (bfobj[i]$posi[2]+bfobj[i]$r>h)
      {
        drop.info.i=rbind(drop.info.i,c(bfobj[i]$ind,bfobj[i]$posi+c(0,-h),bfobj[i]$r))
      }
    }
    
    drop.info=rbind(drop.info,drop.info.i)
  }
  draw.circle(drop.info[,2],drop.info[,3],drop.info[,4],
              border = "blue",
              col = color.scale(drop.info[,2],0,0,1,alpha = 0.5))
  dev.off()
}

######################################################################
#Evolve in between two calculated time
#nstep: total step of evolution
Evol<-function(h,w,nstep,lambda,alpha,rcr,dr,dt,nuclr)
{
  ######################################################################
  #Nucleation (routine)
  search.pool.min=NULL
  search.pool.min=Nucl(t,lambda,h,w,dr,nuclr,search.pool.min)
  #ncoal=ncoal+Neg_dist(t,ncoal,h,w,search.pool)
  ######################################################################
  #Find initial min_dist
  min_dist0=matrix(c(-1,-1,0,0,max(h,w)),nrow=1)
  min_dist=min_dist0
  Min_dist_info<-Min_dist(t,min_dist,search.pool.min)
  min_dist=Min_dist_info[1]
  #search.pool.min=Min_dist_info[2]
  ######################################################################
  bfallobj[names(bfobj)]<<-bfobj
  for (i in 1:nstep)
  {
    ######################################################################
    #Grow (condense) and check if there is coalescence happen between [t,t+dtp]
    coalinfo<-Grow_Chk_coal_depa(alpha,ncoal,t,dt,h,w,min_dist,rcr)
    coalseq<<-c(coalseq,coalinfo[1])
    ncoal<<-ncoal+coalinfo[1]
    dtp=coalinfo[2]
    min_dist=coalinfo[3]
    depart.list<<-c(depart.list,coalinfo[4])
    search.pool.min=coalinfo[5]
    ######################################################################
    #Nucleation (routine)
    t<<-t+dtp
    tseq<<-c(tseq,t)
    search.pool.min=Nucl(t,lambda,h,w,t,dr,nuclr,search.pool.min)
    ######################################################################
    #Find min_dist after nucleation for preparation of next round Grow_Chk_coal
    #only need to go through those droplets in preivous min_dist and newly nucleated
    #but the minimum distance may be higher, so that should start with min_dist0
    min_dist0=matrix(c(-1,-1,0,0,max(h,w)),nrow=1)
    min_dist=min_dist0
    Min_dist_info<-Min_dist(t,min_dist,search.pool.min)
    min_dist=Min_dist_info[1]
    #search.pool.min=Min_dist_info[2]
    ######################################################################
    #store the updtated information into bfallobj
    bfallobj[names(bfobj)]<<-bfobj
    Plot_bf(i,t,h,w)
  }
}


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
