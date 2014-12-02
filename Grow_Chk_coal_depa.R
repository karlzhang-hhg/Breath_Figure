source("Update_coal.R")
source("Neg_dist.R")
source("Min_dist.R")

######################################################################
#Check if there is any coalescence should have happened
#in time interval [t,t+dt]: if does, return a new dt, the
#time elapsed until the first coalescence in this time
#interval, since the last calculated time t.
#First grow; then check coalescence; then check departure. 
#Update the bfobj (growth & coalescence).
Grow_Chk_coal_depa<-function(alpha,ncoal,t,dt,h,w,min_dist,rcr,search.pool.min)
{
  ncoal.old=ncoal
  depart.list=NULL #if, after coalesecence, radius is larger than rcr, the droplet departs from the substrate
  ######################################################################
  #for potential coalescence happening in time interval [t,t+dt]
  print(c("min_dist:",min_dist,is.matrix(min_dist)))
  if (min_dist[1,5]>2*alpha*dt) #there is no coalescence in time interval [t,t+dt]
  {
    print("large separation")
    for (i in 1:length(bfobj))
    {
      #Note that condensation couldn't result in departure
      bfobj[[i]]$r<<-bfobj[[i]]$r+alpha*dt
      #print(c(bfobj[[i]]$ind,bfobj[[i]]$r))
      bfobj[[i]]$log.tr<<-rbind(bfobj[[i]]$log.tr,c(t+dt,bfobj[[i]]$posi,bfobj[[i]]$r))
    }
    for (i in 1:nrow(min_dist))
    {
      min_dist[i,5]=min_dist[i,5]-2*alpha*dt
    }
    search.pool=as.vector(min_dist[,1:2])
    return(list(num_coal=ncoal-ncoal.old,dt,min_dist,depart.list,search.pool))
  }
  else #there is some coalescence in time interval [t,t+dt]; update dt to dtp; update droplet radius to time t+dtp
  {
    #First add the radius
    dtp=min_dist[1,5]/2/alpha
    #print(c("dtp:",dtp))
    for (i in 1:length(bfobj))
    {
      #Note that condensation couldn't result in departure
      bfobj[[i]]$r<<-bfobj[[i]]$r+alpha*dtp
      #print(c(bfobj[[i]]$ind,bfobj[[i]]$r))
      bfobj[[i]]$log.tr<<-rbind(bfobj[[i]]$log.tr,c(t+dtp,bfobj[[i]]$posi,bfobj[[i]]$r))
    }
    for (i in 1:nrow(min_dist))
    {
      min_dist[i,5]=min_dist[i,5]-2*alpha*dtp
    }
    #Then, deal with coalescence
    search.pool.new=NULL #only droplets coalesced need to be checked (put into search.pool.new)
    delete.list=NULL #those absorbed droplets need to be deleted from bfobj
    
    #Because every droplet can only appear once in min_dist, so that I don't need worry about
    #whether previous coalescence would result in the case of unphysical later coalescence
    #where droplet-pair distance is larger than zero
    for (i in 1:nrow(min_dist))
    {      
      ncoal=ncoal+1L
      #The larger drop before coalescence stay. The smaller one is regarded as absorbed
      if (bfobj[[as.character(min_dist[i,1])]]$r<bfobj[[as.character(min_dist[i,2])]]$r)
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
      if (bfobj[[as.character(big)]]$r>rcr)
      {
        depart.list=c(depart.list,big) #included those departed
        bfobj[[as.character(big)]]&jump.t<<-t+dtp
        search.pool.new=c(search.pool.new,big)
        delete.list=c(delete.list,small) #included those absorbed
      }
      else
      {
        search.pool.new=c(search.pool.new,big)
        delete.list=c(delete.list,small)
      }
    }
    stopifnot(!is.null(delete.list))
    bfallobj[names(bfobj)]<<-bfobj
    bfobj<<-bfobj[!(names(bfobj) %in% as.character(delete.list))]
    ######################################################################
    #Check if there is still some coalescence because of previous drop-pair coalescence
    #need to happen using Neg_dist
    search.pool=c(search.pool.new,search.pool.min)
    Neg_dist_info=Neg_dist(t+dtp,ncoal,h,w,depart.list,search.pool)
    ncoal=ncoal+Neg_dist_info[[1]]
    depart.list=Neg_dist_info[[2]]
    ######################################################################
    #Find min_dist after coalescence (a new breath figure)
    #Once coalescence happen, min_dist should be calculated over all existing droplets
    min_dist0=matrix(c(-1,-1,0,0,max(h,w)),nrow=1)
    min_dist=min_dist0
    Min_dist_info<-Min_dist(t+dtp,min_dist,as.integer(names(bfobj)))
    min_dist=Min_dist_info[[1]]
    #search.pool.min is a search pool for finding minimum gap in droplet pattern
    search.pool.min=Min_dist_info[[2]]
#     ######################################################################
#     Neg_dist_info=Neg_dist(t+dtp,ncoal,h,w,depart.list,search.pool.min)
#     ncoal=ncoal+Neg_dist_info[[1]]
#     depart.list=Neg_dist_info[[2]]
#     
#     ######################################################################
#     #Find min_dist after coalescence (a new breath figure)
#     #Once coalescence happen, min_dist should be calculated over all existing droplets
#     min_dist0=matrix(c(-1,-1,0,0,max(h,w)),nrow=1)
#     min_dist=min_dist0
#     Min_dist_info<-Min_dist(t+dtp,min_dist,as.integer(names(bfobj)))
#     min_dist=Min_dist_info[[1]]
#     #search.pool.min is a search pool for finding minimum gap in droplet pattern
#     search.pool.min=Min_dist_info[[2]]
    
    print(c("search.pool.min:",search.pool.min))
    num_coal=ncoal-ncoal.old
    return(list(num_coal,dtp,min_dist,depart.list,search.pool.min))
  }
}