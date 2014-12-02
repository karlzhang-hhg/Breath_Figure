source("Dist.R")
source("Update_coal.R")

######################################################################
#Find neg_dist at time t. neg_dist indicate coalescence
#Update bfobj if there is some coalescence
#Stop when no coalescence happens at time t
Neg_dist<-function(t,ncoal,h,w,depart.list,search.pool)
{
  ncoal.old=ncoal
  repk=0
  repeat{
    repk=repk+1
    print(c("repk:",repk))
    ndrop=length(bfobj)
    #only those coalesced droplets need to be checked in the next loop (search.pool.new)
    search.pool.new=NULL
    neg_dist=NULL
    #check if there is negative distance
    for (i in 1:ndrop)
    {
      for (j in search.pool[as.character(search.pool) %in% names(bfobj)[i+1L:ndrop]])
      {
        #if the two droplets for distance calculation are the same, move to the next iteration
        stopifnot(bfobj[[i]]$ind!=j)
        cal_dist=Dist(bfobj[[i]]$posi,bfobj[[as.character(j)]]$posi,h,w)
        gap=cal_dist[3]-bfobj[[i]]$r-bfobj[[as.character(j)]]$r
        if(gap<=0)
        {
          neg_dist=matrix(rbind(neg_dist,c(bfobj[[i]]$ind,j,cal_dist)),ncol=5)
        }
      }
    }
    
    if (is.null(neg_dist)) #a valid breath figure; nothing needs to change
    {
      break
    }
    else #some coalescence need to happen
    {
      #The search.pool.min should be large enough, in case missing
      search.pool.new=unique(as.vector(neg_dist[,1:2]))
      ######################################################################
      #Sort the neg_dist so that droplets can only appear once in search.pool.neg or neg_dist
      search.pool.neg=as.vector(neg_dist[1,1:2])
      delete.pair=NULL
      if (nrow(neg_dist)>1)
      {
        for (i in 2:nrow(neg_dist))
        {
          if ((neg_dist[i,1] %in% search.pool.neg)|(neg_dist[i,2] %in% search.pool.neg))
          {
            delete.pair=c(delete.pair,i)
          }
          else
          {
            search.pool.neg=c(search.pool.neg,as.vector(neg_dist[i,1:2]))
          }
        }
      }
      if (!is.null(delete.pair))
      {
        neg_dist=matrix(as.vector(neg_dist[-delete.pair,]),ncol=5)
      }
      ######################################################################
      delete.list=NULL
      print(c("neg_dist:",neg_dist,is.matrix(neg_dist)))
      for (i in 1:nrow(neg_dist))
      {
        ncoal=ncoal+1L
        print(c("i",i,neg_dist[i,1],neg_dist[i,2]))
        if (bfobj[[as.character(neg_dist[i,1])]]$r<bfobj[[as.character(neg_dist[i,2])]]$r)
        {
          small=neg_dist[i,1]
          big=neg_dist[i,2]         
        }
        else
        {
          small=neg_dist[i,2]
          big=neg_dist[i,1]
        }
        #print(c("i",i,is.matrix(neg_dist)))
        #update information for both droplets
        Update_coal(ncoal,t,small,big,neg_dist[i,3:4],ds = 0,h,w)
        #Check if the radius after coalescence is larger than rcr
        #If does, the drop jump away
        #Otherwise, it stays
        if (bfobj[[as.character(big)]]$r>rcr)
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
              #stopifnot(big %in% depart.list) #condensation wouldn't result in departure
              depart.list=depart.list[depart.list[]!=small]
            }
            depart.list=unique(c(depart.list,big))
          }
        }
        #search.pool.new=unique(c(search.pool.new,big))
        delete.list=c(delete.list,small) #included those absorbed
      }
      search.pool=subset(search.pool,!(search.pool %in% delete.list))
      stopifnot(!is.null(delete.list))
      bfallobj[names(bfobj)]<<-bfobj
      #Delete droplets absorbed or departing from bfobj
      bfobj<<-bfobj[!(names(bfobj) %in% as.character(delete.list))]
    }
  }
  bfallobj[names(bfobj)]<<-bfobj
  bfobj<<-bfobj[!(names(bfobj) %in% as.character(depart.list))]
  num_coal=ncoal-ncoal.old
  return(list(num_coal,depart.list))
}