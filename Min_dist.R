source("Dist.R")

######################################################################
#Find min_dist after nucleation
Min_dist<-function(t,min_dist,search.pool.min)
{
  ndrop=length(bfobj)
  neg_dist_coll=NULL
  #find the minimum non-negative distance among the droplet pattern
  for (i in 1:ndrop)
  {
    for (j in search.pool.min[as.character(search.pool.min) %in% names(bfobj)[i+1L:ndrop]])
    {
      stopifnot(bfobj[[i]]$ind!=j)
      cal_dist=Dist(bfobj[[i]]$posi,bfobj[[as.character(j)]]$posi,h,w)
      gap=cal_dist[3]-bfobj[[i]]$r-bfobj[[as.character(j)]]$r
      if (gap<=0)
      {
        print(c("!!!!gap<=0:",bfobj[[i]]$ind,j,gap))
        neg_dist_coll=c(neg_dist_coll,bfobj[[i]]$ind,j)
      }
      if (gap>0 & gap<min_dist[1,5])
      {
        min_dist=matrix(c(bfobj[[i]]$ind,j,cal_dist[1:2],gap),nrow=1)
      }
      else
      {
        if (gap==min_dist[1,5])
        {
          min_dist=rbind(min_dist,c(bfobj[[i]]$ind,j,cal_dist[1:2],gap))
        }
      }
    }
  }
  #The search.pool.min should be large enough, in case missing
  search.pool.min.old=unique(c(as.vector(min_dist[,1:2]),neg_dist_coll))
  ######################################################################
  #Sort the min_dist so that droplets can only appear once in search.pool.min or min_dist
  search.pool.min=as.vector(min_dist[1,1:2])
  delete.pair=NULL
  if (nrow(min_dist)>1)
  {
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
  }
  if (!is.null(delete.pair))
  {
    min_dist=as.matrix(as.vector(min_dist[-delete.pair,]),ncol=5)
  }
  #search.pool.min=as.vector(min_dist[,1:2])
#   print(c("search.pool.min.old:",search.pool.min.old))
  return(list(min_dist,search.pool.min.old))
}