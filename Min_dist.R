source("Dist.R")

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
      if (bfobj[[i]]$ind==bfobj[[as.character(j)]]$ind) {next}
      cal_dist=Dist(bfobj[[i]]$posi,bfobj[[as.character(j)]]$posi,h,w)
      gap=cal_dist[3]-bfobj[[i]]$r-bfobj[[as.character(j)]]$r
      stopifnot(gap>=0)
      if (gap<min_dist[1,5])
      {
        min_dist=matrix(c(bfobj[[i]]$ind,j,cal.dist[1:2],gap),nrow=1)
      }
      else
      {
        if (gap==min_dist[1,5])
        {
          min_dist=rbind(min_dist,c(bfobj[[i]]$ind,j,cal.dist[1:2],gap))
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