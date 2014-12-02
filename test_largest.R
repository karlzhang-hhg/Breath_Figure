test_largest<-function(subset)
{
  rmax=0
  imax=0
  ind.range=names(bfobj)[!(names(bfobj) %in% as.character(subset))]
  for (i in 1:length(bfobj[ind.range]))
  {
    if (bfobj[ind.range][[i]]$r>rmax)
    {
      rmax=bfobj[ind.range][[i]]$r
      imax=bfobj[ind.range][[i]]$ind
    }
  }
  return(list(imax,rmax))
}

cal_center_dist<-function(i1,i2)
{
  return(Dist(bfobj[[i1]]$posi,bfobj[[i2]]$posi,h,w))
}
