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

#Find number of coalescence for all droplets in bfallobj
#Check if the length of bfallobj[[i]]$abso.drop and 
#   bfallobj[[i]]$coal match with each other for
#   absorbed droplets and existing droplets
drop_num_coal<-function()
{
  res=NULL
  for (i in 1:length(bfallobj))
  {
    n1=length(bfallobj[[i]]$abso.drop)
    n2=length(bfallobj[[i]]$coal)
    #print(c(n1,n2))
    stopifnot((n1==n2 & bfallobj[[i]]$deat.t<0)|(n1==n2-1 & bfallobj[[i]]$deat.t>0))
    res=c(res,n1)
  }
  return(res)
}

#Find number of coalescence for all droplets in bfobj
#Check if the length of bfallobj[[i]]$abso.drop and 
#   bfallobj[[i]]$coal match with each other for
#   absorbed droplets and existing droplets
drop_num_coal_exis<-function()
{
  res=NULL
  for (i in 1:length(bfobj))
  {
    n1=length(bfobj[[i]]$abso.drop)
    n2=length(bfobj[[i]]$coal)
    #print(c(n1,n2))
    if(!((n1==n2 & bfobj[[i]]$deat.t<0)|(n1==n2-1 & bfobj[[i]]$deat.t>0)))
    {
      print(c(bfobj[[i]]$ind,n1,n2))
      stop("error")
    }
    res=c(res,n1)
  }
  return(res)
}

for (i in 1:length(bfallobj))
{
  if (length(bfallobj[[i]]$coal)!=0)
  {
    for (j in 1:length(bfallobj[[i]]$coal))
    {
      if (is.null(bfallobj[[i]]$coal[[j]]))
      {
        print(bfallobj[[i]]$ind)
        break
      }
    }
  }
}