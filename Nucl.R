library(spatstat)
source("Crea_drop.R")
source("Dist.R")


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
  #plot(nucl)
  if (length(bfobj)==0)
  {
    ndrop=ndrop+1L
    ndropall=ndropall+1L
    search.pool=c(search.pool,ndropall)
    #The name of Crea_drop wouldn't be passed, so that first asign element; then change name
    bfobj[[ndrop]]<<-Crea_drop(ndropall,t,nucl$x[1],nucl$y[1],nuclr)
    names(bfobj)[ndrop]<<-as.character(ndropall)
  }
  for (i in 2:nucl$n)
  {
    ndrop<-length(bfobj)
    flag.thin=0 #thin or not
    for (j in 1:ndrop)
    {
      if (Dist(c(nucl$x[i],nucl$y[i]),bfobj[[j]]$posi,h,w)[3]<=nuclr+dr+bfobj[[j]]$r)
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
      bfobj[[ndrop]]<<-Crea_drop(ndropall,t,nucl$x[i],nucl$y[i],nuclr)
      names(bfobj)[ndrop]<<-as.character(ndropall)
    }
  }
  return(search.pool)
}