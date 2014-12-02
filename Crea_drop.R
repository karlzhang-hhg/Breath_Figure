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
             log.tr=matrix(c(t,x,y,inir),nrow=1),#log of time, position, and radius
             abso.drop=as.vector(as.integer()),#indices of droplets absorbed to generate this droplet
             coal=list()
  )
  #print(node)
  colnames(node$log.tr)<-c('t','x','y','radius')
  return(node)           
}
