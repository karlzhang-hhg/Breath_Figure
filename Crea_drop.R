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
  #print(node)
  return(node)           
}
