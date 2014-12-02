######################################################################
#Update information after one coalescence between two droplets
#Delete the smaller droplet
#Update information for the big droplet in bfobj
#Updata information in bfallobj
Update_coal<-function(ncoal,t,small,big,edge.circ,ds,h,w)
{
#   stopifnot(is.integer(ncoal))
#   stopifnot(is.integer(small))
#   stopifnot(is.integer(big))
  #convert the index number to character
  print(c("ncoal:",ncoal))
  small=as.character(small)
  big=as.character(big)
  #update information for the small droplet after a coalescence
  bfobj[[small]]$deat.t<<-t
  bfobj[[small]]$log.tr<<-rbind(bfobj[[small]]$log.tr,c(t,bfobj[[small]]$posi,bfobj[[small]]$r))
  ######################################################################
  #Store all infromation about coalescence
  coalall[[ncoal]]<<-list(ind=ncoal,t=t,invo.drop=matrix(as.integer(c(small,big)),nrow=1),delt.s=ds)#invo.drop initially is NULL
  names(coalall)[ncoal]<<-as.character(ncoal)
  ######################################################################
  
  
  
  #update information for the big droplet after a coalescence
  bfobj[[big]]$r<<-((bfobj[[small]]$r)^3+(bfobj[[big]]$r)^3)^(1/3)
  #Consider the circlic effect
  if (edge.circ[1]==0 & edge.circ[2]==0)
  {
    #mass weighted average center position as the new position
    bfobj[[big]]$posi<<-bfobj[[big]]$posi+(bfobj[[small]]$posi-bfobj[[big]]$posi)*bfobj[[small]]$r^3/(bfobj[[small]]$r^3+bfobj[[big]]$r^3)
  }
  else
  {
    if (edge.circ[1]==1)
    {
      bfobj[[big]]$posi<<-(bfobj[[big]]$posi+(bfobj[[small]]$posi-bfobj[[big]]$posi-c(w,0)*sign(bfobj[[small]]$posi[1]-bfobj[[big]]$posi[1]))*bfobj[[small]]$r^3/(bfobj[[small]]$r^3+bfobj[[big]]$r^3))%%w
    }
    else #edge.circ[2]==1
    {
      bfobj[[big]]$posi<<-(bfobj[[big]]$posi+(bfobj[[small]]$posi-bfobj[[big]]$posi-c(0,h)*sign(bfobj[[small]]$posi[2]-bfobj[[big]]$posi[2]))*bfobj[[small]]$r^3/(bfobj[[small]]$r^3+bfobj[[big]]$r^3))%%h
    }
  }
  #The information of coalescence from the small absorbed droplet should be included into the big one
  bfobj[[big]]$log.tr<<-rbind(bfobj[[small]]$log.tr,bfobj[[big]]$log.tr,c(t,bfobj[[big]]$posi[1],bfobj[[big]]$posi[2],bfobj[[big]]$r))
  bfobj[[big]]$abso.drop<<-c(bfobj[[big]]$abso.drop,bfobj[[small]]$abso.drop,as.integer(small))
  
  #Update coalescence information
  if (length(bfobj[[small]]$coal)!=0) #even list() is not NULL. it just has 0 length
  {
    bfobj[[big]]$coal<<-as.list(c(bfobj[[big]]$coal,bfobj[[small]]$coal))
  }
#   if (is.null(bfobj[[big]]$coal[[1]]))
#   {
#     bfobj[[big]]$coal[[length(bfobj[[big]]$coal)+1L]]<<-coalall[[ncoal]]
#   }
#   else
#   {
  bfobj[[big]]$coal[[length(bfobj[[big]]$coal)+1L]]<<-coalall[[ncoal]]
#   }
  bfobj[[small]]$coal[[length(bfobj[[small]]$coal)+1L]]<<-coalall[[ncoal]]
  if (is.null(bfobj[[big]]$coal[[1]]))
  {
    print(c(bfobj[[big]],bfobj[[small]]))
    stopifnot(!is.null(bfobj[[big]]$coal[[1]]))
  }
  #   bfallobj[names(bfobj)]<<-bfobj
  #   bfobj<<-bfobj[!(names(bfobj) %in% as.character(delete.list))]
}