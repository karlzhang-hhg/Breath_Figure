library(plotrix)

######################################################################
#Plot the frame of current time t
Plot_bf_online<-function(j,t,h,w,dt)
{
  drop.info=NULL
  for (i in 1:length(bfobj))
  {
    drop.info.i=c(bfobj[[i]]$ind,bfobj[[i]]$posi,bfobj[[i]]$r)
    #Handle cyclic effects from boundaries
    #if droplet exceed left and right boundaries
    if (bfobj[[i]]$posi[1]-bfobj[[i]]$r<0)
    {
      drop.info.i=rbind(drop.info.i,c(bfobj[[i]]$ind,bfobj[[i]]$posi+c(w,0),bfobj[[i]]$r))
    }
    else
    {
      if (bfobj[[i]]$posi[1]+bfobj[[i]]$r>w)
      {
        drop.info.i=rbind(drop.info.i,c(bfobj[[i]]$ind,bfobj[[i]]$posi+c(-w,0),bfobj[[i]]$r))
      }
    }
    #if droplet exceed upper and lower boundaries
    if (bfobj[[i]]$posi[2]-bfobj[[i]]$r<0)
    {
      drop.info.i=rbind(drop.info.i,c(bfobj[[i]]$ind,bfobj[[i]]$posi+c(0,h),bfobj[[i]]$r))
    }
    else
    {
      if (bfobj[[i]]$posi[2]+bfobj[[i]]$r>h)
      {
        drop.info.i=rbind(drop.info.i,c(bfobj[[i]]$ind,bfobj[[i]]$posi+c(0,-h),bfobj[[i]]$r))
      }
    }
    
    drop.info=rbind(drop.info,drop.info.i)
  }
  #print(c("plotdrops",nrow(drop.info)))
#   old.wd<-getwd()
#   setwd(paste0(getwd(),"/images"))
  pdf(paste0(j,".pdf"))
  #pdf(paste0(j,"-",(t-((t/dt*1000)%%1L)*dt/1000),".pdf"))
  plot.new()
  par(mar=c(0,0,0,0))
  print(par("plt"))
  plot.window(xlim=c(0,w),ylim=c(0,h),xaxs="i",yaxs="i",asp=1)
  #rect(0,0,w,h)
  #box()
  #title(paste0(j,"-",(t-((t/dt*1000)%%1L)*dt/1000)))
  draw.circle(drop.info[,2],drop.info[,3],drop.info[,4],
              border = "blue",
              col = color.scale(drop.info[,4],c(1,0),0,c(0,1),alpha = 0.5),
              lty="blank")
  
  dev.off()
#   setwd(old.wd)
}