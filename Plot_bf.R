library(plotrix)

######################################################################
#Plot the frame of current time t
Plot_bf<-function(i,t,h,w)
{
  setwd("C:/My Life Style/Courses of Statistics/STA 561 Probability Machine Learning/Project/R-code/image")
  pdf(paste0(i,"-",t,".pdf"))
  plot.new()
  plot.window(xlim=c(0,w),ylim=c(0,h),xaxs="i",yaxs="i")
  box()
  title(mian=paste0(i,"-",t))
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
  draw.circle(drop.info[,2],drop.info[,3],drop.info[,4],
              border = "blue",
              col = color.scale(drop.info[,2],0,0,1,alpha = 0.5))
  dev.off()
}