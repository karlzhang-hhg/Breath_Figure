load(tseq)
load(bfallobj)

#Evolutionary droplet size distribution
size_dis_t=list()
for (i in 1:length(bfallobj))
{
  k=1
  for (j in 1:length(tseq))
  {
    if ((bfallobj[[i]]$deat.t<0 & bfallobj[[i]]$jump.t<0)|
          (bfallobj[[i]]$deat.t>tseq[j]) | ( bfallobj[[i]]$jump.t>tseq[j]))
    {
      while (k<nrow(bfallobj[[i]]$log.tr))
      {
        if (bfallobj[[i]]$log.tr[k+1,1]<=tseq[j])
        {
          k=k+1
        }
        else
        {
          break
        }
      }
      if (bfallobj[[i]]$log.tr[k,1]==tseq[j])
      {
        if (j>length(size_dis_t))
        {
          size_dis_t[[j]]<-matrix(c(bfallobj[[i]]$ind,bfallobj[[i]]$log.tr[k,4]),ncol=2)
        }
        else
        {
          size_dis_t[[j]]<-rbind(size_dis_t[[j]],c(bfallobj[[i]]$ind,bfallobj[[i]]$log.tr[k,4]))
        }
        k=k+1
      }
    }
    if (k>nrow(bfallobj[[i]]$log.tr)) break
  }
  #print(i)
}

#Find maximum r for range of density plot
rmax=0
for(i in 1:length(bfallobj))
{
  if(rmax<bfallobj[[i]]$r)
  {
    rmax=bfallobj[[i]]$r
  }
}

#plot(density(size_dis_t[[1500]][,2]))
plot(density(size_dis_t[[1]][,2],from=0,to=2*rmax),xlim=c(0,rmax))
lines(density(size_dis_t[[1500]][,2]))

#0
pdf("size_dis-0.pdf")
plot(density(size_dis_t[[1]][,2],from=0),main="Drop size distribution at t=0",xlab="Radius",ylab="Density")
dev.off()

#10
pdf("size_dis-0.064.pdf")
plot(density(size_dis_t[[11]][,2],from=0),main="Drop size distribution at t=0.064",xlab="Radius",ylab="Density")
dev.off()

#100
pdf("size_dis-0.135.pdf")
plot(density(size_dis_t[[101]][,2],from=0),main="Drop size distribution at t=0.135",xlab="Radius",ylab="Density")
dev.off()

#200
pdf("size_dis-0.173.pdf")
plot(density(size_dis_t[[201]][,2],from=0),main="Drop size distribution at t=0.173",xlab="Radius",ylab="Density")
dev.off()

#400
pdf("size_dis-0.261.pdf")
plot(density(size_dis_t[[401]][,2],from=0),main="Drop size distribution at t=0.261",xlab="Radius",ylab="Density")
dev.off()

#600
pdf("size_dis-0.390.pdf")
plot(density(size_dis_t[[101]][,2],from=0),main="Drop size distribution at t=0.390",xlab="Radius",ylab="Density")
dev.off()

#800
pdf("size_dis-0.504.pdf")
plot(density(size_dis_t[[801]][,2],from=0),main="Drop size distribution at t=0.504",xlab="Radius",ylab="Density")
dev.off()

plotc=c(10,100,500,300,seq(200,2000,200))
# plotc=c(500)
for (i in 1:length(plotc))
{
  pdf(paste0("size_dis-",tseq[plotc[i]+1]-((tseq[plotc[i]+1]*1000)%%1)/1000,"-",plotc[i],".pdf"))
  plot(density(size_dis_t[[plotc[i]+1]][,2],from=0),ylim=c(0,50),main=paste0("Drop size distribution at t=",tseq[plotc[i]+1]-((tseq[plotc[i]+1]*1000)%%1)/1000),xlab="Radius",ylab="Density")
  dev.off()
}

plot_drop_distr_t<-function(plotc,size_dis_t,tseq,y_lim,rmax)
{
  lwd<-seq(0.5,2.5,length=length(plotc))
  plot(density(size_dis_t[[plotc[1]+1]][,2],from=0),xlim=c(0,rmax),ylim=y_lim,
       main=paste0("Drop size distribution at Different time"),
       xlab="Radius",ylab="Density",lty=1,lwd=lwd[1])
  legend_name<-c(paste0("t=",tseq[plotc[1]+1]-((tseq[plotc[1]+1]*1000)%%1)/1000))
  if (length(plotc)==1) return()
  for (i in 2:length(plotc))
  {
    legend_name<-c(legend_name,paste0("t=",tseq[plotc[i]+1]-((tseq[plotc[i]+1]*1000)%%1)/1000))
    lines(density(size_dis_t[[plotc[i]+1]][,2],from=0),col=i,lty=i,lwd=lwd[i])
  }
  legend(x = "topright",legend = legend_name,
         bty="n",
         text.col=1:length(plotc),
         lty=1:length(plotc),pch=NA,
         lwd=lwd,
         col=1:length(plotc))
  par("yaxs"="r")
}

plotc<-c(10,50,100,300,400,seq(500,4000,500))

pdf("Asymptotic Distribution.pdf")
plot_drop_distr_t(plotc,size_dis_t,tseq,c(0,80),rmax)
dev.off()

#Jumping size, time, and number of coalescence
jump_info=NULL
for (i in 1:length(depart.list))
{
  const=length(bfallobj[[depart.list[[i]]]]$coal)
  nofcoal=0
  while (bfallobj[[depart.list[[i]]]]$coal[[const-nofcoal]]$t==bfallobj[[depart.list[[i]]]]$jump.t)
  {
    nofcoal=nofcoal+1
    if(nofcoal>=const) break
  }
  jump_info=rbind(jump_info,c(bfallobj[[depart.list[[i]]]]$jump.t,bfallobj[[depart.list[[i]]]]$r,as.character(nofcoal)))
}

hist(as.numeric(jump_info[,2]))
plot(density(as.numeric(jump_info[,2]),from=0))
plot(density(as.numeric(jump_info[,2]),from=0))

plot(as.numeric(jump_info[,1]),as.numeric(jump_info[,2]))

pdf("jump_times-t.pdf")
plot(as.numeric(jump_info[,1]),(1:nrow(jump_info)),pch=1,cex=0.1,main="Times of disappearance before t",xlab="t",ylab="Accumulated times of disappearance")
dev.off()

plot(as.integer(jump_info[,3]))

#Coalescence times v.s. time
coal_times_t=rep(0,length(tseq))
accu_coal_times_t=rep(0,length(tseq))
k=1
for(i in 1:length(coalall))
{
  while (coalall[[i]]$t>tseq[k])
  {
    k=k+1
  }
  if (coalall[[i]]$t==tseq[k])
  {
    coal_times_t[k]=coal_times_t[k]+1
  }
  accu_coal_times_t[k]=sum(coal_times_t[1:k])
}

pdf("coal_times-t.pdf")
plot(tseq,accu_coal_times_t,pch=1,cex=0.1,main="Times of coalescence before t",xlab="t",ylab="Accumulated times of coalescence")
dev.off()

#Sort information in bfallobj
bfallobj_sort<-bfallobj
for(i in 1:length(bfallobj_sort))
{
  delete.row=NULL
  if(nrow(bfallobj_sort[[i]]$log.tr)>1)
  {
    for(j in 1:(nrow(bfallobj_sort[[i]]$log.tr)-1))
    {
      if(bfallobj_sort[[i]]$log.tr[j,1]==bfallobj_sort[[i]]$log.tr[j+1,1])
      {
        delete.row=c(delete.row,j)
      }
    }
    if (!is.null(delete.row))
    {
      bfallobj_sort[[i]]$log.tr<-bfallobj_sort[[i]]$log.tr[-delete.row,]
    }
  }
}

# #Coalsecence times v.s. size
# coal_drop_t=matrix(rep(NA,length(bfallobj_sort)*length(tseq)),ncol=length(tseq))
# for(i in 1:length(bfallobj_all))
# {
#   k=1
#   for(j in 1:length(tseq))
#   {
#     if(bfallobj[[i]]$coal[[k]]$t==tseq[j])
#     {
#       coal_drop_t[i,j]=bfallobj[[i]]$log.tr[]
#       k=k+1
#       if (k>nrow(bfall))
#     }
#   }
# }

#Coalsecence times v.s. size
coal_drop<-function(bfallobj_t,tseq)
{
  coal_drop_t=list()
  for(i in 1:length(bfallobj_t))
  {
    k=1
    for(j in 1:length(tseq))
    {
      if (length(bfallobj_t[[i]]$coal)==0) break
      if(bfallobj_t[[i]]$coal[[k]]$t==tseq[j])
      {
        r<-max(bfallobj_t[[i]]$log.tr[bfallobj_t[[i]]$log.tr[,1]==bfallobj_t[[i]]$coal[[k]]$t,4])
        if (k>length(coal_drop_t))
        {
          coal_drop_t<-c(coal_drop_t, list(matrix(c(bfallobj_t[[i]]$coal[[k]]$t,r),nrow=1)))
        }
        else
        {
          coal_drop_t[[k]]<-rbind(coal_drop_t[[k]],c(bfallobj_t[[i]]$coal[[k]]$t,r))
        }
        k=k+1
        if (k>length(bfallobj_t[[i]]$coal)) break
      }
    }
  }
  return(coal_drop_t)
}

coal_drop_t<-coal_drop(bfallobj_sort,tseq)


#Find the maximum time of coalescence and return the index of the droplets in bfallobj_sort or bfallobj
maxncoal=0
maxncoal_ind=0

for (i in 1:length(bfallobj_sort))
{
  if (length(bfallobj_sort[[i]]$coal)>maxncoal)
  {
    maxncoal=length(bfallobj_sort[[i]]$coal)
    maxncoal_ind=bfallobj_sort[[i]]$ind
  }
}

maxncoal_drop<-bfallobj_sort[[maxncoal_ind]]
coal_t<-unlist(lapply(maxncoal_drop$coal,function(x) x$t))


plot(1:length(maxncoal_drop$coal),maxncoal_drop$log.tr[maxncoal_drop$log.tr[,1] %in% coal_t,4],
     pch=1,cex=0.3,
     main="Times of coalescence v.s. size",xlab="# of coalescence",
     ylab="Size of droplets with /n the most coalescence times")


CEX=0.7

pdf("coal_times-size.pdf")
plot(1:length(coal_drop_t),unlist(lapply(coal_drop_t,function(x) mean(x[,2]))),
     pch=1,cex=CEX,main="Times of coalescence v.s. size",xlab="# of coalescence",ylab="Size of droplets")
points(1:length(coal_drop_t),unlist(lapply(coal_drop_t,function(x) mean(x[,2]))) +sqrt(unlist(lapply(coal_drop_t,function(x) var(x[,2])))),
       pch=1,cex=CEX,col="red")
points(1:length(coal_drop_t),unlist(lapply(coal_drop_t,function(x) mean(x[,2]))) -sqrt(unlist(lapply(coal_drop_t,function(x) var(x[,2])))),
       pch=1,cex=CEX,col="red")
points(1:length(maxncoal_drop$coal),maxncoal_drop$log.tr[maxncoal_drop$log.tr[,1] %in% coal_t,4],
       pch=3,cex=CEX,col="blue")
legend(x = "bottomright",legend = c("Mean radius","1 STD upper","1 STD lower", "The case with largest coalescence time"),
       bty="n",
       text.col=c("black","red","red","blue"),
       lty=NA,pch=c(1,1,1,3),
       lwd=1,
       col=c("black","red","red","blue"))
par("yaxs"="r")
dev.off()

#At every moment, the size distribution of existing droplets with different times of coalescence

drop_gener_t<-function(bfallobj_t,tseq,t,y_lim,nsample)
{
  for (i in 1:length(tseq))
  {
    if (t<tseq[i])
    {
      nt<-i-1
      break
    }
  }
  drop_list<-NULL
  for (i in 1:length(bfallobj_t))
  {
    if (bfallobj_t[[i]]$deat.t>=t | bfallobj_t[[i]]$jump.t>=t)
    {
      drop_list<-c(drop_list,i)
    }
  }
  
  plot_drop_gener_t<-function(drop_gener_t,t)
  {
    plot(density(drop_gener_t1[[1]][,2],from=0,to=2*rmax),xlim=c(0,rmax),ylim=y_lim,
         main=paste("Distribution of different \n generation (coalescence) at t=",t),
         xlab="Drop size",
         ylab="Density")
    legend_name<-c(paste0("g1","-",length(drop_gener_t[[1]])))
    plotted_one<-c(1)
    for (i in 2:length(drop_gener_t))
    {
      if (length(drop_gener_t[[i]])>nsample)
      {
        lines(density(drop_gener_t1[[i]][,2]),col=i)
        legend_name<-c(legend_name,paste0("g",i,"-",length(drop_gener_t[[i]])))
        plotted_one<-c(plotted_one,i)
      }
    }
    legend(x = "topright",legend = legend_name,
           bty="n",
           text.col=1:length(drop_gener_t)[plotted_one],
           lty=1,pch=NA,
           lwd=1,
           col=1:length(drop_gener_t)[plotted_one])
    par("yaxs"="r")
  }
  
  drop_gener_t1<-coal_drop(bfallobj_t[drop_list],tseq[1:nt])
  plot_drop_gener_t(drop_gener_t1,t)
  
  return(drop_gener_t1)
}

for (i in 1:23)
{
  pdf(paste0("Distr-gener-",i/10,".pdf"))
  drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,i/10,c(0,100),10)
  dev.off()
}

drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,0.5,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,0.8,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,0.9,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.0,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.1,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.2,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.3,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.4,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.5,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.6,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.7,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.8,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,1.9,c(0,100),10)
drop_gener_t1<-drop_gener_t(bfallobj_sort,tseq,2.0,c(0,100),10)








