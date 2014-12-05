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
  print(i)
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
plot(tseq,accu_coal_times_t,pch=1,cex=0.1,main="Times of coalecence before t",xlab="t",ylab="Accumulated times of coalescence")
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

#Coalsecence times v.s. size
coal_drop_t=matrix(rep(NA,length(bfallobj)*length(tseq)),ncol=length(tseq))
for(i in 1:length(bfallobj))
{
  k=1
  for(j in 1:length(tseq))
  {
    if(bfallobj[[i]]$coal[[k]]$t==tseq[j])
    {
      coal_drop_t[i,j]=bfallobj[[i]]$log.tr[]
      k=k+1
      if (k>nrow(bfall))
    }
  }
}


#At every moment, the size distribution of existing droplets with different times of coalescence



