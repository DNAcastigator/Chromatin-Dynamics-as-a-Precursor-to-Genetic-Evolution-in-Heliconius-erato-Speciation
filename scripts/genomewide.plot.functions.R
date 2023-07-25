
library(cowplot)
library(ggforce)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggpubr)
library(pheatmap)
library(RUVSeq)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(Rcpp)
library(reshape2)
library(ggrastr)
library(GenomicRanges)
library(rtracklayer)

###### data
geneslist=read.table("../total/genes.start.pan.txt",h=F,sep="\t")
Vvl=data.frame(V1=410278646,V2=410279821)
chrom.numb=data.frame(a=c(seq(1,20),"Z"),b=(chrom.SE$end+chrom.SE$start)/2)
###### nota
sweep.not=read.table("sweeps/sweep.pan.nota.final.1l.txt")
sweep.not=sweep.not[sweep.not$V2< 688640245,]

####dem
sweep.dem=read.table("sweeps/sweep.pan.dem.final.1l.txt")
sweep.dem=sweep.dem[sweep.dem$V2< 688640245,]

#### etylus
sweep.ety=read.table("sweeps/sweep.pan.ety.final.1l.txt")
sweep.ety=sweep.ety[sweep.ety$V2< 688640245,]


### hyd
sweep.hyd=read.table("sweeps/sweep.pan.hyd.final.1l.txt")
sweep.hyd=sweep.hyd[sweep.hyd$V2< 688640245,]

###favorinus
sweep.favo=read.table("sweeps/favo.pan.total.1l.txt")
sweep.favo=sweep.favo[sweep.favo$V2< 688640245,]

###chestertonii
sweep.chest=read.table("sweeps/chest.pan.total.1l.txt")
sweep.chest=sweep.chest[sweep.chest$V2< 688640245,]



#######notaxety 1k fst

notaxety.1k.fst=read.table("fst/notavsety.1k.FST.2022.1l.txt",h=F,sep=" ")
colnames(notaxety.1k.fst)=c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst")
notaxety.1k.fst$position=notaxety.1k.fst$start
notaxety.1k.fst$Fst[notaxety.1k.fst$Fst<0]<-0
notaxety.1k.fst=notaxety.1k.fst[notaxety.1k.fst$start< 688640245,]

notaxety.10k.fst=read.table("fst/notaxety.10k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
#colnames(notaxety.10k.fst)=c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst")
notaxety.10k.fst$position=notaxety.10k.fst$start
notaxety.10k.fst$Fst[notaxety.10k.fst$Fst<0]<-0
notaxety.10k.fst=notaxety.10k.fst[notaxety.10k.fst$start< 688640245,]



notaxety.5k.fst=read.table("fst/notaxety.5k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
#colnames(notaxety.5k.fst)=c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst")
notaxety.5k.fst$position=notaxety.5k.fst$start
notaxety.5k.fst$Fst[notaxety.5k.fst$Fst<0]<-0
notaxety.5k.fst=notaxety.5k.fst[notaxety.5k.fst$start< 688640245,]

########### FSt 50k
notaxety.fst=read.table("fst/notavsety.50k.FST.2022.1l.txt",h=F,sep=" ")
colnames(notaxety.fst)=c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst")
notaxety.fst$Fst[notaxety.fst$Fst<0]<-0
notaxety.fst$position=notaxety.fst$start
notaxety.fst=notaxety.fst[notaxety.fst$start< 688640245,]



####dem hyd
demxhd.fst=read.table("fst/demoxhyda.50k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxhd.fst$Fst[demxhd.fst$Fst<0]<-0
demxhd.fst$position=demxhd.fst$start
demxhd.fst=demxhd.fst[demxhd.fst$start< 688640245,]

demxhd.1k.fst=read.table("fst/demoxhyda.1k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxhd.1k.fst$Fst[demxhd.1k.fst$Fst<0]<-0
demxhd.1k.fst$position=demxhd.1k.fst$start
demxhd.1k.fst=demxhd.1k.fst[demxhd.1k.fst$start< 688640245,]

demxhd.10k.fst=read.table("fst/demoxhyda.10k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxhd.10k.fst$Fst[demxhd.10k.fst$Fst<0]<-0
demxhd.10k.fst$position=demxhd.10k.fst$start
demxhd.10k.fst=demxhd.10k.fst[demxhd.10k.fst$start< 688640245,]

demxhd.5k.fst=read.table("fst/demoxhyda.5k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxhd.5k.fst$Fst[demxhd.5k.fst$Fst<0]<-0
demxhd.5k.fst$position=demxhd.5k.fst$start
demxhd.5k.fst=demxhd.5k.fst[demxhd.5k.fst$start< 688640245,]

############### dem chest
demxches.fst=read.table("fst/demvschest.50k.FST.1l.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxches.fst$Fst[demxches.fst$Fst<0]<-0
demxches.fst$position=demxches.fst$start
demxches.fst=demxches.fst[demxches.fst$start< 688640245,]

demxches.1k.fst=read.table("fst/demvschest.1k.FST.1l.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxches.1k.fst$Fst[demxches.1k.fst$Fst<0]<-0
demxches.1k.fst$position=demxches.1k.fst$start
demxches.1k.fst=demxches.1k.fst[demxches.1k.fst$start< 688640245,]

demxches.5k.fst=read.table("fst/demvschest.5k.FST.1l.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxches.5k.fst$Fst[demxches.5k.fst$Fst<0]<-0
demxches.5k.fst$position=demxches.5k.fst$start
demxches.5k.fst=demxches.5k.fst[demxches.5k.fst$start< 688640245,]

#################demo fav

demxfav.fst=read.table("fst/demvsfav.50k.FST.1l.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxfav.fst$Fst[demxfav.fst$Fst<0]<-0
demxfav.fst$position=demxfav.fst$start
demxfav.fst=demxfav.fst[demxfav.fst$start< 688640245,]

demxfav.1k.fst=read.table("fst/demvsfav.1k.FST.1l.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxfav.1k.fst$Fst[demxfav.1k.fst$Fst<0]<-0
demxfav.1k.fst$position=demxfav.1k.fst$start
demxfav.1k.fst=demxfav.1k.fst[demxfav.1k.fst$start< 688640245,]

demxfav.5k.fst=read.table("fst/demvsfav.5k.FST.1l.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxfav.5k.fst$Fst[demxfav.5k.fst$Fst<0]<-0
demxfav.5k.fst$position=demxfav.5k.fst$start
demxfav.5k.fst=demxfav.5k.fst[demxfav.5k.fst$start< 688640245,]


############# demo ety
demxety.fst=read.table("fst/demoxety.50k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxety.fst$Fst[demxety.fst$Fst<0]<-0
demxety.fst$position=demxety.fst$start
demxety.fst=demxety.fst[demxety.fst$start< 688640245,]

demxety.1k.fst=read.table("fst/demoxety.1k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxety.1k.fst$Fst[demxety.1k.fst$Fst<0]<-0
demxety.1k.fst$position=demxety.1k.fst$start
demxety.1k.fst=demxety.1k.fst[demxety.1k.fst$start< 688640245,]

demxety.5k.fst=read.table("fst/demoxety.5k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxety.5k.fst$Fst[demxety.5k.fst$Fst<0]<-0
demxety.5k.fst$position=demxety.5k.fst$start
demxety.5k.fst=demxety.5k.fst[demxety.5k.fst$start< 688640245,]



###########funtions
##### ATAC



ATAC_genomewide_plot_final<-function(dataf,original,unique1,unique2,factor,color1,color2,gene,log,spleen,limit,chest){  
  dataf$color=ifelse(dataf$log2FoldChange >0,color1,color2)
  dataf=dataf[dataf$padj<0.05,]
  
  species1=original[rownames(original) %in% rownames(dataf[dataf$color==color1,]),][,c(1,2,3)]
  rownames(species1)=str_split_fixed(rownames(species1),"-",2)[,1]
  species1=rbind(species1,unique1)
  #species1=rbind(species1,noHomo1)
  species1=species1/factor[1:3]
  species1$aver=rowMeans(species1)
  species1$start=as.numeric(str_split_fixed(rownames(species1),"_",3)[,2])
  species1$end=as.numeric(str_split_fixed(rownames(species1),"_",3)[,3])
  
  
  
  if(chest=="yes"){species2=original[rownames(original) %in% rownames(dataf[dataf$color==color2,]),][,c(4,5)]
  }else{species2=original[rownames(original) %in% rownames(dataf[dataf$color==color2,]),][,c(4,5,6)]}
  rownames(species2)=str_split_fixed(rownames(species2),"-",2)[,2]
  species2=rbind(species2,unique2)
  #species2=rbind(species2,noHomo2)
  
  if(chest=="yes"){species2=species2/factor[4:5]
  }else{species2=species2/factor[4:6]}
  species2$aver=rowMeans(species2)
  species2$start=as.numeric(str_split_fixed(rownames(species2),"_",3)[,2])
  species2$end=as.numeric(str_split_fixed(rownames(species2),"_",3)[,3])
  
  
  
  stepp=1000000
  seqe=seq(1,688640245,by=stepp)
  df=data.frame(pos=seqe)
  df$pos2=df$pos+stepp-1
  #df$pos2=df$pos+(stepp*2)-1
  
  countAVE<-function(file,start,end){
    return(mean(file[(file[5]>=start & file[6]<=end),4],na.rm=TRUE))
  }
  countAVE_ches<-function(file,start,end){
    return(mean(file[(file[4]>=start & file[5]<=end),3],na.rm=TRUE))
  }
  
  if(chest=="yes"){
    df$AVEX=apply(df,1,function(x) countAVE(species1,x[1],x[2]))
    df$AVEY=apply(df,1,function(x) countAVE_ches(species2,x[1],x[2]))
    df[is.na(df)]<-0
  }else{  
    df$AVEX=apply(df,1,function(x) countAVE(species1,x[1],x[2]))
    df$AVEY=apply(df,1,function(x) countAVE(species2,x[1],x[2]))
    df[is.na(df)]<-0
  }
if (log=="yes"){
  df$AVEX=log(df$AVEX)
  df$AVEY=log(df$AVEY)
  df$AVEX[which(!is.finite(df$AVEX))]<-0
  df$AVEY[which(!is.finite(df$AVEY))]<-0
}
  
  df$diff=df$AVEX-df$AVEY
  df$col=ifelse(df$diff>0,color1,color2)
  
  
  
  
  spline_list=list()
  
  for (i in seq(1,nrow(chrom.SE))){
    df.c=df[df$pos>chrom.SE$start[i] & df$pos<chrom.SE$end[i],]
    #spline_list[[i]]<-as.data.frame(spline(df.c$pos[df.c$diff!=0], df.c$diff[df.c$diff!=0],n=15))
    spline_list[[i]]<-as.data.frame(spline(df.c$pos,df.c$diff ,n=15))
  }   
  
  #ys=max(abs(df$diff))+2
  ys=limit
  plot=ggplot(df)+
    geom_rect(data=chrom.SE, mapping=aes(xmin=start, xmax=end, ymin=-ys, ymax=ys),fill=chrom.SE$col,alpha=0.5)+
    geom_bar(aes(x=pos,y=diff),fill=df$col,alpha=1,size=0.5,stat="identity", position="dodge")+
    #geom_bar(aes(x=pos,y=AVEX),fill=color1,alpha=1,size=0.5,stat="identity", position="dodge")+
    #geom_bar(aes(x=pos,y=-AVEY),fill=color2,alpha=1,size=0.5,stat="identity", position="dodge")+
    {if (gene=="yes")geom_point(data=geneslist,aes(x=V5,y=ys-((ys/10)*3)))}+
    {if (gene=="yes")geom_text(data=geneslist,aes(x=V5,y=ys-((ys/10)*7)),label=geneslist$V1,hjust=0, vjust=-2)}+
    {if (gene=="yes")geom_text(data=chrom.numb,aes(x=b,y=ys-((ys/10)*1)),label=chrom.numb$a,size=5,alpha=0.5)}+
    geom_segment(data=chrom.SE,aes(x = start, y = 0, xend = end, yend = 0))+
    scale_y_continuous(name="Average(reads count)",expand = c(0, 0),label=abs)+scale_x_continuous(name="",expand = c(0, 0))+
    coord_cartesian(ylim=c(-limit, limit))+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank())
  
  
  
  if (spleen=="yes"){
  for (i in seq(1,length(spline_list))){
    plot=plot+geom_link2(data = spline_list[[i]], aes(x = x, y = y,colour = after_stat(y < 0)),size=0.5)
  }
  plot=plot+scale_colour_manual( values = c(color1, color2))+ theme(legend.position = "none")
  }
  
  return(plot)
}





#formula = y ~ poly(x, 25)

Fst_genomewide_plot<-function(datase,datase1k){
plot=ggplot(datase)+
  geom_rect(data=chrom.SE, mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1),fill=chrom.SE$col,alpha=0.5)+
  #geom_smooth(data=datase1k,aes(x=position,y=Fst), method = "glm", formula = y ~ splines::bs(x, 15),se = FALSE,color="gray60")+
  geom_point(data=datase1k,aes(x=position,y=Fst),col="gray60",size=0.5)+
  geom_point(aes(x=position,y=Fst),col="black",size=0.5)+
  geom_segment(data=chrom.SE,aes(x = start, y = 0, xend = end, yend = 0))+
  scale_y_continuous(name="Fst",expand = c(0, 0),limits = c(0,1),label=abs)+scale_x_continuous(name="",expand = c(0, 0))+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#spline_list=list()

#for (i in seq(1,nrow(chrom.SE))){
#  df.c=datase1k[datase1k$position>chrom.SE$start[i] & datase1k$position<chrom.SE$end[i],]
#  spline_list[[i]]<-as.data.frame(spline(df.c$position,df.c$Fst ,n=10))
#}

#for (i in seq(1,length(spline_list))){
#  plot=plot+geom_link2(data = spline_list[[i]], aes(x = x, y = y,colour = after_stat(y < 0)),col="gray60")
#}

return(plot)
}

sweep_genomewide_plot<-function(sweep,sweep2,color1,color2,axis){
  plot=ggplot(sweep)+
    geom_rect(data=chrom.SE, mapping=aes(xmin=start, xmax=end, ymin=-axis, ymax=axis),fill=chrom.SE$col,alpha=0.5)+
    geom_point(aes(x=V2,y=V3),col=color1,size=0.5)+
    geom_point(data=sweep2,aes(x=V2,y=-V3),col=color2,size=0.5)+
    geom_hline(yintercept=quantile(sweep$V3,0.999),linetype = "dashed")+
    geom_hline(yintercept=-quantile(sweep2$V3,0.999),linetype = "dashed")+
    geom_segment(data=chrom.SE,aes(x = start, y = 0, xend = end, yend = 0))+
    scale_y_continuous(name="CLR", expand = c(0, 0),limits = c(-axis,axis),label=abs)+scale_x_continuous(name="",expand = c(0, 0))+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  return(plot)
}



genename<- function(bk,gene,chrom){
ys=0  
g=ggplot(geneslist)+
  {if(bk=="yes")geom_rect(data=chrom.SE, mapping=aes(xmin=start, xmax=end, ymin=-0, ymax=0.7),fill=chrom.SE$col,alpha=0.5)}+
   {if(gene=="yes") geom_point(data=geneslist,aes(x=V5,y=ys+0.4))}+
  #{if(gene=="yes") geom_text(data=geneslist,aes(x=V5,y=ys+0.4),label=geneslist$V1,hjust=0, vjust=0)} +
  {if(chrom=="yes")geom_text(data=chrom.numb,aes(x=b,y=ys+0.4),label=chrom.numb$a,size=5,alpha=1)}+
  
  scale_y_continuous(name="genes", expand = c(0, 0),limits=c(0,0.7))+scale_x_continuous(name="",expand = c(0, 0),limits = c(0,688640245))+
  theme_void()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
return(g)  
}
##########call function


################################################################
########### notabilis vs etylus
chrom=genename("yes","no","no")
gene=genename("no","yes","no")
GW.ATAC.FW=ATAC_genomewide_plot_final(notavsety.FW3of3[[1]],notavety3of3,notabilis.FW.unique,etylus.FW.unique,notabilis.FW.3of3.NoHomo_unique.count,etylus.FW.3of3.NoHomo_unique.count,notavsety.FW3of3[[2]],"#009e73","#d55e00","no","no","no",500)
GW.ATAC.HW=ATAC_genomewide_plot_final(notavsety.HW3of3[[1]],notavety3of3.HW,notabilis.HW.unique,etylus.HW.unique,notabilis.HW.3of3.NoHomo_unique.count,etylus.HW.3of3.NoHomo_unique.count,notavsety.HW3of3[[2]],"#009e73","#d55e00","no","no","no",500)
GW.sweep=sweep_genomewide_plot(sweep.not,sweep.ety,"#009e73","#d55e00",300)
GW_fst=Fst_genomewide_plot(notaxety.fst,notaxety.5k.fst)
#GW.plot=plot_grid(chrom,GW.ATAC.FW,GW.ATAC.HW,GW.sweep,GW_fst,gene,ncol = 1, align = "v",rel_heights= c(0.2,1,1,1,1,0.2))
#ggsave(GW.plot,file="newplotffolder_sonobodycomplains/notabilisVSeylus.GW.Feb2023.ultimate.text4.png",device="png",dpi=300, width = 12, height = 8,bg="white")
#ggsave(GW.plot,file="newplotffolder_sonobodycomplains/notabilisVSeylus.GW.Feb2023.ultimate.pdf",device="pdf",dpi=300, width = 12, height = 8)


########### demophoon vs hydara
GW.ATAC.FW.2=ATAC_genomewide_plot_final(demoxhyda.FW3of3[[1]],demoxhyda3of3,demophoon.FW.unique,hydara.FW.unique,demophoon.FW.3of3.NoHomo_unique.count,hydara.FW.3of3.NoHomo_unique.count,demoxhyda.FW3of3[[2]],"#0072b2","#cc79a7","yes","no","no",500,"no")
GW.ATAC.HW.2=ATAC_genomewide_plot_final(demoxhyda.HW3of3[[1]],demoxhyda3of3.HW,demophoon.HW.unique,hydara.HW.unique,demophoon.HW.3of3.NoHomo_unique.count,hydara.HW.3of3.NoHomo_unique.count,demoxhyda.HW3of3[[2]],"#0072b2","#cc79a7","no","no","no",500)
GW.sweep.2=sweep_genomewide_plot(sweep.dem,sweep.hyd,"#0072b2","#cc79a7",200)
GW_fst.2=Fst_genomewide_plot(demxhd.fst,demxhd.5k.fst)
#GW.plot.2=plot_grid(chrom,GW.ATAC.FW,GW.ATAC.HW,GW.sweep,GW_fst,gene,ncol = 1, align = "v",rel_heights= c(0.2,1,1,1,1,0.2))
#ggsave(GW.plot,file="newplotffolder_sonobodycomplains/demophoonVShydara.GW.Feb2023.ultimate.png",device="png",dpi=300, width = 12, height = 8,bg="white")
#ggsave(GW.plot,file="newplotffolder_sonobodycomplains/demophoonVShydara.GW.Feb2023.ultimate.pdf",device="pdf",dpi=300, width = 12, height = 8,bg="white")






GW.ATAC.FW.2=GW.ATAC.FW.2+theme(plot.margin = unit(c(0.2,0,0,0), "cm"),axis.text.y=element_text(size=5))
GW.ATAC.HW.2=GW.ATAC.HW.2+theme(plot.margin = unit(c(0.2,0,0,0), "cm"),axis.text.y=element_text(size=5))
GW.sweep.2=GW.sweep.2+theme(plot.margin = unit(c(0.2,0,0,0), "cm"),axis.text.y=element_text(size=5))
GW_fst.2=GW_fst.2+theme(plot.margin = unit(c(0.2,0,0,0), "cm"),axis.text.y=element_text(size=5))
GW.ATAC.FW=GW.ATAC.FW+theme(plot.margin = unit(c(0.2,0,0,0), "cm"),axis.text.y=element_text(size=5))
GW.ATAC.HW=GW.ATAC.HW+theme(plot.margin = unit(c(0.2,0,0,0), "cm"),axis.text.y=element_text(size=5))
GW.sweep=GW.sweep+theme(plot.margin = unit(c(0.2,0,0,0), "cm"),axis.text.y=element_text(size=5))
GW_fst=GW_fst+theme(plot.margin = unit(c(0.2,0,0,0), "cm"),axis.text.y=element_text(size=5))
gene=gene+theme(plot.margin = unit(c(0,0,0,0), "cm"))

super_megaplot=plot_grid(chrom,GW.ATAC.FW.2,GW.ATAC.HW.2,GW.sweep.2,GW_fst.2,GW.ATAC.FW,GW.ATAC.HW,GW.sweep,GW_fst,gene,ncol = 1, align = "v",rel_heights= c(0.2,rep(1,8),0.2))


ggsave(super_megaplot,file="newplotffolder_sonobodycomplains/supermega.GW.Feb2023.ultimate.png",device="png",dpi=300, width = 12, height = 8,bg="white")







########### demophoon vs chest



GW_fst_chest=Fst_genomewide_plot(demxches.fst,demxches.5k.fst)




########### demophoon vs favorinus


GW.ATAC.FW_fav=ATAC_genomewide_plot_final(demoxfavo.FW3of3[[1]],demvsfavo.FW,demophoonxfav.FW.unique,favorinus.FW.unique,demoxfavo.FW3of3[[2]],"#0072b2","#661100","no","no","no",500,"no")

GW.ATAC.HW_fav=ATAC_genomewide_plot_final(demoxfavo.HW3of3[[1]],demvsfavo.HW,demophoonxfav.HW.unique,favorinus.HW.unique,demoxfavo.HW3of3[[2]],"#0072b2","#332288","no","no","no",500,"no")

GW_fst_fav=Fst_genomewide_plot(demxfav.fst,demxfav.5k.fst)

GW.sweep.fav=sweep_genomewide_plot(sweep.dem,sweep.favo,"#0072b2","red",200)

GW.plot.fav=plot_grid(GW.ATAC.FW_fav,GW.ATAC.HW_fav,GW.sweep.fav,GW_fst_fav,ncol = 1, align = "v")

ggsave(GW.plot.fav,file="newplotffolder_sonobodycomplains/demophoonVSfavor.GW.Feb2023.ultimate.png",device="png",dpi=300, width = 12, height = 8,bg="white")




############# demophoon vs chestertonii

GW.ATAC.FW_chest=ATAC_genomewide_plot_final(demoxchest.FW3of3[[1]],demvschesto.FW,demophoonxchest.FW.unique,chesteronii.FW.unique,demoxchest.FW3of3[[2]],"#0072b2","#E69F00","no","no","no",500,"yes")

GW.ATAC.HW_chest=ATAC_genomewide_plot_final(demoxchest.HW3of3[[1]],demvschesto.HW,demophoonxchest.HW.unique,chesteronii.HW.unique,demoxchest.HW3of3[[2]],"#0072b2","#E69F00","no","no","no",500,"yes")

GW_fst_chest=Fst_genomewide_plot(demxches.fst,demxches.5k.fst)

GW.sweep.chest=sweep_genomewide_plot(sweep.dem,sweep.chest,"#0072b2","goldenrod1",500)

GW.plot.chest=plot_grid(GW.ATAC.FW_chest,GW.ATAC.HW_chest,GW.sweep.chest,GW_fst_chest,ncol = 1, align = "v")

ggsave(GW.plot.chest,file="newplotffolder_sonobodycomplains/demophoonVSchest.GW.Feb2023.ultimate.png",device="png",dpi=300, width = 12, height = 8,bg="white")


############## demophoon etylus



GW.ATAC.FW_ety=ATAC_genomewide_plot_final(demoxety.FW3of3[[1]],demvsety.FW,demophoonxety.FW.unique,etylus.FW.unique,demoxety.FW3of3[[2]],"#0072b2","#d55e00","no","no","no",500,"no")

GW.ATAC.HW_ety=ATAC_genomewide_plot_final(demoxety.HW3of3[[1]],demvsety.HW,demophoonxety.HW.unique,etylus.HW.unique,demoxety.HW3of3[[2]],"#0072b2","#d55e00","no","no","no",500,"no")

GW_fst_ety=Fst_genomewide_plot(demxety.fst,demxety.5k.fst)

GW.sweep.ety=sweep_genomewide_plot(sweep.dem,sweep.ety,"#0072b2","#d55e00",200)

GW.plot.ety=plot_grid(GW.ATAC.FW_ety,GW.ATAC.HW_ety,GW.sweep.ety,GW_fst_ety,ncol = 1, align = "v")


ggsave(GW.plot.ety,file="newplotffolder_sonobodycomplains/demophoonVSetylus.GW.Feb2023.ultimate.png",device="png",dpi=300, width = 12, height = 8,bg="white")



########################################
GW.ATAC.FW_hyd=ATAC_genomewide_plot_final(demoxhyda.FW3of3[[1]],demoxhyda3of3,demophoon.FW.unique,hydara.FW.unique,demoxhyda.FW3of3[[2]],"#0072b2","#cc79a7","no","no","no",500,"no")
GW_fst.2=GW_fst.2+theme(axis.title.y=element_blank())
GW_fst_fav=GW_fst_fav+theme(axis.title.y=element_blank())
GW_fst_ety=GW_fst_ety+theme(axis.title.y=element_blank())
GW_fst_chest=GW_fst_chest+theme(axis.title.y=element_blank())

otherfigure=plot_grid(chrom,GW.ATAC.FW_hyd,GW.ATAC.FW_fav,GW.ATAC.FW_ety,GW.ATAC.FW_chest,GW_fst.2,GW_fst_fav,GW_fst_ety,GW_fst_chest,gene,ncol = 1, align = "v",rel_heights= c(0.2,rep(1,8),0.2))
ggsave(otherfigure,file="newplotffolder_sonobodycomplains/correfigure.part1-2.png",device="png",dpi=300, width = 5, height = 7,bg="white")



#######################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
#ATAC_genomewide_plot_final<-function(dataf,original,unique1,unique2,factor,color1,color2,gene,log,spleen,limit,chest){

  dataf=demoxchest.FW3of3[[1]]
    original=demvschesto.FW
    unique1=demophoonxchest.FW.unique
    unique2=chesteronii.FW.unique
    factor=demoxchest.FW3of3[[2]]
    color1="#0072b2"
    color2="red"
    gene="yes"
    log="no"
    spleen="no"
    limit=500
    chest="yes"
    
  dataf$color=ifelse(dataf$log2FoldChange >0,color1,color2)
  dataf=dataf[dataf$padj<0.05,]
  
  species1=original[rownames(original) %in% rownames(dataf[dataf$color==color1,]),][,c(1,2,3)]
  rownames(species1)=str_split_fixed(rownames(species1),"-",2)[,1]
  species1=rbind(species1,unique1)
  #species1=rbind(species1,noHomo1)
  species1=species1/factor[1:3]
  species1$aver=rowMeans(species1)
  species1$start=as.numeric(str_split_fixed(rownames(species1),"_",3)[,2])
  species1$end=as.numeric(str_split_fixed(rownames(species1),"_",3)[,3])
  
  
  
  if(chest=="yes"){species2=original[rownames(original) %in% rownames(dataf[dataf$color==color2,]),][,c(4,5)]
  }else{species2=original[rownames(original) %in% rownames(dataf[dataf$color==color2,]),][,c(4,5,6)]}
  rownames(species2)=str_split_fixed(rownames(species2),"-",2)[,2]
  species2=rbind(species2,unique2)
  #species2=rbind(species2,noHomo2)
  
  if(chest=="yes"){species2=species2/factor[4:5]
  }else{species2=species2/factor[4:6]}
  species2$aver=rowMeans(species2)
  species2$start=as.numeric(str_split_fixed(rownames(species2),"_",3)[,2])
  species2$end=as.numeric(str_split_fixed(rownames(species2),"_",3)[,3])
  
  
  
  stepp=1000000
  seqe=seq(1,688640245,by=stepp)
  df=data.frame(pos=seqe)
  df$pos2=df$pos+stepp-1
  #df$pos2=df$pos+(stepp*2)-1
  
  countAVE<-function(file,start,end){
    return(mean(file[(file[5]>=start & file[6]<=end),4],na.rm=TRUE))
  }
  countAVE_ches<-function(file,start,end){
    return(mean(file[(file[4]>=start & file[5]<=end),3],na.rm=TRUE))
  }
  
  if(chest=="yes"){
    df$AVEX=apply(df,1,function(x) countAVE(species1,x[1],x[2]))
    df$AVEY=apply(df,1,function(x) countAVE_ches(species2,x[1],x[2]))
    df[is.na(df)]<-0
  }else{  
  df$AVEX=apply(df,1,function(x) countAVE(species1,x[1],x[2]))
  df$AVEY=apply(df,1,function(x) countAVE(species2,x[1],x[2]))
  df[is.na(df)]<-0
  }
  if (log=="yes"){
    df$AVEX=log(df$AVEX)
    df$AVEY=log(df$AVEY)
    df$AVEX[which(!is.finite(df$AVEX))]<-0
    df$AVEY[which(!is.finite(df$AVEY))]<-0
  }
  
  df$diff=df$AVEX-df$AVEY
  df$col=ifelse(df$diff>0,color1,color2)
  
  
  
  
  spline_list=list()
  
  for (i in seq(1,nrow(chrom.SE))){
    df.c=df[df$pos>chrom.SE$start[i] & df$pos<chrom.SE$end[i],]
    #spline_list[[i]]<-as.data.frame(spline(df.c$pos[df.c$diff!=0], df.c$diff[df.c$diff!=0],n=15))
    spline_list[[i]]<-as.data.frame(spline(df.c$pos,df.c$diff ,n=15))
  }   
  
  #ys=max(abs(df$diff))+2
  ys=limit
  plot=ggplot(df)+
    geom_rect(data=chrom.SE, mapping=aes(xmin=start, xmax=end, ymin=-ys, ymax=ys),fill=chrom.SE$col,alpha=0.5)+
    geom_bar(aes(x=pos,y=diff),fill=df$col,alpha=1,size=0.5,stat="identity", position="dodge")+
    #geom_bar(aes(x=pos,y=AVEX),fill=color1,alpha=1,size=0.5,stat="identity", position="dodge")+
    #geom_bar(aes(x=pos,y=-AVEY),fill=color2,alpha=1,size=0.5,stat="identity", position="dodge")+
    {if (gene=="yes")geom_point(data=geneslist,aes(x=V5,y=ys-((ys/10)*3)))}+
    {if (gene=="yes")geom_text(data=geneslist,aes(x=V5,y=ys-((ys/10)*7)),label=geneslist$V1,hjust=0, vjust=-2)}+
    {if (gene=="yes")geom_text(data=chrom.numb,aes(x=b,y=ys-((ys/10)*1)),label=chrom.numb$a,size=5,alpha=0.5)}+
    geom_segment(data=chrom.SE,aes(x = start, y = 0, xend = end, yend = 0))+
    scale_y_continuous(name="Average(reads count)",expand = c(0, 0),label=abs)+scale_x_continuous(name="",expand = c(0, 0))+
    coord_cartesian(ylim=c(-limit, limit))+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank())
  
  
  
  if (spleen=="yes"){
    for (i in seq(1,length(spline_list))){
      plot=plot+geom_link2(data = spline_list[[i]], aes(x = x, y = y,colour = after_stat(y < 0)),size=0.5)
    }
    plot=plot+scale_colour_manual( values = c(color1, color2))+ theme(legend.position = "none")
  }
  
plot
