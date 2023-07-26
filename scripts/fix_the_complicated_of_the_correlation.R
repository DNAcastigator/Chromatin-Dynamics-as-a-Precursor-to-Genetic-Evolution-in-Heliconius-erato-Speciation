library(tidyverse)
library(cowplot)
###########################################################################################################################################################

freq_fun<-function(dataf,position,dat){
  if(dat=="Fst"){
    freq.df=dataf[,c(2,3,9)]
  }else{
    freq.df=dataf[,c(2,1,3)]
    freq.df$V1=dataf$V2+999
    colnames(freq.df)=c("start","end","Fst")
    
  }
  freq.df$frq=apply(freq.df,1,function(x){
    len=length(position[position>x[1] & position<x[2]])
    #/(x[2]-x[1])
    return(len)  
  })
  return(freq.df)
}
#########################

Fst.dist=function(data,limit,sigla){
    uno=na.omit(data[data$Fst>=limit[1] & data$Fst<limit[2],])
  due=data.frame(chrom="pan",start=uno$start,end=uno$end)
  due$sigla=sigla
  return(due)
}


#################

Fst.dist.loop=function(data){
 
  fst.dist.1=Fst.dist(data,c(0,0.25),"0-0.25")
  fst.dist.2=Fst.dist(data,c(0.25,0.50),"0.25-0.50")
  fst.dist.3=Fst.dist(data,c(0.50,0.75),"0.50-0.75")
  fst.dist.4=Fst.dist(data,c(0.75,1),"0.75-1")
  
  fst.dist.total=do.call("rbind", list(fst.dist.1,fst.dist.2,fst.dist.3,fst.dist.4))
  fst.dist.total$chrom="pan"
  return(fst.dist.total)
}

Fst.dist.loop=function(data){
  fst.dist.1=Fst.dist(data,c(0,0.1),"0-0.1")
  fst.dist.2=Fst.dist(data,c(0.1,0.2),"0.1-0.2")
  fst.dist.3=Fst.dist(data,c(0.2,0.3),"0.2-0.3")
  fst.dist.4=Fst.dist(data,c(0.3,0.4),"0.3-0.4")
  fst.dist.5=Fst.dist(data,c(0.4,0.5),"0.4-0.5")
  fst.dist.6=Fst.dist(data,c(0.5,0.6),"0.5-0.6")
  fst.dist.7=Fst.dist(data,c(0.6,0.7),"0.6-0.7")
  fst.dist.8=Fst.dist(data,c(0.7,0.8),"0.7-0.8")
  fst.dist.9=Fst.dist(data,c(0.8,0.9),"0.8-0.9")
  fst.dist.10=Fst.dist(data,c(0.9,1),"0.9-1")
  
  
  fst.dist.total=do.call("rbind", list(fst.dist.1,fst.dist.2,fst.dist.3,fst.dist.4,fst.dist.5,fst.dist.6,fst.dist.7,fst.dist.8,fst.dist.9,fst.dist.10))
  fst.dist.total$chrom="pan"
  return(fst.dist.total)
}




Fst.stat=function(fsts,position,total_pos,n,totplot){
  
position.range=makeGRangesFromDataFrame(DataFrame(chrom="pan",start=position,end=position+1))
position.total=makeGRangesFromDataFrame(DataFrame(chrom="pan",start=total_pos,end=total_pos+1))
test=makeGRangesFromDataFrame(fsts,keep.extra.columns=T)

test$over=countOverlaps(test,position.range)
test$over_total=countOverlaps(test,position.total)
resu=data.frame(test) %>% group_by(sigla) %>% summarise(final_count = sum(over),total=sum(over_total)) %>% data.frame()


random=data.frame(base::transform(test,sigla=sample(test$sigla)))%>% group_by(sigla) %>% summarise(final_count = sum(over)/sum(over_total)) %>% data.frame()
#random2=data.frame(random) %>% group_by(sigla) %>% summarise(final_count = sum(over)) %>% data.frame()
for (i in seq(1:n)){
ran=data.frame(base::transform(test,sigla=sample(test$sigla)))%>% group_by(sigla) %>% summarise(final_count = sum(over)/sum(over_total)) %>% data.frame()
random=rbind(random,ran)
}

resu_random=na.omit(random) %>% group_by(sigla) %>% summarise(random_count = mean(final_count),se=sqrt(sd(final_count)/n())) %>% data.frame()
#resu_random=na.omit(random) %>% group_by(sigla) %>% summarise(random_count = mean(final_count),se=sd(final_count)) %>% data.frame()
################

confronto=cbind(resu,resu_random[c(2,3)])
confronto$random_count[!is.finite(confronto$random_count)] <- 0

confronto=confronto[confronto$total>0,]
statistic=apply(confronto,1,function(x){
extest=binom.test(as.numeric(x[2]), as.numeric(x[3]), p = as.numeric(x[4]),
           alternative =  "g",
           conf.level = 0.95)

#output=data.frame(p_value=extest$p.value,con_int_low=extest$conf.int[1],con_int_hig=extest$conf.int[2])


return(extest$p.value)
})
#confronto=cbind(confronto,do.call("rbind",statistic))
confronto=cbind(confronto,statistic)
#confronto$statistic=confronto$statistic*10

print(confronto)
g=ggplot(confronto)+
  geom_violin(data=random, aes(x=sigla, y=final_count)) + 
  geom_point(aes(x=sigla,y=final_count/total),col="red",shape=ifelse(p.adjust(confronto$statisti,method = "holm")<0.05,8,16))+
  #geom_point(aes(x=sigla,y=random_count),col="blue")+
  #geom_errorbar(aes(x=sigla,ymin = random_count-se, ymax = random_count+se),col="blue")+
  {if(totplot=="yes") ylim(c(0,1)) }+
  theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1),axis.title.y=element_blank())
return(list(g,confronto))
}



################################################### rank plot
rank_plot<-function(dataf,position,freq.df,num,dat){
  if (dat=="sweep"){
    end=tail(dataf$V2,1)+999
    dataf=dataf[order(dataf$V2),]
    dataf=base::transform(dataf, V5 = c(dataf$V2[-1]-1, end))
    subset=lapply(position,function(x){
      res=dataf$V3[x>dataf$V2 & x<dataf$V5] 
      
      return(res)
      
    })
  }else{
  subset=lapply(position,function(x){
    res=dataf$Fst[x>dataf$start & x<dataf$end] 
    
    return(res)
    
  })
  }
  subset=do.call(rbind.data.frame, subset)
  subset=na.omit(subset)
  colnames(subset)="Fst"
  subset$sigla="test"
  subset=subset[order(subset$Fst),]
  subset$rank=seq(1:nrow(subset))
  
  n=nrow(subset)
  
  
  
  #random=data.frame(Fst=sample(freq.df$Fst,size=n,prob=freq.df$frq,replace=TRUE))
  #random$sigla="random"
  ############
  freq.df=freq.df[!is.na(freq.df$Fst),]
  total=data.frame(Fst=sample(freq.df$Fst,size=n,prob=freq.df$frq,replace=TRUE),rank=0)
  total=total[order(total$Fst),]
  total$rank=seq(1:nrow(total))
  
  for (i in seq(1:num)){
      single=data.frame(Fst=sample(freq.df$Fst,size=n,prob=freq.df$frq,replace=TRUE),rank=0)
      single=single[order(single$Fst),]
      single$rank=seq(1:nrow(single))
      total=rbind(total,single)
    
  }
  rankandom=data.frame(total %>% group_by(rank) %>% summarise(Fst = mean(Fst)),sigla="random")
  random=rankandom[,c(2,3,1)]
  #random$sigla="random"
  
  ###############
  final=rbind(subset,random)
  
  
  ranksweep2=subset %>% slice(-seq(1, (n()/10)*9.5))
  
  rankandom2=rankandom %>% slice(-seq(1, (n()/10)*9.5))
  
  print(ks.test(subset$Fst,rankandom$Fst,alternative="less",simulate.p.value =TRUE))
  print(nrow(rankandom))
  print(nrow(subset))
  #f=ggplot(rbind(ranksweep2,rankandom2))+geom_line(aes(x=as.numeric(rank),y=Fst,col=sigla))+geom_point(aes(x=as.numeric(rank),y=Fst,col=sigla))+theme_classic()
  f=ggplot(rbind(subset,rankandom))+geom_line(aes(x=as.numeric(rank),y=Fst,col=sigla))+geom_point(aes(x=as.numeric(rank),y=Fst,col=sigla))+theme_classic()
  if(dat=="sweep") {f=f+ylab("CLR")}
  oupput=full_join(subset,rankandom,by="rank")
  return(list(f,oupput))
}





Fst.dist.sweep=function(data,limit,sigla){
  uno=na.omit(data[data$V3>=limit[1] & data$V3<limit[2],])
  due=data.frame(chrom="pan",start=uno$V2,end=uno$V5)
  due$sigla=sigla
  return(due)
}


#################

Fst.dist.loop_swep=function(data){
  data$V3=data$V3/max(data$V3)
  data=data[order(data$V2),]
  end=tail(data$V2,1)+999
  print(head(data))
  data=base::transform(data, V5 = c(data$V2[-1]-1, end))
  print(head(data))
  fst.dist.1=Fst.dist.sweep(data,c(0,0.1),"0-0.1")
  fst.dist.2=Fst.dist.sweep(data,c(0.1,0.2),"0.1-0.2")
  fst.dist.3=Fst.dist.sweep(data,c(0.2,0.3),"0.2-0.3")
  fst.dist.4=Fst.dist.sweep(data,c(0.3,0.4),"0.3-0.4")
  fst.dist.5=Fst.dist.sweep(data,c(0.4,0.5),"0.4-0.5")
  fst.dist.6=Fst.dist.sweep(data,c(0.5,0.6),"0.5-0.6")
  fst.dist.7=Fst.dist.sweep(data,c(0.6,0.7),"0.6-0.7")
  fst.dist.8=Fst.dist.sweep(data,c(0.7,0.8),"0.7-0.8")
  fst.dist.9=Fst.dist.sweep(data,c(0.8,0.9),"0.8-0.9")
  fst.dist.10=Fst.dist.sweep(data,c(0.9,1),"0.9-1")
  
  
  fst.dist.total=do.call("rbind", list(fst.dist.1,fst.dist.2,fst.dist.3,fst.dist.4,fst.dist.5,fst.dist.6,fst.dist.7,fst.dist.8,fst.dist.9,fst.dist.10))
  fst.dist.total$chrom="pan"
  return(fst.dist.total)
}









######################################################
############################################################ generate data

########complete nota ety
notabilisxetylus.complete=data.frame(str_split_fixed(c(rownames(notabilisFWall),rownames(notabilisHWall),rownames(etylusFWall),rownames(etylusHWall)),"_",3))
colnames(notabilisxetylus.complete)=c("chrom","start","end")
notabilisxetylus.complete$start=as.integer(as.character(notabilisxetylus.complete$start))
notabilisxetylus.complete$end=as.integer(as.character(notabilisxetylus.complete$end))

notxety_total_pos=notabilisxetylus.complete$start+((notabilisxetylus.complete$end-notabilisxetylus.complete$start)/2)


###################only nota FW
data=notabilisFW.3of3

data$start=as.numeric(str_split_fixed(rownames(data),"_",3)[,2])
data$end=as.numeric(str_split_fixed(rownames(data),"_",3)[,3])
position=data$start+((data$end-data$start)/2)

#position=notxety_total_pos
position=data.frame(test[test$sigla=="0.75-1" & test$over_total>0]) %>% summarise(pos=start+((end-start)/2));position=position$pos


#################real DE
position1=notavsety.FW3of3[[1]]$pos[notavsety.FW3of3[[1]]$padj <0.05]
position2=notavsety.FW3of3[[1]]$pos2[notavsety.FW3of3[[1]]$padj <0.05]
position3=notavsety.HW3of3[[1]]$pos[notavsety.HW3of3[[1]]$padj <0.05]
position4=notavsety.HW3of3[[1]]$pos2[notavsety.HW3of3[[1]]$padj <0.05]
position_DE=c(position1+((position2-position1)/2),position3+((position4-position3)/2))

#################### unique

position1=notabilis.HW.unique.bed$start+((notabilis.HW.unique.bed$end-notabilis.HW.unique.bed$end)/2)
position2=notabilis.FW.unique.bed$start+((notabilis.FW.unique.bed$end-notabilis.FW.unique.bed$end)/2)
position3=etylus.HW.unique.bed$start+((etylus.HW.unique.bed$end-etylus.HW.unique.bed$end)/2)
position4=etylus.FW.unique.bed$start+((etylus.FW.unique.bed$end-etylus.FW.unique.bed$end)/2)
#position_unique=c(position1,position2,position3,position4)
position_Deunique_notaxety=c(position1,position2,position3,position4,position_DE)
#position=c(position1,position2)
###########################################################################test the functions
#freq.df_notaetyFst_1k=freq_fun(notaxety.1k.fst,notxety_total_pos,"Fst")

rank1=rank_plot(notaxety.1k.fst,position_Deunique_notaxety,freq.df_notaetyFst_1k,1000,"Fst")

fst.dist.notacety=Fst.dist.loop(notaxety.1k.fst)
split1=Fst.stat(fst.dist.notacety,position_Deunique_notaxety,notxety_total_pos,1000,"no")
comb1=plot_grid(rank1,split1[[1]],ncol=1)
ggsave(comb1,file="newplotffolder_sonobodycomplains/correlation/notabilisVSeylus.corr.DEunique.zoom.png",device="png",dpi=300,bg="white")








############################################### demophoon hydara

position1=demoxhyda.FW3of3[[1]]$pos[demoxhyda.FW3of3[[1]]$padj <0.05]
position2=demoxhyda.FW3of3[[1]]$pos2[demoxhyda.FW3of3[[1]]$padj <0.05]
position3=demoxhyda.HW3of3[[1]]$pos[demoxhyda.HW3of3[[1]]$padj <0.05]
position4=demoxhyda.HW3of3[[1]]$pos2[demoxhyda.HW3of3[[1]]$padj <0.05]
position_DE=c(position1+((position2-position1)/2),position3+((position4-position3)/2))

position1=demophoon.HW.unique.bed$start+((demophoon.HW.unique.bed$end-demophoon.HW.unique.bed$end)/2)
position2=demophoon.FW.unique.bed$start+((demophoon.FW.unique.bed$end-demophoon.FW.unique.bed$end)/2)
position3=hydara.HW.unique.bed$start+((hydara.HW.unique.bed$end-hydara.HW.unique.bed$end)/2)
position4=hydara.FW.unique.bed$start+((hydara.FW.unique.bed$end-hydara.FW.unique.bed$end)/2)
#position_unique=c(position1,position2,position3,position4)
position_Deunique_demohyda=c(position1,position2,position3,position4,position_DE)
#position=c(position1,position2)




demophoonxhydara.complete=data.frame(str_split_fixed(c(rownames(demophoonFWall),rownames(demophoonHWall),rownames(hydaraFWall),rownames(hydaraHWall)),"_",3))
colnames(demophoonxhydara.complete)=c("chrom","start","end")
demophoonxhydara.complete$start=as.integer(as.character(demophoonxhydara.complete$start))
demophoonxhydara.complete$end=as.integer(as.character(demophoonxhydara.complete$end))

demhyd_total_pos=demophoonxhydara.complete$start+((demophoonxhydara.complete$end-demophoonxhydara.complete$start)/2)
########################################################################### test the functions
freq.df_demohydFst_1k=freq_fun(demxhd.1k.fst,demhyd_total_pos,"Fst")

rank2=rank_plot(demxhd.1k.fst,position_Deunique_demohyda,freq.df_demohydFst_1k,1000,"Fst")

fst.dist.demohyda=Fst.dist.loop(demxhd.1k.fst)
split2=Fst.stat(fst.dist.demohyda,position_Deunique_demohyda,notxety_total_pos,1000,"no")
comb2=plot_grid(rank2,split2[[1]],ncol=1)
ggsave(comb2,file="newplotffolder_sonobodycomplains/correlation/demophoonvshydara.corr.DEunique1.zoom.png",device="png",dpi=300,bg="white")


rank11=rank1+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none")
rank21=rank2+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none")

split11=split1[[1]]+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none")
split21=split2[[1]]+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none")


comb1=plot_grid(rank11,split11,ncol=1)
comb2=plot_grid(rank21,split21,ncol=1)
supercomb=plot_grid(comb1,comb2,ncol=2)
ggsave(supercomb,file="newplotffolder_sonobodycomplains/correlation/HZs.correlation.png",device="png",dpi=300,bg="white")

################################################################
##################################################################
#######################################################
########################################### sweep

notabilisxetylus.pos=data.frame(str_split_fixed(c(rownames(notabilisFWall),rownames(notabilisHWall)),"_",3))
colnames(notabilisxetylus.pos)=c("chrom","start","end")
notabilisxetylus.pos$start=as.integer(as.character(notabilisxetylus.pos$start))
notabilisxetylus.pos$end=as.integer(as.character(notabilisxetylus.pos$end))

notxety_nota_pos=notabilisxetylus.pos$start+((notabilisxetylus.pos$end-notabilisxetylus.pos$start)/2)





freq.df_nota=freq_fun(sweep.not,notxety_total_pos,"sweep")
freq.df_ety=freq_fun(sweep.ety,notxety_total_pos,"sweep")
freq.df_demo=freq_fun(sweep.dem,demhyd_total_pos,"sweep")
freq.df_hyd=freq_fun(sweep.hyd,demhyd_total_pos,"sweep")


rank_nota=rank_plot(sweep.not,position_Deunique_notaxety,freq.df_nota,1000,"sweep")
testo_nota=Fst.dist.loop_swep(sweep.not)
split_nota=Fst.stat(testo_nota,position_Deunique_notaxety,notxety_total_pos,1000,"no")

rank_ety=rank_plot(sweep.ety,position_Deunique_notaxety,freq.df_ety,1000,"sweep")
testo_ety=Fst.dist.loop_swep(sweep.ety)
split_ety=Fst.stat(testo_ety,position_Deunique_notaxety,notxety_total_pos,1000,"no")

rank_demo=rank_plot(sweep.dem,position_Deunique_demohyda,freq.df_demo,1000,"sweep")
testo_demo=Fst.dist.loop_swep(sweep.dem)
split_demo=Fst.stat(testo_demo,position_Deunique_demohyda,demhyd_total_pos,1000,"no")

rank_hyd=rank_plot(sweep.hyd,position_Deunique_demohyda,freq.df_demo,1000,"sweep")
testo_hyd=Fst.dist.loop_swep(sweep.hyd)
split_hyd=Fst.stat(testo_hyd,position_Deunique_demohyda,demhyd_total_pos,1000,"no")

rank_nota2=rank_nota+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1),legend.position = "none")
rank_ety2=rank_ety+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1),legend.position = "none")
rank_demo2=rank_demo+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1),legend.position = "none")
rank_hyd2=rank_hyd+theme(axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1),legend.position = "none")

comb_rank=plot_grid(rank_nota2,rank_ety2,rank_demo2,rank_hyd2,nrow=1)
comb_split=plot_grid(split_nota[[1]],split_ety[[1]],split_demo[[1]],split_hyd[[1]],nrow=1)
comb=plot_grid(comb_rank,comb_split,nrow=2,rel_heights = c(1,0.5))
finalcomb=plot_grid(supercomb,comb,nrow=2)

ggsave(finalcomb,file="newplotffolder_sonobodycomplains/correlation/definitve.correlation.pdf",device="pdf",dpi=300,bg="white")




####################################### demohppon favorinus
position.demfav.1=demophoonxfav.HW.unique.bed$start+((demophoonxfav.HW.unique.bed$end-demophoonxfav.HW.unique.bed$end)/2)
position.demfav.2=demophoonxfav.FW.unique.bed$start+((demophoonxfav.FW.unique.bed$end-demophoonxfav.FW.unique.bed$end)/2)
position.demfav.3=favorinus.HW.unique.bed$start+((favorinus.HW.unique.bed$end-favorinus.HW.unique.bed$end)/2)
position.demfav.4=favorinus.FW.unique.bed$start+((favorinus.FW.unique.bed$end-favorinus.FW.unique.bed$end)/2)

position_Deunique_demofava=c(position.demfav.1,position.demfav.2,position.demfav.3,position.demfav.4)


demophoonxfavorinus.complete=data.frame(str_split_fixed(c(rownames(demophoonFWall),rownames(demophoonHWall),rownames(favorinusFWall),rownames(favorinusHWall)),"_",3))
colnames(demophoonxfavorinus.complete)=c("chrom","start","end")
demophoonxfavorinus.complete$start=as.integer(as.character(demophoonxfavorinus.complete$start))
demophoonxfavorinus.complete$end=as.integer(as.character(demophoonxfavorinus.complete$end))

demfav_total_pos=demophoonxfavorinus.complete$start+((demophoonxfavorinus.complete$end-demophoonxfavorinus.complete$start)/2)
########################################################################### test the functions

freq.df_demofavFst_1k=freq_fun(demxfav.1k.fst,demfav_total_pos,"Fst")

rank2_fav=rank_plot(demxfav.1k.fst,position_Deunique_demofava,freq.df_demofavFst_1k,1000,"Fst")

fst.dist.demofava=Fst.dist.loop(demxfav.1k.fst)
split2_fav=Fst.stat(fst.dist.demofava,position_Deunique_demofava,demfav_total_pos,1000,"no")
comb2_fav=plot_grid(rank2_fav,split2_fav[[1]],ncol=1)
ggsave(comb2_fav,file="newplotffolder_sonobodycomplains/correlation/demophoonvsfavorinus.corr.DEunique1.zoom.png",device="png",dpi=300,bg="white")




################################ demophoon vs chest

position.demchest.1=demophoonxchest.HW.unique.bed$start+((demophoonxchest.HW.unique.bed$end-demophoonxchest.HW.unique.bed$end)/2)
position.demchest.2=demophoonxchest.FW.unique.bed$start+((demophoonxchest.FW.unique.bed$end-demophoonxchest.FW.unique.bed$end)/2)
position.demchest.3=chestertonii.HW.unique.bed$start+((chestertonii.HW.unique.bed$end-chestertonii.HW.unique.bed$end)/2)
position.demchest.4=chestertonii.FW.unique.bed$start+((chestertonii.FW.unique.bed$end-chestertonii.FW.unique.bed$end)/2)

position_Deunique_demochesta=c(position.demchest.1,position.demchest.2,position.demchest.3,position.demchest.4)


demophoonxchestertonii.complete=data.frame(str_split_fixed(c(rownames(demophoonFWall),rownames(demophoonHWall),rownames(chestertoniiFWall),rownames(chestertoniiHWall)),"_",3))
colnames(demophoonxchestertonii.complete)=c("chrom","start","end")
demophoonxchestertonii.complete$start=as.integer(as.character(demophoonxchestertonii.complete$start))
demophoonxchestertonii.complete$end=as.integer(as.character(demophoonxchestertonii.complete$end))

demchest_total_pos=demophoonxchestertonii.complete$start+((demophoonxchestertonii.complete$end-demophoonxchestertonii.complete$start)/2)
########################################################################### test the functions

freq.df_demochestFst_1k=freq_fun(demxches.1k.fst,demchest_total_pos,"Fst")

rank2_chest=rank_plot(demxches.1k.fst,position_Deunique_demochesta,freq.df_demochestFst_1k,1000,"Fst")

fst.dist.demochesta=Fst.dist.loop(demxches.1k.fst)
split2_chest=Fst.stat(fst.dist.demochesta,position_Deunique_demochesta,demchest_total_pos,1000,"no")
comb2_chest=plot_grid(rank2_chest,split2_chest[[1]],ncol=1)
ggsave(comb2_chest,file="newplotffolder_sonobodycomplains/correlation/demophoonvschestertonni.corr.DEunique1.zoom.png",device="png",dpi=300,bg="white")

####################### demophoon etylus
position.demety.1=demophoonxety.HW.unique.bed$start+((demophoonxety.HW.unique.bed$end-demophoonxety.HW.unique.bed$end)/2)
position.demety.2=demophoonxety.FW.unique.bed$start+((demophoonxety.FW.unique.bed$end-demophoonxety.FW.unique.bed$end)/2)
position.demety.3=etylusxdem.HW.unique.bed$start+((etylusxdem.HW.unique.bed$end-etylusxdem.HW.unique.bed$end)/2)
position.demety.4=etylusxdem.FW.unique.bed$start+((etylusxdem.FW.unique.bed$end-etylusxdem.FW.unique.bed$end)/2)

position_Deunique_demoetya=c(position.demety.1,position.demety.2,position.demety.3,position.demety.4)


demophoonxetylus.complete=data.frame(str_split_fixed(c(rownames(demophoonFWall),rownames(demophoonHWall),rownames(etylusFWall),rownames(etylusHWall)),"_",3))
colnames(demophoonxetylus.complete)=c("chrom","start","end")
demophoonxetylus.complete$start=as.integer(as.character(demophoonxetylus.complete$start))
demophoonxetylus.complete$end=as.integer(as.character(demophoonxetylus.complete$end))

demety_total_pos=demophoonxetylus.complete$start+((demophoonxetylus.complete$end-demophoonxetylus.complete$start)/2)
########################################################################### test the functions

freq.df_demoetyFst_1k=freq_fun(demxety.1k.fst,demety_total_pos,"Fst")

rank2_ety=rank_plot(demxety.1k.fst,position_Deunique_demoetya,freq.df_demoetyFst_1k,1000,"Fst")

fst.dist.demoetya=Fst.dist.loop(demxety.1k.fst)
split2_ety=Fst.stat(fst.dist.demoetya,position_Deunique_demoetya,demety_total_pos,1000,"no")
comb2_ety=plot_grid(rank2_ety,split2_ety[[1]],ncol=1)
ggsave(comb2_ety,file="newplotffolder_sonobodycomplains/correlation/demophoonvsetylus.corr.DEunique1.zoom.png",device="png",dpi=300,bg="white")



###################summary
primo=split2[[1]]+geom_point(data=split2[[2]],aes(x=sigla,y=final_count/total),col="#cc79a7",shape=ifelse(split2[[2]]$statistic*10<0.05,8,16),size=2)#+ggtitle("dem vs hyd")
secondo=split2_fav[[1]]+geom_point(data=split2_fav[[2]],aes(x=sigla,y=final_count/total),col="#661100",shape=ifelse(split2_fav[[2]]$statistic*10<0.05,8,16),size=2)#+ggtitle("dem vs fav")
terzo=split2_ety[[1]]+geom_point(data=split2_ety[[2]],aes(x=sigla,y=final_count/total),col="#d55e00",shape=ifelse(split2_ety[[2]]$statistic*10<0.05,8,16),size=2)#+ggtitle("dem vs ety")
quarto=split2_chest[[1]]+geom_point(data=split2_chest[[2]],aes(x=sigla,y=final_count/total),col="#E69F00",shape=ifelse(split2_chest[[2]]$statistic*10<0.05,8,16),size=2)#+ggtitle("dem vs ety")



null=ggplot()
full_corr_adj=plot_grid(null,primo,secondo,terzo,quarto,ncol=1,align="v",rel_heights = c(1,rep(0.8,4)))
ggsave(full_corr_adj,file="newplotffolder_sonobodycomplains/correlation/full.corr.adj.pdf",device="pdf",dpi=300,bg="white",width=5, height=5)
library(Cairo)
ggsave(full_corr_adj,file="newplotffolder_sonobodycomplains/correlation/full_corr_adj.svg",device=svg,dpi=300,bg="white",width = 5, height = 5)



full_corr=plot_grid(split2[[1]]+ggtitle("dem vs hyd"),split2_fav[[1]]+ggtitle("dem vs fav"),split2_ety[[1]]+ggtitle("dem vs ety"),split2_chest[[1]]+ggtitle("dem vs chest"),ncol=1)
ggsave(full_corr,file="newplotffolder_sonobodycomplains/correlation/full.corr.pdf",device="pdf",dpi=300,bg="white")



ranksummary=data.frame(pop=c("hydara","favorinus","etylus","chestertonii"),distance=c(0.059805,0.110330,0.125850,0.016899),pvalue=c(0.076410,0.000157,0.000001,0.459600))
ranksummary$pop=factor(ranksummary$pop, levels = ranksummary$pop)
ranksummary$fst=c(0.06716223,0.1006078,0.1160202,0.2708378)
###Interpopulation nucleotide diversity  mean(dat$Dxy-(dat$pi_1+dat$pi_2)/2,na.rm=T)
ranksummary$Vxy=c(0.001449584,0.00474223,0.005429103,0.01330291) 
plot_summ=ggplot(ranksummary,aes(x=fst,y=distance,col=pop))+geom_point(shape=ifelse(ranksummary$pvalue<0.05,18,20),size=4)+
    scale_color_manual(values = c("#cc79a7", "red", "#d55e00","goldenrod1")) +theme_classic()+ theme(legend.position = "none")
ggplot(ranksummary)+geom_bar(aes(x=pop,y=distance,fill=pop),col=ifelse(ranksummary$pvalue<0.05,"black",NA),size=1.5,stat="identity", position="dodge",width = 0.5)+scale_fill_manual(values = c("#cc79a7", "red", "#d55e00","goldenrod1")) +theme_classic()                                                                                                             
########## demo vs notab
notabilis=data.frame(pop ="notabilis",distance=0.10919,pvalue=0.0004998,fst=0.1060611)

ranksummary_comp=rbind(ranksummary,notabilis)
#plot_sum_comp=
ggplot(ranksummary_comp,aes(x=fst,y=distance,col=pop))+geom_point(shape=ifelse(ranksummary_comp$pvalue<0.05,18,20),size=4)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("#cc79a7", "red", "#d55e00","goldenrod1","#009e73")) +theme_classic()+ theme(legend.position = "none")


ggsave(plot_summ,file="newplotffolder_sonobodycomplains/correlation/summary.corr.pdf",device="pdf",dpi=300,bg="white",width = 5, height = 5)
ggsave(plot_sum_comp,file="newplotffolder_sonobodycomplains/correlation/summary.not.corr.pdf",device="pdf",dpi=300,bg="white")

########################
######################sweeps
freq.df_fav=freq_fun(sweep.favo,demfav_total_pos,"sweep")

rank_fav=rank_plot(sweep.favo,position_Deunique_demofava,freq.df_fav,1000,"sweep")
testo_fav=Fst.dist.loop_swep(sweep.favo)
split_fav=Fst.stat(testo_fav,position_Deunique_demofava,demfav_total_pos,1000,"no")


freq.df_chest=freq_fun(sweep.chest,demchest_total_pos,"sweep")

rank_chest=rank_plot(sweep.chest,position_Deunique_demochesta,freq.df_chest,1000,"sweep")
testo_chest=Fst.dist.loop_swep(sweep.chest)
split_chest=Fst.stat(testo_chest,position_Deunique_demochesta,demchest_total_pos,1000,"no")


freq.df_ety2=freq_fun(sweep.ety,demety_total_pos,"sweep")

rank_ety2=rank_plot(sweep.ety,position_Deunique_demoetya,freq.df_ety2,1000,"sweep")
testo_ety2=Fst.dist.loop_swep(sweep.ety)
split_ety2=Fst.stat(testo_ety2,position_Deunique_demoetya,demety_total_pos,1000,"no")



ggsave(plot_grid(rank_ety,split_ety[[1]],ncol=1),file="newplotffolder_sonobodycomplains/correlation/etylysvsdem.sweep..corr.DEunique1.zoom.png",device="png",dpi=300,bg="white")
ggsave(plot_grid(rank_fav,split_fav[[1]],ncol=1),file="newplotffolder_sonobodycomplains/correlation/favorinusvsdem.sweep..corr.DEunique1.zoom.png",device="png",dpi=300,bg="white")
ggsave(plot_grid(rank_chest,split_chest[[1]],ncol=1),file="newplotffolder_sonobodycomplains/correlation/chestertoniivsdem.sweep..corr.DEunique1.zoom.png",device="png",dpi=300,bg="white")





input=split1[[2]]
forexcel=function(input){
pvalue=format(p.adjust(input$statistic,method = "holm"), decimal.mark = ',')
print(data.frame(cbind(input$final_count,pvalue)))
}

###################################   fst nota

demxnota.1k.fst=read.table("fst/demoxnota.1k.FST.2022.txt",h=F,sep=" ",col.names = c("scaffold","start","end","mid","sites","pi_1","pi_2","Dxy","Fst"))
demxnota.1k.fst$Fst[demxnota.1k.fst$Fst<0]<-0
demxnota.1k.fst$position=demxnota.1k.fst$start
demxnota.1k.fst=demxnota.1k.fst[demxnota.1k.fst$start< 688640245,]


notabilisxdem.HW.unique.bed=read.table("demohpoon vs notabilis/DOWN/notabilisxdem.HW.unique.bed",col.names=c("chrom","start","end"))
notabilisxdem.FW.unique.bed=read.table("demohpoon vs notabilis/DOWN/notabilisxdem.FW.unique.bed",col.names=c("chrom","start","end"))

demophoonxnota.HW.unique.bed=read.table("demohpoon vs notabilis/DOWN/demophoonxnotabilis.HW.unique2.bed",col.names=c("chrom","start","end"))
demophoonxnota.FW.unique.bed=read.table("demohpoon vs notabilis/DOWN/demophoonxnotabilis.FW.unique.bed",col.names=c("chrom","start","end"))

demophoonFW.nota.all=read.table("demohpoon vs notabilis/DOWN/demophoon.strict.DOWN.FW.all.bed",col.names=c("chrom","start","end"))
demophoonHW.nota.all=read.table("demohpoon vs notabilis/DOWN/demophoon.strict.DOWN.HW.all.bed",col.names=c("chrom","start","end"))


position1=demophoonxnota.HW.unique.bed$start+((demophoonxnota.HW.unique.bed$end-demophoonxnota.HW.unique.bed$end)/2)
position2=demophoonxnota.FW.unique.bed$start+((demophoonxnota.FW.unique.bed$end-demophoonxnota.FW.unique.bed$end)/2)
position3=notabilisxdem.HW.unique.bed$start+((notabilisxdem.HW.unique.bed$end-notabilisxdem.HW.unique.bed$end)/2)
position4=notabilisxdem.FW.unique.bed$start+((notabilisxdem.FW.unique.bed$end-notabilisxdem.FW.unique.bed$end)/2)
position_Deunique_demonota=c(position1,position2,position3,position4)
position_Deunique_demonota=c(position1,position3,position4)
position_Deunique_demonota=c(position3,position4)

position1=head(position1,length(position1)/2)
position2=head(position2,length(position2)/2)

demophoonxnotabilis.complete=data.frame(str_split_fixed(c(rownames(notabilisFWall),rownames(notabilisHWall)),"_",3))

colnames(demophoonxnotabilis.complete)=c("chrom","start","end")
demophoonxnotabilis.complete=rbind(demophoonxnotabilis.complete,rbind(demophoonxnotabilis.complete,(demophoonFW.nota.all)),demophoonHW.nota.all)
demophoonxnotabilis.complete$start=as.integer(as.character(demophoonxnotabilis.complete$start))
demophoonxnotabilis.complete$end=as.integer(as.character(demophoonxnotabilis.complete$end))

demnota_total_pos=demophoonxnotabilis.complete$start+((demophoonxnotabilis.complete$end-demophoonxnotabilis.complete$start)/2)


freq.df_demonotaFst_1k=freq_fun(demxnota.1k.fst,demnota_total_pos,"Fst")

rank2_nota=rank_plot(demxnota.1k.fst,position_Deunique_demonota,freq.df_demonotaFst_1k,1000,"Fst")
#print(ks.test(rank2_nota[[2]]$Fst.x,rank2_nota[[2]]$Fst.y,alternative="less",simulate.p.value =TRUE))

fst.dist.demonotaa=Fst.dist.loop(demxnota.1k.fst)
split2_nota=Fst.stat(fst.dist.demonotaa,position_Deunique_demonota,demnota_total_pos,1000,"no")
comb2_nota=plot_grid(rank2_nota[[1]],split2_nota[[1]],ncol=1)


###########################################################
##############################################################


