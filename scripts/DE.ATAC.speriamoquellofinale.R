library(RUVSeq)
library(tidyverse)
library(DESeq2)

library(cowplot)
library(ggforce)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

library(pheatmap)



library(RColorBrewer)
library(Rcpp)
library(reshape2)
library(ggrastr)
library(GenomicRanges)
library(rtracklayer)

library(ggpubr)


###################################################################################################   FUNCTIONS



########## RUV parameters
RUVparam<-function(subset,subsetheader,plot,fact,K,con){
  
  subset=select(subset,c(subsetheader$sample))
    filter <- apply(subset, 1, function(x) length(x[x>5])>5)
  filtered <- subset[filter,]
  
  #filtered=subset
  
  set <- newSeqExpressionSet(as.matrix(filtered),
                             phenoData = data.frame(subsetheader, row.names=colnames(filtered)))
  set <- betweenLaneNormalization(set, which="upper")
  
  #plotPCA(as.matrix(filtered), col=colors[as.factor(subsetheader$species)], cex=1.2,main="filtered")
  #plotPCA(set, col=colors[as.factor(subsetheader$species)], cex=1.2,main="set")
  
  design <- model.matrix(~subsetheader$species, data=pData(set))
  
  y <- DGEList(counts=counts(set), group=subsetheader$species)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  res <- residuals(fit, type="deviance")
  genes <- rownames(filtered)
  ###############################################################################################################
  
  differences <- makeGroups(subsetheader$species)
  ############################################################################################################
  
  
  
  colors <- brewer.pal(12, "Set3")
  ############
  if (plot == "yes"){
    par(mfrow = c(2,2))
    plotPCA(set, col=colors[as.factor(subsetheader$species)], cex=1.2,main="raw")
    
    for(k in 1:3) {
      set.x=RUVr(set, genes, k=k, res)
      #set.x= RUVs(set, genes, k=k, differences)
      plotPCA(set.x, col=colors[as.factor(subsetheader$species)], cex=1.2,main=k)
    }
  }
  if (fact == "no"){
  n <- readline(prompt="Enter desired design (in W_): ")  
  }else{
    n=K
  }
  set.x=RUVr(set, genes, k=4, res)
  reads=data.frame(counts(set.x))
  print(pData(set.x))
  print(head(reads))
  
  ddsXxYlFW <- DESeqDataSetFromMatrix(countData = reads,
                                      colData = pData(set.x),
                                      design = formula(paste("~",n,"+species",sep="")))
  
  
  atacDDSXxYlFW <- DESeq(ddsXxYlFW,test="Wald")
  
  keep <- rowSums(counts(atacDDSXxYlFW) >= 20) >= 3
  atacDDSXxYlFW <- atacDDSXxYlFW[keep,]
  #test=assays(atacDDSXxYlFW)[["cooks"]]
  #test=counts(atacDDSXxYlFW,normalized=T)
  test=sizeFactors(atacDDSXxYlFW)
  #plotPCA(as.matrix(reads), col=colors[as.factor(header.notaxety.FW$species)], cex=1.2,main="reads")
  #plotPCA(as.matrix(test), col=colors[as.factor(header.notaxety.FW$species)], cex=1.2,main="test")
  
  #print(reads)
  #atacDDSXxYlFW <- DESeq(ddsXxYlFW,test="Wald")
  #con=c("species","nota","ety")
  #con=c("species","dem","ety")
  #con=c("species","dem","hyd")
  res_FW <- results(atacDDSXxYlFW,contrast=con,lfcThreshold=0.5,alpha=0.05)
  
  FW.vsall.filt=na.omit(data.frame(res_FW))
  print(summary(res_FW))
  #ten=quantile(na.omit(FW.vsall.filt[FW.vsall.filt$padj<0.05,]$baseMean), prob = 1 - 75/100)
  #FW.vsall.filt=as.data.frame(FW.vsall.filt[FW.vsall.filt$baseMean>ten,])
  FW.vsall.filt$pan=ifelse(FW.vsall.filt$log2FoldChange > 0,str_split_fixed(rownames(FW.vsall.filt),"-",2)[,1],str_split_fixed(rownames(FW.vsall.filt),"-",2)[,2])
  FW.vsall.filt$pos=as.numeric(str_split_fixed(FW.vsall.filt$pan,"_",3)[,2])
  FW.vsall.filt$pos2=as.numeric(str_split_fixed(FW.vsall.filt$pan,"_",3)[,3])
  FW.vsall.filt=FW.vsall.filt[FW.vsall.filt$pos2<688640245,]
   # print(nrow(FW.vsall.filt[(FW.vsall.filt$padj<0.05) & (FW.vsall.filt$log2FoldChange>0),]))
  #print(nrow(FW.vsall.filt[(FW.vsall.filt$padj<0.01) & (FW.vsall.filt$log2FoldChange<0),]))
  
  ggplot(FW.vsall.filt[(FW.vsall.filt$padj<0.05),])+geom_point(aes(x=pos,y=log2FoldChange))
  result=list(FW.vsall.filt,test)
  return(result)
  
}

normalie<-function(evenbefore1,evenbefore2,factor){
#ratio=before/after
#factor=ratio[1,]
#print(ratio)
before.norm1=as.data.frame(mapply('/', evenbefore1, factor[1:3]))
before.norm2=as.data.frame(mapply('/', evenbefore2, factor[4:6]))
rownames(before.norm1)=rownames(evenbefore1)
rownames(before.norm2)=rownames(evenbefore2)
lista=list(before.norm1,before.norm2)
return(lista)
}


############################################################################################## load variables


##########chromsomes
chrom.SE=read.table("chrom.SE")

######headers
header.nota=read.table("../allDE//notabilis.header.txt",col.names =c("sample","stage","tissue"),sep="\t")
header.nota$species="nota"

header.ety=read.table("../allDE/etylus.header.txt",col.names =c("sample","stage","tissue"),sep="\t")
header.ety$species="ety"

header.dem=read.table("../allDE/demophoon.header.txt",col.names =c("sample","stage","tissue"),sep="\t")
header.dem$species="dem"

header.hyd=read.table("../allDE/hydara.header.txt",col.names =c("sample","stage","tissue"),sep="\t")
header.hyd$species="hyd"

header.fav=read.table("../allDE/favorinus.header.txt",col.names =c("sample","stage","tissue"),sep="\t")
header.fav$species="fav"

header.ches=read.table("../allDE/chest.header.txt",col.names =c("sample","stage","tissue"),sep="\t")
header.ches$species="chest"



header.total=rbind(header.nota,header.ety) 
header.total=rbind(header.total,header.dem) 
header.total=rbind(header.total,header.hyd)
header.total=rbind(header.total,header.fav) 
header.total=rbind(header.total,header.ches) 

#######################################################################################populations

########################################### notabilis vs etylus


####################### FW


########3 out of 3
notabilisFWall=read.table("notabilis vs etylus/strict/notabilis.strict.FW.all.count.txt",h=F,col.names=c("chrom","start","end","E_Not3_FW","Not1_FW","Not2_FW"))
rownames(notabilisFWall)=paste(notabilisFWall$chrom,notabilisFWall$start,notabilisFWall$end,sep="_")
notabilisFWall=notabilisFWall[c(4,5,6)] ################check the order!!!!


etylusFWall=read.table("notabilis vs etylus/strict/etylus.strict.FW.all.count.txt",h=F,col.names=c("chrom","start","end","E_ety2_FW","Et4_FW","Et5_FW"))
rownames(etylusFWall)=paste(etylusFWall$chrom,etylusFWall$start,etylusFWall$end,sep="_")
etylusFWall=etylusFWall[c(4,5,6)]


notavety3of3=read.table("notabilis vs etylus/strict/notabilis.strictxetylus.strict.FW.3of3.final.txt",row.names=NULL)
#notavety3of3$row.names=NULL
rownames(notavety3of3)=paste(paste(notavety3of3$V1,notavety3of3$V2,notavety3of3$V3,sep="_"),paste(notavety3of3$V7,notavety3of3$V8,notavety3of3$V9,sep="_"),sep="-")
notavety3of3=notavety3of3[c(4,5,6,10,11,12)]

names=c(colnames(notabilisFWall),colnames(etylusFWall))
colnames(notavety3of3)=names


######## uniqu
notabilis.FW.unique.bed=read.table("notabilis vs etylus/strict/notabilis.strict.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(notabilis.FW.unique.bed)=paste("pan",notabilis.FW.unique.bed$start,notabilis.FW.unique.bed$end,sep="_")

etylus.FW.unique.bed=read.table("notabilis vs etylus/strict/etylus.strict.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(etylus.FW.unique.bed)=paste("pan",etylus.FW.unique.bed$start,etylus.FW.unique.bed$end,sep="_")


####### 2of 3
notabilisFW.20f3=read.table("notabilis vs etylus/strict/notabilis.strict.FW.2of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(notabilisFW.20f3)=paste("pan",notabilisFW.20f3$start,notabilisFW.20f3$end,sep="_")

etylusFW.20f3=read.table("notabilis vs etylus/strict/etylus.strict.FW.2of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(etylusFW.20f3)=paste("pan",etylusFW.20f3$start,etylusFW.20f3$end,sep="_")

###### 3 of 3

notabilisFW.3of3=read.table("notabilis vs etylus/strict/notabilis.strict.FW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(notabilisFW.3of3)=paste("pan",notabilisFW.3of3$start,notabilisFW.3of3$end,sep="_")

etylusFW.3of3=read.table("notabilis vs etylus/strict/etylus.strict.FW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(etylusFW.3of3)=paste("pan",etylusFW.3of3$start,etylusFW.3of3$end,sep="_")


####HW


notabilisHWall=read.table("notabilis vs etylus/strict/notabilis.strict.HW.all.count.txt",h=F,col.names=c("chrom","start","end","E_Not3_HW","Not1_HW","Not2_HW"))
rownames(notabilisHWall)=paste(notabilisHWall$chrom,notabilisHWall$start,notabilisHWall$end,sep="_")
notabilisHWall=notabilisHWall[c(4,5,6)] ################check the order!!!!


etylusHWall=read.table("notabilis vs etylus/strict/etylus.strict.HW.all.count.txt",h=F,col.names=c("chrom","start","end","E_ety2_HW","Et4_HW","Et5_HW"))
rownames(etylusHWall)=paste(etylusHWall$chrom,etylusHWall$start,etylusHWall$end,sep="_")
etylusHWall=etylusHWall[c(4,5,6)]


notavety3of3.HW=read.table("notabilis vs etylus/strict/notabilis.strictxetylus.strict.HW.3of3.final.txt",row.names=NULL)
#notavety3of3.HW$row.names=NULL
rownames(notavety3of3.HW)=paste(paste(notavety3of3.HW$V1,notavety3of3.HW$V2,notavety3of3.HW$V3,sep="_"),paste(notavety3of3.HW$V7,notavety3of3.HW$V8,notavety3of3.HW$V9,sep="_"),sep="-")
notavety3of3.HW=notavety3of3.HW[c(4,5,6,10,11,12)]

names=c(colnames(notabilisHWall),colnames(etylusHWall))
colnames(notavety3of3.HW)=names

notabilis.HW.unique.bed=read.table("notabilis vs etylus/strict/notabilis.strict.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(notabilis.HW.unique.bed)=paste("pan",notabilis.HW.unique.bed$start,notabilis.HW.unique.bed$end,sep="_")

etylus.HW.unique.bed=read.table("notabilis vs etylus/strict/etylus.strict.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(etylus.HW.unique.bed)=paste("pan",etylus.HW.unique.bed$start,etylus.HW.unique.bed$end,sep="_")



####### 20f3
notabilisHW.20f3=read.table("notabilis vs etylus/strict/notabilis.strict.HW.2of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(notabilisHW.20f3)=paste("pan",notabilisHW.20f3$start,notabilisHW.20f3$end,sep="_")

etylusHW.20f3=read.table("notabilis vs etylus/strict/etylus.strict.HW.2of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(etylusHW.20f3)=paste("pan",etylusHW.20f3$start,etylusHW.20f3$end,sep="_")

###### 3 of 3

notabilisHW.3of3=read.table("notabilis vs etylus/strict/notabilis.strict.HW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(notabilisHW.3of3)=paste("pan",notabilisHW.3of3$start,notabilisHW.3of3$end,sep="_")

etylusHW.3of3=read.table("notabilis vs etylus/strict/etylus.strict.HW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(etylusHW.3of3)=paste("pan",etylusHW.3of3$start,etylusHW.3of3$end,sep="_")





##################################################### demophoon hydara
###### FW

###############
demophoonFWall=read.table("demophoon x hydara//demophoon.strict.FW.all.count.txt",h=F,col.names=c("chrom","start","end","LI7_demophoon_FW","FW-pboy","LB_41"))
rownames(demophoonFWall)=paste(demophoonFWall$chrom,demophoonFWall$start,demophoonFWall$end,sep="_")
demophoonFWall=demophoonFWall[c(5,6,4)] ################check the order!!!!


hydaraFWall=read.table("demophoon x hydara/hydara.strict.FW.all.count.txt",h=F,col.names=c("chrom","start","end","LB_20","LB_27","LB_29"))
rownames(hydaraFWall)=paste(hydaraFWall$chrom,hydaraFWall$start,hydaraFWall$end,sep="_")
hydaraFWall=hydaraFWall[c(4,5,6)]


demoxhyda3of3=read.table("demophoon x hydara/demoxhydara.strict.FW.3of3.final.txt",row.names=NULL)
#demxhyda3of3$row.names=NULL
rownames(demoxhyda3of3)=paste(paste(demoxhyda3of3$V1,demoxhyda3of3$V2,demoxhyda3of3$V3,sep="_"),paste(demoxhyda3of3$V7,demoxhyda3of3$V8,demoxhyda3of3$V9,sep="_"),sep="-")
demoxhyda3of3=demoxhyda3of3[c(5,6,4,10,11,12)]

names=c(colnames(demophoonFWall),colnames(hydaraFWall))
colnames(demoxhyda3of3)=names


demophoon.FW.unique.bed=read.table("demophoon x hydara/demophoon.strict.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(demophoon.FW.unique.bed)=paste("pan",demophoon.FW.unique.bed$start,demophoon.FW.unique.bed$end,sep="_")

hydara.FW.unique.bed=read.table("demophoon x hydara/hydara.strict.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(hydara.FW.unique.bed)=paste("pan",hydara.FW.unique.bed$start,hydara.FW.unique.bed$end,sep="_")

####### 20f3
demophoonFW.20f3=read.table("demophoon x hydara/demophoon.strict.FW.2of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(demophoonFW.20f3)=paste("pan",demophoonFW.20f3$start,demophoonFW.20f3$end,sep="_")

hydaraFW.20f3=read.table("demophoon x hydara/hydara.strict.FW.2of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(hydaraFW.20f3)=paste("pan",hydaraFW.20f3$start,hydaraFW.20f3$end,sep="_")

###### 3 of 3

demophoonFW.3of3=read.table("demophoon x hydara/demophoon.strict.FW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(demophoonFW.3of3)=paste("pan",demophoonFW.3of3$start,demophoonFW.3of3$end,sep="_")

hydaraFW.3of3=read.table("demophoon x hydara/hydara.strict.FW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(hydaraFW.3of3)=paste("pan",hydaraFW.3of3$start,hydaraFW.3of3$end,sep="_")

#HW
demophoonHWall=read.table("demophoon x hydara//demophoon.strict.HW.all.count.txt",h=F,col.names=c("chrom","start","end","LI7_demophoon_HW","E3_HW","LB_42"))
rownames(demophoonHWall)=paste(demophoonHWall$chrom,demophoonHWall$start,demophoonHWall$end,sep="_")
demophoonHWall=demophoonHWall[c(5,6,4)] ################check the order!!!!


hydaraHWall=read.table("demophoon x hydara/hydara.strict.HW.all.count.txt",h=F,col.names=c("chrom","start","end","LB_21","LB_28","LB_30"))
rownames(hydaraHWall)=paste(hydaraHWall$chrom,hydaraHWall$start,hydaraHWall$end,sep="_")
hydaraHWall=hydaraHWall[c(4,5,6)]


demoxhyda3of3.HW=read.table("demophoon x hydara/demoxhydara.strict.HW.3of3.final.txt",row.names=NULL)
#demxhyda3of3$row.names=NULL
rownames(demoxhyda3of3.HW)=paste(paste(demoxhyda3of3.HW$V1,demoxhyda3of3.HW$V2,demoxhyda3of3.HW$V3,sep="_"),paste(demoxhyda3of3.HW$V7,demoxhyda3of3.HW$V8,demoxhyda3of3.HW$V9,sep="_"),sep="-")
demoxhyda3of3.HW=demoxhyda3of3.HW[c(5,6,4,10,11,12)]

names=c(colnames(demophoonHWall),colnames(hydaraHWall))
colnames(demoxhyda3of3.HW)=names

demophoon.HW.unique.bed=read.table("demophoon x hydara/demophoon.strict.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(demophoon.HW.unique.bed)=paste("pan",demophoon.HW.unique.bed$start,demophoon.HW.unique.bed$end,sep="_")

hydara.HW.unique.bed=read.table("demophoon x hydara/hydara.strict.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(hydara.HW.unique.bed)=paste("pan",hydara.HW.unique.bed$start,hydara.HW.unique.bed$end,sep="_")


####### 20f3
demophoonHW.20f3=read.table("demophoon x hydara/demophoon.strict.HW.2of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(demophoonHW.20f3)=paste("pan",demophoonHW.20f3$start,demophoonHW.20f3$end,sep="_")

hydaraHW.20f3=read.table("demophoon x hydara/hydara.strict.HW.2of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(hydaraHW.20f3)=paste("pan",hydaraHW.20f3$start,hydaraHW.20f3$end,sep="_")

###### 3 of 3

demophoonHW.3of3=read.table("demophoon x hydara/demophoon.strict.HW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(demophoonHW.3of3)=paste("pan",demophoonHW.3of3$start,demophoonHW.3of3$end,sep="_")

hydaraHW.3of3=read.table("demophoon x hydara/hydara.strict.HW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(hydaraHW.3of3)=paste("pan",hydaraHW.3of3$start,hydaraHW.3of3$end,sep="_")


#########################################################################################################################################################
######################################## demophoon facvorinus

favorinusFWall=read.table("demophoon vs favorinus/favorinus.FW.all.count.txt",h=F,col.names=c("chrom","start","end","LB_16","LB_18","LB_6"))
rownames(favorinusFWall)=paste(favorinusFWall$chrom,favorinusFWall$start,favorinusFWall$end,sep="_")
favorinusFWall=favorinusFWall[c(4,5,6)] ################check the order!!!!


demvsfavo.FW=read.table("demophoon vs favorinus/demophoonxfavorinus.FW.3of3.final.txt",row.names=NULL)
#demvsfavo.FW$row.names=NULL
rownames(demvsfavo.FW)=paste(paste(demvsfavo.FW$V1,demvsfavo.FW$V2,demvsfavo.FW$V3,sep="_"),paste(demvsfavo.FW$V7,demvsfavo.FW$V8,demvsfavo.FW$V9,sep="_"),sep="-")
demvsfavo.FW=demvsfavo.FW[c(4,5,6,10,11,12)]

names=c(colnames(demophoonFWall),colnames(favorinusFWall))
colnames(demvsfavo.FW)=names


######## uniqu
favorinus.FW.unique.bed=read.table("demophoon vs favorinus/favorinus.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(favorinus.FW.unique.bed)=paste("pan",favorinus.FW.unique.bed$start,favorinus.FW.unique.bed$end,sep="_")

demophoonxfav.FW.unique.bed=read.table("demophoon vs favorinus/demophoonxfav.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(demophoonxfav.FW.unique.bed)=paste("pan",demophoonxfav.FW.unique.bed$start,demophoonxfav.FW.unique.bed$end,sep="_")

###### 3 of 3

favorinusFW.3of3=read.table("demophoon vs favorinus/favorinus.FW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(favorinusFW.3of3)=paste("pan",favorinusFW.3of3$start,favorinusFW.3of3$end,sep="_")


###HW
favorinusHWall=read.table("demophoon vs favorinus/favorinus.HW.all.count.txt",h=F,col.names=c("chrom","start","end","LB_17","LB_19","LB_7"))
rownames(favorinusHWall)=paste(favorinusHWall$chrom,favorinusHWall$start,favorinusHWall$end,sep="_")
favorinusHWall=favorinusHWall[c(4,5,6)] ################check the order!!!!


demvsfavo.HW=read.table("demophoon vs favorinus/demophoonxfavorinus.HW.3of3.final.txt",row.names=NULL)
#demvsfavo.HW$row.names=NULL
rownames(demvsfavo.HW)=paste(paste(demvsfavo.HW$V1,demvsfavo.HW$V2,demvsfavo.HW$V3,sep="_"),paste(demvsfavo.HW$V7,demvsfavo.HW$V8,demvsfavo.HW$V9,sep="_"),sep="-")
demvsfavo.HW=demvsfavo.HW[c(4,5,6,10,11,12)]

names=c(colnames(demophoonHWall),colnames(favorinusHWall))
colnames(demvsfavo.HW)=names


######## uniqu
favorinus.HW.unique.bed=read.table("demophoon vs favorinus/favorinus.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(favorinus.HW.unique.bed)=paste("pan",favorinus.HW.unique.bed$start,favorinus.HW.unique.bed$end,sep="_")

demophoonxfav.HW.unique.bed=read.table("demophoon vs favorinus/demophoonxfav.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(demophoonxfav.HW.unique.bed)=paste("pan",demophoonxfav.HW.unique.bed$start,demophoonxfav.HW.unique.bed$end,sep="_")

###### 3 of 3

favorinusHW.3of3=read.table("demophoon vs favorinus/favorinus.HW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(favorinusHW.3of3)=paste("pan",favorinusHW.3of3$start,favorinusHW.3of3$end,sep="_")



###############################################3 demophoon vs chestertonii

chestertoniiFWall=read.table("demophoon vs chestertonii/chestertonii.FW.all.count.txt",h=F,col.names=c("chrom","start","end","Chest2_FW","E_Ches1_FW"))
rownames(chestertoniiFWall)=paste(chestertoniiFWall$chrom,chestertoniiFWall$start,chestertoniiFWall$end,sep="_")
chestertoniiFWall=chestertoniiFWall[c(4,5)] ################check the order!!!!


demvschesto.FW=read.table("demophoon vs chestertonii/demophoonxchestertonii.FW.3of3.final.txt",row.names=NULL)
#demvschesto.FW$row.names=NULL
rownames(demvschesto.FW)=paste(paste(demvschesto.FW$V1,demvschesto.FW$V2,demvschesto.FW$V3,sep="_"),paste(demvschesto.FW$V7,demvschesto.FW$V8,demvschesto.FW$V9,sep="_"),sep="-")
demvschesto.FW=demvschesto.FW[c(4,5,6,10,11,12)]

names=c(colnames(demophoonFWall),colnames(chestertoniiFWall),"NA")
colnames(demvschesto.FW)=names


######## uniqu
chestertonii.FW.unique.bed=read.table("demophoon vs chestertonii/chestertonii.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(chestertonii.FW.unique.bed)=paste("pan",chestertonii.FW.unique.bed$start,chestertonii.FW.unique.bed$end,sep="_")

demophoonxchest.FW.unique.bed=read.table("demophoon vs chestertonii/demophoonxches.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(demophoonxchest.FW.unique.bed)=paste("pan",demophoonxchest.FW.unique.bed$start,demophoonxchest.FW.unique.bed$end,sep="_")

###### 3 of 3

chestertoniiFW.3of3=read.table("demophoon vs chestertonii/chestertonii.FW.2of2.bed",h=F,col.names=c("chrom","start","end"))
rownames(chestertoniiFW.3of3)=paste("pan",chestertoniiFW.3of3$start,chestertoniiFW.3of3$end,sep="_")


###HW
chestertoniiHWall=read.table("demophoon vs chestertonii/chestertonii.HW.all.count.txt",h=F,col.names=c("chrom","start","end","Chest2_HW","E_Ches1_Hw"))
rownames(chestertoniiHWall)=paste(chestertoniiHWall$chrom,chestertoniiHWall$start,chestertoniiHWall$end,sep="_")
chestertoniiHWall=chestertoniiHWall[c(4,5)] ################check the order!!!!


demvschesto.HW=read.table("demophoon vs chestertonii/demophoonxchestertonii.HW.3of3.final.txt",row.names=NULL)
#demvschesto.HW$row.names=NULL
rownames(demvschesto.HW)=paste(paste(demvschesto.HW$V1,demvschesto.HW$V2,demvschesto.HW$V3,sep="_"),paste(demvschesto.HW$V7,demvschesto.HW$V8,demvschesto.HW$V9,sep="_"),sep="-")
demvschesto.HW=demvschesto.HW[c(4,5,6,10,11,12)]

names=c(colnames(demophoonHWall),colnames(chestertoniiHWall))
colnames(demvschesto.HW)=names


######## uniqu
chestertonii.HW.unique.bed=read.table("demophoon vs chestertonii/chestertonii.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(chestertonii.HW.unique.bed)=paste("pan",chestertonii.HW.unique.bed$start,chestertonii.HW.unique.bed$end,sep="_")

demophoonxchest.HW.unique.bed=read.table("demophoon vs chestertonii/demophoonxches.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(demophoonxchest.HW.unique.bed)=paste("pan",demophoonxchest.HW.unique.bed$start,demophoonxchest.HW.unique.bed$end,sep="_")

###### 3 of 3

chestertoniiHW.3of3=read.table("demophoon vs chestertonii/chestertonii.HW.2of2.bed",h=F,col.names=c("chrom","start","end"))
rownames(chestertoniiHW.3of3)=paste("pan",chestertoniiHW.3of3$start,chestertoniiHW.3of3$end,sep="_")

################################### demophoon vs etylus


etylusxdemFWall=read.table("demo x ety/etylus.FW.all.count.txt",h=F,col.names=c("chrom","start","end","E_ety2_FW","Et4_FW","Et5_FW"))
rownames(etylusxdemFWall)=paste(etylusxdemFWall$chrom,etylusxdemFWall$start,etylusxdemFWall$end,sep="_")
etylusxdemFWall=etylusxdemFWall[c(4,5,6)] ################check the order!!!!


demvsety.FW=read.table("demo x ety//demophoonxetylus.FW.3of3.final.txt",row.names=NULL)
#demvsety.FW$row.names=NULL
rownames(demvsety.FW)=paste(paste(demvsety.FW$V1,demvsety.FW$V2,demvsety.FW$V3,sep="_"),paste(demvsety.FW$V7,demvsety.FW$V8,demvsety.FW$V9,sep="_"),sep="-")
demvsety.FW=demvsety.FW[c(4,5,6,10,11,12)]

names=c(colnames(demophoonFWall),colnames(etylusxdemFWall))
colnames(demvsety.FW)=names


######## uniqu
  etylusxdem.FW.unique.bed=read.table("demo x ety/etylusxdem.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(etylusxdem.FW.unique.bed)=paste("pan",etylusxdem.FW.unique.bed$start,etylusxdem.FW.unique.bed$end,sep="_")

demophoonxety.FW.unique.bed=read.table("demo x ety/demophoonxetylus.FW.unique.bed",col.names=c("chrom","start","end"))
rownames(demophoonxety.FW.unique.bed)=paste("pan",demophoonxety.FW.unique.bed$start,demophoonxety.FW.unique.bed$end,sep="_")

###### 3 of 3

etylusxdemFW.3of3=read.table("demo x ety//etylus.FW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(etylusxdemFW.3of3)=paste("pan",etylusxdemFW.3of3$start,etylusxdemFW.3of3$end,sep="_")


###HW

etylusxdemHWall=read.table("demo x ety/etylus.HW.all.count.txt",h=F,col.names=c("chrom","start","end","E_ety2_HW","Et4_HW","Et5_HW"))
rownames(etylusxdemHWall)=paste(etylusxdemHWall$chrom,etylusxdemHWall$start,etylusxdemHWall$end,sep="_")
etylusxdemHWall=etylusxdemHWall[c(4,5,6)] ################check the order!!!!


demvsety.HW=read.table("demo x ety//demophoonxetylus.HW.3of3.final.txt",row.names=NULL)
#demvsety.HW$row.names=NULL
rownames(demvsety.HW)=paste(paste(demvsety.HW$V1,demvsety.HW$V2,demvsety.HW$V3,sep="_"),paste(demvsety.HW$V7,demvsety.HW$V8,demvsety.HW$V9,sep="_"),sep="-")
demvsety.HW=demvsety.HW[c(4,5,6,10,11,12)]

names=c(colnames(demophoonHWall),colnames(etylusxdemHWall))
colnames(demvsety.HW)=names


######## uniqu
etylusxdem.HW.unique.bed=read.table("demo x ety/etylusxdem.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(etylusxdem.HW.unique.bed)=paste("pan",etylusxdem.HW.unique.bed$start,etylusxdem.HW.unique.bed$end,sep="_")

demophoonxety.HW.unique.bed=read.table("demo x ety/demophoonxetylus.HW.unique.bed",col.names=c("chrom","start","end"))
rownames(demophoonxety.HW.unique.bed)=paste("pan",demophoonxety.HW.unique.bed$start,demophoonxety.HW.unique.bed$end,sep="_")

###### 3 of 3

etylusxdemHW.3of3=read.table("demo x ety//etylus.HW.3of3.bed",h=F,col.names=c("chrom","start","end"))
rownames(etylusxdemHW.3of3)=paste("pan",etylusxdemHW.3of3$start,etylusxdemHW.3of3$end,sep="_")


############################################################################# EXECUTION
################################notabilis vs etylus
#################### FW

#header.notaxety.FW=header.total[(header.total$tissue=="forewing") & (header.total$species=="nota" | header.total$species=="ety") ,]
#ind=match(colnames(notaxety),header.notaxety.FW$sample) 
#header.notaxety.FW=header.notaxety.FW[ind,]
#notavsety.FW=RUVparam(notavety,header.notaxety.FW,"no","","") # W_1
#notavsety.norm=normalie(notabilisFW,etylusFW,notavety,notavsety.FW[[2]])


####3 out of 3
header.notaxety.FW=header.total[(header.total$tissue=="forewing") & (header.total$species=="nota" | header.total$species=="ety") ,]

notavsety.FW3of3=RUVparam(notavety3of3,header.notaxety.FW,"no","yes","W_1+W_2",c("species","nota","ety")) 
notavsety.norm.all=normalie(notabilisFWall,etylusFWall,notavsety.FW3of3[[2]])

notabilis.FW.unique=notavsety.norm.all[[1]][rownames(notavsety.norm.all[[1]]) %in% rownames(notabilis.FW.unique.bed),]
notabilis.FW.unique=notabilis.FW.unique[rowMeans(notabilis.FW.unique) > quantile(rowMeans(notabilis.FW.unique),0.25),]

etylus.FW.unique=notavsety.norm.all[[2]][rownames(notavsety.norm.all[[2]]) %in% rownames(etylus.FW.unique.bed),]
etylus.FW.unique=etylus.FW.unique[rowMeans(etylus.FW.unique) > quantile(rowMeans(etylus.FW.unique),0.25),]




#####################   HW
header.notaxety.HW=header.total[(header.total$tissue=="hindwing") & (header.total$species=="nota" | header.total$species=="ety") ,]

notavsety.HW3of3=RUVparam(notavety3of3.HW,header.notaxety.HW,"no","yes","W_1+W_2",c("species","nota","ety")) 
notavsety.norm.all.HW=normalie(notabilisHWall,etylusHWall,notavsety.HW3of3[[2]])

notabilis.HW.unique=notavsety.norm.all.HW[[1]][rownames(notavsety.norm.all.HW[[1]]) %in% rownames(notabilis.HW.unique.bed),]
#notabilis.HW.unique=notabilis.HW.unique[rowMeans(notabilis.HW.unique) > quantile(rowMeans(notabilis.HW.unique),0.25),]

etylus.HW.unique=notavsety.norm.all.HW[[2]][rownames(notavsety.norm.all.HW[[2]]) %in% rownames(etylus.HW.unique.bed),]
#etylus.HW.unique=etylus.HW.unique[rowMeans(etylus.HW.unique) > quantile(rowMeans(etylus.HW.unique),0.25),]

  

 
###################################################### demophoon x hydara
#######################################FW
header.demoxhyda.FW=header.total[(header.total$tissue=="forewing") & (header.total$species=="dem" | header.total$species=="hyd") ,]
header.demoxhyda.FW=header.demoxhyda.FW[header.demoxhyda.FW$sample != "E3_FW",]
levels(header.demoxhyda.FW$sample)[match("FW-pboy",levels(header.demoxhyda.FW$sample))] <- "FW.pboy"




demoxhyda.FW3of3=RUVparam(demoxhyda3of3,header.demoxhyda.FW,"yes","no","W_1+W_2",c("species","dem","hyd")) 
demxhyda.norm.all=normalie(demophoonFWall,hydaraFWall,demoxhyda.FW3of3[[2]])

demophoon.FW.unique=demxhyda.norm.all[[1]][rownames(demxhyda.norm.all[[1]]) %in% rownames(demophoon.FW.unique.bed),]
demophoon.FW.unique=demophoon.FW.unique[rowMeans(demophoon.FW.unique) > quantile(rowMeans(demophoon.FW.unique),0.25),]

hydara.FW.unique=demxhyda.norm.all[[2]][rownames(demxhyda.norm.all[[2]]) %in% rownames(hydara.FW.unique.bed),]
hydara.FW.unique=hydara.FW.unique[rowMeans(hydara.FW.unique) > quantile(rowMeans(hydara.FW.unique),0.25),]







#####################################HW


header.demoxhyda.HW=header.total[(header.total$tissue=="hindwing") & (header.total$species=="dem" | header.total$species=="hyd") ,]
demoxhyda.HW3of3=RUVparam(demoxhyda3of3.HW,header.demoxhyda.HW,"yes","yes","W_1+W_2",c("species","dem","hyd")) 
demxhyda.norm.all.HW=normalie(demophoonHWall,hydaraHWall,demoxhyda.HW3of3[[2]])

demophoon.HW.unique=demxhyda.norm.all.HW[[1]][rownames(demxhyda.norm.all.HW[[1]]) %in% rownames(demophoon.HW.unique.bed),]
demophoon.HW.unique=demophoon.HW.unique[rowMeans(demophoon.HW.unique) > quantile(rowMeans(demophoon.HW.unique),0.25),]

hydara.HW.unique=demxhyda.norm.all.HW[[2]][rownames(demxhyda.norm.all.HW[[2]]) %in% rownames(hydara.HW.unique.bed),]
hydara.HW.unique=hydara.HW.unique[rowMeans(hydara.HW.unique) > quantile(rowMeans(hydara.HW.unique),0.25),]





################################################################### demophhon vs favorinus
header.demoxfavo.FW=header.total[(header.total$tissue=="forewing") & (header.total$species=="dem" | header.total$species=="fav") ,]
header.demoxfavo.FW=header.demoxfavo.FW[header.demoxfavo.FW$sample != "E3_FW",]
levels(header.demoxfavo.FW$sample)[match("FW-pboy",levels(header.demoxfavo.FW$sample))] <- "FW.pboy"

demoxfavo.FW3of3=RUVparam(demvsfavo.FW,header.demoxfavo.FW,"no","yes","W_1+W_2",c("species","dem","fav")) 
demoxfavo.norm.all=normalie(demophoonFWall,favorinusFWall,demoxfavo.FW3of3[[2]])


demophoonxfav.FW.unique=demoxfavo.norm.all[[1]][rownames(demoxfavo.norm.all[[1]]) %in% rownames(demophoonxfav.FW.unique.bed),]
demophoonxfav.FW.unique=demophoonxfav.FW.unique[rowMeans(demophoonxfav.FW.unique) > quantile(rowMeans(demophoonxfav.FW.unique),0.25),]

favorinus.FW.unique=demoxfavo.norm.all[[2]][rownames(demoxfavo.norm.all[[2]]) %in% rownames(favorinus.FW.unique.bed),]
favorinus.FW.unique=favorinus.FW.unique[rowMeans(favorinus.FW.unique) > quantile(rowMeans(favorinus.FW.unique),0.25),]



#### hw
header.demoxfavo.HW=header.total[(header.total$tissue=="hindwing") & (header.total$species=="dem" | header.total$species=="fav") ,]

demoxfavo.HW3of3=RUVparam(demvsfavo.HW,header.demoxfavo.HW,"no","yes","W_1+W_2",c("species","dem","fav")) 
demoxfavo.norm.all.HW=normalie(demophoonHWall,favorinusHWall,demoxfavo.HW3of3[[2]])


demophoonxfav.HW.unique=demoxfavo.norm.all[[1]][rownames(demoxfavo.norm.all[[1]]) %in% rownames(demophoonxfav.HW.unique.bed),]
demophoonxfav.HW.unique=demophoonxfav.HW.unique[rowMeans(demophoonxfav.HW.unique) > quantile(rowMeans(demophoonxfav.HW.unique),0.25),]

favorinus.HW.unique=demoxfavo.norm.all[[2]][rownames(demoxfavo.norm.all[[2]]) %in% rownames(favorinus.HW.unique.bed),]
favorinus.HW.unique=favorinus.HW.unique[rowMeans(favorinus.HW.unique) > quantile(rowMeans(favorinus.HW.unique),0.25),]


######################################## demophoon chesteronii


header.demoxchest.FW=header.total[(header.total$tissue=="forewing") & (header.total$species=="dem" | header.total$species=="chest") ,]
header.demoxchest.FW=header.demoxchest.FW[header.demoxchest.FW$sample != "E3_FW",]
levels(header.demoxchest.FW$sample)[match("FW-pboy",levels(header.demoxchest.FW$sample))] <- "FW.pboy"

demvschesto.FW=demvschesto.FW[-6]
demoxchest.FW3of3=RUVparam_chest(demvschesto.FW,header.demoxchest.FW,"no","yes","W_1",c("species","dem","chest")) 
demoxchest.norm.all=normalie_chest(demophoonFWall,chestertoniiFWall,demoxchest.FW3of3[[2]])


demophoonxchest.FW.unique=demoxchest.norm.all[[1]][rownames(demoxchest.norm.all[[1]]) %in% rownames(demophoonxchest.FW.unique.bed),]
demophoonxchest.FW.unique=demophoonxchest.FW.unique[rowMeans(demophoonxchest.FW.unique) > quantile(rowMeans(demophoonxchest.FW.unique),0.25),]

chesteronii.FW.unique=demoxchest.norm.all[[2]][rownames(demoxchest.norm.all[[2]]) %in% rownames(chestertonii.FW.unique.bed),]
chesteronii.FW.unique=chesteronii.FW.unique[rowMeans(chesteronii.FW.unique) > quantile(rowMeans(chesteronii.FW.unique),0.25),]


############

header.demoxchest.HW=header.total[(header.total$tissue=="hindwing") & (header.total$species=="dem" | header.total$species=="chest") ,]

demvschesto.HW=demvschesto.HW[-6]
demoxchest.HW3of3=RUVparam_chest(demvschesto.HW,header.demoxchest.HW,"no","yes","W_1",c("species","dem","chest")) 
demoxchest.norm.all.HW=normalie_chest(demophoonHWall,chestertoniiHWall,demoxchest.HW3of3[[2]])


demophoonxchest.HW.unique=demoxchest.norm.all[[1]][rownames(demoxchest.norm.all[[1]]) %in% rownames(demophoonxchest.HW.unique.bed),]
demophoonxchest.HW.unique=demophoonxchest.HW.unique[rowMeans(demophoonxchest.HW.unique) > quantile(rowMeans(demophoonxchest.HW.unique),0.25),]

chesteronii.HW.unique=demoxchest.norm.all[[2]][rownames(demoxchest.norm.all[[2]]) %in% rownames(chestertonii.HW.unique.bed),]
chesteronii.HW.unique=chesteronii.HW.unique[rowMeans(chesteronii.HW.unique) > quantile(rowMeans(chesteronii.HW.unique),0.25),]











#############3 demo x ety 
## FW

header.demoxety.FW=header.total[(header.total$tissue=="forewing") & (header.total$species=="dem" | header.total$species=="ety") ,]
header.demoxety.FW=header.demoxety.FW[header.demoxety.FW$sample != "E3_FW",]
levels(header.demoxety.FW$sample)[match("FW-pboy",levels(header.demoxety.FW$sample))] <- "FW.pboy"

demoxety.FW3of3=RUVparam(demvsety.FW,header.demoxety.FW,"yes","no","W_1+W_2",c("species","dem","ety")) 
demoxety.norm.all=normalie(demophoonFWall,etylusFWall,demoxety.FW3of3[[2]])


demophoonxety.FW.unique=demoxety.norm.all[[1]][rownames(demoxety.norm.all[[1]]) %in% rownames(demophoonxety.FW.unique.bed),]
demophoonxety.FW.unique=demophoonxety.FW.unique[rowMeans(demophoonxety.FW.unique) > quantile(rowMeans(demophoonxety.FW.unique),0.25),]

etylusxdem.FW.unique=demoxety.norm.all[[2]][rownames(demoxety.norm.all[[2]]) %in% rownames(etylusxdem.FW.unique.bed),]
etylusxdem.FW.unique=etylusxdem.FW.unique[rowMeans(etylusxdem.FW.unique) > quantile(rowMeans(etylusxdem.FW.unique),0.25),]

####HW

header.demoxety.HW=header.total[(header.total$tissue=="hindwing") & (header.total$species=="dem" | header.total$species=="ety") ,]


demoxety.HW3of3=RUVparam(demvsety.HW,header.demoxety.HW,"yes","no","W_1+W_2",c("species","dem","ety")) 
demoxety.HW.norm.all=normalie(demophoonHWall,etylusHWall,demoxety.HW3of3[[2]])


demophoonxety.HW.unique=demoxety.HW.norm.all[[1]][rownames(demoxety.HW.norm.all[[1]]) %in% rownames(demophoonxety.HW.unique.bed),]
demophoonxety.HW.unique=demophoonxety.HW.unique[rowMeans(demophoonxety.HW.unique) > quantile(rowMeans(demophoonxety.HW.unique),0.25),]

etylusxdem.HW.unique=demoxety.HW.norm.all[[2]][rownames(demoxety.HW.norm.all[[2]]) %in% rownames(etylusxdem.HW.unique.bed),]
etylusxdem.HW.unique=etylusxdem.HW.unique[rowMeans(etylusxdem.HW.unique) > quantile(rowMeans(etylusxdem.HW.unique),0.25),]








################################
#################################

uniquelist=c(
rownames(notabilis.FW.unique[rowMeans(notabilis.FW.unique)>250,]),
rownames(notabilis.HW.unique[rowMeans(notabilis.HW.unique)>250,]),
rownames(etylus.FW.unique[rowMeans(etylus.FW.unique)>250,]),
rownames(etylus.HW.unique[rowMeans(etylus.HW.unique)>250,]),
notavsety.FW3of3[[1]][notavsety.FW3of3[[1]]$padj<0.05,]$pan,
notavsety.HW3of3[[1]][notavsety.HW3of3[[1]]$padj<0.05,]$pan)

uniquelist=data.frame(str_split_fixed(uniquelist,"_",3))
colnames(uniquelist)=c("chrom","start","end")

write.table(uniquelist,"notaxety.totalistforFST.bed",quote=F,row.names = F,col.names = F)


#############################################
RUVparam_chest<-function(subset,subsetheader,plot,fact,K,con){
  
  subset=select(subset,c(subsetheader$sample))
  filter <- apply(subset, 1, function(x) length(x[x>5])>4)
  filtered <- subset[filter,]
  
  #filtered=subset
  
  set <- newSeqExpressionSet(as.matrix(filtered),
                             phenoData = data.frame(subsetheader, row.names=colnames(filtered)))
  set <- betweenLaneNormalization(set, which="upper")
  
  #plotPCA(as.matrix(filtered), col=colors[as.factor(subsetheader$species)], cex=1.2,main="filtered")
  #plotPCA(set, col=colors[as.factor(subsetheader$species)], cex=1.2,main="set")
  
  design <- model.matrix(~subsetheader$species, data=pData(set))
  
  y <- DGEList(counts=counts(set), group=subsetheader$species)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  res <- residuals(fit, type="deviance")
  genes <- rownames(filtered)
  ###############################################################################################################
  
  differences <- makeGroups(subsetheader$species)
  ############################################################################################################
  
  
  
  colors <- brewer.pal(12, "Set3")
  ############
  if (plot == "yes"){
    par(mfrow = c(2,2))
    plotPCA(set, col=colors[as.factor(subsetheader$species)], cex=1.2,main="raw")
    
    for(k in 1:3) {
      set.x=RUVr(set, genes, k=k, res)
      #set.x= RUVs(set, genes, k=k, differences)
      plotPCA(set.x, col=colors[as.factor(subsetheader$species)], cex=1.2,main=k)
    }
  }
  if (fact == "no"){
    n <- readline(prompt="Enter desired design (in W_): ")  
  }else{
    n=K
  }
  set.x=RUVr(set, genes, k=4, res)
  reads=data.frame(counts(set.x))
  print(pData(set.x))
  print(head(reads))
  
  ddsXxYlFW <- DESeqDataSetFromMatrix(countData = reads,
                                      colData = pData(set.x),
                                      design = formula(paste("~",n,"+species",sep="")))
  
  
  atacDDSXxYlFW <- DESeq(ddsXxYlFW,test="Wald")
  
  keep <- rowSums(counts(atacDDSXxYlFW) >= 20) >= 3
  atacDDSXxYlFW <- atacDDSXxYlFW[keep,]
  #test=assays(atacDDSXxYlFW)[["cooks"]]
  #test=counts(atacDDSXxYlFW,normalized=T)
  test=sizeFactors(atacDDSXxYlFW)
  #plotPCA(as.matrix(reads), col=colors[as.factor(header.notaxety.FW$species)], cex=1.2,main="reads")
  #plotPCA(as.matrix(test), col=colors[as.factor(header.notaxety.FW$species)], cex=1.2,main="test")
  
  #print(reads)
  #atacDDSXxYlFW <- DESeq(ddsXxYlFW,test="Wald")
  #con=c("species","nota","ety")
  #con=c("species","dem","ety")
  #con=c("species","dem","hyd")
  res_FW <- results(atacDDSXxYlFW,contrast=con,lfcThreshold=0.5,alpha=0.05)
  
  FW.vsall.filt=na.omit(data.frame(res_FW))
  print(summary(res_FW))
  #ten=quantile(na.omit(FW.vsall.filt[FW.vsall.filt$padj<0.05,]$baseMean), prob = 1 - 75/100)
  #FW.vsall.filt=as.data.frame(FW.vsall.filt[FW.vsall.filt$baseMean>ten,])
  FW.vsall.filt$pan=ifelse(FW.vsall.filt$log2FoldChange > 0,str_split_fixed(rownames(FW.vsall.filt),"-",2)[,1],str_split_fixed(rownames(FW.vsall.filt),"-",2)[,2])
  FW.vsall.filt$pos=as.numeric(str_split_fixed(FW.vsall.filt$pan,"_",3)[,2])
  FW.vsall.filt$pos2=as.numeric(str_split_fixed(FW.vsall.filt$pan,"_",3)[,3])
  FW.vsall.filt=FW.vsall.filt[FW.vsall.filt$pos2<688640245,]
  # print(nrow(FW.vsall.filt[(FW.vsall.filt$padj<0.05) & (FW.vsall.filt$log2FoldChange>0),]))
  #print(nrow(FW.vsall.filt[(FW.vsall.filt$padj<0.01) & (FW.vsall.filt$log2FoldChange<0),]))
  
  ggplot(FW.vsall.filt[(FW.vsall.filt$padj<0.05),])+geom_point(aes(x=pos,y=log2FoldChange))
  result=list(FW.vsall.filt,test)
  return(result)
  
}


normalie_chest<-function(evenbefore1,evenbefore2,factor){
  #ratio=before/after
  #factor=ratio[1,]
  #print(ratio)
  before.norm1=as.data.frame(mapply('/', evenbefore1, factor[1:3]))
  before.norm2=as.data.frame(mapply('/', evenbefore2, factor[4:5]))
  rownames(before.norm1)=rownames(evenbefore1)
  rownames(before.norm2)=rownames(evenbefore2)
  lista=list(before.norm1,before.norm2)
  return(lista)
}
