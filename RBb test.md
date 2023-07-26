# Ruggieri-Bellin binomial test
The RBb test splits the genome into 10 Fst intervals (increasing by 0.1), then calculates an expected probability to find a unique ATAC peak in each interval (using the average of 1000 random sampling), and then performs a binomial test for each interval using the empirical number of unique ATAC-peak as “number of successes”, the total number of ATAC peaks in that Fst interval as “number of trials” and the expected probability as “hypothesized probability of success”. The test outputs a p value for each test (null hypothesis “greater”) that then are adjusted using holm correction.
the following picture tries to summarize the previous, rather wordy explanation.
![image](https://github.com/DNAcastigator/summer-project/assets/47642926/a83005c1-ac3f-4231-872d-f297ce2b8404)

the function that calculates it takes as input 1) the genomewide fst intervals distribution for the two populations studied (calculated with another custom function, reported later); 2) the position of the unique ATAC-seq for the two population studied; 3) the position of all the ATAC peaks identified in the two populations:
```R
Fst.stat=function(fsts,position,total_pos,n,totplot){
  
position.range=makeGRangesFromDataFrame(DataFrame(chrom="pan",start=position,end=position+1))
position.total=makeGRangesFromDataFrame(DataFrame(chrom="pan",start=total_pos,end=total_pos+1))
test=makeGRangesFromDataFrame(fsts,keep.extra.columns=T)

test$over=countOverlaps(test,position.range)
test$over_total=countOverlaps(test,position.total)
resu=data.frame(test) %>% group_by(sigla) %>% summarise(final_count = sum(over),total=sum(over_total)) %>% data.frame()


random=data.frame(base::transform(test,sigla=sample(test$sigla)))%>% group_by(sigla) %>% summarise(final_count = sum(over)/sum(over_total)) %>% data.frame()
for (i in seq(1:n)){
ran=data.frame(base::transform(test,sigla=sample(test$sigla)))%>% group_by(sigla) %>% summarise(final_count = sum(over)/sum(over_total)) %>% data.frame()
random=rbind(random,ran)
}

resu_random=na.omit(random) %>% group_by(sigla) %>% summarise(random_count = mean(final_count),se=sqrt(sd(final_count)/n())) %>% data.frame()
################

confronto=cbind(resu,resu_random[c(2,3)])
confronto$random_count[!is.finite(confronto$random_count)] <- 0

confronto=confronto[confronto$total>0,]
statistic=apply(confronto,1,function(x){
extest=binom.test(as.numeric(x[2]), as.numeric(x[3]), p = as.numeric(x[4]),
           alternative =  "g",
           conf.level = 0.95)

return(extest$p.value)
})
```
The function then produces some violin plot with the random sampling distributions and adds a colored point that represents the empirical probability to find a unique ATAC-peak in the specific fst intervals (x-axis), the asterisks show where the binomial test was significant (the p.values are adjusted with holm), an example:
![demophoonvsfavorinus corr DEunique1 zoom crop](https://github.com/DNAcastigator/summer-project/assets/47642926/e65994ec-6676-4e3f-b2b5-497de8b67dc1)
```R

confronto=cbind(confronto,statistic)
g=ggplot(confronto)+
  geom_violin(data=random, aes(x=sigla, y=final_count)) + 
  geom_point(aes(x=sigla,y=final_count/total),col="red",shape=ifelse(p.adjust(confronto$statisti,method = "holm")<0.05,8,16))+
  #geom_point(aes(x=sigla,y=random_count),col="blue")+
  #geom_errorbar(aes(x=sigla,ymin = random_count-se, ymax = random_count+se),col="blue")+
  {if(totplot=="yes") ylim(c(0,1)) }+
  theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1),axis.title.y=element_blank())
return(list(g,confronto))
}
```
two functions needed to process the data before they can be used for RBb test and calculate the genomewide fst intervals distribution are:
```R
Fst.dist=function(data,limit,sigla){
    uno=na.omit(data[data$Fst>=limit[1] & data$Fst<limit[2],])
  due=data.frame(chrom="pan",start=uno$start,end=uno$end)
  due$sigla=sigla
  return(due)
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

```
So, an example of how to run this test is this:
```R
position1=demophoon.HW.unique.bed$start+((demophoon.HW.unique.bed$end-demophoon.HW.unique.bed$end)/2)
position2=demophoon.FW.unique.bed$start+((demophoon.FW.unique.bed$end-demophoon.FW.unique.bed$end)/2)
position3=hydara.HW.unique.bed$start+((hydara.HW.unique.bed$end-hydara.HW.unique.bed$end)/2)
position4=hydara.FW.unique.bed$start+((hydara.FW.unique.bed$end-hydara.FW.unique.bed$end)/2)
position_unique=c(position1,position2,position3,position4)


demophoonxhydara.complete=data.frame(str_split_fixed(c(rownames(demophoonFWall),rownames(demophoonHWall),rownames(hydaraFWall),rownames(hydaraHWall)),"_",3))
colnames(demophoonxhydara.complete)=c("chrom","start","end")
demophoonxhydara.complete$start=as.integer(as.character(demophoonxhydara.complete$start))
demophoonxhydara.complete$end=as.integer(as.character(demophoonxhydara.complete$end))

demhyd_total_pos=demophoonxhydara.complete$start+((demophoonxhydara.complete$end-demophoonxhydara.complete$start)/2)
########################################################################### test the functions
freq.df_demohydFst_1k=freq_fun(demxhd.1k.fst,demhyd_total_pos,"Fst")

fst.dist.demohyda=Fst.dist.loop(demxhd.1k.fst)
split2=Fst.stat(fst.dist.demohyda,position_Deunique_demohyda,notxety_total_pos,1000,"no")
```
the full script used in this work (that also performs [Ks.test](https://github.com/DNAcastigator/summer-project/blob/main/Kolmogorov%20Smirnov%20test.md)) is [here](https://github.com/DNAcastigator/summer-project/blob/main/scripts/correlation_statistic.R)
