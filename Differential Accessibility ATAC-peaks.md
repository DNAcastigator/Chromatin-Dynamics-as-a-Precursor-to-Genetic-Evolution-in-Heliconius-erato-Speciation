# Differential Accessibility pipeline
The DE analysis is divided into two main parts: Removal of unwanted variation (RUVseq) and actual Differential accessibility.
I used a single custom function to perform all the steps together but I will break it down here.

The full script for this step can be found [here](https://github.com/DNAcastigator/summer-project/blob/main/scripts/DE.ATAC.R)

## RUVseq (R)
The R package 'Removal of Unwanted Variation' eliminates various sources of noise in the data, including batch, library preparation, and other nuisance effects. It employs between-sample normalization methods as proposed in Risso et al. (2014).

in this work, I used specifically RUVr, which removes variability based on the residuals (e.g., deviance residuals) obtained from a preliminary GLM regression of the counts against the covariates of interest. Prior to this step, you may consider applying normalization using a method like upper-quartile normalization.

In this first part of the function, we proceed with upper-quartile normalization and then with the calculation of residuals; note that `subset` is the gene count matrix and `subsetheader` contains the sample metadata.

```R
RUVparam<-function(subset,subsetheader,plot,fact,K,con){
  
  subset=select(subset,c(subsetheader$sample))
    filter <- apply(subset, 1, function(x) length(x[x>5])>5)
  filtered <- subset[filter,]
  
  set <- newSeqExpressionSet(as.matrix(filtered),
                             phenoData = data.frame(subsetheader, row.names=colnames(filtered)))
  set <- betweenLaneNormalization(set, which="upper")
 
  design <- model.matrix(~subsetheader$species, data=pData(set))
  
  y <- DGEList(counts=counts(set), group=subsetheader$species)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  res <- residuals(fit, type="deviance")
  genes <- rownames(filtered)
```
In the following step, a series of scatterplots are plotted to show the effect of different factors of unwanted variation, the user may then input the desired number of factos.
```R
  
  colors <- brewer.pal(12, "Set3")
  ############
  if (plot == "yes"){
    par(mfrow = c(2,2))
    plotPCA(set, col=colors[as.factor(subsetheader$species)], cex=1.2,main="raw")
    
    for(k in 1:3) {
      set.x=RUVr(set, genes, k=k, res)
      plotPCA(set.x, col=colors[as.factor(subsetheader$species)], cex=1.2,main=k)
    }
  }
  if (fact == "no"){
  n <- readline(prompt="Enter desired design (in W_): ")  
  }else{
    n=K
  }
```

We then add the desired factors in the DEseq2 design formula and proceed with the DE analysis
```R
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
  test=sizeFactors(atacDDSXxYlFW)
  res_FW <- results(atacDDSXxYlFW,contrast=con,lfcThreshold=0.5,alpha=0.05)
```
after some filtering and cleaning, the function returns the DA peaks output from DEseq2 as well as the factor of normalization obtained from the analysis for each sample.
```R
FW.vsall.filt=na.omit(data.frame(res_FW))
  print(summary(res_FW))
  FW.vsall.filt$pan=ifelse(FW.vsall.filt$log2FoldChange > 0,str_split_fixed(rownames(FW.vsall.filt),"-",2)[,1],str_split_fixed(rownames(FW.vsall.filt),"-",2)[,2])
  FW.vsall.filt$pos=as.numeric(str_split_fixed(FW.vsall.filt$pan,"_",3)[,2])
  FW.vsall.filt$pos2=as.numeric(str_split_fixed(FW.vsall.filt$pan,"_",3)[,3])
  FW.vsall.filt=FW.vsall.filt[FW.vsall.filt$pos2<688640245,]
  result=list(FW.vsall.filt,test)
  return(result)
  
}
```
a real example of this function would look something like this (using H. e. demophoon and H. e. hydara):
```R
demoxhyda.FW3of3=RUVparam(demoxhyda3of3,header.demoxhyda.FW,"yes","no","W_1+W_2",c("species","dem","hyd"))
```
Finally, we can now use the normalization factor from the previous function to normalize the totality of the peaks and not only the ones that are shared among the two populations, this will be useful for plotting later
```R
normalie<-function(evenbefore1,evenbefore2,factor){
before.norm1=as.data.frame(mapply('/', evenbefore1, factor[1:3]))
before.norm2=as.data.frame(mapply('/', evenbefore2, factor[4:6]))
rownames(before.norm1)=rownames(evenbefore1)
rownames(before.norm2)=rownames(evenbefore2)
lista=list(before.norm1,before.norm2)
return(lista)
}


demxhyda.norm.all=normalie(demophoonFWall,hydaraFWall,demoxhyda.FW3of3[[2]])
```
Now is possible to use this data to generate plots like the one in Figure 1 of the main paper, the script is available [here](https://github.com/DNAcastigator/summer-project/blob/main/scripts/genomewide.plot.functions.R)
