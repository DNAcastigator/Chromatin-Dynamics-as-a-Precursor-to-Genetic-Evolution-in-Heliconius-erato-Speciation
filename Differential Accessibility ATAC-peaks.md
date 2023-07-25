# Differential Accessibility pipeline
The DE analysis is divided into two main parts: Removal of unwanted variation and actual Differential accessibility. I used a single custom function to perform all the steps together but I will break it down here.
## RUV-seq
description----------


```
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
```
