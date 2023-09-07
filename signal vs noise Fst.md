
The Fst analysis is often performed on wider window sizes (50-10k bp) than those used in this paper. In our work, we aimed to identify signals that may be too narrow for the wider window sizes, so we employed 1k bp windows. However, this approach results in an increase in background noise since there are more data points that may show a higher signal due to technical issues such as genotyping errors, which may not reflect actual biological differences. To address this issue in ATAC-seq analysis and Fst, we employed a set of statistical tests that you can find in our previous project:

- [K-S test](https://github.com/DNAcastigator/summer-project/blob/main/Kolmogorov%20Smirnov%20test.md)
- [Ruggieri-Bellin binomial test](https://github.com/DNAcastigator/summer-project/blob/main/RBb%20test.md)

However, to identify broader signals that may pertain to genes rather than regulatory elements, we developed a different strategy.
Initially, the Fst data was segmented into windows of 10 points, wherein we calculated both the number of points exceeding the average Fst value and the maximum Fst value within each window.
To determine the authenticity of a signal, we established specific criteria. 

-The signal must demonstrate a higher point density above the average (density exceeding the 0.99 quantile of the overall density calculation).

-It must include at least one data point surpassing the 0.99 percentile of the entire Fst distribution.

in R the function looks something liek this:

```R
lil_function_fst_density=function(data){
  data=na.omit(data)
data$group <- as.numeric(cut_number(data$position, nrow(data)/10))
#trashold=quantile(data$Fst,0.95)
trashold=mean(data$Fst)
 testfst=data %>% group_by(group) %>% summarise(max=max(Fst),density=sum(Fst>trashold))
 if (quantile(testfst$density,0.99)==max(testfst$density)){
   testpos=testfst[testfst$density>=quantile(testfst$density,0.99) & testfst$max>quantile(data$Fst,0.99),] 
    } else {
 testpos=testfst[testfst$density>quantile(testfst$density,0.99) & testfst$max>quantile(data$Fst,0.99),]
    }
 newdata=data[data$group %in% testpos$group,]
 position=newdata %>% group_by(group) %>% summarise(start=min(position),end=max(position)) %>% mutate(group="pan") %>% data.frame()
  return(position)
}

```

Here's an example of a genuine signal identified by the function (multiple windows combined):
![real_fst_signal](https://github.com/DNAcastigator/summer-project/assets/47642926/15a011d2-2bd5-4f54-9bc7-51a7e041bbda)

And here's what a "fake" signal looks like (note that the fake signal still contains some points above the max Fst threshold, but lacks a higher point density above the average):
![fake_fst_signal](https://github.com/DNAcastigator/summer-project/assets/47642926/2fa570a7-2314-419b-b975-e96debf84a9d)
