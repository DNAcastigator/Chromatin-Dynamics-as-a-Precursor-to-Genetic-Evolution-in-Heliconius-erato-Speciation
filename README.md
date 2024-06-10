# summer-project pipeline
Hello there

this is a brief summary of the pipelines used in the Ruggieri et al 2024 paper. I included tutorials and full scripts of those parts that may be new:
## Genomes assemblies and pan-genome generation
These processes have already been described in  Ruggieri et al 2022, and you can find information in [Dr. Steven van Belleghem github](https://github.com/StevenVB12/Genomics) 
## Resequencing data and population genomics
I will not indulge in these basic processes since they have already been discussed in previous papers.
### Resequencing
-BWA (aligning to pan-genome)

-Picard (remove PCR duplicates)

-GATK (genotype calling and quality filtering)
### pop genetics
-SweepFinder2 (selective sweep)

-Steve Martin's [popgenWindows.py](https://github.com/simonhmartin/genomics_general)) (Fst)

## ATAC-seq data preprocessing

-trimmomatic (cleaning)

-bowtie2 (aligning to reference genome)

-Picard (remove PCR duplicate)

-MACS2 (peak calling)


## TUTORIALS

After the peaks have been identified for each sample, the first step is to select the right peaks subset, count the reads, and identify the peaks that are unique and the ones that
are shared among two pair of populations -> [ATAC-seq analysis](https://github.com/DNAcastigator/summer-project/blob/main/ATAC-seq%20analysis.md) (tutorial)

The following step is to perform Differential accessibility analysis ->[ATAC-peaks DA](https://github.com/DNAcastigator/summer-project/blob/main/Differential%20Accessibility%20ATAC-peaks.md) (tutorial)

It is now possible to create some plots as in Figure 1 format the main paper <-[genomewide plots](https://github.com/DNAcastigator/summer-project/blob/main/scripts/genomewide.plot.functions.R) (script)

Alternatively, there is the possibility to perform some statistical analysis to capture a possible connection between the values of Fst and the ATAC-seq position

-[K-S test](https://github.com/DNAcastigator/summer-project/blob/main/Kolmogorov%20Smirnov%20test.md) (tutorial)

-[Custom Binomial Test](https://github.com/DNAcastigator/summer-project/blob/main/Custom%20binomial%20test.md) (tutorial)

### Extra analysis
Other analyses you may find interesting concern the identification of gene-wise Fst signal in the Fst data calculated in windows of 1k.-> [noise vs signal for fst](https://github.com/DNAcastigator/summer-project/blob/main/signal%20vs%20noise%20Fst.md)
