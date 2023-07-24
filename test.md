# ATAC-seq post peak calling processing
## reads counting (BASH)
The actual script for this step can be found [here](https://github.com/DNAcastigator/summer-project/blob/main/ATAC-peak%20post%20peak%20calling%20processing.sh)

After the basic reads alignment and peak calling, the first step is to combine the ATAC peaks identified per population (pooling the three samples for each population). this job has to be done for Forewing and Hindwing for each population, for the peaks present in 3 samples out of 3, and for the one found in at least 2 samples out of 3.
the following example shows one population (H e. demophoon), forewing, 3 out of 3 samples:
```
cat LI7_demophoon_FW.MACSQ.strict_peaks.narrowPeak FW-pboy.MACSQ.strict_peaks.narrowPeak LB_41.MACSQ.strict_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1>2' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - | bedtools merge -i - > final_plan/demophoon.strict.FW.3of3.bed
```
to find the list of all the peaks identified in each population, the code is a lot simpler:
```
cat LI7_demophoon_FW.MACSQ.strict_peaks.narrowPeak FW-pboy.MACSQ.strict_peaks.narrowPeak LB_41.MACSQ.strict_peaks.narrowPeak |sort -k2,2n -k3,3n | bedtools merge -i - >final_plan/demophoon.strict.FW.all.bed
```
finally, to identify the subset of peaks unique in a population, we need to subtract from the 3 out of 3 peaks set of our population of interest, the totality of the peaks identified in the other population we are comparing, In this example H. e. demophoon and H. e. hydara:
```
bedtools subtract -A -a demophoon.strict.FW.3of3.bed -b hydara.strict.FW.all.bed >demophoon.strict.FW.unique.bed
```
Now that we have all the peak subsets that we need, we can proceed to count the reads that align in the specific intervals of the ATAC peaks. The count is performed on the total list of peaks identified in a population, they will be subsetted later, so that the counting proceddure happens only once.
In this example the analysis is run on a Linux cluster system , regulated with slurm:
```
#!/bin/bash
#SBATCH --mem-per-cpu=100gb
#SBATCH --time=72:00:00
#SBATCH --job-name=ATAC.pan
#SBATCH --error=results/ATAC.pan.count.2022-%a.err
#SBATCH --output=results/ATAC.pan.count.2022-%a.out
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

module load bedtools
module load samtools

####
ID=$((SLURM_ARRAY_TASK_ID -1))

#####input data directory
DIR=/work/rpapa/aruggieri/eratopop/ATAC/pangenome

#####   demophoon FW

samples=(LI7_demophoon_FW FW-pboy LB_41)

bedtools coverage -counts -b "$DIR"/BAM/$(echo "${samples[ID]}").cut.trim.filtered.sorted.rm.bam -a "$DIR"/MACS2_sub/final_plan/demophoon.strict.FW.all.bed > "$DIR"/counts_sub/final_counts/$(echo "${samples[ID]}").all.count.out
```
We then combine the results from each sample of a population together and we subset it to select only the list of 3 out of 3 peaks:
```
paste LI7_demophoon_FW.all.count.out FW-pboy.all.count.out LB_41.all.count.out | cut -f1,2,3,4,8,12 |sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g'>all_count/demophoon.strict.FW.all.count.txt
sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g' ../../../MACS2_sub/final_plan/demophoon.strict.FW.3of3.bed | bedtools intersect -wa -a demophoon.strict.FW.all.count.txt -b - > demophoon.strict.FW.3of3.count.txt
```

OK, this final step is a little rough, so stay with me.
The last step allows us to find those peaks that are shared between two different populations, the main problem is that some bigger peaks in one population aligned with multiple smaller peaks in the other population, so the smaller peaks were merged, and the reads were summed
So, First, we identified the shared peaks that don't give any problem, using always h. e.  demophoon and h. e.  hydara as example:
```
bedmap --echo --echo-map --fraction-either 0.5 --count demophoon.strict.FW.3of3.count.txt  hydara.strict.FW.3of3.count.txt |  awk -F"|" '$3==1' | cut -f1,2 -d"|" | tr "|" "\t"   > demoxhydara.strict.FW.3of3.single.txt
```
then we found those more complicated cases:
```
bedmap --echo --echo-map --fraction-either 0.5 --count demophoon.strict.FW.3of3.count.txt  hydara.strict.FW.3of3.count.txt | grep ";"  > demoxhydara.strict.FW.3of3.multy.txt
```
In this last file `demoxhydara.strict.FW.3of3.multy.txt` the pipe `|` separates the two population and semilcolon `;` separate the multiple peaks from the second population that align with one of the first population. to remove this problem we use:
```
while read line;
do
  uno=$(echo $line | cut -f1 -d"|");
  echo $uno >>demoxhyda.strict.FW.3of3.1.txt;
  echo $line | cut -f2 -d"|" | awk -F";" '{for (i=1;i<=NF;i++) printf("%s\n",$i)}'|tr " " "\t" | bedtools merge -i - -d 100000000 -c 4,5,6 -o sum,sum,sum  >>demoxhyda.strict.FW.3of3.2.txt;
done<demoxhydara.strict.FW.3of3.multy.txt

paste demoxhyda.strict.FW.3of3.1.txt demoxhyda.strict.FW.3of3.2.txt | tr " " "\t" >demoxhydara.strict.FW.3of3.multy.fixed.txt
cat demoxhydara.strict.FW.3of3.multy.fixed.txt demoxhydara.strict.FW.3of3.single.txt | sort -k2,2n -k3,3n >demoxhydara.strict.FW.3of3.almostfinal.txt
```
Finally, we want to solve the problem of the multiple peaks from the first population that align with a single bigger one from the second (that is, the opposite situation described before)
```
cat demoxhydara.strict.FW.3of3.almostfinal.txt | bedtools merge -i - -c 4,5,6,7,8,9,10,11,12 -o sum,sum,sum,distinct,min,max,sum,sum,sum | awk 'BEGIN {FS="\t"; OFS="\t"} {print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}' | bedtools merge -i - -c 4,5,6,7,8,9,10,11,12 -o sum,sum,sum,distinct,min,max,sum,sum,sum | awk 'BEGIN {FS="\t"; OFS="\t"} {print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}'>demoxhydara.strict.FW.3of3.final.txt
```

we are now ready to proceed with the differential accessibility analysis 

