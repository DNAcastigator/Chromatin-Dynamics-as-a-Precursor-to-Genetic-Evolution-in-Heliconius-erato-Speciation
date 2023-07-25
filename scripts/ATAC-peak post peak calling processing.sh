


sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g' ../../../MACS2_sub/final_plan/demophoon.strict.HW.unique.bed | bedtools intersect -wa -a demophoon.strict.HW.all.count.txt -b - > demophoon.strict.HW.unique.count.txt

sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g' ../../../MACS2_sub/final_plan/hydara.strict.HW.unique.bed | bedtools intersect -wa -a hydara.strict.HW.all.count.txt -b - > hydara.strict.HW.unique.count.txt



##################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

################################################################################################################################################################## Forewing and hindwing


cat LI7_demophoon_FW.MACSQ.strict_peaks.narrowPeak FW-pboy.MACSQ.strict_peaks.narrowPeak LB_41.MACSQ.strict_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1>2' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - | bedtools merge -i - > final_plan/demophoon.strict.FW.3of3.bed


cat LB_20.MACSQ.strict_peaks.narrowPeak LB_27.MACSQ.strict_peaks.narrowPeak LB_29.MACSQ.strict_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1>2' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - | bedtools merge -i - > final_plan/hydara.strict.FW.3of3.bed




cat LI7_demophoon_HW.MACSQ.strict_peaks.narrowPeak E3_HW.MACSQ.strict_peaks.narrowPeak LB_42.MACSQ.strict_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1>2' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - | bedtools merge -i - > final_plan/demophoon.strict.HW.3of3.bed


cat LB_21.MACSQ.strict_peaks.narrowPeak LB_28.MACSQ.strict_peaks.narrowPeak LB_30.MACSQ.strict_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1>2' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - | bedtools merge -i - > final_plan/hydara.strict.HW.3of3.bed

############# 2 out of 3


cat LI7_demophoon_FW.MACSQ.strict_peaks.narrowPeak FW-pboy.MACSQ.strict_peaks.narrowPeak LB_41.MACSQ.strict_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1==2' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - | bedtools merge -i - > final_plan/demophoon.strict.FW.2of3.bed

cat LB_20.MACSQ.strict_peaks.narrowPeak LB_27.MACSQ.strict_peaks.narrowPeak LB_29.MACSQ.strict_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1==2' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - | bedtools merge -i - > final_plan/hydara.strict.FW.2of3.bed



cat LI7_demophoon_HW.MACSQ.strict_peaks.narrowPeak E3_HW.MACSQ.strict_peaks.narrowPeak LB_42.MACSQ.strict_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1==2' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - | bedtools merge -i - > final_plan/demophoon.strict.HW.2of3.bed

cat LB_21.MACSQ.strict_peaks.narrowPeak LB_28.MACSQ.strict_peaks.narrowPeak LB_30.MACSQ.strict_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1==2' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - | bedtools merge -i - > final_plan/hydara.strict.HW.2of3.bed


################## all
cat LI7_demophoon_FW.MACSQ.strict_peaks.narrowPeak FW-pboy.MACSQ.strict_peaks.narrowPeak LB_41.MACSQ.strict_peaks.narrowPeak |sort -k2,2n -k3,3n | bedtools merge -i - >final_plan/demophoon.strict.FW.all.bed

cat LB_20.MACSQ.strict_peaks.narrowPeak LB_27.MACSQ.strict_peaks.narrowPeak LB_29.MACSQ.strict_peaks.narrowPeak |sort -k2,2n -k3,3n | bedtools merge -i - >final_plan/hydara.strict.FW.all.bed


cat LI7_demophoon_HW.MACSQ.strict_peaks.narrowPeak E3_HW.MACSQ.strict_peaks.narrowPeak LB_42.MACSQ.strict_peaks.narrowPeak |sort -k2,2n -k3,3n | bedtools merge -i - >final_plan/demophoon.strict.HW.all.bed

cat LB_21.MACSQ.strict_peaks.narrowPeak LB_28.MACSQ.strict_peaks.narrowPeak LB_30.MACSQ.strict_peaks.narrowPeak |sort -k2,2n -k3,3n | bedtools merge -i - >final_plan/hydara.strict.HW.all.bed


##########unique

bedtools subtract -A -a demophoon.strict.FW.3of3.bed -b hydara.strict.FW.all.bed >demophoon.strict.FW.unique.bed
bedtools subtract -A -a hydara.strict.FW.3of3.bed -b demophoon.strict.FW.all.bed > hydara.strict.FW.unique.bed

bedtools subtract -A -a demophoon.strict.HW.3of3.bed -b hydara.strict.HW.all.bed >demophoon.strict.HW.unique.bed
bedtools subtract -A -a hydara.strict.HW.3of3.bed -b demophoon.strict.HW.all.bed > hydara.strict.HW.unique.bed

##############################################################



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


###### cat Not1_FW_peaks.narrowPeak Not2_FW_peaks.narrowPeak E_Not3_FW_peaks.narrowPeak | sort -k2,2n -k3,3n | bedmap --count --echo-map-range --fraction-both 0.01 --delim '\t' - | awk '$1>1' - | cut -f2- -| sort -k2,2n -k3,3n | uniq - >pan.notabilis.FW.bed
####
ID=$((SLURM_ARRAY_TASK_ID -1))

#####panome
DIR=/work/rpapa/aruggieri/eratopop/ATAC/pangenome


#####   demophoon HW

samples=(LI7_demophoon_FW FW-pboy LB_41)


bedtools coverage -counts -b "$DIR"/BAM/$(echo "${samples[ID]}").cut.trim.filtered.sorted.rm.bam -a "$DIR"/MACS2_sub/final_plan/demophoon.strict.FW.all.bed > "$DIR"/counts_sub/final_counts/$(echo "${samples[ID]}").all.count.out



######### hydara


samples=(LB_20 LB_27 LB_29)


bedtools coverage -counts -b "$DIR"/BAM/$(echo "${samples[ID]}").cut.trim.filtered.sorted.rm.bam -a "$DIR"/MACS2_sub/final_plan/hydara.strict.FW.all.bed > "$DIR"/counts_sub/final_counts/$(echo "${samples[ID]}").all.count.out



samples=(LI7_demophoon_HW E3_HW LB_42))


bedtools coverage -counts -b "$DIR"/BAM/$(echo "${samples[ID]}").cut.trim.filtered.sorted.rm.bam -a "$DIR"/MACS2_sub/final_plan/demophoon.strict.HW.all.bed > "$DIR"/counts_sub/final_counts/$(echo "${samples[ID]}").all.count.out



######### hydara


samples=(LB_21 LB_28 LB_30)


bedtools coverage -counts -b "$DIR"/BAM/$(echo "${samples[ID]}").cut.trim.filtered.sorted.rm.bam -a "$DIR"/MACS2_sub/final_plan/hydara.strict.HW.all.bed > "$DIR"/counts_sub/final_counts/$(echo "${samples[ID]}").all.count.out




######################## afetr couting with ALL

paste LI7_demophoon_FW.all.count.out FW-pboy.all.count.out LB_41.all.count.out | cut -f1,2,3,4,8,12 |sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g'>all_count/demophoon.strict.FW.all.count.txt


paste LB_20.all.count.out LB_27.all.count.out LB_29.all.count.out | cut -f1,2,3,4,8,12 |sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g' >all_count/hydara.strict.FW.all.count.txt



paste LI7_demophoon_HW.all.count.out E3_HW.all.count.out LB_42.all.count.out | cut -f1,2,3,4,8,12 |sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g'>all_count/demophoon.strict.HW.all.count.txt


paste LB_21.all.count.out LB_28.all.count.out LB_30.all.count.out | cut -f1,2,3,4,8,12 |sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g' >all_count/hydara.strict.HW.all.count.txt


######################## FW



sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g' ../../../MACS2_sub/final_plan/demophoon.strict.FW.3of3.bed | bedtools intersect -wa -a demophoon.strict.FW.all.count.txt -b - > demophoon.strict.FW.3of3.count.txt


sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g' ../../../MACS2_sub/final_plan/hydara.strict.FW.3of3.bed | bedtools intersect -wa -a hydara.strict.FW.all.count.txt -b - > hydara.strict.FW.3of3.count.txt



bedmap --echo --echo-map --fraction-either 0.5 --count demophoon.strict.FW.3of3.count.txt  hydara.strict.FW.3of3.count.txt |  awk -F"|" '$3==1' | cut -f1,2 -d"|" | tr "|" "\t"   > demoxhydara.strict.FW.3of3.single.txt

bedmap --echo --echo-map --fraction-either 0.5 --count demophoon.strict.FW.3of3.count.txt  hydara.strict.FW.3of3.count.txt | grep ";"  > demoxhydara.strict.FW.3of3.multy.txt

while read line; do  uno=$(echo $line | cut -f1 -d"|"); echo $uno >>demoxhyda.strict.FW.3of3.1.txt;  echo $line | cut -f2 -d"|" | awk -F";" '{for (i=1;i<=NF;i++) printf("%s\n",$i)}'|tr " " "\t" | bedtools merge -i - -d 100000000 -c 4,5,6 -o sum,sum,sum  >>demoxhyda.strict.FW.3of3.2.txt; done<demoxhydara.strict.FW.3of3.multy.txt
paste demoxhyda.strict.FW.3of3.1.txt demoxhyda.strict.FW.3of3.2.txt | tr " " "\t" >demoxhydara.strict.FW.3of3.multy.fixed.txt

cat demoxhydara.strict.FW.3of3.multy.fixed.txt demoxhydara.strict.FW.3of3.single.txt | sort -k2,2n -k3,3n >demoxhydara.strict.FW.3of3.almostfinal.txt

cat demoxhydara.strict.FW.3of3.almostfinal.txt | bedtools merge -i - -c 4,5,6,7,8,9,10,11,12 -o sum,sum,sum,distinct,min,max,sum,sum,sum | awk 'BEGIN {FS="\t"; OFS="\t"} {print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}' | bedtools merge -i - -c 4,5,6,7,8,9,10,11,12 -o sum,sum,sum,distinct,min,max,sum,sum,sum | awk 'BEGIN {FS="\t"; OFS="\t"} {print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}'>demoxhydara.strict.FW.3of3.final.txt




##################### HW

sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g' ../../../MACS2_sub/final_plan/demophoon.strict.HW.3of3.bed | bedtools intersect -wa -a demophoon.strict.HW.all.count.txt -b - > demophoon.strict.HW.3of3.count.txt


sed 's/SeqSeqPan;0|SeqSeqPan_dem_hyd_fav_ety_not_che.xmfa/pan/g' ../../../MACS2_sub/final_plan/hydara.strict.HW.3of3.bed | bedtools intersect -wa -a hydara.strict.HW.all.count.txt -b - > hydara.strict.HW.3of3.count.txt




bedmap --echo --echo-map --fraction-either 0.5 --count demophoon.strict.HW.3of3.count.txt  hydara.strict.HW.3of3.count.txt |  awk -F"|" '$3==1' | cut -f1,2 -d"|" | tr "|" "\t"   > demoxhydara.strict.HW.3of3.single.txt

bedmap --echo --echo-map --fraction-either 0.5 --count demophoon.strict.HW.3of3.count.txt  hydara.strict.HW.3of3.count.txt | grep ";"  > demoxhydara.strict.HW.3of3.multy.txt

while read line; do  uno=$(echo $line | cut -f1 -d"|"); echo $uno >>demoxhyda.strict.HW.3of3.1.txt;  echo $line | cut -f2 -d"|" | awk -F";" '{for (i=1;i<=NF;i++) printf("%s\n",$i)}'|tr " " "\t" | bedtools merge -i - -d 100000000 -c 4,5,6 -o sum,sum,sum  >>demoxhyda.strict.HW.3of3.2.txt; done<demoxhydara.strict.HW.3of3.multy.txt
paste demoxhyda.strict.HW.3of3.1.txt demoxhyda.strict.HW.3of3.2.txt | tr " " "\t" >demoxhydara.strict.HW.3of3.multy.fixed.txt

cat demoxhydara.strict.HW.3of3.multy.fixed.txt demoxhydara.strict.HW.3of3.single.txt | sort -k2,2n -k3,3n >demoxhydara.strict.HW.3of3.almostfinal.txt

cat demoxhydara.strict.HW.3of3.almostfinal.txt | bedtools merge -i - -c 4,5,6,7,8,9,10,11,12 -o sum,sum,sum,distinct,min,max,sum,sum,sum | awk 'BEGIN {FS="\t"; OFS="\t"} {print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}' | bedtools merge -i - -c 4,5,6,7,8,9,10,11,12 -o sum,sum,sum,distinct,min,max,sum,sum,sum | awk 'BEGIN {FS="\t"; OFS="\t"} {print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}'>demoxhydara.strict.HW.3of3.final.txt
















