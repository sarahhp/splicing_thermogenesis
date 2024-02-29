
## Generate not_sig promoter windows
bedtools slop -s -l 2000 -r 500 -i ../histone_profile/lc_3db_not_sig_TSSes.bed -g ../../annotations/hg38.chrom.sizes |\
 bedtools sort -i stdin | bedtools merge -s -i stdin -c 4,5,6,7 -o distinct,absmax,distinct,absmin > lc_3db_not_sig_promoters.merged.bed

## Get white-beige intersecting promoters
bedtools intersect -a lc_3db_beige_promoters.bed -b lc_3db_white_promoters.bed  -wa  > white_and_beige_promoters.beige.bed
bedtools intersect -a lc_3db_white_promoters.bed -b lc_3db_beige_promoters.bed  -wa  > white_and_beige_promoters.white.bed
cat white_and_beige_promoters.white.bed white_and_beige_promoters.beige.bed |\
 bedtools sort -i stdin | bedtools merge -s -i stdin -c 4,5,6,7 -o distinct,absmax,distinct,absmin > white_and_beige_promoters.merged.bed

## Subtract intersects and merge overlapping promoters white and beige
bedtools subtract -a lc_3db_beige_promoters.bed -b white_and_beige_promoters.merged.bed |\
 bedtools merge -s -i stdin -c 4,5,6,7 -o distinct,absmax,distinct,absmin > lc_3db_beige_promoters.uniq.bed
bedtools subtract -a lc_3db_white_promoters.bed -b white_and_beige_promoters.merged.bed |\
 bedtools merge -s -i stdin -c 4,5,6,7 -o distinct,absmax,distinct,absmin > lc_3db_white_promoters.uniq.bed


