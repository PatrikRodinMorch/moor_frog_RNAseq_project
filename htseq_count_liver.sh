#!/bin/bash -l


module load bioinfo-tools
module load htseq

bam_path=/crex/proj/snic2019-30-54/reference_guided/liver
gtf=/crex/proj/snic2019-30-54/reference_guided/GCF_905171775.1_aRanTem1.1_genomic.gtf

htseq-count -f bam -r name --strand=reverse \
$bam_path/Sample_SL-2404-P1sample1_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample9_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample13_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample20_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample26_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample27_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample29_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample36_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample37_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample38_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample41_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample45_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample48_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample58_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample59_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample61_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample62_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample71_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample77_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample78_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample79_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample86_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample87_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample89_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample96_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample98_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample99_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample113_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample114_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample120_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample121_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample122_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample125_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample134_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample135_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample136_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample137_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample138_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample155_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample156_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample161_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample167_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample170_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample174_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample175_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample180_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample187_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample191_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample5_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample6_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample14_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample17_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample18_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample21_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample22_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample30_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample31_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample42_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample49_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample50_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample54_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample56_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample63_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample64_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample69_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample73_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample75_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample80_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample84_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P1sample93_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample100_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample104_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample105_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample106_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample110_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample116_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample117_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample128_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample129_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample132_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample139_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample145_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample146_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample149_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample151_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample152_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample157_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample162_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample163_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample171_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample176_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample181_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample185_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample188_Aligned.sortedByCoord.out.bam \
$bam_path/Sample_SL-2404-P2sample192_Aligned.sortedByCoord.out.bam \
$gtf > gene_counts_htseq_liver.txt
