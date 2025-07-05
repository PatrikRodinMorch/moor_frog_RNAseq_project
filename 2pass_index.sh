#!/bin/bash -l


module load bioinfo-tools
module load star/2.7.9a

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir ./2pass_index/ \
--genomeFastaFiles ../GCF_905171775.1_aRanTem1.1_genomic.fna \
--sjdbGTFfile ../GCF_905171775.1_aRanTem1.1_genomic.gtf \
--sjdbFileChrStartEnd 2pass_SJ_out_filtered.tab \
--sjdbOverhang 149 \
--outFileNamePrefix GCF_905171775.1_aRanTem1.1


