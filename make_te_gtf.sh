#!/bin/bash -l


module load bioinfo-tools
module load star/2.7.9a

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir . \
--genomeFastaFiles /crex/proj/snic2019-30-54/reference_guided/GCF_905171775.1_aRanTem1.1_genomic.fna \
--sjdbGTFfile /crex/proj/snic2019-30-54/reference_guided/GCF_905171775.1_aRanTem1.1_genomic.gtf \
--sjdbOverhang 149 \
--outFileNamePrefix GCF_905171775.1_aRanTem1.1
