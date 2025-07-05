#!/bin/bash -l


module load bioinfo-tools
module load RepeatMasker


awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' GCF_905171775.1_aRanTem1.1_genomic.fna > GCF_905171775.1_aRanTem1.1_genomic_upper.fna
RepeatMasker -pa 20 -s -species Silurana GCF_905171775.1_aRanTem1.1_genomic_upper.fna
