#!/bin/bash -l


ml bioinfo-tools
ml AGAT

agat_sp_add_introns.pl --gff GCF_905171775.1_aRanTem1.1_genomic.gff --out GCF_905171775.1_aRanTem1.1_genomic_with_INTRONs.gff


agat_sp_manage_introns.pl --gff GCF_905171775.1_aRanTem1.1_genomic_with_INTRONs.gff --out intron_lengths

