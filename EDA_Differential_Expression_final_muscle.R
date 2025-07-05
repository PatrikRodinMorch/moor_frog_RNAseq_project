#Load packages
library(DESeq2)
library(BiocParallel)
library(edgeR)
library(PCAtools)
library(tidyverse)
library(ggfortify)
library(EnhancedVolcano)
library(ggplot2)
library(ggbeeswarm)
library(pheatmap)
library(patchwork)
library(glmGamPoi)
register(MulticoreParam(12))

### Testing differential expression as a response of temperature and local pond characteristics in moor frog populations.

setwd("/home/prm/Desktop/RNA_seq_stuff/Differential_expression_arvalis_genes_TE/raw_liver/")


# Load data and pre-process data
rawcounts_muscle <- read.table("TE_count_matrix_final.txt", header = T, 
                               row.names=1, com='', check.names=F)
metadata_muscle <- read.csv("muscle_metadata2.txt", header = TRUE, sep = ",", 
                            row.names = 1)

#Rename two individuals with wrong label
colnames(rawcounts_muscle)[47] <- "15_P1_F6_I3"
colnames(rawcounts_muscle)[23] <- "15_P3_F4_I3"

for (i in seq(from=1, to=14)){
  metadata_muscle[,i] <- as.factor(metadata_muscle[,i])
}



all(colnames(rawcounts_muscle) == rownames(metadata_muscle))


## Exploratory data analysis

## I filter individuals based on gosner stage. I remove all entries with a gosner stage that does not make sense i.e. 42, 0.

# Filter out gosner stages 31, 32 ,33, 34, 42. We do however have unbalanced populations (see table below). This is also true for muscle.
muscle_sub <-metadata_muscle %>%
  filter(Gosner %in% c(31, 32, 33, 34))



# Filter expression matrix based on gosner stage 31,32,33
rawcounts_muscle_sub <- rawcounts_muscle[, colnames(rawcounts_muscle) %in% rownames(muscle_sub)]
metadata_muscle_sub <- metadata_muscle[rownames(metadata_muscle) %in% rownames(muscle_sub), ]
rownames(rawcounts_muscle_sub) <- gsub('gene-', '', rownames(rawcounts_muscle_sub))



## Here we can see how many individuals there are for each population and temperature.

# Summarize number of individuals in each gosner stage per population and temperature for muscle
muscle_sub %>% 
  group_by(Pop, Temp) %>%
  summarize(n())

muscle_sub %>% 
  group_by(Temp) %>%
  summarize(n())


## Here we see how the unnormalized count distribution looks like among samples. We also see how normalization improves the distribution.

statusCol_muscle <- as.numeric(factor(metadata_muscle_sub$Temp)) + 1

(rawcounts_muscle <- boxplot(rawcounts_muscle_sub, 
                             xlab="", 
                             ylab="raw counts (muscle)",
                             las=2,
                             col=statusCol_muscle,outline=FALSE, cex.lab=1.2, main= "Raw counts (muscle)"))


muscle_sub_deseqobj <- DESeqDataSetFromMatrix(countData = round(rawcounts_muscle_sub), 
                                              colData = metadata_muscle_sub, design = ~Gosner + Temp)
muscle_sub_deseqobj_filt <- muscle_sub_deseqobj[rowSums(counts(muscle_sub_deseqobj) > 10) >= 43, ]

muscle_sf_sub <- DESeq2::estimateSizeFactors(muscle_sub_deseqobj_filt)
muscle_norm_sub <- counts(muscle_sf_sub, normalized = TRUE)



(normalized_counts_muscle <- boxplot(muscle_norm_sub, 
                                     xlab="", 
                                     ylab="Normalized counts",
                                     las=2,
                                     col=statusCol_muscle,outline=FALSE, cex.lab=1.2, main= "Normalized expression (muscle)"))


###################################
dds <- DESeq2::estimateDispersions(muscle_sf_sub)
plotDispEsts(dds, ylim=c(0.001, 6))
###################################

## Here we have PCA plots. As for grouping, I used temperature and i decided to include population 1 (fast development) against all other (slow developer), as well as temporary and permanent pond. The latter is perfectly correlated with canopy cover. Populations 1(P1) and 2(P4) both have low canopy cover 0% vs 10% and both are temporary ponds. The other three have high canopy cover and are permanent ponds.


## muscle - Temperature
## PC1-PC2
## It is still quite weird that temperature by itself does not have a stronger effect.
muscle_sub_deseqobj_filt <- muscle_sub_deseqobj[rowSums(counts(muscle_sub_deseqobj) > 10) >= 43, ]

muscle_sub_deseqobj_filt <- estimateSizeFactors(muscle_sub_deseqobj_filt)
vst_muscle <- vst(muscle_sub_deseqobj_filt, blind=F)
vst_muscle_mat <- assay(vst_muscle)
pcdat_muscle_sub <- prcomp(t(vst_muscle_mat))

PC1 <- autoplot(pcdat_muscle_sub, data = metadata_muscle_sub, colour = 'Temp', shape = 'Gosner', x = 1, y = 2, frame = TRUE, frame.type = 'norm')


## PC2-PC3
PC2 <- autoplot(pcdat_muscle_sub, data = metadata_muscle_sub, colour = 'Temp', shape = 'Gosner', x = 2, y = 3, frame = TRUE, frame.type = 'norm')


## PC3-PC4
PC3 <- autoplot(pcdat_muscle_sub, data = metadata_muscle_sub, colour = 'Temp', shape = 'Gosner', x = 3, y = 4, frame = TRUE, frame.type = 'norm')

PC1+PC2/PC3

## Differential expression - Temp
muscle_sub_deseqobj <- DESeqDataSetFromMatrix(countData = round(rawcounts_muscle_sub), 
                                              colData = metadata_muscle_sub, design = ~Gosner +Temp)
muscle_sub_deseqobj_filt <- muscle_sub_deseqobj[rowSums(counts(muscle_sub_deseqobj) > 10) >= 43, ]
muscle_temp_muscle_DE <- DESeq(muscle_sub_deseqobj_filt, parallel = T)
#DE_results_temp_muscle <- results(muscle_temp_muscle_DE, alpha = 0.1, contrast=c("Temp", "20", "15"))
DE_results_temp_muscleLFC <- lfcShrink(muscle_temp_muscle_DE, coef = "Temp_20_vs_15")
DE_results_temp_muscleLFC$padj[is.na(DE_results_temp_muscleLFC$padj)] <- 1
summary(DE_results_temp_muscleLFC)

EnhancedVolcano(DE_results_temp_muscleLFC,
                lab = rownames(DE_results_temp_muscleLFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = '20C vs 15C',
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


## Extract DEGs for GO
degs_muscle_temp_up <- DE_results_temp_muscleLFC[(DE_results_temp_muscleLFC$log2FoldChange >= 1) & (DE_results_temp_muscleLFC$padj <= 0.05),]
degs_muscle_temp_down <- DE_results_temp_muscleLFC[(DE_results_temp_muscleLFC$log2FoldChange <= -1) & (DE_results_temp_muscleLFC$padj <= 0.05),]
degs_muscle <- rbind(degs_muscle_temp_up, degs_muscle_temp_down)



write.table(rownames(degs_muscle_temp_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_temp_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_muscle_temp_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_temp_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_muscle), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_temp_DEG.txt", row.names = F, col.names = F, quote = F)


##### Plasticity ########
## Pond permanency
## Splitting pond perm and just compare temporary_15 vs temporary_20 and permanent_15 vs permanent_20.
## Permanent ponds
## PC1-PC2
muscle_sub_deseqobj_hydro <- DESeqDataSetFromMatrix(countData = round(rawcounts_muscle_sub), 
                                                    colData = metadata_muscle_sub, design = ~Gosner + Hydro_temp)


muscle_sub_deseqobj_filt_perm <- muscle_sub_deseqobj_hydro[, muscle_sub_deseqobj_hydro$Hydro_temp %in% c("Permanent_15", "Permanent_20")]
metadata_muscle_sub_perm <- metadata_muscle_sub[metadata_muscle_sub$Hydro_temp %in% c("Permanent_15", "Permanent_20"), ]
muscle_sub_deseqobj_filt_perm_filt <- muscle_sub_deseqobj_filt_perm[rowSums(counts(muscle_sub_deseqobj_filt_perm) > 10) >= 24, ]
muscle_sub_deseqobj_filt_perm_filt$Hydro_temp <- droplevels(muscle_sub_deseqobj_filt_perm_filt$Hydro_temp)


muscle_sub_deseqobj_filt_perm_filt <- estimateSizeFactors(muscle_sub_deseqobj_filt_perm_filt)
vst_muscle_perm <- vst(muscle_sub_deseqobj_filt_perm_filt, blind=F)
vst_muscle_perm_mat <- assay(vst_muscle_perm)

pcdat_muscle_perm <- prcomp(t(vst_muscle_perm_mat))

PC1 <- autoplot(pcdat_muscle_perm, data = metadata_muscle_sub_perm, colour = 'Hydro_temp', x = 1, y = 2, frame = TRUE, frame.type = 'norm')


## PC2-PC3

PC2 <- autoplot(pcdat_muscle_perm, data = metadata_muscle_sub_perm, colour = 'Hydro_temp', x = 2, y = 3, frame = TRUE, frame.type = 'norm')


## PC3-PC4
PC3 <- autoplot(pcdat_muscle_perm, data = metadata_muscle_sub_perm, colour = 'Hydro_temp', x = 3, y = 4, frame = TRUE, frame.type = 'norm')

PC1+PC2/PC3

## Differential expression - permanent ponds
muscle_hydro_temp_muscle_DE <- DESeq(muscle_sub_deseqobj_filt_perm_filt, parallel = T)
#DE_results_hydro_temp_muscle <- results(muscle_hydro_temp_muscle_DE, alpha = 0.05, contrast=c("Hydro_temp", "Permanent_15", "Permanent_20"))
DE_results_hydro_perm_muscle_LFC <- lfcShrink(muscle_hydro_temp_muscle_DE, coef = "Hydro_temp_Permanent_20_vs_Permanent_15")
DE_results_hydro_perm_muscle_LFC$padj[is.na(DE_results_hydro_perm_muscle_LFC$padj)] <- 1
summary(DE_results_hydro_perm_muscle_LFC)


EnhancedVolcano(DE_results_hydro_perm_muscle_LFC,
                lab = rownames(DE_results_hydro_perm_muscle_LFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Permannent: 20C vs 15C',
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


## Extract DEGs for GO
degs_hydro_perm_muscle_up <- DE_results_hydro_perm_muscle_LFC[(DE_results_hydro_perm_muscle_LFC$log2FoldChange >= 1) & (DE_results_hydro_perm_muscle_LFC$padj <= 0.05),]
degs_hydro_perm_muscle_down <- DE_results_hydro_perm_muscle_LFC[(DE_results_hydro_perm_muscle_LFC$log2FoldChange <= -1) & (DE_results_hydro_perm_muscle_LFC$padj <= 0.05),]
hydro_perm_muscle_DEG <- rbind(degs_hydro_perm_muscle_up, degs_hydro_perm_muscle_down)

write.table(rownames(degs_hydro_perm_muscle_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_hydro_temp_permanent_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_hydro_perm_muscle_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGOgenes_muscle_hydro_temp_permanent_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(hydro_perm_muscle_DEG), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_hydro_temp_permanent_DEG.txt", row.names = F, col.names = F, quote = F)


## Temporary ponds
## PC1-PC2
muscle_sub_deseqobj_filt_temp <- muscle_sub_deseqobj_hydro[, muscle_sub_deseqobj_hydro$Hydro_temp %in% c("Temporary_15", "Temporary_20")]
metadata_muscle_sub_temp <- metadata_muscle_sub[metadata_muscle_sub$Hydro_temp %in% c("Temporary_15", "Temporary_20"), ]
muscle_sub_deseqobj_filt_temp_filt <- muscle_sub_deseqobj_filt_temp[rowSums(counts(muscle_sub_deseqobj_filt_temp) > 10) >= 16, ]
muscle_sub_deseqobj_filt_temp_filt$Hydro_temp <- droplevels(muscle_sub_deseqobj_filt_temp_filt$Hydro_temp)


muscle_sub_deseqobj_filt_temp_filt <- estimateSizeFactors(muscle_sub_deseqobj_filt_temp_filt)
vst_muscle_temp <- vst(muscle_sub_deseqobj_filt_temp_filt, blind=F)
vst_muscle_temp_mat <- assay(vst_muscle_temp)

pcdat_muscle_temp <- prcomp(t(vst_muscle_temp_mat))

PC1 <- autoplot(pcdat_muscle_temp, data = metadata_muscle_sub_temp, colour = 'Hydro_temp', x = 1, y = 2, frame = TRUE, frame.type = 'norm')

## PC2-PC3
PC2 <- autoplot(pcdat_muscle_temp, data = metadata_muscle_sub_temp, colour = 'Hydro_temp', x = 2, y = 3, frame = TRUE, frame.type = 'norm')

## PC3-PC4
PC3 <- autoplot(pcdat_muscle_temp, data = metadata_muscle_sub_temp, colour = 'Hydro_temp', x = 3, y = 4, frame = TRUE, frame.type = 'norm')

PC1+PC2/PC3

## Differential expression
muscle_sub_deseqobj_filt_temp_filt$Hydro_temp <- droplevels(muscle_sub_deseqobj_filt_temp_filt$Hydro_temp)
muscle_hydro_temp_muscle_DE <- DESeq(muscle_sub_deseqobj_filt_temp_filt, parallel = T)
#DE_results_hydro_temp_muscle <- results(muscle_hydro_temp_muscle_DE, alpha = 0.05, contrast=c("Hydro_temp", "Temporary_15", "Temporary_20"))
DE_results_hydro_temp_muscle_LFC <- lfcShrink(muscle_hydro_temp_muscle_DE, coef = "Hydro_temp_Temporary_20_vs_Temporary_15")
DE_results_hydro_temp_muscle_LFC$padj[is.na(DE_results_hydro_temp_muscle_LFC$padj)] <- 1
summary(DE_results_hydro_temp_muscle_LFC)



EnhancedVolcano(DE_results_hydro_temp_muscle_LFC,
                lab = rownames(DE_results_hydro_temp_muscle_LFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Temporary: 20C vs 15C',
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


## Extract DEGs for GO
degs_muscle_hydro_temp_temp_up <- DE_results_hydro_temp_muscle_LFC[(DE_results_hydro_temp_muscle_LFC$log2FoldChange >= 1) & (DE_results_hydro_temp_muscle_LFC$padj <= 0.05),]
degs_muscle_hydro_temp_temp_down <- DE_results_hydro_temp_muscle_LFC[(DE_results_hydro_temp_muscle_LFC$log2FoldChange <= -1) & (DE_results_hydro_temp_muscle_LFC$padj <= 0.05),]
degs_muscle_hydro_temp_temporary_DEG <- rbind(degs_muscle_hydro_temp_temp_up, degs_muscle_hydro_temp_temp_down)

write.table(rownames(degs_muscle_hydro_temp_temp_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_hydro_temp_temporary_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_muscle_hydro_temp_temp_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_hydro_temp_temporary_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_muscle_hydro_temp_temporary_DEG), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_hydro_temp_temporary_DEG.txt", row.names = F, col.names = F, quote = F)


##### Divergence ########
## Pond permanency
## Splitting pond perm and just compare temporary_15 vs temporary_20 and permanent_15 vs permanent_20.
## Permanent ponds
## PC1-PC2
muscle_sub_deseqobj_hydro <- DESeqDataSetFromMatrix(countData = round(rawcounts_muscle_sub), 
                                                    colData = metadata_muscle_sub, design = ~Gosner + Hydro_temp)


muscle_sub_deseqobj_filt_15 <- muscle_sub_deseqobj_hydro[, muscle_sub_deseqobj_hydro$Hydro_temp %in% c("Temporary_15", "Permanent_15")]
metadata_muscle_sub_15 <- metadata_muscle_sub[metadata_muscle_sub$Hydro_temp %in% c("Temporary_15", "Permanent_15"), ]
muscle_sub_deseqobj_filt_15_filt <- muscle_sub_deseqobj_filt_15[rowSums(counts(muscle_sub_deseqobj_filt_15) > 10) >= 24, ]
muscle_sub_deseqobj_filt_15$Hydro_temp <- droplevels(muscle_sub_deseqobj_filt_15$Hydro_temp)


muscle_sub_deseqobj_filt_15 <- estimateSizeFactors(muscle_sub_deseqobj_filt_15)
vst_muscle_15 <- vst(muscle_sub_deseqobj_filt_15, blind=F)
vst_muscle_15_mat <- assay(vst_muscle_15)

pcdat_muscle_15 <- prcomp(t(vst_muscle_15_mat))

PC1 <- autoplot(pcdat_muscle_15, data = metadata_muscle_sub_15, colour = 'Hydro_temp', x = 1, y = 2, frame = TRUE, frame.type = 'norm')


## PC2-PC3

PC2 <- autoplot(pcdat_muscle_15, data = metadata_muscle_sub_15, colour = 'Hydro_temp', x = 2, y = 3, frame = TRUE, frame.type = 'norm')


## PC3-PC4
PC3 <- autoplot(pcdat_muscle_15, data = metadata_muscle_sub_15, colour = 'Hydro_temp', x = 3, y = 4, frame = TRUE, frame.type = 'norm')

PC1+PC2/PC3

## Differential expression
muscle_hydro_15_muscle_DE <- DESeq(muscle_sub_deseqobj_filt_15, parallel = T)
#DE_results_hydro_temp_muscle <- results(muscle_hydro_temp_muscle_DE, alpha = 0.05, contrast=c("Hydro_temp", "Permanent_15", "Permanent_20"))
DE_results_hydro_15_muscle_LFC <- lfcShrink(muscle_hydro_15_muscle_DE, coef = "Hydro_temp_Temporary_15_vs_Permanent_15")
DE_results_hydro_15_muscle_LFC$padj[is.na(DE_results_hydro_15_muscle_LFC$padj)] <- 1
summary(DE_results_hydro_15_muscle_LFC)


EnhancedVolcano(DE_results_hydro_15_muscle_LFC,
                lab = rownames(DE_results_hydro_15_muscle_LFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Temporary 15C vs permanent 15C',
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


## Extract DEGs for GO
degs_hydro_15_muscle_up <- DE_results_hydro_15_muscle_LFC[(DE_results_hydro_15_muscle_LFC$log2FoldChange >= 1) & (DE_results_hydro_15_muscle_LFC$padj <= 0.05),]
degs_hydro_15_muscle_down <- DE_results_hydro_15_muscle_LFC[(DE_results_hydro_15_muscle_LFC$log2FoldChange <= -1) & (DE_results_hydro_15_muscle_LFC$padj <= 0.05),]
hydro_15_muscle_DEG <- rbind(degs_hydro_15_muscle_up, degs_hydro_15_muscle_down)

write.table(rownames(degs_hydro_15_muscle_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_hydro_temp_15_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_hydro_15_muscle_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGOgenes_muscle_hydro_temp_15_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(hydro_15_muscle_DEG), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_hydro_temp_15_DEG.txt", row.names = F, col.names = F, quote = F)


## 20
muscle_sub_deseqobj_hydro <- DESeqDataSetFromMatrix(countData = round(rawcounts_muscle_sub), 
                                                    colData = metadata_muscle_sub, design = ~Hydro_temp)


muscle_sub_deseqobj_filt_20 <- muscle_sub_deseqobj_hydro[, muscle_sub_deseqobj_hydro$Hydro_temp %in% c("Temporary_20", "Permanent_20")]
metadata_muscle_sub_20 <- metadata_muscle_sub[metadata_muscle_sub$Hydro_temp %in% c("Temporary_20", "Permanent_20"), ]
muscle_sub_deseqobj_filt_20_filt <- muscle_sub_deseqobj_filt_20[rowSums(counts(muscle_sub_deseqobj_filt_20) > 10) >= 24, ]
muscle_sub_deseqobj_filt_20$Hydro_temp <- droplevels(muscle_sub_deseqobj_filt_20$Hydro_temp)


muscle_sub_deseqobj_filt_20 <- estimateSizeFactors(muscle_sub_deseqobj_filt_20)
vst_muscle_20 <- vst(muscle_sub_deseqobj_filt_20, blind=F)
vst_muscle_20_mat <- assay(vst_muscle_20)

pcdat_muscle_20 <- prcomp(t(vst_muscle_20_mat))

PC1 <- autoplot(pcdat_muscle_20, data = metadata_muscle_sub_20, colour = 'Hydro_temp', x = 1, y = 2, frame = TRUE, frame.type = 'norm')


## PC2-PC3

PC2 <- autoplot(pcdat_muscle_20, data = metadata_muscle_sub_20, colour = 'Hydro_temp', x = 2, y = 3, frame = TRUE, frame.type = 'norm')


## PC3-PC4
PC3 <- autoplot(pcdat_muscle_20, data = metadata_muscle_sub_20, colour = 'Hydro_temp', x = 3, y = 4, frame = TRUE, frame.type = 'norm')

PC1+PC2/PC3

## Differential expression
muscle_hydro_20_muscle_DE <- DESeq(muscle_sub_deseqobj_filt_20, parallel = T)
#DE_results_hydro_temp_muscle <- results(muscle_hydro_temp_muscle_DE, alpha = 0.05, contrast=c("Hydro_temp", "Permanent_20", "Permanent_20"))
DE_results_hydro_20_muscle_LFC <- lfcShrink(muscle_hydro_20_muscle_DE, coef = "Hydro_temp_Temporary_20_vs_Permanent_20")
DE_results_hydro_20_muscle_LFC$padj[is.na(DE_results_hydro_20_muscle_LFC$padj)] <- 1
summary(DE_results_hydro_20_muscle_LFC)


EnhancedVolcano(DE_results_hydro_20_muscle_LFC,
                lab = rownames(DE_results_hydro_20_muscle_LFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Temporary 20C vs permanent 20C',
                FCcutoff = 1,
                pCutoff = 0.05,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


## Extract DEGs for GO
degs_hydro_20_muscle_up <- DE_results_hydro_20_muscle_LFC[(DE_results_hydro_20_muscle_LFC$log2FoldChange >= 1) & (DE_results_hydro_20_muscle_LFC$padj <= 0.05),]
degs_hydro_20_muscle_down <- DE_results_hydro_20_muscle_LFC[(DE_results_hydro_20_muscle_LFC$log2FoldChange <= -1) & (DE_results_hydro_20_muscle_LFC$padj <= 0.05),]
hydro_20_muscle_DEG <- rbind(degs_hydro_20_muscle_up, degs_hydro_20_muscle_down)

write.table(rownames(degs_hydro_20_muscle_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_hydro_temp_20_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_hydro_20_muscle_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGOgenes_muscle_hydro_temp_20_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(hydro_20_muscle_DEG), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/muscle/ForGO/genes_muscle_hydro_temp_20_DEG.txt", row.names = F, col.names = F, quote = F)
