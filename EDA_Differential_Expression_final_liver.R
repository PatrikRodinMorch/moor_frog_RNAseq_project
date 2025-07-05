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
library(apeglm)
register(MulticoreParam(12))


### Testing differential expression as a response of temperature and local pond characteristics in moor frog populations.

setwd("/home/prm/Desktop/RNA_seq_stuff/Differential_expression_arvalis_genes_TE/raw_liver/")


# Load data and pre-process data
rawcounts_liver <- read.table("TE_count_matrix_final_filtered.txt", header = T, 
                              row.names=1, com='', check.names=F)
metadata_liver <- read.csv("liver_metadata2.txt", header = TRUE, sep = ",", 
                           row.names = 1)


for (i in seq(from=1, to=14)){
  metadata_liver[,i] <- as.factor(metadata_liver[,i])
}



all(colnames(rawcounts_liver) == rownames(metadata_liver))


## Exploratory data analysis

## I filter individuals based on gosner stage. I remove all entries with a gosner stage that does not make sense i.e. 42, 0.

# Filter out gosner stages 31, 32 ,33, 34, 42. We do however have unbalanced populations (see table below). This is also true for muscle.
liver_sub <-metadata_liver %>%
  filter(!(rownames(metadata_liver) %in% c("20_P1_F9_I1")))



# Filter expression matrix based on gosner stage 31,32,33
rawcounts_liver_sub <- rawcounts_liver[, colnames(rawcounts_liver) %in% rownames(liver_sub)]
metadata_liver_sub <- metadata_liver[rownames(metadata_liver) %in% rownames(liver_sub), ]
rownames(rawcounts_liver_sub) <- gsub('gene-', '', rownames(rawcounts_liver_sub))



## Here we can see how many individuals there are for each population and temperature.

# Summarize number of individuals in each gosner stage per population and temperature for liver
liver_sub %>% 
  group_by(Pop, Temp) %>%
  summarize(n())

liver_sub %>% 
  group_by(Temp) %>%
  summarize(n())


## Here we see how the unnormalized count distribution looks like among samples. We also see how normalization improves the distribution.

statusCol_liver <- as.numeric(factor(metadata_liver_sub$Temp)) + 1

(rawcounts_liver <- boxplot(rawcounts_liver_sub, 
        xlab="", 
        ylab="raw counts (liver)",
        las=2,
        col=statusCol_liver,outline=FALSE, cex.lab=1.2, main= "Raw counts (Liver)"))


liver_sub_deseqobj <- DESeqDataSetFromMatrix(countData = round(rawcounts_liver_sub), 
                                                  colData = metadata_liver_sub, design = ~Gosner + Temp)
liver_sub_deseqobj_filt <- liver_sub_deseqobj[rowSums(counts(liver_sub_deseqobj) > 5) >= 46, ]

liver_sf_sub <- DESeq2::estimateSizeFactors(liver_sub_deseqobj_filt)
liver_norm_sub <- counts(liver_sf_sub, normalized = TRUE)



(normalized_counts_liver <- boxplot(liver_norm_sub, 
        xlab="", 
        ylab="Normalized counts",
        las=2,
        col=statusCol_liver,outline=FALSE, cex.lab=1.2, main= "Normalized expression (Liver)"))


###################################
dds <- DESeq2::estimateDispersions(liver_sf_sub)
plotDispEsts(dds, ylim=c(0.001, 6))
###################################

## Here we have PCA plots. As for grouping, I used temperature and i decided to include population 1 (fast development) against all other (slow developer), as well as temporary and permanent pond. The latter is perfectly correlated with canopy cover. Populations 1(P1) and 2(P4) both have low canopy cover 0% vs 10% and both are temporary ponds. The other three have high canopy cover and are permanent ponds.


## Liver - Temperature
## PC1-PC2
## It is still quite weird that temperature by itself does not have a stronger effect.
liver_sub_deseqobj_filt <- liver_sub_deseqobj[rowSums(counts(liver_sub_deseqobj) > 5) >= 46, ]

liver_sub_deseqobj_filt <- estimateSizeFactors(liver_sub_deseqobj_filt)
vst_liver <- vst(liver_sub_deseqobj_filt, blind=F)
vst_liver_mat <- assay(vst_liver)
pcdat_liver_sub <- prcomp(t(vst_liver_mat))

PC1 <- autoplot(pcdat_liver_sub, data = metadata_liver_sub, colour = 'Temp', shape = 'Gosner', x = 1, y = 2, frame = TRUE, frame.type = 'norm')
PC1

## PC2-PC3
PC2 <- autoplot(pcdat_liver_sub, data = metadata_liver_sub, colour = 'Temp', shape = 'Gosner', x = 2, y = 3, frame = TRUE, frame.type = 'norm')
PC2

## PC3-PC4
PC3 <- autoplot(pcdat_liver_sub, data = metadata_liver_sub, colour = 'Temp', shape = 'Gosner', x = 3, y = 4, frame = TRUE, frame.type = 'norm')
PC3

PC1+PC2/PC3

## Differential expression - Temp
liver_sub_deseqobj <- DESeqDataSetFromMatrix(countData = round(rawcounts_liver_sub), 
                                                  colData = metadata_liver_sub, design = ~Gosner +Temp)
liver_sub_deseqobj_filt <- liver_sub_deseqobj[rowSums(counts(liver_sub_deseqobj) > 5) >= 46, ]
liver_temp_liver_DE <- DESeq(liver_sub_deseqobj_filt, parallel = T)
#DE_results_temp_liver <- results(liver_temp_liver_DE, alpha = 0.1, contrast=c("Temp", "20", "15"))
DE_results_temp_liverLFC <- lfcShrink(liver_temp_liver_DE, coef = "Temp_20_vs_15")
DE_results_temp_liverLFC$padj[is.na(DE_results_temp_liverLFC$padj)] <- 1
summary(DE_results_temp_liverLFC)

EnhancedVolcano(DE_results_temp_liverLFC,
    lab = rownames(DE_results_temp_liverLFC),
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
degs_liver_temp_up <- DE_results_temp_liverLFC[(DE_results_temp_liverLFC$log2FoldChange >= 1) & (DE_results_temp_liverLFC$padj <= 0.05),]
degs_liver_temp_down <- DE_results_temp_liverLFC[(DE_results_temp_liverLFC$log2FoldChange <= -1) & (DE_results_temp_liverLFC$padj <= 0.05),]
degs_liver <- rbind(degs_liver_temp_up, degs_liver_temp_down)



write.table(rownames(degs_liver_temp_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_temp_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_liver_temp_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_temp_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_liver), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_temp_DEG.txt", row.names = F, col.names = F, quote = F)


##### Plasticity ########
## Pond permanency
## Splitting pond perm and just compare temporary_15 vs temporary_20 and permanent_15 vs permanent_20.
## Permanent ponds
## PC1-PC2
liver_sub_deseqobj_hydro <- DESeqDataSetFromMatrix(countData = round(rawcounts_liver_sub), 
                                                  colData = metadata_liver_sub, design = ~Gosner + Hydro_temp)


liver_sub_deseqobj_filt_perm <- liver_sub_deseqobj_hydro[, liver_sub_deseqobj_hydro$Hydro_temp %in% c("Permanent_15", "Permanent_20")]
metadata_liver_sub_perm <- metadata_liver_sub[metadata_liver_sub$Hydro_temp %in% c("Permanent_15", "Permanent_20"), ]
liver_sub_deseqobj_filt_perm_filt <- liver_sub_deseqobj_filt_perm[rowSums(counts(liver_sub_deseqobj_filt_perm) > 10) >= 24, ]
liver_sub_deseqobj_filt_perm_filt$Hydro_temp <- droplevels(liver_sub_deseqobj_filt_perm_filt$Hydro_temp)


liver_sub_deseqobj_filt_perm_filt <- estimateSizeFactors(liver_sub_deseqobj_filt_perm_filt)
vst_liver_perm <- vst(liver_sub_deseqobj_filt_perm_filt, blind=F)
vst_liver_perm_mat <- assay(vst_liver_perm)

pcdat_liver_perm <- prcomp(t(vst_liver_perm_mat))

PC1 <- autoplot(pcdat_liver_perm, data = metadata_liver_sub_perm, colour = 'Hydro_temp', x = 1, y = 2, frame = TRUE, frame.type = 'norm')


## PC2-PC3

PC2 <- autoplot(pcdat_liver_perm, data = metadata_liver_sub_perm, colour = 'Hydro_temp', x = 2, y = 3, frame = TRUE, frame.type = 'norm')


## PC3-PC4
PC3 <- autoplot(pcdat_liver_perm, data = metadata_liver_sub_perm, colour = 'Hydro_temp', x = 3, y = 4, frame = TRUE, frame.type = 'norm')

PC1+PC2/PC3

## Differential expression - permanent ponds
liver_hydro_temp_liver_DE <- DESeq(liver_sub_deseqobj_filt_perm_filt, parallel = T)
#DE_results_hydro_temp_liver <- results(liver_hydro_temp_liver_DE, alpha = 0.05, contrast=c("Hydro_temp", "Permanent_15", "Permanent_20"))
DE_results_hydro_perm_liver_LFC <- lfcShrink(liver_hydro_temp_liver_DE, coef = "Hydro_temp_Permanent_20_vs_Permanent_15")
DE_results_hydro_perm_liver_LFC$padj[is.na(DE_results_hydro_perm_liver_LFC$padj)] <- 1
summary(DE_results_hydro_perm_liver_LFC)


EnhancedVolcano(DE_results_hydro_perm_liver_LFC,
    lab = rownames(DE_results_hydro_perm_liver_LFC),
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
degs_hydro_perm_liver_up <- DE_results_hydro_perm_liver_LFC[(DE_results_hydro_perm_liver_LFC$log2FoldChange >= 1) & (DE_results_hydro_perm_liver_LFC$padj <= 0.05),]
degs_hydro_perm_liver_down <- DE_results_hydro_perm_liver_LFC[(DE_results_hydro_perm_liver_LFC$log2FoldChange <= -1) & (DE_results_hydro_perm_liver_LFC$padj <= 0.05),]
hydro_perm_liver_DEG <- rbind(degs_hydro_perm_liver_up, degs_hydro_perm_liver_down)

write.table(rownames(degs_hydro_perm_liver_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_hydro_temp_permanent_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_hydro_perm_liver_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGOgenes_liver_hydro_temp_permanent_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(hydro_perm_liver_DEG), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_hydro_temp_permanent_DEG.txt", row.names = F, col.names = F, quote = F)


## Temporary ponds
## PC1-PC2
liver_sub_deseqobj_filt_temp <- liver_sub_deseqobj_hydro[, liver_sub_deseqobj_hydro$Hydro_temp %in% c("Temporary_15", "Temporary_20")]
metadata_liver_sub_temp <- metadata_liver_sub[metadata_liver_sub$Hydro_temp %in% c("Temporary_15", "Temporary_20"), ]
liver_sub_deseqobj_filt_temp_filt <- liver_sub_deseqobj_filt_temp[rowSums(counts(liver_sub_deseqobj_filt_temp) > 10) >= 16, ]
liver_sub_deseqobj_filt_temp_filt$Hydro_temp <- droplevels(liver_sub_deseqobj_filt_temp_filt$Hydro_temp)


liver_sub_deseqobj_filt_temp_filt <- estimateSizeFactors(liver_sub_deseqobj_filt_temp_filt)
vst_liver_temp <- vst(liver_sub_deseqobj_filt_temp_filt, blind=F)
vst_liver_temp_mat <- assay(vst_liver_temp)

pcdat_liver_temp <- prcomp(t(vst_liver_temp_mat))

PC1 <- autoplot(pcdat_liver_temp, data = metadata_liver_sub_temp, colour = 'Hydro_temp', x = 1, y = 2, frame = TRUE, frame.type = 'norm')

## PC2-PC3
PC2 <- autoplot(pcdat_liver_temp, data = metadata_liver_sub_temp, colour = 'Hydro_temp', x = 2, y = 3, frame = TRUE, frame.type = 'norm')

## PC3-PC4
PC3 <- autoplot(pcdat_liver_temp, data = metadata_liver_sub_temp, colour = 'Hydro_temp', x = 3, y = 4, frame = TRUE, frame.type = 'norm')

PC1+PC2/PC3

## Differential expression
liver_sub_deseqobj_filt_temp_filt$Hydro_temp <- droplevels(liver_sub_deseqobj_filt_temp_filt$Hydro_temp)
liver_hydro_temp_liver_DE <- DESeq(liver_sub_deseqobj_filt_temp_filt, parallel = T)
#DE_results_hydro_temp_liver <- results(liver_hydro_temp_liver_DE, alpha = 0.05, contrast=c("Hydro_temp", "Temporary_15", "Temporary_20"))
DE_results_hydro_temp_liver_LFC <- lfcShrink(liver_hydro_temp_liver_DE, coef = "Hydro_temp_Temporary_20_vs_Temporary_15")
DE_results_hydro_temp_liver_LFC$padj[is.na(DE_results_hydro_temp_liver_LFC$padj)] <- 1
summary(DE_results_hydro_temp_liver_LFC)



  EnhancedVolcano(DE_results_hydro_temp_liver_LFC,
    lab = rownames(DE_results_hydro_temp_liver_LFC),
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
degs_liver_hydro_temp_temp_up <- DE_results_hydro_temp_liver_LFC[(DE_results_hydro_temp_liver_LFC$log2FoldChange >= 1) & (DE_results_hydro_temp_liver_LFC$padj <= 0.05),]
degs_liver_hydro_temp_temp_down <- DE_results_hydro_temp_liver_LFC[(DE_results_hydro_temp_liver_LFC$log2FoldChange <= -1) & (DE_results_hydro_temp_liver_LFC$padj <= 0.05),]
degs_liver_hydro_temp_temporary_DEG <- rbind(degs_liver_hydro_temp_temp_up, degs_liver_hydro_temp_temp_down)

write.table(rownames(degs_liver_hydro_temp_temp_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_hydro_temp_temporary_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_liver_hydro_temp_temp_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_hydro_temp_temporary_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_liver_hydro_temp_temporary_DEG), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_hydro_temp_temporary_DEG.txt", row.names = F, col.names = F, quote = F)


##### Divergence ########
## Pond permanency
## Splitting pond perm and just compare temporary_15 vs temporary_20 and permanent_15 vs permanent_20.
## Permanent ponds
## PC1-PC2
liver_sub_deseqobj_hydro <- DESeqDataSetFromMatrix(countData = round(rawcounts_liver_sub), 
                                                   colData = metadata_liver_sub, design = ~Gosner + Hydro_temp)


liver_sub_deseqobj_filt_15 <- liver_sub_deseqobj_hydro[, liver_sub_deseqobj_hydro$Hydro_temp %in% c("Temporary_15", "Permanent_15")]
metadata_liver_sub_15 <- metadata_liver_sub[metadata_liver_sub$Hydro_temp %in% c("Temporary_15", "Permanent_15"), ]
liver_sub_deseqobj_filt_15_filt <- liver_sub_deseqobj_filt_15[rowSums(counts(liver_sub_deseqobj_filt_15) > 10) >= 24, ]
liver_sub_deseqobj_filt_15$Hydro_temp <- droplevels(liver_sub_deseqobj_filt_15$Hydro_temp)


liver_sub_deseqobj_filt_15 <- estimateSizeFactors(liver_sub_deseqobj_filt_15)
vst_liver_15 <- vst(liver_sub_deseqobj_filt_15, blind=F)
vst_liver_15_mat <- assay(vst_liver_15)

pcdat_liver_15 <- prcomp(t(vst_liver_15_mat))

PC1 <- autoplot(pcdat_liver_15, data = metadata_liver_sub_15, colour = 'Hydro_temp', x = 1, y = 2, frame = TRUE, frame.type = 'norm')


## PC2-PC3

PC2 <- autoplot(pcdat_liver_15, data = metadata_liver_sub_15, colour = 'Hydro_temp', x = 2, y = 3, frame = TRUE, frame.type = 'norm')


## PC3-PC4
PC3 <- autoplot(pcdat_liver_15, data = metadata_liver_sub_15, colour = 'Hydro_temp', x = 3, y = 4, frame = TRUE, frame.type = 'norm')

PC1+PC2/PC3

## Differential expression
liver_hydro_15_liver_DE <- DESeq(liver_sub_deseqobj_filt_15, parallel = T)
#DE_results_hydro_temp_liver <- results(liver_hydro_temp_liver_DE, alpha = 0.05, contrast=c("Hydro_temp", "Permanent_15", "Permanent_20"))
DE_results_hydro_15_liver_LFC <- lfcShrink(liver_hydro_15_liver_DE, coef = "Hydro_temp_Temporary_15_vs_Permanent_15")
DE_results_hydro_15_liver_LFC$padj[is.na(DE_results_hydro_15_liver_LFC$padj)] <- 1
summary(DE_results_hydro_15_liver_LFC)


EnhancedVolcano(DE_results_hydro_15_liver_LFC,
                lab = rownames(DE_results_hydro_15_liver_LFC),
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
degs_hydro_15_liver_up <- DE_results_hydro_15_liver_LFC[(DE_results_hydro_15_liver_LFC$log2FoldChange >= 1) & (DE_results_hydro_15_liver_LFC$padj <= 0.05),]
degs_hydro_15_liver_down <- DE_results_hydro_15_liver_LFC[(DE_results_hydro_15_liver_LFC$log2FoldChange <= -1) & (DE_results_hydro_15_liver_LFC$padj <= 0.05),]
hydro_15_liver_DEG <- rbind(degs_hydro_15_liver_up, degs_hydro_15_liver_down)

write.table(rownames(degs_hydro_15_liver_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_hydro_temp_15_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_hydro_15_liver_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGOgenes_liver_hydro_temp_15_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(hydro_15_liver_DEG), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_hydro_temp_15_DEG.txt", row.names = F, col.names = F, quote = F)


## 20
liver_sub_deseqobj_hydro <- DESeqDataSetFromMatrix(countData = round(rawcounts_liver_sub), 
                                                   colData = metadata_liver_sub, design = ~Hydro_temp)


liver_sub_deseqobj_filt_20 <- liver_sub_deseqobj_hydro[, liver_sub_deseqobj_hydro$Hydro_temp %in% c("Temporary_20", "Permanent_20")]
metadata_liver_sub_20 <- metadata_liver_sub[metadata_liver_sub$Hydro_temp %in% c("Temporary_20", "Permanent_20"), ]
liver_sub_deseqobj_filt_20_filt <- liver_sub_deseqobj_filt_20[rowSums(counts(liver_sub_deseqobj_filt_20) > 10) >= 24, ]
liver_sub_deseqobj_filt_20$Hydro_temp <- droplevels(liver_sub_deseqobj_filt_20$Hydro_temp)


liver_sub_deseqobj_filt_20 <- estimateSizeFactors(liver_sub_deseqobj_filt_20)
vst_liver_20 <- vst(liver_sub_deseqobj_filt_20, blind=F)
vst_liver_20_mat <- assay(vst_liver_20)

pcdat_liver_20 <- prcomp(t(vst_liver_20_mat))

PC1 <- autoplot(pcdat_liver_20, data = metadata_liver_sub_20, colour = 'Hydro_temp', x = 1, y = 2, frame = TRUE, frame.type = 'norm')


## PC2-PC3

PC2 <- autoplot(pcdat_liver_20, data = metadata_liver_sub_20, colour = 'Hydro_temp', x = 2, y = 3, frame = TRUE, frame.type = 'norm')


## PC3-PC4
PC3 <- autoplot(pcdat_liver_20, data = metadata_liver_sub_20, colour = 'Hydro_temp', x = 3, y = 4, frame = TRUE, frame.type = 'norm')

PC1+PC2/PC3

## Differential expression
liver_hydro_20_liver_DE <- DESeq(liver_sub_deseqobj_filt_20, parallel = T)
#DE_results_hydro_temp_liver <- results(liver_hydro_temp_liver_DE, alpha = 0.05, contrast=c("Hydro_temp", "Permanent_20", "Permanent_20"))
DE_results_hydro_20_liver_LFC <- lfcShrink(liver_hydro_20_liver_DE, coef = "Hydro_temp_Temporary_20_vs_Permanent_20")
DE_results_hydro_20_liver_LFC$padj[is.na(DE_results_hydro_20_liver_LFC$padj)] <- 1
summary(DE_results_hydro_20_liver_LFC)


EnhancedVolcano(DE_results_hydro_20_liver_LFC,
                lab = rownames(DE_results_hydro_20_liver_LFC),
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
degs_hydro_20_liver_up <- DE_results_hydro_20_liver_LFC[(DE_results_hydro_20_liver_LFC$log2FoldChange >= 1) & (DE_results_hydro_20_liver_LFC$padj <= 0.05),]
degs_hydro_20_liver_down <- DE_results_hydro_20_liver_LFC[(DE_results_hydro_20_liver_LFC$log2FoldChange <= -1) & (DE_results_hydro_20_liver_LFC$padj <= 0.05),]
hydro_20_liver_DEG <- rbind(degs_hydro_20_liver_up, degs_hydro_20_liver_down)

write.table(rownames(degs_hydro_20_liver_up), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_hydro_temp_20_up.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(degs_hydro_20_liver_down), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGOgenes_liver_hydro_temp_20_down.txt", row.names = F, col.names = F, quote = F)
write.table(rownames(hydro_20_liver_DEG), file = "/home/prm/Skrivbord/RNAseq_ref_mapping/expr_matrices/DESeq2/Liver/ForGO/genes_liver_hydro_temp_20_DEG.txt", row.names = F, col.names = F, quote = F)
