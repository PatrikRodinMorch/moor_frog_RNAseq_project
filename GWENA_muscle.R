library("GWENA")
library("DESeq2")
library("tidyverse")


### Exploring gene co-expression networks.

setwd("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/Differential_expression/expr_matrices/DESeq2")


threads_to_use <- 8

rawcounts_muscle <- read.table("muscle/gene_counts_htseq_muscle_analysis.txt", header = T, 
                               row.names=1, com='', check.names=F)
metadata_muscle <- read.csv("muscle/muscle_metadata2.txt", header = TRUE, sep = ",", 
                            row.names = 1) 
colnames(rawcounts_muscle)[47] <- "15_P1_F6_I3"
colnames(rawcounts_muscle)[23] <- "15_P3_F4_I3"


all(colnames(rawcounts_muscle) == rownames(metadata_muscle))

muscle_sub <-metadata_muscle %>%
  filter(Gosner %in% c(31, 32, 33, 34))


# Filter expression matrix based on gosner stage 31,32,33
rawcounts_muscle_sub <- rawcounts_muscle[, colnames(rawcounts_muscle) %in% rownames(muscle_sub)]
metadata_muscle_sub <- metadata_muscle[rownames(metadata_muscle) %in% rownames(muscle_sub), ]
rownames(rawcounts_muscle_sub) <- gsub('gene-', '', rownames(rawcounts_muscle_sub))


###### Normalize expression counts
muscle_sub_deseqobj <- DESeqDataSetFromMatrix(countData = round(rawcounts_muscle_sub), 
                                              colData = metadata_muscle_sub, design = ~ 1)
muscle_sf_sub <- DESeq2::estimateSizeFactors(muscle_sub_deseqobj)

###1
muscle_norm_sub <- counts(muscle_sf_sub, normalized = TRUE)
muscle_norm_sub_t <- t(muscle_norm_sub)
ncol(muscle_norm_sub_t)


#### GWENA

### Gene filtering
filtered_muscle <- filter_RNA_seq(muscle_norm_sub_t, min_count = 5, method = "at least one")
ncol(filtered_muscle)
ncol(na.omit(filtered_muscle))
##### Network building
net <- build_net(filtered_muscle, cor_func = "spearman", 
                 n_threads = threads_to_use, fit_cut_off = 0.83)
net$metadata$power
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]

#### Modules detection
modules <- detect_modules(filtered_muscle, 
                          net$network, 
                          detailled_result = TRUE,
                          merge_threshold = 0.25)

# Number of modules before merging :
length(unique(modules$modules_premerge))

# Number of modules after merging: 
length(unique(modules$modules))

#### Plot modules
layout_mod_merge <- plot_modules_merge(
  modules_premerge = modules$modules_premerge, 
  modules_merged = modules$modules)

###Number of genes
ggplot2::ggplot(data.frame(modules$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")

####
#plot_expression_profiles(filtered_muscle, modules$modules)


#### Phenotypic association
phenotype_association <- associate_phenotype(modules$modules_eigengenes, 
                                             metadata_muscle_sub %>% dplyr::select(Temp, LP, MAM))
plot_modules_phenotype(phenotype_association)


for (n in 1:length(unique(modules$modules))){
  module_number <- n-1
  write.table(modules$modules[[n]], file = paste0("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/Differential_expression/expr_matrices/GWENA/muscle/module_",module_number,"_genes.txt"), row.names = F, col.names = F, quote = F)
  
}



