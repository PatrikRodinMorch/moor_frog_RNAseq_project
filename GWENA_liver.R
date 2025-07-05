library("GWENA")
library("DESeq2")
library("tidyverse")

### Exploring gene co-expression networks.

setwd("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/Differential_expression/expr_matrices/DESeq2")


threads_to_use <- 8

rawcounts_liver <- read.table("liver/gene_counts_htseq_liver_analysis.txt", header = T, 
                              row.names=1, com='', check.names=F)
metadata_liver <- read.csv("liver/liver_metadata2.txt", header = TRUE, sep = ",", 
                           row.names = 1) 


all(colnames(rawcounts_liver) == rownames(metadata_liver))

liver_sub <-metadata_liver %>%
  filter(Gosner %in% c(31, 32, 33, 34))


# Filter expression matrix based on gosner stage 31,32,33
rawcounts_liver_sub <- rawcounts_liver[, colnames(rawcounts_liver) %in% rownames(liver_sub)]
metadata_liver_sub <- metadata_liver[rownames(metadata_liver) %in% rownames(liver_sub), ]
rownames(rawcounts_liver_sub) <- gsub('gene-', '', rownames(rawcounts_liver_sub))


###### Normalize expression counts
liver_sub_deseqobj <- DESeqDataSetFromMatrix(countData = round(rawcounts_liver_sub), 
                                             colData = metadata_liver_sub, design = ~1)
liver_sf_sub <- DESeq2::estimateSizeFactors(liver_sub_deseqobj)

###1
liver_norm_sub <- counts(liver_sf_sub, normalized = TRUE)
liver_norm_sub_t <- t(liver_norm_sub)
ncol(liver_norm_sub_t)


#### GWENA

### Gene filtering
filtered_liver <- filter_RNA_seq(liver_norm_sub_t, min_count = 5, method = "at least one")
ncol(filtered_liver)
ncol(na.omit(filtered_liver))
##### Network building
net <- build_net(filtered_liver, cor_func = "spearman", 
                 n_threads = threads_to_use, fit_cut_off = 0.94)
net$metadata$power
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]

#### Modules detection
modules <- detect_modules(filtered_liver, 
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
#plot_expression_profiles(filtered_liver, modules$modules)


#### Phenotypic association
phenotype_association <- associate_phenotype(modules$modules_eigengenes, 
                                             metadata_liver_sub %>% dplyr::select(Temp, LP, MAM))
plot_modules_phenotype(phenotype_association)


for (n in 1:length(unique(modules$modules))){
  module_number <- n-1
  write.table(modules$modules[[n]], file = paste0("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/Differential_expression/expr_matrices/GWENA/liver/module_",module_number,"_genes.txt"), row.names = F, col.names = F, quote = F)
  
}


