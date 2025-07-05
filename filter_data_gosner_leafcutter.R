library(tidyverse)


#### Data wrangling to prepare for differential splicing analysis.

setwd("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/leafcutter/")



##### Liver
liver <- read.table(gzfile("./liver/intron_clustering_perind_numers.counts.gz"))
colnames(liver) <- gsub("_INTRON_Aligned.out.bam", "", colnames(liver))

metadata_liver <- read.csv("./liver_metadata2.txt", header = TRUE, sep = ",",row.names = 1)
rownames(metadata_liver) <- gsub("-", ".", rownames(metadata_liver))
rownames(metadata_liver) <- gsub("_noquotes.txt", "", rownames(metadata_liver))


liver_sub <-metadata_liver %>%
  filter(Gosner %in% c(31, 32, 33, 34))


liver_filtered <- liver[,colnames(liver) %in% rownames(liver_sub)]

liver_sub_match <- liver_sub[match(colnames(liver_filtered), rownames(liver_sub)), ]

all(colnames(liver_filtered) == rownames(liver_sub_match))


write.table(liver_filtered[,-13], file = "./liver/intron_clustering_perind_numers.counts.FINAL", row.names = T, col.names = T, quote = F)
write.table(liver_sub_match[-13,], file = "./liver/liver_metadata_final.txt", row.names = T, col.names = F, quote = F)
system("gzip -f ./liver/intron_clustering_perind_numers.counts.FINAL")


###### Population separately
### Balance the design

population <- c("P1", "P2", "P3", "P4", "P5")

liver_sub_match2 <- liver_sub_match[-13,]
liver_filtered2 <- liver_filtered[,-13]

all(colnames(liver_filtered2) == rownames(liver_sub_match2))

for (p in population) {

  muscle_sub_p <- muscle_sub_match2[muscle_sub_match2$Pop %in% c(p), ]
  counts_p <- muscle_filtered2[, colnames(muscle_filtered2) %in% rownames(muscle_sub_p)]

  l1 <- length(muscle_sub_p[muscle_sub_p$Temp == "15",]$Temp)
  l2 <- length(muscle_sub_p[muscle_sub_p$Temp == "20",]$Temp)
  
  
  if (l1 > l2){
    metadata_sample <- slice(muscle_sub_p[muscle_sub_p$Temp == "15",], 1:min(l1,l2))
    balanced_metadata <- rbind(metadata_sample, muscle_sub_p[muscle_sub_p$Temp == "20",])
    rawcounts_muscle_sub <- counts_p[, colnames(counts_p) %in% rownames(balanced_metadata)]
    balanced_metadata_match <- balanced_metadata[match(colnames(rawcounts_muscle_sub), rownames(balanced_metadata)), ]
    
    print(paste0("l1>l2",all(colnames(rawcounts_muscle_sub) == rownames(balanced_metadata_match))))
    
  } else if (l1 < l2) {
    metadata_sample <- slice(muscle_sub_p[muscle_sub_p$Temp == "20",], 1:min(l1,l2))
    balanced_metadata <- rbind(muscle_sub_p[muscle_sub_p$Temp == "15",], metadata_sample)
    rawcounts_muscle_sub <- counts_p[, colnames(counts_p) %in% rownames(balanced_metadata)]
    balanced_metadata_match <- balanced_metadata[match(colnames(rawcounts_muscle_sub), rownames(balanced_metadata)), ]
    
    print(paste0("l1<l2",all(colnames(rawcounts_muscle_sub) == rownames(balanced_metadata_match))))
    
  } else if (l1 == l2) {
    balanced_metadata <- muscle_sub_p
    rawcounts_muscle_sub <- counts_p[, colnames(counts_p) %in% rownames(balanced_metadata)]
    balanced_metadata_match <- balanced_metadata[match(colnames(rawcounts_muscle_sub), rownames(balanced_metadata)), ]
    
    print(paste0("l1==l2",all(colnames(rawcounts_muscle_sub) == rownames(balanced_metadata_match))))
  }
  
  write.table(rawcounts_muscle_sub, file = paste0("./muscle/", p, "_intron_clustering_perind_numers.counts.FINAL"), row.names = T, col.names = T, quote = F)
  write.table(balanced_metadata_match[-c(3)], file = paste0("./muscle/",p,"_muscle_metadata_final.txt"), row.names = T, col.names = F, quote = F)

}



################ Muscle

muscle <- read.table(gzfile("./muscle/intron_clustering_perind_numers.counts.gz"))
colnames(muscle) <- gsub("_INTRON_Aligned.out.bam", "", colnames(muscle))

metadata_muscle <- read.csv("./muscle_metadata2.txt", header = TRUE, sep = ",",row.names = 1)
rownames(metadata_muscle) <- gsub("-", ".", rownames(metadata_muscle))
rownames(metadata_muscle) <- gsub("_noquotes.txt", "", rownames(metadata_muscle))


muscle_sub <-metadata_muscle %>%
  filter(Gosner %in% c(31, 32, 33, 34))


muscle_filtered <- muscle[,colnames(muscle) %in% rownames(muscle_sub)]

muscle_sub_match <- muscle_sub[match(colnames(muscle_filtered), rownames(muscle_sub)), ]

all(colnames(muscle_filtered) == rownames(muscle_sub_match))


write.table(muscle_filtered[,-c(12,17,18)], file = "./muscle/intron_clustering_perind_numers.counts.FINAL", row.names = T, col.names = T, quote = F)
write.table(muscle_sub_match[-c(12,17,18),], file = "./muscle/muscle_metadata_final.txt", row.names = T, col.names = F, quote = F)
system("gzip -f ./muscle/intron_clustering_perind_numers.counts.FINAL")


###### Population separately
### Balance the design

population <- c("P1", "P2", "P3", "P4", "P5")
muscle_sub_match2 <- muscle_sub_match[,-c(12,17,18)]
muscle_filtered2 <- muscle_filtered[-c(12,17,18),]

all(colnames(muscle_filtered2) == rownames(muscle_sub_match2))

for (p in population) {
  
  muscle_sub_p <- muscle_sub_match2[muscle_sub_match2$Pop %in% c(p), ]
  counts_p <- muscle_filtered2[, colnames(muscle_filtered2) %in% rownames(muscle_sub_p)]
  
  l1 <- length(muscle_sub_p[muscle_sub_p$Temp == "15",]$Temp)
  l2 <- length(muscle_sub_p[muscle_sub_p$Temp == "20",]$Temp)
  
  
  if (l1 > l2){
    metadata_sample <- slice(muscle_sub_p[muscle_sub_p$Temp == "15",], 1:min(l1,l2))
    balanced_metadata <- rbind(metadata_sample, muscle_sub_p[muscle_sub_p$Temp == "20",])
    rawcounts_muscle_sub <- counts_p[, colnames(counts_p) %in% rownames(balanced_metadata)]
    balanced_metadata_match <- balanced_metadata[match(colnames(rawcounts_muscle_sub), rownames(balanced_metadata)), ]
    
    print(paste0("l1>l2",all(colnames(rawcounts_muscle_sub) == rownames(balanced_metadata_match))))
    
  } else if (l1 < l2) {
    metadata_sample <- slice(muscle_sub_p[muscle_sub_p$Temp == "20",], 1:min(l1,l2))
    balanced_metadata <- rbind(muscle_sub_p[muscle_sub_p$Temp == "15",], metadata_sample)
    rawcounts_muscle_sub <- counts_p[, colnames(counts_p) %in% rownames(balanced_metadata)]
    balanced_metadata_match <- balanced_metadata[match(colnames(rawcounts_muscle_sub), rownames(balanced_metadata)), ]
    
    print(paste0("l1<l2",all(colnames(rawcounts_muscle_sub) == rownames(balanced_metadata_match))))
    
  } else if (l1 == l2) {
    balanced_metadata <- muscle_sub_p
    rawcounts_muscle_sub <- counts_p[, colnames(counts_p) %in% rownames(balanced_metadata)]
    balanced_metadata_match <- balanced_metadata[match(colnames(rawcounts_muscle_sub), rownames(balanced_metadata)), ]
    
    print(paste0("l1==l2",all(colnames(rawcounts_muscle_sub) == rownames(balanced_metadata_match))))
  }
  
  write.table(rawcounts_muscle_sub, file = paste0("./muscle/", p, "_intron_clustering_perind_numers.counts.FINAL"), row.names = T, col.names = T, quote = F)
  write.table(balanced_metadata_match[-c(2)], file = paste0("./muscle/",p,"_muscle_metadata_final.txt"), row.names = T, col.names = F, quote = F)
  
}
