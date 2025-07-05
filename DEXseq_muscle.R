library(DEXSeq)
library(gtools)
library(EnhancedVolcano)
register(MulticoreParam(20))

#### Script to test differential exon usage as a response to temperature in moor frog tail muscle tissue.

setwd("/proj/naiss2023-23-494/patrik/reference_guided/liver/htseq/dexseq/analysis/muscle/")

### Loading, cleaning and normalizing data
 
metadata_muscle <- read.csv("/proj/naiss2023-23-494/patrik/reference_guided/liver/htseq/dexseq/analysis/muscle/metadata/muscle_metadata2.txt", header = TRUE, sep = ",",
                             row.names = 1)
 
countFiles = list.files("/proj/naiss2023-23-494/patrik/reference_guided/liver/htseq/dexseq/analysis/muscle/", pattern=".txt$", full.names=TRUE)
basename(mixedsort(countFiles))
write.csv(basename(mixedsort(countFiles)), file="new_sample_names.txt", quote = F, row.names = F, col.names = F)
 
flattenedFile = list.files("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/DEXseq/", pattern="gff$", full.names=TRUE)
basename(flattenedFile)
basename(mixedsort(countFiles))
 
rownames(counts_muscle) <- gsub('gene-', '', rownames(counts_muscle))
 
for (i in seq(from=1, to=14)){
   metadata_muscle[,i] <- as.factor(metadata_muscle[,i])
}
 
metadata_muscle_sub <- metadata_muscle[, c("Temp", "Gosner", "Pop")]
metadata_muscle$Temp <- factor(metadata_muscle$Temp, levels=c(20,15))
 
 
## Differential exon usage temp
dxd_temp = DEXSeqDataSetFromHTSeq(
   mixedsort(countFiles),
   sampleData=metadata_muscle_sub,
   design= ~ sample + exon + Temp:exon,
   flattenedfile=flattenedFile )
 
 
formulaFullModel    =  ~ sample + exon + Gosner:exon + Pop:exon + Temp:exon
formulaReducedModel =  ~ sample + exon + Gosner:exon + Pop:exon 
 
dxd_temp = estimateSizeFactors(dxd_temp)
dxd_temp = estimateDispersions(dxd_temp,BPPARAM=MulticoreParam(20), formula = formulaFullModel)
pdf(file="dispersion_muscle.pdf")
plotDispEsts(dxd_temp)
dev.off()
 
 
dxd_temp = testForDEU(dxd_temp, reducedModel = formulaReducedModel, 
                       fullModel = formulaFullModel,
                       BPPARAM=MulticoreParam(20))
 
dxd_temp=estimateExonFoldChanges(dxd_temp,fitExpToVar="Temp",BPPARAM=MulticoreParam(20))

dxr1_temp=DEXSeqResults(dxd_temp)

dxr1_temp@listData$padj[is.na(dxr1_temp@listData$padj)]<-1


de_exons_muscle_temp_up<-dxr1_temp[(dxr1_temp$log2fold_20_15>=0.58)&(dxr1_temp$padj<=0.05),]
de_exons_muscle_temp_down<-dxr1_temp[(dxr1_temp$log2fold_20_15<=-0.58)&(dxr1_temp$padj<=0.05),]
de_exons_muscle<-rbind(de_exons_muscle_temp_up,de_exons_muscle_temp_down)

sign.genes_temp_up<-unique(de_exons_muscle_temp_up$groupID)
sign.genes_temp_up<-gsub('gene-','',sign.genes_temp_up)

sign.genes_temp_down<-unique(de_exons_muscle_temp_down$groupID)
sign.genes_temp_down<-gsub('gene-','',sign.genes_temp_down)

sign.genes_temp<-unique(de_exons_muscle$groupID)
sign.genes_temp<-gsub('gene-','',sign.genes_temp)

write.table(sign.genes_temp_up,file="/home/prm/Skrivbord/Final_RNAseq_small_scale_results/DEXseq/muscle/exons_muscle_temp_up.txt",row.names=F,col.names=F,quote=F)
write.table(sign.genes_temp_down,file="/home/prm/Skrivbord/Final_RNAseq_small_scale_results/DEXseq/muscle/exons_muscle_temp_down.txt",row.names=F,col.names=F,quote=F)
write.table(sign.genes_temp,file="/home/prm/Skrivbord/Final_RNAseq_small_scale_results/DEXseq/muscle/exons_muscle_temp_total.txt",row.names=F,col.names=F,quote=F)

plotMA(dxr1_temp, cex=0.8 )
 
plotDEXSeq(dxr1_temp, "gene-ILF3", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, fitExpToVar = "Temp")
 
 
 
 
 
save.image(file="DEXSeq_muscle.RData")
DEXSeqHTML( dxr1_dev_Fast, fitExpToVar="Development_temp", FDR=0.05, color=c("#FF000080", "#0000FF80") )

load("DEXSeq_muscle.RData")

##### Per population
population <- c("P1", "P2", "P3", "P4", "P5")

pop <- c()
up_regulated <- c()
down_regulated <- c()
total_DEGs <- c()


for (p in population) {
  metadata_muscle_sub_pop <- metadata_muscle_sub[metadata_muscle_sub$Pop %in% c(p), ]
  l1 <- length(metadata_muscle_sub_pop[metadata_muscle_sub_pop$Temp == "15",]$Temp)
  l2 <- length(metadata_muscle_sub_pop[metadata_muscle_sub_pop$Temp == "20",]$Temp)
  
  
  if (l1 > l2){
    metadata_sample <- dplyr::slice(metadata_muscle_sub_pop[metadata_muscle_sub_pop$Temp == "15",], 1:min(l1,l2))
    balanced_metadata <- rbind(metadata_sample, metadata_muscle_sub_pop[metadata_muscle_sub_pop$Temp == "20",])
    
  } else if (l1 < l2) {
    metadata_sample <- dplyr::slice(metadata_muscle_sub_pop[metadata_muscle_sub_pop$Temp == "20",], 1:min(l1,l2))
    balanced_metadata <- rbind(metadata_muscle_sub_pop[metadata_muscle_sub_pop$Temp == "15",], metadata_sample)
    
    
  } else if (l1 == l2) {
    balanced_metadata <- metadata_muscle_sub_pop
    
  }
  
  
  countFiles = list.files("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/DEXseq/muscle_files", pattern=".txt$", full.names=TRUE)
  
  
  countFiles_sub = paste0("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/DEXseq/muscle_files/",basename(mixedsort(countFiles))[basename(mixedsort(countFiles)) %in% rownames(balanced_metadata)])
  
  
  dxd_temp_p = DEXSeqDataSetFromHTSeq(
    mixedsort(countFiles_sub),
    sampleData=balanced_metadata,
    design= ~ sample + exon + Temp:exon,
    flattenedfile=flattenedFile )
  
  
  formulaFullModel    =  ~ sample + exon + Gosner:exon + Temp:exon
  formulaReducedModel =  ~ sample + exon + Gosner:exon
  
  dxd_temp_p = estimateSizeFactors(dxd_temp_p)
  dxd_temp_p = estimateDispersions(dxd_temp_p, BPPARAM=MulticoreParam(20), formula = formulaFullModel)
  pdf(file = paste0("dispersion_muscle","_",p,".pdf"))
  plotDispEsts(dxd_temp_p)
  dev.off()
  
  
  dxd_temp_p = testForDEU(dxd_temp_p, reducedModel = formulaReducedModel, 
                          fullModel = formulaFullModel,
                          BPPARAM=MulticoreParam(20))
  
  dxd_temp_p = estimateExonFoldChanges(dxd_temp_p, fitExpToVar="Temp",BPPARAM=MulticoreParam(20))
  
  dxr1_temp_p = DEXSeqResults(dxd_temp_p)
  
  dxr1_temp_p@listData$padj[is.na(dxr1_temp_p@listData$padj)] <- 1
  
  
  de_exons_muscle_temp_up <- dxr1_temp_p[(dxr1_temp_p$log2fold_20_15 >= 0.58) & (dxr1_temp_p$padj <= 0.05),]
  de_exons_muscle_temp_down <- dxr1_temp_p[(dxr1_temp_p$log2fold_20_15 <= -0.58) & (dxr1_temp_p$padj <= 0.05),]
  de_exons_muscle <- rbind(de_exons_muscle_temp_up, de_exons_muscle_temp_down)
  
  sign.genes_temp_up <- unique(de_exons_muscle_temp_up$groupID)
  sign.genes_temp_up <- gsub('gene-', '', sign.genes_temp_up)
  
  sign.genes_temp_down <- unique(de_exons_muscle_temp_down$groupID)
  sign.genes_temp_down <- gsub('gene-', '', sign.genes_temp_down)
  
  sign.genes_temp <- unique(de_exons_muscle$groupID)
  sign.genes_temp <- gsub('gene-', '', sign.genes_temp)
  
  up_regulated <- c(up_regulated, length(sign.genes_temp_up))
  down_regulated <- c(down_regulated, length(sign.genes_temp_down))
  pop <- c(pop, p)
  total_DEGs <- c(total_DEGs, c(length(sign.genes_temp_up) + length(sign.genes_temp_down)))
  
  write.table(sign.genes_temp_up, file = paste0("exons_muscle_temp_up_",p,".txt"), row.names = F, col.names = F, quote = F)
  write.table(sign.genes_temp_down, file = paste0("exons_muscle_temp_down_",p,".txt"), row.names = F, col.names = F, quote = F)
  write.table(sign.genes_temp, file = paste0("exons_muscle_temp_total_",p,".txt"), row.names = F, col.names = F, quote = F)
  
  
}

pop_DEG <- data.frame(pop, up_regulated, down_regulated, total_DEGs) |> 
  mutate(Pop = recode(pop, "P4"="P18", "P2"="P4", "P3"="P26", "P5"="P14", "P1" = "P1")) 

colnames(pop_DEG) <- c("pop", "up_regulated", "down_regulated", "total_DEGs", "Pop")
pop_DEG2 <- gather(pop_DEG[,c("up_regulated", "down_regulated", "Pop")], "up_regulated", "degs", -Pop) 


pdf(file = paste0("muscle_within_population_DEGs.pdf"))  


ggplot(pop_DEG2, aes(fill=up_regulated, y=degs, x=Pop)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw() + 
  #scale_fill_discrete(guide = guide_legend(reverse=TRUE))
  scale_fill_discrete(name = "", labels = c("Down-regulated", "Up-regulated")) + 
  xlab("") + 
  ylab("Number of genes with differential exon usage")

dev.off()
