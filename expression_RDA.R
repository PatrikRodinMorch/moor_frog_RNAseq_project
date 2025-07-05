library("vegan")
library("psych")


### Using a redundancy analysis to explore phenotypic seperation in gene expression.

setwd("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/RDA/")


DE_liver <- t(read.table("liver_vst_DE_for_pheno_corr.txt", header = TRUE, row.names = 1))
DE_liver_metadata <- read.table("liver_vst_metadata_for_pheno_corr.txt", header = TRUE)
DE_muscle <- t(read.table("muscle_vst_DE_for_pheno_corr.txt", header = TRUE, row.names = 1))
DE_muscle_metadata <- read.table("muscle_vst_metadata_for_pheno_corr.txt", header = TRUE)

introns_liver <- t(read.table("Introns_liver_vst_for_pheno_corr.txt", header = TRUE, row.names = 1))
introns_liver_metadata <- read.table("intron_liver_metadata_final.txt", header = TRUE)
introns_muscle <- t(read.table("Introns_muscle_vst_for_pheno_corr.txt", header = TRUE, row.names = 1))
introns_muscle_metadata <- read.table("intron_muscle_metadata_final.txt", header = TRUE)


for (i in seq(from=2, to=2)){
  DE_liver_metadata[,i] <- as.factor(DE_liver_metadata[,i])
  DE_muscle_metadata[,i] <- as.factor(DE_muscle_metadata[,i])
}

for (i in seq(from=4, to=4)){
  introns_liver_metadata[,i] <- as.factor(introns_liver_metadata[,i])
  introns_muscle_metadata[,3] <- as.factor(introns_muscle_metadata[,3])
}


DE_liver.rda <- rda(introns_muscle ~ LP + MAM + Temp, data=introns_muscle_metadata, scale=T)
summary(DE_liver.rda)
RsquareAdj(DE_liver.rda)
screeplot(DE_liver.rda)
vif.cca(DE_liver.rda)

DE_liver.rda.signif.full <- anova.cca(DE_liver.rda, parallel=getOption("mc.cores"))
DE_liver.rda.signif.full

signif.axis <- anova.cca(wolf.rda, by="axis", parallel=getOption("mc.cores"))


#### Plot

eco <- introns_muscle_metadata$Pop
bg <- c("#ff7f00","#1f78b4","#ffff33","#33a02c","#e31a1c") # 6 nice colors for our ecotypes

# axes 1 & 2
plot(DE_liver.rda, type="n", scaling=3)
points(DE_liver.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(DE_liver.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the wolves
text(DE_liver.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
