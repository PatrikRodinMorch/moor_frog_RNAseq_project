library("o2plsda")

setwd("/home/prm/Skrivbord/Final_RNAseq_small_scale_results/RDA/")


DE_liver <- t(read.table("liver_vst_DE_for_pheno_corr.txt", header = TRUE, row.names = 1))
DE_liver_metadata <- read.table("liver_vst_metadata_for_pheno_corr.txt", header = TRUE)
DE_liver_pred <- subset(DE_liver_metadata, select=c("LP", "MAM", "Temp"))
DE_liver_metadata$Temp <- as.factor(DE_liver_metadata$Temp)

DE_muscle <- t(read.table("muscle_vst_DE_for_pheno_corr.txt", header = TRUE, row.names = 1))
DE_muscle_metadata <- read.table("muscle_vst_metadata_for_pheno_corr.txt", header = TRUE)
DE_muscle_pred <- subset(DE_muscle_metadata, select=c("LP", "MAM", "Temp"))
DE_muscle_metadata$Temp <- as.factor(DE_muscle_metadata$Temp)


introns_liver <- t(read.table("Introns_liver_vst_for_pheno_corr.txt", header = TRUE, row.names = 1))
introns_liver_metadata <- read.table("intron_liver_metadata_final.txt", header = TRUE)
introns_liver_pred <- subset(introns_liver_metadata, select=c("LP", "MAM", "Temp"))
introns_liver_metadata$Temp <- as.factor(introns_liver_metadata$Temp)

introns_muscle <- t(read.table("Introns_muscle_vst_for_pheno_corr.txt", header = TRUE, row.names = 1))
introns_muscle_metadata <- read.table("intron_muscle_metadata_final.txt", header = TRUE)
introns_muscle_pred <- subset(introns_muscle_metadata, select=c("LP", "MAM", "Temp"))
introns_muscle_metadata$Temp <- as.factor(introns_muscle_metadata$Temp)


X = scale(introns_muscle, scale = TRUE)
Y = scale(introns_muscle_pred, scale = TRUE)

set.seed(123)
## nr_folds : cross validation k-fold (suggest 10)
## ncores : parallel paramaters for large datasets
cv <- o2cv(X, Y, 1:1,1:1,1:1, group = introns_muscle_metadata$Temp, nr_folds = 10)
fit <- o2pls(X,Y,1,1,1)
summary(fit)

Xl <- loadings(fit,loading="Xjoint")
Xs <- scores(fit,score="Xjoint")
plot(fit,type="score",var="Xjoint", introns_muscle_metadata$Temp, label = FALSE)

plot(fit,type="loading",var="Xjoint", group=introns_muscle_metadata$Temp,repel=F,rotation=TRUE)


res <- oplsda(fit,introns_muscle_metadata$Temp, nc=1)
plot(res,type="score", group=introns_muscle_metadata$Temp,repel=TRUE, label = FALSE)

vip <- vip(res)
plot(res,type="vip", group=introns_muscle_metadata$Temp, repel = FALSE,order=TRUE)
