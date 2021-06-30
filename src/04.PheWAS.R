setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/08.OAS1/")
library(PheWAS)
phenotypes <- readRDS("ICD.rds")
final <- readRDS("WES.dosage1.rds")
final <- final %>% mutate(status = case_when(ENSG00000174840.CADD.dosage == 2 ~ 1,
                                          TRUE ~ 0))
final <- final %>% mutate(status = ENSG00000174840.CADD.dosage)

genotypes <- data.frame(final$IID, final$status)
colnames(genotypes) <- c("id", "status")
covariates <- data.frame(final$IID, final$AGE, final$SEX)
colnames(covariates) <- c("id", "age", "sex")

geno1 <- data.frame(genotypes$id)
colnames(geno1) <- "id"
a <- setdiff(geno1$id, phenotypes$id)

tmp <- data.frame(matrix(0, length(a), length(colnames(phenotypes))))
colnames(tmp) <- colnames(phenotypes)
tmp$id <- a
tmp[,2:dim(tmp)[2]] <- FALSE
pheno <- rbind(phenotypes, tmp)

results <- phewas(phenotypes=pheno,genotypes=genotypes, covariates = covariates)
saveRDS(results, file="phewasPDE12_additive.rds")

results <- readRDS("phewasPDE12_additive.rds")
results2 <- results[is.na(results$p) == FALSE,]
results2$q <- p.adjust(results2$p, method = "BH")
results2 <- results2[,c(1,4,5,6,7,10,11,16)]
results2$OR_L <- exp(results2$beta + qnorm(0.025) * results2$SE)
results2$OR_U <- exp(results2$beta + qnorm(0.975) * results2$SE)
results2 <- results2[,c(1,4,9,10,5,8,6,7)]
#write.table(results2, file="phewasPDE12.txt", sep="\t", row.names = F, col.names = T, quote=F)
head(results2[order(results2$p, decreasing = F),])
png("phewasPDE12_additive.png", width = 1000, height = 600)
phewasManhattan(results2, annotate.level = 5e-3, size.x.labels = 20, annotate.size=5, max.x = 1400, annotate.angle=0, title= paste0("homo/comphet carrieres (N=34185), non-carriers (N=145882)"))
dev.off()
