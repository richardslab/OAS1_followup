setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/08.OAS1/")

final <- readRDS("WES.dosage.rds")

#withdrawal
w <- fread("/scratch/richards/tomoko.nakanishi/09.COVID19/data/01.UKBB/w27449_20210201.csv")
final <- final %>% filter(!(FID %in% w$V1))

pc <- fread("~/09.COVID19/scratch/data/01.UKBB/EMC-White-British.20pc.txt")
colnames(pc) <- c("FID", paste0("PC",1:20))

final <- final %>% merge(pc, by="FID")

infile <- "~/scratch/09.COVID19/data/01.UKBB/ukb27449_20688_fetch_data.tab.gz"
t <- fread(infile, sep="\t")
d <- data.frame(t$f.eid, t$f.31.0.0, t$f.54.0.0, t$f.53.0.0, t$f.21003.0.0, t$f.22000.0.0)
d$GTARRAY <- NA
d$GTARRAY[d$t.f.22000.0.0 < 0] <- "UKBiLEVE"
d$GTARRAY[d$t.f.22000.0.0 > 0] <- "UKAxiom"
d$GTARRAY <- ifelse(d$GTARRAY == "UKAxiom", 1, 0)
d1 <- d[,c(1,2,5)]
colnames(d1)[1:3] <- c("FID", "SEX", "AGE")

final <- final %>% merge(d1, by ="FID")
saveRDS(final, "WES.dosage1.rds")

##PheWAS
hesin <- fread("~/09.COVID19/data/01.UKBB/hesin_diag_20210530.txt.gz")
icd10 <- hesin[is.na(hesin$diag_icd10) == FALSE,]
icd10 <- icd10[,c(1,7)]
names(icd10) <- c("ID", "code")
icd10$vocabulary_id <- "ICD10CM"
icd9 <- hesin[is.na(hesin$diag_icd9) == FALSE,]
icd9 <- icd9 %>% filter(diag_icd9 != "")
icd9 <- icd9[,c(1,5)]
names(icd9) <- c("ID", "code")
icd9$vocabulary_id <- "ICD9CM"
icd <- data.frame(rbind(icd9, icd10))
icd$count <- 1


data <- merge(final, icd, by.x="FID", by.y="ID",all.x=T)
data$code[is.na(data$code)] <- 0
data$vocabulary_id[is.na(data$vocabulary_id)] <- 0

icd10 <- read.table("~/my_project/data/PheWAS/UKB_ICD10_2phecode.tsv", sep="\t", header=T, quote="")
tmp1 <- data[data$vocabulary_id == "ICD10CM",]
colnames(icd10)[1] <- "code"
tmp2 <- merge(tmp1, icd10, by="code", all.x = T, sorted =F)
tmp2 <- tmp2[is.na(tmp2$phecode) == FALSE,]

icd9 <- read.table("~/my_project/data/PheWAS/UKB_ICD9_2phecoderev.tsv", sep="\t", header=T, quote="")
icd9$coding = as.character(icd9$coding)
tmp3 <- data[data$vocabulary_id == "ICD9CM",]
colnames(icd9)[1] <- "code"
tmp4 <- merge(tmp3, icd9, by="code", all.x = T, sorted =F)
tmp4 <- tmp4[is.na(tmp4$phecode) == FALSE,]

tmp5 <- data.frame(rbind(tmp2, tmp4))
icdrev <- data.frame(tmp5$IID, tmp5$vocabulary_id, tmp5$phecode, tmp5$count)
colnames(icdrev) <- c("id", "vocabulary_id", "code", "count")
#icdrev1 <- icdrev[is.na(icdrev$code) == FALSE,]
icdrev <- icdrev %>% mutate(code = as.character(code))

id.sex <- data.frame(final$IID, final$SEX)
colnames(id.sex) <- c("ID", "sex")
id.sex <- ifelse(id.sex$sex == 1, "M", "F")

icdrev$vocabulary_id = "phecode"
library(PheWAS)
phenotypes <- createPhenotypes(icdrev, aggregate.fun=sum, id.sex=id.sex, min.code.count = 1, translate = F, add.phecode.exclusions = T)
saveRDS(phenotypes, file="ICD.rds")
