setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/08.OAS1")
library(SomaDataIO)
library(tidyr)
library(dplyr)
library(lubridate)
library(data.table)
library(ggplot2)
#sample.adat <- read.adat("/home/richards/tomoko.nakanishi/09.COVID19/data/Somalogic/SS-200150_v4_ACDPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.adat")
sample.adat <- read.adat("/home/richards/tomoko.nakanishi/09.COVID19/data/Somalogic/SS-200150_v4_ACDPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.medNormRefSMP.adat")
meta.vec <- getMeta(sample.adat)
tbl <- getFeatureData(sample.adat)
ex_set <- adat2eSet(sample.adat)
adat_long <- meltExpressionSet(ex_set)
basic <- readRDS("/project/richards/tomoko.nakanishi/09.COVID19/scratch/05.BQC/BQC_phenotype/basic_JGH_CHUM_20210428.rds")
vap <- readRDS("/project/richards/tomoko.nakanishi/09.COVID19/scratch/05.BQC/BQC_phenotype/VAPinformation_JGH_CHUM_20210418.rds")

data <- vap %>% inner_join(basic, by="anonymized_patient_id")
data <- data %>% mutate(time_from_first_sx = as.numeric(as.Date(date_sampling) - as.Date(covid19_first_symptoms_date)))

adat_wide_select <- adat_long %>% filter(Organism == "Human") %>% select(c("SubjectID", "AptName", "value")) %>% drop_na(SubjectID) %>% 
  tidyr::spread(key = AptName, value = value)

dat_all <- adat_wide_select %>% inner_join(data, by=c("SubjectID"="VAPcode")) %>% mutate(hospital = ifelse(grepl("CHUM", anonymized_patient_id), "CHUM","JGH"))
final <- dat_all

final <- final[order(final$time_from_first_sx, decreasing = F),]

geno1 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/scratch/05.BQC/01.genotypeQC/v3.1/raws/12.raw") %>% select(c("FID","chr12:112919388:G:A_G"))
geno2 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/scratch/05.BQC/01.genotypeQC/v3.1/raws/21.raw") %>% select(c("FID","chr21:33242905:T:C_T"))
geno <- inner_join(geno1, geno2, by="FID")
map <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v3.1/01.batch/sampleid.studyid.sex.map", header=F)
geno <- geno %>% inner_join(map, by=c("FID"="V1"))

pc1 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v3.1/09.PC/EUR.sample.PCs") 
#pc1 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v3.1/09.PC/AFR.sample.PCs") 

geno <- geno %>% merge(pc1, by="FID", all.x=T)

final <- final %>% merge(geno, by.x="anonymized_patient_id", by.y="V2", all.x=T)
cor.test(final$AKAP7.18399.1, final$OAS1.10361.25)

final <- final[order(final$time_from_first_sx, decreasing=F),]
final <- final %>% mutate(case = case_when(A2 == 1 ~ "03Critical",
                                           covid19_test == 1 ~ "02Mild",
                                           TRUE ~ "01Non-COVID"))

final <- final %>% mutate(snp1 = round(`chr12:112919388:G:A_G`),
                          snp2 = round(`chr21:33242905:T:C_T`))

###non-infectious
dat1 <- final %>% filter((covid19_test == 1 & time_from_first_sx > 30) | (covid19_test == 0))  %>% group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == last(OAS1.10361.25)) %>% ungroup() %>%
  distinct(anonymized_patient_id, .keep_all=TRUE)

somamarlist <- colnames(dat1)[3:5034]
out <- data.frame(matrix(0, length(somamarlist), 4))
colnames(out) <- c("SOMAmar", "beta", "se", "pvalue")
for(i in seq(1,length(somamarlist))){
  out$SOMAmar[i] <- somamarlist[i]
  tmp1 <- dat1 %>% rename(value = somamarlist[i]) %>%
    mutate(value = scale(log(value)))
  LM <- glm(A2 ~ value + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=tmp1, family="binomial")
  sum <- summary(LM)$coefficients
  out[i,2:4] <- sum[2,c(1,2,4)]
}


##acute infectious
dat1 <- final %>% filter((covid19_test == 1 & time_from_first_sx <= 14) | (covid19_test == 0))  %>% group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == first(OAS1.10361.25)) %>% ungroup() %>%
  distinct(anonymized_patient_id, .keep_all=TRUE)

somamarlist <- colnames(dat1)[3:5034]
out1 <- data.frame(matrix(0, length(somamarlist), 4))
colnames(out1) <- c("SOMAmar", "beta", "se", "pvalue")
for(i in seq(1,length(somamarlist))){
  out1$SOMAmar[i] <- somamarlist[i]
  tmp1 <- dat1 %>% rename(value = somamarlist[i]) %>%
    mutate(value = scale(log(value)))
  LM <- glm(A2 ~ value + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=tmp1, family="binomial")
  sum <- summary(LM)$coefficients
  out1[i,2:4] <- sum[2,c(1,2,4)]
}

hist(out1$pvalue)

out1_rev <- out1 %>% filter(pvalue < 0.05 & beta > 0)
out1_rev <- out1_rev[order(out1_rev$pvalue, decreasing = F),]
out_rev <- out %>% filter(SOMAmar %in% out1_rev$SOMAmar)

hist(out_rev$pvalue)

out_rev1 <- out_rev %>% filter(pvalue < 0.05)

ggplot(out_rev1, aes(x=beta, y=-log10(pvalue))) + geom_point()
