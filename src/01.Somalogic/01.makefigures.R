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

#AKAP7.18399.1
#"IRF9.12439.67" "IRF1.17462.19" "IRF2.12801.33" "IRF6.9999.1"   "IRF4.19564.61" "IRF3.17151.84"
#"PDE5A.5256.86"  "PDE4A.18918.86" "PDE4D.5255.22"  "PDE6D.13491.40" "PDE7A.5178.5"   "PDE9A.5201.50" 
#[7] "PDE1A.5253.1"   "PDE2A.5246.64"  "PDE11A.5252.33" "PDE3A.5254.69"  "SPDEF.10012.5"  "PDE5A.16805.5" 

data <- vap %>% inner_join(basic, by="anonymized_patient_id")
data <- data %>% mutate(time_from_first_sx = as.numeric(as.Date(date_sampling) - as.Date(covid19_first_symptoms_date)))

adat_long_select <- adat_long %>% filter(AptName %in% c("OAS1.10361.25", "AKAP7.18399.1","IRF3.17151.84","IFNB1.14127.240"))
adat_wide_select <- adat_long_select %>% select(c("SubjectID", "AptName", "value")) %>% drop_na(SubjectID) %>% 
  tidyr::spread(key = AptName, value = value)
cor(adat_wide_select[,-1])
dat_all <- adat_wide_select %>% inner_join(data, by=c("SubjectID"="VAPcode")) %>% mutate(hospital = ifelse(grepl("CHUM", anonymized_patient_id), "CHUM","JGH"))

dat_luminex <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/BQC_phenotype/Luminex_vRNA_JGH_CHUM.rds") %>% select(-c("anonymized_patient_id", "VAP", "date_of_sampling"))
dat_luminex <- dat_luminex[!(duplicated(dat_luminex)),]

final <- merge(dat_all, dat_luminex, by.x="SubjectID", by.y="VAPcode", all.x=T)
final <- final[order(final$time_from_first_sx, decreasing = F),]
cor(final[,c("Ncopy","AKAP7.18399.1")], use = "complete.obs")
geno1 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/scratch/05.BQC/01.genotypeQC/v3.1/raws/12.raw") %>% select(c("FID","chr12:112919388:G:A_G"))
geno2 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/scratch/05.BQC/01.genotypeQC/v3.1/raws/21.raw") %>% select(c("FID","chr21:33242905:T:C_T"))
geno <- inner_join(geno1, geno2, by="FID")
map <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v3.1/01.batch/sampleid.studyid.sex.map", header=F)
geno <- geno %>% inner_join(map, by=c("FID"="V1"))

pc1 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v3.1/03.Ancestry/all.sample.ancestry") 
#pc1 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v3.1/09.PC/AFR.sample.PCs") 

geno <- geno %>% merge(pc1, by.x="FID", by.y="ID", all.x=T)

final <- final %>% merge(geno, by.x="anonymized_patient_id", by.y="V2.x", all.x=T)
cor.test(final$AKAP7.18399.1, final$OAS1.10361.25)

final <- final[order(final$time_from_first_sx, decreasing=F),]
final <- final %>% mutate(case = case_when(A2 == 1 ~ "03Critical",
                                           covid19_test == 1 ~ "02Mild",
                                           TRUE ~ "01Non-COVID"))

final <- final %>% mutate(snp1 = round(`chr12:112919388:G:A_G`),
                          snp2 = round(`chr21:33242905:T:C_T`))

##withdrawal
w <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/BQC_phenotype/JGH/20210615/withdrawal.list", header=F)
w <- w %>% mutate(id = paste0("JGH_",V1))
final <- final %>% filter(!(anonymized_patient_id %in% w$id))
final <- final %>% mutate(AKAP7 = log(AKAP7.18399.1),
                        OAS1 = log(OAS1.10361.25),
                        IRF3 = log(IRF3.17151.84))

##longitudinal 

dat1 <- final %>% group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == last(OAS1.10361.25)) %>% ungroup() %>%
  distinct(anonymized_patient_id, .keep_all=TRUE)
table(dat1$C2)

######### ask Mount Sinai to replicate Figure 1A Figure 1B Figure 2A ##################
#final > all longitudinal sample
final <- final %>% mutate(tmp = 1)
tmp <- final %>% nest(data = -tmp) %>%
  mutate(fit = map(data, ~ lm(OAS1 ~ age_at_diagnosis + sex + time_to_processing, data = .x)),
         augmented = map(fit, augment)) %>%
  unnest(augmented)
final$OAS1_rev <- tmp$.std.resid
png("/home/richards/tomoko.nakanishi/my_project/repo/OAS1/results/Figure1A.png", width=700, height = 400)
ggplot(final, aes(x=time_from_first_sx, y=OAS1_rev, col=as.factor(C2))) + geom_line(alpha=0.2,size=0.4,aes(group=as.factor(anonymized_patient_id),col=as.factor(C2)))+ 
  theme_bw() + geom_smooth(se=T,size=1, aes(group=as.factor(C2),col=as.factor(C2))) + xlim(0,30) +
  xlab("Days from covid symptoms") + 
  ylab("standardised OAS1") +
  theme(legend.title = element_text( size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20)) +
  scale_colour_brewer(palette="Set1", labels = c("Negative (N=108)", "Positive (N=395)"), name="COVID-19 status") +
  guides(colour = guide_legend(reverse=TRUE))
dev.off()

tmp <- final %>% filter(C2 == 1) %>% filter(snp1 %in% c(0,2))
dat1 <- tmp %>% group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == last(OAS1.10361.25)) %>% ungroup() %>%
  distinct(anonymized_patient_id, .keep_all=TRUE)
table(dat1$snp1)

png("/home/richards/tomoko.nakanishi/my_project/repo/OAS1/results/Figure1B.png", width=700, height = 400)
ggplot(tmp, aes(x=time_from_first_sx, y=OAS1_rev, col=as.factor(snp1))) + geom_line(alpha=0.2,size=0.4,aes(group=as.factor(anonymized_patient_id),col=as.factor(snp1)))+ 
  theme_bw() + geom_smooth(se=T,size=1, aes(group=as.factor(snp1),col=as.factor(snp1))) + xlim(0,30) +
  xlab("Days from covid symptoms") + 
  ylab("standardised OAS1") +
  theme(legend.title = element_text( size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20)) +
  scale_colour_brewer(palette="Dark2", labels = c("AA (N=150)", "GG (N=73)"), name="rs10774671") +
  guides(colour = guide_legend(reverse=TRUE))
dev.off()

###non-infectious
dat1 <- final %>% filter((covid19_test == 1 & time_from_first_sx >= 30) | (covid19_test == 0))  %>% group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == last(OAS1.10361.25)) %>% ungroup() %>%
  distinct(anonymized_patient_id, .keep_all=TRUE)
#dat1 <- dat1 %>% filter(UMAP_POP == "EUR")
dat1 <- dat1 %>% mutate(tmp = 1)
tmp <- dat1 %>% nest(data = -tmp) %>%
  mutate(fit = map(data, ~ lm(OAS1 ~ age_at_diagnosis + sex + time_to_processing, data = .x)),
         augmented = map(fit, augment)) %>%
  unnest(augmented)
dat1$OAS1_rev <- tmp$.std.resid


#dat1 <- dat1 %>% filter(time_to_processing < 80)


#dat1 <- dat1 %>% mutate(AKAP7 = log(AKAP7.18399.1),
#                        OAS1 = log(OAS1.10361.25),
#                        IRF3 = log(IRF3.17151.84))

genotype_names <- c(
  '0'="AA",
  '1'="GA",
  '2'="GG"
)

#dat1 <- dat1 %>% filter(!is.na(PC1))
png("/home/richards/tomoko.nakanishi/my_project/repo/OAS1/results/Figure2A.png", width=800, height = 500)
ggplot(dat1, aes(x=case, y=OAS1_rev, fill=case)) + geom_jitter(aes(color=case)) +
  geom_violin(trim = F) + geom_boxplot(outlier.colour = NA) + facet_wrap(~snp1, labeller=as_labeller(genotype_names)) + theme_classic() + 
  xlab("") +
  ylab("standardised plasma OAS1") + theme(legend.title = element_text( size = 20),
                                            legend.text = element_text(size = 20),
                                            axis.title.x = element_blank(),
                                            axis.text.x = element_blank(),
                                            axis.text.y = element_text(size=20),
                                            axis.title.y = element_text(size=20),
                                           strip.text.x = element_text(size = 20)) +
  guides(colour=FALSE) + scale_fill_brewer(palette="Set2", labels = c("Negative", "Mild", "Critical"), name="COVID-19 status") + 
  scale_colour_brewer(palette="Set2") 
dev.off()

######### ask Mount Sinai to replicate ##################
table(dat1$snp1, dat1$case)

dat1 <- final %>% filter((covid19_test == 1 & time_from_first_sx >= 30) | (covid19_test == 0))  %>% group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == last(OAS1.10361.25)) %>% ungroup() %>%
  distinct(anonymized_patient_id, .keep_all=TRUE)
#dat1 <- dat1 %>% filter(UMAP_POP == "EUR")
dat1 <- dat1 %>% mutate(tmp = 1)
tmp <- dat1 %>% nest(data = -tmp) %>%
  mutate(fit = map(data, ~ lm(OAS1 ~ time_to_processing, data = .x)),
         augmented = map(fit, augment)) %>%
  unnest(augmented)
dat1$OAS1_rev <- tmp$.std.resid

out <- data.frame(matrix(0, 6, 7))
colnames(out) <- c("phenotype", "adj", "beta", "se", "pval", "Ncase", "Ncontrol")
LM <- glm(A2 ~ OAS1_rev + snp1 + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family="binomial")
out$phenotype[1] <- "Critical vs Mild+Negative"
out$adj[1] <- "SNP"
out[1,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[1,6] <- sum(dat1$A2 == 1)
out[1,7] <- sum(dat1$A2 == 0)
LM <- glm(A2 ~ OAS1_rev + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family="binomial")
summary(LM)
out$phenotype[2] <- "Critical vs Mild+Negative"
out$adj[2] <- "None"
out[2,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[2,6] <- sum(dat1$A2 == 1)
out[2,7] <- sum(dat1$A2 == 0)
LM <- glm(C2 ~ OAS1_rev + snp1 + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family="binomial")
summary(LM)
out$phenotype[3] <- "Positive vs Negative"
out$adj[3] <- "SNP"
out[3,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[3,6] <- sum(dat1$C2 == 1)
out[3,7] <- sum(dat1$C2 == 0)
LM <- glm(C2 ~ OAS1_rev + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family="binomial")
summary(LM)
out$phenotype[4] <- "Positive vs Negative"
out$adj[4] <- "None"
out[4,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[4,6] <- sum(dat1$C2 == 1)
out[4,7] <- sum(dat1$C2 == 0)
dat2 <- dat1 %>% filter(case %in% c("02Mild", "03Critical"))
LM <- glm(A2 ~ OAS1_rev + snp1 + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat2, family="binomial")
summary(LM)
out$phenotype[5] <- "Critical vs Mild"
out$adj[5] <- "SNP"
out[5,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[5,6] <- sum(dat2$A2 == 1)
out[5,7] <- sum(dat2$A2 == 0)
LM <- glm(A2 ~ OAS1_rev + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat2, family="binomial")
summary(LM)
out$phenotype[6] <- "Critical vs Mild"
out$adj[6] <- "None"
out[6,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[6,6] <- sum(dat2$A2 == 1)
out[6,7] <- sum(dat2$A2 == 0)

out <- out %>% mutate(OR = exp(beta), UL = exp(beta + qnorm(0.975)*se),
                      LL = exp(beta - qnorm(0.975)*se))
out <- out %>% mutate(phenotype = case_when(phenotype == "Critical vs Mild+Negative" ~ "Critical (32) vs Mild+Negative (195)",
                                            phenotype == "Critical vs Mild" ~ "Critical (32) vs Mild (87)",
                                            phenotype == "Positive vs Negative" ~ "Positive (119) vs Negative (108)"))
out$phenotype <- factor(out$phenotype, levels = c("Critical (32) vs Mild+Negative (195)", "Critical (32) vs Mild (87)", "Positive (119) vs Negative (108)"))

png("/home/richards/tomoko.nakanishi/my_project/repo/OAS1/results/Figure2B.png", width=600, height = 500)
p1 <- ggplot(out, aes(x=adj, y=OR, ymin=LL, ymax=UL, color=adj)) +
  geom_pointrange(aes(col=adj), lwd=0.8) + geom_hline(aes(fill=adj), yintercept =1, linetype=2) +
  xlab("") + ylab("Odds Ratio (95% Confidence Interval) per 1SD of OAS1") + geom_text(aes(label=round(OR,1), y=OR, col=adj, hjust = 0.5, vjust = 1.3), size=5) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=adj), width=0.1, cex=1) + ylim(0.01,1.5)+
  facet_wrap(~phenotype,  strip.position = 'left', nrow = 10) + theme_minimal() +
  scale_color_manual(values = c("#CC6600", "#CC0000"), labels = c("basic + rs10774671", "basic"), name="adjustment")+ labs(color="Group") +
  theme(plot.title=element_text(size=15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15,face="bold"), 
        legend.title = element_text(size=15,face="bold"), 
        strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=15,face="bold"))+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))
p1
dev.off()

###AKAP7
dat1 <- final %>% group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == last(OAS1.10361.25)) %>% ungroup() %>%
  distinct(anonymized_patient_id, .keep_all=TRUE)
table(dat1$C2)

final <- final %>% mutate(tmp = 1)
tmp <- final %>% nest(data = -tmp) %>%
  mutate(fit = map(data, ~ lm(AKAP7 ~ time_to_processing, data = .x)),
         augmented = map(fit, augment)) %>%
  unnest(augmented)
final$AKAP7_rev <- tmp$.std.resid
png("/home/richards/tomoko.nakanishi/my_project/repo/OAS1/results/Figure3A.png", width=700, height = 400)
ggplot(final, aes(x=time_from_first_sx, y=AKAP7_rev, col=as.factor(C2))) + geom_line(alpha=0.2,size=0.4,aes(group=as.factor(anonymized_patient_id),col=as.factor(C2)))+ 
  theme_bw() + geom_smooth(se=T,size=1, aes(group=as.factor(C2),col=as.factor(C2))) + xlim(0,30) +
  xlab("Days from covid symptoms") + 
  ylab("standardised AKAP7") +
  theme(legend.title = element_text( size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20)) +
  scale_colour_brewer(palette="Set1", labels = c("Negative (N=108)", "Positive (N=395)"), name="COVID-19 status") +
  guides(colour = guide_legend(reverse=TRUE))
dev.off()

dat1 <- final %>% group_by(anonymized_patient_id) %>% 
  filter(OAS1.10361.25 == first(OAS1.10361.25)) %>% ungroup() %>%
  distinct(anonymized_patient_id, .keep_all=TRUE)

out <- data.frame(matrix(0, 6, 7))
colnames(out) <- c("phenotype", "adj", "beta", "se", "pval", "Ncase", "Ncontrol")
LM <- glm(A2 ~ AKAP7_rev + snp1 + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family="binomial")
out$phenotype[1] <- "Critical vs Mild+Negative"
out$adj[1] <- "SNP"
out[1,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[1,6] <- sum(dat1$A2 == 1)
out[1,7] <- sum(dat1$A2 == 0)
LM <- glm(A2 ~ AKAP7_rev + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family="binomial")
summary(LM)
out$phenotype[2] <- "Critical vs Mild+Negative"
out$adj[2] <- "None"
out[2,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[2,6] <- sum(dat1$A2 == 1)
out[2,7] <- sum(dat1$A2 == 0)
LM <- glm(C2 ~ AKAP7_rev + snp1 + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family="binomial")
summary(LM)
out$phenotype[3] <- "Positive vs Negative"
out$adj[3] <- "SNP"
out[3,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[3,6] <- sum(dat1$C2 == 1)
out[3,7] <- sum(dat1$C2 == 0)
LM <- glm(C2 ~ AKAP7_rev + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat1, family="binomial")
summary(LM)
out$phenotype[4] <- "Positive vs Negative"
out$adj[4] <- "None"
out[4,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[4,6] <- sum(dat1$C2 == 1)
out[4,7] <- sum(dat1$C2 == 0)
dat2 <- dat1 %>% filter(case %in% c("02Mild", "03Critical"))
LM <- glm(A2 ~ AKAP7_rev + snp1 + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat2, family="binomial")
summary(LM)
out$phenotype[5] <- "Critical vs Mild"
out$adj[5] <- "SNP"
out[5,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[5,6] <- sum(dat2$A2 == 1)
out[5,7] <- sum(dat2$A2 == 0)
LM <- glm(A2 ~ AKAP7_rev + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5, data=dat2, family="binomial")
summary(LM)
out$phenotype[6] <- "Critical vs Mild"
out$adj[6] <- "None"
out[6,3:5] <- summary(LM)$coefficients[2,c(1,2,4)]
out[6,6] <- sum(dat2$A2 == 1)
out[6,7] <- sum(dat2$A2 == 0)

out <- out %>% mutate(OR = exp(beta), UL = exp(beta + qnorm(0.975)*se),
                      LL = exp(beta - qnorm(0.975)*se))
out <- out %>% mutate(phenotype = case_when(phenotype == "Critical vs Mild+Negative" ~ "Critical (140) vs Mild+Negative (363)",
                                            phenotype == "Critical vs Mild" ~ "Critical (395) vs Mild (108)",
                                            phenotype == "Positive vs Negative" ~ "Positive (140) vs Negative (255)"))
out$phenotype <- factor(out$phenotype, levels = c("Critical (140) vs Mild+Negative (363)", "Critical (395) vs Mild (108)", "Positive (140) vs Negative (255)"))


png("/home/richards/tomoko.nakanishi/my_project/repo/OAS1/results/Figure3B.png", width=600, height = 500)
p1 <- ggplot(out, aes(x=adj, y=OR, ymin=LL, ymax=UL, color=adj)) +
  geom_pointrange(aes(col=adj), lwd=0.8) + geom_hline(aes(fill=adj), yintercept =1, linetype=2) +
  xlab("") + ylab("Odds Ratio (95% Confidence Interval) per 1SD of AKAP7") + geom_text(aes(label=round(OR,1), y=OR, col=adj, hjust = 0.5, vjust = 1.3), size=5) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=adj), width=0.1, cex=1) + ylim(0.9,2.5)+
  facet_wrap(~phenotype,  strip.position = 'left', nrow = 10) + theme_minimal() +
  scale_color_manual(values = c("#CC6600", "#CC0000"), labels = c("basic + rs10774671", "basic"), name="adjustment")+ labs(color="Group") +
  theme(plot.title=element_text(size=15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15,face="bold"), 
        legend.title = element_text(size=15,face="bold"), 
        strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=15,face="bold"))+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))
p1
dev.off()
