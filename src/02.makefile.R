setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/08.OAS1/")

#ENSG00000174840 PDE12 pLI/LOEUF = 0.04/0.6
ENSG00000174840.LOF <- fread("ENSG00000174840.LoF.raw")
ENSG00000174840.LOF <- ENSG00000174840.LOF %>% mutate_at(.vars = vars(c(`chr3:57556384:G:A_G`:`chr3:57557688:G:A_G`)), 
                                                         .funs = funs(2 - .))

ENSG00000174840.LOF <- ENSG00000174840.LOF %>% mutate(dosage = select(., c(`chr3:57556384:G:A_G`:`chr3:57557688:G:A_G`)) %>%
                                                        rowSums(na.rm = T))


#ENSG00000174840 PDE12 pLI/LOEUF = 0.04/0.6
ENSG00000174840.CADD <- fread("ENSG00000174840.CADD.raw")
ENSG00000174840.CADD <- ENSG00000174840.CADD %>% mutate_at(.vars = vars(c(`chr3:57556382:G:C_G`:`chr3:57559981:G:C_G`)), 
                                                         .funs = funs(2 - .))

ENSG00000174840.CADD <- ENSG00000174840.CADD %>% mutate(dosage = select(., c(`chr3:57556382:G:C_G`:`chr3:57559981:G:C_G`)) %>%
                                                        rowSums(na.rm = T))

ENSG00000174840.CADD <- ENSG00000174840.CADD %>% mutate(dosage = ifelse(dosage > 2, 2, dosage))

#ENSG00000118507 AKAP7 pLI/LOUEF = 0/1.22
ENSG00000118507.LOF <- fread("ENSG00000118507.LoF.raw")
ENSG00000118507.LOF <- ENSG00000118507.LOF %>% mutate_at(.vars = vars(c(`chr6:131165218:G:A_G`:`chr6:131281529:G:A_G`)), 
                                                         .funs = funs(2 - .))

ENSG00000118507.LOF <- ENSG00000118507.LOF %>% mutate(dosage = select(., c(`chr6:131165218:G:A_G`:`chr6:131281529:G:A_G`)) %>%
                                                        rowSums(na.rm = T))

ENSG00000118507.CADD <- fread("ENSG00000118507.CADD.raw")
ENSG00000118507.CADD <- ENSG00000118507.CADD %>% mutate_at(.vars = vars(c(`chr6:131135766:G:C_G`:`chr6:131281720:G:C_G`)), 
                                                         .funs = funs(2 - .))

ENSG00000118507.CADD <- ENSG00000118507.CADD %>% mutate(dosage = select(., c(`chr6:131135766:G:C_G`:`chr6:131281720:G:C_G`)) %>%
                                                        rowSums(na.rm = T))
ENSG00000118507.CADD <- ENSG00000118507.CADD %>% mutate(dosage = ifelse(dosage > 2, 2, dosage))

#
ENSG00000174840.CADD <- ENSG00000174840.CADD %>% select(c(FID, IID, dosage)) 
ENSG00000174840.CADD <- ENSG00000174840.CADD %>% rename(ENSG00000174840.CADD.dosage = dosage)
ENSG00000118507.CADD <- ENSG00000118507.CADD %>% select(c(FID, IID, dosage)) 
ENSG00000118507.CADD <- ENSG00000118507.CADD %>% rename(ENSG00000118507.CADD.dosage = dosage)

final <- merge(ENSG00000174840.CADD, ENSG00000118507.CADD, by=c("FID","IID"))
final %>% saveRDS("WES.dosage.rds")
