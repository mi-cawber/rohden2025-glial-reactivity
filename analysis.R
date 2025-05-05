library(arm)
library(MASS)
library(merTools)
library(lme4)
library(broom)
library(ggplot2)
library(tidyr)
library(sjPlot)
library(sjmisc)
library(dplyr)
library(lmerTest)
library(wesanderson)
library(effects)
library(MuMIn)
library(lsmeans)
library(glmmTMB)
library(optimx)
library(readxl)
library(dplyr)
library(table1)


#### ADNI ####

Ng <- read.csv  ("C://BLENNOW_CSF_NG.csv")

GFAP <- read.csv ("C://CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv")

GAP <- read.csv ("C://BLENNOW_LAB_CSF_GAP.csv")

Abeta <- read.csv("C://UPENNBIOMK_ROCHE_ELECSYS.csv")

adnimerge <- read.csv ("C://adnimerge.csv")

#Filter by base line
adnimerge_bl = adnimerge %>%
  filter(VISCODE == "bl") 

GFAP_bl = GFAP %>%
  filter(VISCODE2 == "bl") 

GFAP_bl_sel = GFAP_bl %>%
  select(RID,  X3034.1) 

Ng_bl = Ng %>%
  filter(VISCODE == "bl") 

Abeta_bl = Abeta %>%
  filter(VISCODE2 == "bl") 

GAP_bl = GAP %>%
  filter(VISCODE2 == "bl") 


# Perform the merge using the correct column names
adni_merge1 <- merge(adnimerge_bl, GFAP_bl_sel, by = c("RID"), all.x = TRUE)

adni_merge2 <- merge(adni_merge1, Ng_bl, by = c("RID"), all.x = TRUE)

adni_merge3 <- merge(adni_merge2, Abeta_bl_sel, by = c("RID"), all.x = TRUE)

adni_merge4 <- merge(adni_merge3, GAP_bl, by = c("RID"), all.x = TRUE)



#### TRIAD ####

TRIAD_cohort <- read.csv("C://TRIAD_cohort.csv")

TOTAL  <-  TRIAD_cohort

#Filter by base line or visit zero
VM00 <- filter(TOTAL, visit == "VM00")


## To Both Cohorts ##

#Filter of interesse biomarkers (GAP43, NG, GFAP, sTREM2)
VM00 = VM00 %>%
  filter(!is.na(biomarker))
VM00 = VM00 %>% 
  filter(!is.na(biomarker)& biomarker != "NA" )
VM00$biomarker <- as.numeric(VM00$biomarker)

#outliers
mediana <- median(VM00$biomarker)
IQR <- IQR(VM00$biomarker)

Q1 <- quantile(VM00$biomarker, 0.25)
Q3 <- quantile(VM00$biomarker, 0.75)
IQR <- Q3 - Q1
lim_sup <- Q3 + 3*IQR
lim_inf <- Q1 - 3*IQR
outliers <- VM00$colun_name[which(VM00$biomarker < lim_inf | VM00$biomarker > lim_sup)]
outliers <- subset(VM00, biomarker < lim_inf | biomarker > lim_sup)
outliers


VM00 = VM00 %>%
  filter(!is.na(pTau181_CSF))
VM00 = VM00 %>% 
  filter(!is.na(pTau181_CSF_Lumip)& pTau181_CSF != "NA" )
VM00$pTau181_CSF_Lumip <- as.numeric(VM00$pTau181_CSF)

VM00 = VM00 %>%
  filter(!is.na(AB142_CSF))
VM00 = VM00 %>% 
  filter(!is.na(AB142_CSF)& AB142_CSF != "NA" )
VM00$AB142_CSF <- as.numeric(VM00$AB142_CSF)
VM00 <- VM00[VM00$AB142_CSF <= 1700, ]


#Filter demografics#
VM00 = VM00 %>%
  filter(!is.na(sex))
VM00$sex<- factor(VM00$sex, levels = c("F", "M"))

VM00 = VM00 %>%
  filter(!is.na(AGE))
VM00 = VM00 %>% 
  filter(!is.na(AGE)& AGE != "NA" )
VM00$AGE <- as.numeric(VM00$AGE)

VM00 = VM00 %>%
  filter(!is.na(DX))
VM00$DX[VM00$DX == "MCI"] <- "CI"
VM00$DX[VM00$DX == "AD"] <- "CI"
VM00$DX[VM00$DX == "CN"] <- "CU"
VM00 <- VM00 %>% drop_na(DX)

VM00 = VM00 %>%
  filter(!is.na(MMSE))
VM00$MMSE <- as.numeric(VM00$MMSE)

VM00 = VM00 %>%
  filter(!is.na(years_of_education))
VM00 = VM00 %>% 
  filter(!is.na(years_of_education)& years_of_education != "NA" )
VM00$years_of_education <- as.numeric(VM00$years_of_education)

VM00 = VM00 %>%
  filter(!is.na(apoee4_alleles))
VM00$APOE4<- factor(VM00$APOE4, levels = c("1", "0"))

#z-score
VM00$AGEz <- scale(VM00$AGE, center = TRUE, scale = TRUE)

#Groups
Years50 <- filter(VM00, AGE > 50) 

CU = Years50 %>%
  filter(DX == "CU") 

CI = Years50 %>%
  filter(DX ==  "CI")

# Abeta positive and negative to ADNI #

CU_AB_neg = CU %>%
  filter (ABETA_bl >= 977)

CU_AB_pos = CU %>%
  filter (ABETA_bl < 977) 

CI_AB_pos = CI %>%
  filter (ABETA_bl < 977)


# Abeta positive and negative to TRIAD #

CU_AB_neg = CU %>%
  filter (AB42_40 >= 0.068) 

CU_AB_pos = CU %>%
  filter (AB42_40 < 0.068) 

CI_AB_neg = CI %>%
  filter (AB42_40 >= 0.068)

CI_AB_pos = CI %>%
  filter (AB42_40 < 0.068)


#ANCOVA
a_aov <- aov(GFAP_zscore ~ DX, data = ADNI)
emm <- emmeans(a_aov, pairwise ~ DX, adjust = "tukey", pbkrtest.limit = 6000)
summary(emm)

#LM
Model <- lm(Synaptic_Biomarker ~ Glial_cell_Biomarker + sex + AGEz, data = Group) 
summary(Model)
conf_intervals <- confint(Model)
conf_intervals <- confint(Model, level = 0.95)
conf_intervals

#Figures 1 and 2
ggplot(Group, aes (x = Glial_cell_Biomarker, y = Synaptic_Biomarker)) +
  geom_point(fill = "#00fa9a" , shape = 21 , size = 10 )+
  xlab("Biomarker (z-score)")+
  ylab("Biomarker (z-score)") +
  theme_classic()+
  theme(text = element_text(size = 24)) +
  geom_smooth (aes (x = Biomarker, y = Biomarker), method = "lm", se = TRUE, color = "black",
               size = 1 )

#LOWESS and figure 3
GAP43z_norm <-lm(GAP43_normalized ~  age_at_mri + sex, data = Years50)
variable_GAP43z <- resid(GAP43z_norm)
mean_GAP43z <- mean(Years50$GAP43_normalized)
Years50$GAP43z_norm  <-  variable_GAP43z + mean_GAP43z

cutoff_ab <- quantile(Years50$AB42_40, 0.068)
subset_df_ab <- subset(Years50, AB42_40 > cutoff_ab & DX_cat_num == 0)
mean_GAP43 <- mean(subset_df_ab$GAP43z_norm)
sd_GAP43 <- sd(subset_df_ab$GAP43z_norm)
Years50$GAP_ztert <- (Years50$GAP43z_norm - mean_GAP43)/sd_GAP43

NG_norm <-lm(Ng ~  AGEz + sex, data = Years50)
variable_NG <- resid(NG_norm)
mean_NG <- mean(Years50$Ng)
Years50$NG_norm <-  variable_NG + mean_NG

cutoff_ab <- quantile(Years50$AB42_40, 0.068)
subset_df_ab <- subset(Years50, AB42_40 > cutoff_ab & DX_cat_num == 0)
mean_ngs <- mean(subset_df_ab$NG_norm)
sd_ngs <- sd(subset_df_ab$NG_norm)
Years50$ng_ztert <- (Years50$NG_norm - mean_ngs)/sd_ngs

a<-subset(CU_AB_neg)
b<-subset(CU_AB_pos)
c<-subset(CI_AB_pos)


plot(lowess(a$sTREMz, a$GAP_z ,f=.8), col="white", xaxt="none", yaxt="none", xlab="",ylab="", cex=0.1, xlim = c(-1.5, 2.5), ylim = c(-1, 2.5))
plot(lowess(a$sTREMz, a$ng_z ,f=.8), col="white", xaxt="none", yaxt="none", xlab="",ylab="", cex=0.1, xlim = c(-1.5, 2.5), ylim = c(-1, 2.5))

plot(lowess(b$sTREMz, b$GAP_z ,f=.8), col="white", xaxt="none", yaxt="none", xlab="",ylab="", cex=0.1, xlim = c(-1.5, 2.5), ylim = c(-1, 2.5))
plot(lowess(b$sTREMz, b$ng_z ,f=.8), col="white", xaxt="none", yaxt="none", xlab="",ylab="", cex=0.1, xlim = c(-1.5, 2.5), ylim = c(-1, 2.5))

plot(lowess(c$sTREMz, c$GAP_z ,f=.8), col="white", xaxt="none", yaxt="none", xlab="",ylab="", cex=0.1, xlim = c(-1.5, 2.5), ylim = c(-1, 2.5))
plot(lowess(c$sTREMz, c$ng_z ,f=.8), col="white", xaxt="none", yaxt="none", xlab="",ylab="", cex=0.1, xlim = c(-1.5, 2.5), ylim = c(-1, 2.5))


title(xlab="CSF sTREM2 (z-score)", cex.lab= 2)  
title(ylab="Synaptic Markers (z-score)", cex.lab= 2)

axis(1, las=1, font=1,cex.axis=1.2)
axis(2, las=2, font=1, cex.axis=1.2)


lines(lowess(a$sTREMz, a$GAP_z,f=0.95, iter= 1500), col="#00ccff", lwd=c(8))
lines(lowess(a$sTREMz, a$ng_z,f=0.95, iter=1500), col="#0000ff", lwd=c(8))

lines(lowess(b$sTREMz, b$GAP_z,f=0.95, iter= 1500), col="#32cd32", lwd=c(8))
lines(lowess(b$sTREMz, b$ng_z,f=0.95, iter=1500), col="#006400", lwd=c(8))

lines(lowess(c$sTREMz, c$GAP_z,f=0.95, iter= 1500), col="#ff69b4", lwd=c(8))
lines(lowess(c$sTREMz, c$ng_z,f=0.95, iter=1500), col="#c21e56", lwd=c(8))



#Mediation

model.0 <- lm(Synaptic_Biomarker ~ CSFpTau181 + AGEz + SEX , data= CU_AB_neg)
summary(model.0)

model.M <- lm(Glial_cell_Biomarker ~ CSFpTau181 + AGEz + SEX, data= CU_AB_neg)
summary(model.M)

model.Y <- lm(Synaptic_Biomarker ~ CSFpTau181 + Glial_cell_Biomarker + AGEz + SEX , data= CU_AB_neg)
summary(model.Y)

library(mediation)
results <- mediate(model.M, model.Y, treat='Glial cell Biomarker', mediator='CSFpTAU',
                   boot=TRUE, sims=500)

