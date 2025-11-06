## Libraries
library(haven)
library(sas7bdat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(GGally)
library(corrplot)
library(nlme)
library(lmtest)
library(DescTools)

## Import and fix the data
alz <- read_sas("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.sas7bdat")
#alz <- read.csv("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.csv")

head(alz)
summary(alz)

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)

## Create baseline values
alz$ab_base <- alz$abpet0
alz$tau_base <- alz$taupet0
alz$cdrsb_base <- alz$cdrsb0

summary(alz)

## Create longitudinal dataset
alz_df <- data.frame(alz)

alz_long <- alz_df %>%
  pivot_longer(
    # 1️⃣ Seleziona tutte le colonne che iniziano con cdrsb, abpet o taupet,
    #     seguite da un numero (0–6)
    cols = matches("^(bprs|cdrsb|abpet|taupet)\\d+$"),
    
    # 2️⃣ Crea due nuove colonne:
    #     - ".value" → il prefisso (cdrsb, abpet, taupet)
    #     - "year"   → il numero alla fine
    names_to = c(".value", "year"),
    
    # 3️⃣ Divide i nomi in due parti: (prefisso, numero)
    names_pattern = "(bprs|cdrsb|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                            # converte in numerico
    sample = factor(rep(1:nrow(alz_df), each = 7)) # ID per ogni paziente
  )

## Discretize variables

## Maybe any 5 years?
alz_long$age_disc <- (alz_long$age %/% 5) * 5
alz_long$age_disc <- as.factor(alz_long$age_disc)

## bmi any 4
alz_long$bmi_disc <- (alz_long$bmi %/% 4) * 4
alz_long$bmi_disc <- as.factor(alz_long$bmi_disc)

## inkomen any 500
alz_long$inkomen_disc <- (alz_long$inkomen %/% 500) * 500
alz_long$inkomen_disc <- as.factor(alz_long$inkomen_disc)

## adl any 5
alz_long$adl_disc <- (as.numeric(alz_long$adl) %/% 5) * 5
alz_long$adl_disc <- as.factor(alz_long$adl_disc)

## cdrsb any 5
alz_long$cdrsb_disc <- (alz_long$cdrsb %/% 5) * 5
alz_long$cdrsb_disc <- as.factor(alz_long$cdrsb_disc)

## abpet any 0.2
alz_long$abpet_disc <- (alz_long$abpet %/% 0.2) * 0.2
alz_long$abpet_disc <- as.factor(alz_long$abpet_disc)

## taupet any 0.2
alz_long$taupet_disc <- (alz_long$taupet %/% 0.2) * 0.2
alz_long$taupet_disc <- as.factor(alz_long$taupet_disc)

## adl any 5
alz$adl_disc <- (as.numeric(alz$adl) %/% 5) * 5
alz$adl_disc <- as.factor(alz$adl_disc)

## year discrete
alz_long$year_seq <- ave(alz_long$year, alz_long$sample, FUN = function(x) as.integer(factor(x)))



#### LINEAR MIXED MODEL ####

#lme_model_3 <- lme(bprs ~ age + adl_disc + wzc + ab_base + 
#                             (sex + inkomen + adl_disc + edu + bmi + cdrsb_base)*year,
#                   data = alz_long,
#                   random = ~ year | sample,
#                   correlation = corAR1(form = ~ year | sample),
#                   weights = varIdent(form = ~ 1 | year),
#                   na.action = na.omit
#)

## No stavolta sono molto diversi...
## Secondo me possiamo evitare di tenere i pesi nell'approccio lm...

## Start by zero for the linear mixed model
## A good starting point is to keep any variable and an unstructured covariance structure

lme_model_1_full <- lme(bprs ~ (sex + age + edu + bmi + inkomen + job +
                                  adl_num + wzc + cdrsb_base +
                                  ab_base + tau_base) * year,
                        data = alz_long,
                        random = ~ year_seq | sample,
                        correlation = corSymm(form = ~ year_seq | sample),
                        weights = varIdent(form = ~ 1 | year_seq),
                        na.action = na.omit)

## MODELLO TROPPO PESANTE DA ESEGUIRE: non arriva a convergenza
## Come fare per ridurre il tutto?
## Quale potrebbe essere una buona idea?

## Secondo me una buona idea è partire da struttura di covarianza AR(1)

lme_model_1_full_AR <- lme(bprs ~ (sex + age + edu + bmi + inkomen + job +
                                  adl_num + wzc + cdrsb_base +
                                  ab_base + tau_base) * year,
                        data = alz_long,
                        random = ~ year_seq | sample,
                        correlation = corAR1(form = ~ year_seq | sample),
                        weights = varIdent(form = ~ 1 | year_seq),
                        na.action = na.omit)

## Questo ha raggiunto convergenza

summary(lme_model_1_full_AR)

## A livello di baseline: toglierei (forse sex), edu, inkomen, ab (da vedere),
## tau

## A livello di intercetta: toglierei forse sex, edu, bmi, wzc, ab e tau

## Però sono cose che vanno testate?


