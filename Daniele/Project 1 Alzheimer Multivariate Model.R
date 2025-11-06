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
library(effects)
library(ggeffects)
library(MASS)

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



#### MULTIVARIATE MODEL ####

## Start by fitting one of the most general model we can think of

## Maybe try to use heterogeneous AR(1) as a starting point

mult_model_1 <- gls(
  bprs ~ (age + edu + bmi + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
    sex + job) * year ,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

mult_model_1$coefficients

## This is our benchmark model

## Now we would like to reduce the mean structure

summary(mult_model_1)

### ALTERNATIVA: computationally long
#library(MASS)
#stepAIC(mult_model_1)

mult_model_final <- gls(
  bprs ~ age + edu + bmi + inkomen + adl_disc + wzc + cdrsb_base + 
    ab_base + tau_base + sex + job + year + age:year + edu:year +
    inkomen:year + wzc:year + cdrsb_base:year + ab_base:year +
    sex:year,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

## DIFFERENCE COVARIANCE STRUCTURE ##

mult_model_final_exp <- gls(
  bprs ~ age + edu + bmi + inkomen + adl_disc + wzc + cdrsb_base + 
    ab_base + tau_base + sex + job + year + age:year + edu:year +
    inkomen:year + wzc:year + cdrsb_base:year + ab_base:year +
    sex:year,
  correlation = corExp(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

## Confronto

anova(mult_model_final, mult_model_final_exp)

## I can pick the one I want so I will stick to my choice of AR1

### ALTRA PROVA ###

## What does it happen if we stick to the model I inferred from the EDA?

model_nostro <- gls(
  bprs ~ age + inkomen + bmi + tau_base +
    (bmi + job + adl_disc + wzc + sex + ab_base)*year - bmi,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

## LR Test to compare the results

lrtest(mult_model_final, model_nostro)
anova(mult_model_final, model_nostro)

## No we should keep the model we found


#### VECCHIA COPIA DI BACKUP ####

## For example remove sex from intercept

mult_model_2 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_1, mult_model_2)

## High p-value, so remove it. Now try to remove edu*year

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - edu*year,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## No, keep it, try to remove tau*year

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - tau_base*year,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try to remove cdrsb

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - cdrsb_base*year,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try to remove bmi from slope

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - bmi*year,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try to remove age

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - age*year,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try with inkomen

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - inkomen*year,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try with adl_disc

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - adl_disc*year,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try with ab_base

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - ab_base*year,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try with the baselines now

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - age,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try bmi

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - bmi,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try edu

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - inkomen,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try adl_disc

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - adl_disc,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try with wzc

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - wzc,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it, try with cdrbs

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - cdrsb_base,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it try with ab_base

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - ab_base,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it try with tau

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - tau_base,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## keep it try with job

mult_model_3 <- gls(
  bprs ~ (age + bmi + edu + inkomen + adl_disc + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year -
    sex - job,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_2, mult_model_3)

## Ok so this is our final model we did not reduce particularly the mean strcuture, why??
