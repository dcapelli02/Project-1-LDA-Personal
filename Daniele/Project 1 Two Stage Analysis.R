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

#### 2 STAGE MODEL ####

### STAGE 1 MODEL ###

## Start by fitting any time a linear regression

count <- 0
coeff_stage1 <- matrix(nrow = length(alz$patid), ncol = 2)
sigma_stage_1 <- matrix(nrow = length(alz$patid), ncol = 1)

# We assume basically that bprs_i = beta_0i + beta_1i * year_i + eps_i
# We should define a structure for the errors
# Usually it is reasonable to consider esp_i ~ N(0, Sigma)
# and Sigma = sigma^2 * I

for (i in 1:length(alz$patid)) {
  mod_prova <- lm(bprs ~ year, 
                  data = alz_long[c((7*count + 1) : (7*count + 7)), ])
  coeff_stage1[i, ] <- mod_prova$coefficients
  sigma_prov <- sqrt(sum(residuals(mod_prova)^2) / df.residual(mod_prova))
  sigma_stage_1[i] <- sigma_prov
  count = count + 1
}

# Here the stage 1 model is over

### STAGE 2 MODEL ###

## First of all store the results
beta0 <- coeff_stage1[, 1]
beta1 <- coeff_stage1[, 2]

## Then we can perform a regression

## Here usually we consider the matrix of errors D to be unstructured
## so we need to use gls and not lm

## In order to use it we need to construct a new dataframe
## (otherwise gls does not work properly)

alz$beta0 <- beta0
alz$beta1 <- beta1

alz_long_beta <- alz %>%
  mutate(id = row_number()) %>%
  select(id, sex, job, age, inkomen, adl_disc, wzc,
         ab_base, tau_base, edu, bmi, cdrsb_base, beta0, beta1) %>%
  pivot_longer(cols = c(beta0, beta1),
               names_to = "param",
               values_to = "beta")

alz_long_beta$param <- as.factor(alz_long_beta$param)

model_2stage <- gls(beta ~ param + sex + job + age + inkomen +
                      adl_disc + wzc + ab_base + tau_base +
                      edu + bmi + cdrsb_base,
                    data = alz_long_beta,
                    correlation = corSymm(form = ~ 1 | id),
                    weights     = varIdent(form = ~ 1 | param),
                    method      = "ML",
                    na.action   = na.exclude)

model_2stage$coefficients

## Then I should reduce the model?
## But if I do like this I do not take into consideration that I would like
## to exclude different covariates for beta0 and beta1...

## Actually, if we proceed in this way we need to take into consideration
## the same regressors both for beta0 and beta1

## So I prefer to run a lm regression for both the variables, maybe by simply
## considering a different weight for the estimations carried on by the 
## fact that we are considering estimates from a model

## Fix the NA values
# sostituisci NA/0 con un piccolo valore, ad esempio
sigma_stage_1_clean <- ifelse(is.na(sigma_stage_1) | sigma_stage_1 == 0,
                              min(sigma_stage_1[sigma_stage_1 > 0], na.rm=TRUE),
                              sigma_stage_1)

weights_lm <- 1 / sigma_stage_1_clean^2

model_beta0 <- lm(beta0 ~ sex + job + age + inkomen + 
                    adl_disc + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base,
                  #weights = weights_lm,
                  data = alz)

model_beta1 <- lm(beta1 ~ sex + job + age + inkomen + 
                    adl_disc + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base,
                  #weights = weights_lm,
                  data = alz)

## So these are the starting points
## We want to reduce the models

step_beta_0 <- step(model_beta0)

model_beta0_final <- step_beta_0$call

## Se servisse
# lm(formula = beta0 ~ job + age + inkomen + adl_disc + wzc + ab_base + 
# tau_base + edu + bmi + cdrsb_base, data = alz)

step_beta_1 <- step(model_beta1)

model_beta1_final <- step_beta_1$call

## Se servisse
# lm(formula = beta1 ~ adl_disc + wzc + edu + cdrsb_base, data = alz)

coeff_beta0_final <- model_beta0_final$coefficients
coeff_beta1_final <- model_beta1_final$coefficients



#### PROVA CON TRIAL ####

## Stage 1 the same

### STAGE 2 MODEL ###

model_beta0_trial <- lm(beta0 ~ sex + job + age + inkomen + 
                    adl_num + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base + trial,
                  #weights = weights_lm,
                  data = alz)

model_beta1_trial <- lm(beta1 ~ sex + job + age + inkomen + 
                    adl_num + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base + trial,
                  #weights = weights_lm,
                  data = alz)

step0_trial <- step(model_beta0_trial)
step1_trial <- step(model_beta1_trial)

model_beta0_trial_final <- step0_trial$call
model_beta1_trial_final <- step1_trial$call



### E CON ADL CONTINUA? ###

model_beta0_num <- lm(beta0 ~ sex + job + age + inkomen + 
                          adl_num + wzc + ab_base + tau_base +
                          edu + bmi + cdrsb_base,
                        #weights = weights_lm,
                        data = alz)

model_beta1_num <- lm(beta1 ~ sex + job + age + inkomen + 
                          adl_num + wzc + ab_base + tau_base +
                          edu + bmi + cdrsb_base,
                        #weights = weights_lm,
                        data = alz)

step0_num <- step(model_beta0_num)
step1_num <- step(model_beta1_num)

model_beta0_num <- step0_num$call
model_beta1_num <- step1_num$call

#> step0_num$call
#lm(formula = beta0 ~ sex + job + age + adl_num + wzc + ab_base + 
#     bmi + cdrsb_base, data = alz)
#> step1_num$call
#lm(formula = beta1 ~ wzc + edu + cdrsb_base, data = alz)

