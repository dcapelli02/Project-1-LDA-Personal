#### PROVA LETTURA FILE SAS ####

library(haven)
library(sas7bdat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(GGally)
library(broom)
library(glmnet)
alz <- read_sas("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.sas7bdat")
#alz <- read.csv("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.csv")

head(alz)
boxplot(alz$cdrsb0, alz$age)
summary(alz)

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)

summary(alz)

## Now we can create the longitudinal data frame

alz_df <- data.frame(alz)

alz_long <- alz_df %>%
  pivot_longer(
    # 1️⃣ Seleziona tutte le colonne che iniziano con cdrsb, abpet o taupet,
    #     seguite da un numero (0–6)
    cols = matches("^(cdrsb|abpet|taupet)\\d+$"),
    
    # 2️⃣ Crea due nuove colonne:
    #     - ".value" → il prefisso (cdrsb, abpet, taupet)
    #     - "year"   → il numero alla fine
    names_to = c(".value", "year"),
    
    # 3️⃣ Divide i nomi in due parti: (prefisso, numero)
    names_pattern = "(cdrsb|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                            # converte in numerico
    sample = factor(rep(1:nrow(alz_df), each = 7)) # ID per ogni paziente
  )


### NOW CREATE A RESTRRICTED VERSION TO WORK ON ###

casual <- sample(1:length(alz$patid), 50)


## Now we can start by looking at random values for the mean and see
## if we can work on the mean and so on

alz_rist <- alz[casual, ]
cdrsb_rist <- alz_rist[, c(11:17)]

## Now I need to create the data frame in order to work with ggplot
alz_rist_df <- data.frame(alz_rist)
alz_rist_long_cdrsb <- pivot_longer(alz_rist_df, 
                                    cols = c("cdrsb0",  "cdrsb1",  "cdrsb2", 
                                             "cdrsb3",  "cdrsb4",  "cdrsb5",  "cdrsb6"),
                                    names_to = "year",
                                    values_to = "cdrsb")
alz_rist_long_cdrsb$year <- rep(0:6, length(alz_rist_df[,1]))
alz_rist_long_cdrsb$sample <- rep(1:length(alz_rist_df[,1]), each = 7)
alz_rist_long_cdrsb$sample <- as.factor(alz_rist_long_cdrsb$sample)



ggplot(alz_rist_long_cdrsb, aes(x = year, y = cdrsb, group = patid, 
                                color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 0.5, show.legend = FALSE) +
  theme_bw()



## Now we can create the longitudinal data frame

alz_rist_long <- alz_rist_df %>%
  pivot_longer(
    # 1️⃣ Seleziona tutte le colonne che iniziano con cdrsb, abpet o taupet,
    #     seguite da un numero (0–6)
    cols = matches("^(cdrsb|abpet|taupet)\\d+$"),
    
    # 2️⃣ Crea due nuove colonne:
    #     - ".value" → il prefisso (cdrsb, abpet, taupet)
    #     - "year"   → il numero alla fine
    names_to = c(".value", "year"),
    
    # 3️⃣ Divide i nomi in due parti: (prefisso, numero)
    names_pattern = "(cdrsb|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                            # converte in numerico
    sample = factor(rep(1:nrow(alz_rist_df), each = 7)) # ID per ogni paziente
  )


## Now we can go to the visual part

ggplot(alz_rist_long, aes(x = year, y = cdrsb, group = patid, 
                                color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 0.5, show.legend = FALSE) +
  theme_bw() +
  labs(title = "CDRSB against year of measuring")


ggplot(alz_rist_long, aes(x = year, y = abpet, group = patid, 
                                color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 0.5, show.legend = FALSE) +
  theme_bw() +
  labs(title = "ABPET against year of measuring")

ggplot(alz_rist_long, aes(x = year, y = taupet, group = patid, 
                          color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 0.5, show.legend = FALSE) +
  theme_bw() +
  labs(title = "TAUPET against year of measuring")


## Now different combinations

ggplot(alz_rist_long, aes(x = year, y = cdrsb, group = patid, 
                          color = trial, show.legend = FALSE)) + 
  geom_line(alpha = 0.5, show.legend = FALSE) +
  theme_bw() +
  labs(title = "CDRSB against year of measuring")

ggplot(alz_rist_long, aes(x = abpet, y = cdrsb, group = patid, 
                          color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 0.5, show.legend = FALSE) +
  theme_bw() +
  labs(title = "CDRSB against year of measuring")




ggplot(alz_rist_long, aes(x = year, y = cdrsb, group = sex, color = sex)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal()

ggplot(alz_rist_long, aes(x = year, y = cdrsb, group = wzc, color = wzc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal()

ggplot(alz_rist_long, aes(x = year, y = cdrsb)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal()

ggplot(alz_long, aes(x = year, y = cdrsb)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal()

## Qua sembra esserci vagamente una varianza costante, quindi posso usare il 
## semivariogram approach (ma è quella la variance?)

### DA FARE: SEMIVARIOGRAMMA ###ù

### CORRELATION STRUCTURE ###




### VISUALIZZAZIONE DATI ###
ggpairs(alz_rist_long, columns = c(3:10, 19:21))
pairs(alz_rist_long[, c(3:10, 19:21)])
