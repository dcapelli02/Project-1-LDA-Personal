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
alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))

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

## Voglio un attimo studiare i random effects
random_eff_1 <- ranef(lme_model_1_full_AR)
plot(random_eff_1[, 2])

summary(lme_model_1_full_AR)

## A livello di baseline: toglierei (forse sex), edu, inkomen, ab (da vedere),
## tau

## A livello di intercetta: toglierei forse sex, edu, bmi, wzc, ab e tau

## Però sono cose che vanno testate?

lme_model_1_full_naive <- lme(bprs ~ (sex + age + edu + bmi + inkomen + job +
                                        adl_num + wzc + cdrsb_base +
                                        ab_base + tau_base) * year,
                              data = alz_long,
                              random = ~ year_seq | sample,
                              correlation = NULL,
                              weights = varIdent(form = ~ 1 | year_seq),
                              na.action = na.omit)

AIC(lme_model_1_full_AR, lme_model_1_full_naive)

acf(residuals(lme_model_1_full_naive), type = "response")

anova(lme_model_1_full_AR, lme_model_1_full_naive)




#### PROVA A CASO ####

## Cioe tipo parto cosi
## Da quanto ho capito non è fondamentale specificare perfettamente la matrice
## di covarianza qua (ma cerca di capire meglio)

mod_base_lm <- lm(bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
                      adl_num + wzc + cdrsb_base +
                      ab_base + tau_base) * year,
            data = alz_long,
            na.action = na.exclude)

## Poi studio i residui

alz_long$primi_residui <- residuals(mod_base_lm)

## Forse meglio visualizzarne solo alcuni?
ggplot(alz_long[c(1:(7*15)), ], aes(x = year, y = primi_residui, color = sample)) +
  geom_line(show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2)
  #stat_summary(fun = mean, geom = "line", size = 1)

ggplot(alz_long, aes(x = year, y = primi_residui, color = sample)) +
  geom_line(show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2)

## Sembra che la varianza dei residui aumenti lungo il tempo
## Proviamo a vedere se riusciamo a farne il plot

ggplot(alz_long, aes(x = year, y = primi_residui)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "Residuals linked to lm")
  

hist(alz_long$primi_residui[alz_long$n_obs_data == 1])
hist(alz_long$primi_residui[alz_long$n_obs_data == 2])
hist(alz_long$primi_residui[alz_long$n_obs_data == 3])
hist(alz_long$primi_residui[alz_long$n_obs_data == 4])
hist(alz_long$primi_residui[alz_long$n_obs_data == 5])
hist(alz_long$primi_residui[alz_long$n_obs_data == 6])
hist(alz_long$primi_residui[alz_long$n_obs_data == 7])

## Si anche da qua sembra che la varianza aumenti nel tempo
## Ma cosa mi vuol dire questo?

acf(alz_long$primi_residui, na.action = na.exclude)

## Studiamo anche i quadrati dei residui
alz_long$residui_squared <- (alz_long$primi_residui)^2

## Forse meglio visualizzarne solo alcuni?
ggplot(alz_long[c(1:(7*15)), ], aes(x = year, y = residui_squared, color = sample)) +
  geom_line(show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2)
#stat_summary(fun = mean, geom = "line", size = 1)

ggplot(alz_long, aes(x = year, y = residui_squared, color = sample)) +
  geom_line(show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2)

## Vediamo che crescono nel tempo --> sufficiente per dire che dobbbiamo inserire
## random slope nel tempo?


## Poi fitto lme con i random effects

## Primo modello con matrice di covarianza diagonale

lme_primo_modello <- lme(bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
                                   adl_num + wzc + cdrsb_base +
                                   ab_base + tau_base) * year,
                         data = alz_long,
                         random = ~ year_seq | sample,
                         #correlation = corAR1(form = ~ year_seq | sample),
                         #weights = varIdent(form = ~ 1 | year_seq),
                         na.action = na.exclude)

## Studiamo un attimo i residui
alz_long$residui_primo_modello <- residuals(lme_primo_modello, type = "response")

## Forse meglio visualizzarne solo alcuni?
ggplot(alz_long[c(1:(7*15)), ], aes(x = year, y = residui_primo_modello, color = sample)) +
  geom_line(show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2)
#stat_summary(fun = mean, geom = "line", size = 1)

ggplot(alz_long, aes(x = year, y = residui_primo_modello, color = sample)) +
  geom_line(show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2)

## Varianza nel tempo
# Varianza degli effetti casuali
vc <- VarCorr(lme_primo_modello)

# Estrai le varianze e covarianza
var_intercept <- as.numeric(vc[1, "Variance"])
var_slope <- as.numeric(vc[2, "Variance"])
cov_int_slope <- as.numeric(vc[2, "Corr"]) * sqrt(var_intercept * var_slope)

# Varianza residua
var_resid <- as.numeric(vc[nrow(vc), "Variance"])
# Crea una sequenza di valori di year_seq
anni <- sort(unique(alz_long$year_seq))

# Calcola la varianza marginale per ciascun anno
var_marginale <- sapply(anni, function(t) {
  var_intercept + 2 * t * cov_int_slope + t^2 * var_slope + var_resid
})

# Crea un data frame per il plot
df_var <- data.frame(year_seq = anni, varianza = var_marginale)

ggplot(df_var, aes(x = year_seq, y = varianza)) +
  geom_line(color = "steelblue", size = 1.2) +
  labs(title = "Varianza marginale predetta di BPRS nel tempo",
       x = "Anno (year_seq)",
       y = "Varianza marginale") +
  theme_minimal()


## Secondo modello --> troppo pesante questo qua

lme_secondo_modello <- lme(bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
                                   adl_num + wzc + cdrsb_base +
                                   ab_base + tau_base) * year,
                         data = alz_long,
                         random = ~ year_seq | sample,
                         correlation = corAR1(form = ~ year_seq | sample),
                         #weights = varIdent(form = ~ 1 | year_seq),
                         na.action = na.exclude)

## Studiamo un attimo i residui
alz_long$residui_secondo_modello <- residuals(lme_secondo_modello, type = "response")

## Forse meglio visualizzarne solo alcuni?
ggplot(alz_long[c(1:(7*15)), ], aes(x = year, y = residui_secondo_modello, color = sample)) +
  geom_line(show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2)
#stat_summary(fun = mean, geom = "line", size = 1)

ggplot(alz_long, aes(x = year, y = residui_secondo_modello, color = sample)) +
  geom_line(show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2)

anova(lme_primo_modello, lme_secondo_modello)

## Non sembra cambiare molto, però forse è ancora meglio il primo modello qua...
## In fin dei conti non è così insensata come idea quella di non avere serial
## correlation forse

## Poi studio possibile matrice Sigma

## Poi studio come ridurre covarianza

## Poi studio come ridurre la media











#### DO NOT LOOK HERE ####






lme_terzo_modello <- lme(bprs ~ (sex + age + edu + bmi + inkomen + job +
                                   adl_num + wzc + cdrsb_base +
                                   ab_base + tau_base) * year,
                         data = alz_long,
                         random = ~ year_seq | sample,
                         #correlation = corAR1(form = ~ year_seq | sample),
                         weights = varExp(form = ~ year_seq),
                         na.action = na.exclude)


# Estrai parametri di varianza
vc <- VarCorr(lme_terzo_modello)

var_intercept <- as.numeric(vc[1, "Variance"])
var_slope <- as.numeric(vc[2, "Variance"])
cov_int_slope <- as.numeric(vc[2, "Corr"]) * sqrt(var_intercept * var_slope)

# Varianza di base e parametro varExp
base_var <- as.numeric(vc[nrow(vc), "Variance"])
delta <- coef(lme_terzo_modello$modelStruct$varStruct, unconstrained = FALSE)  # parametro varExp

# Calcola la varianza residua per ogni anno in modo continuo
anni <- sort(unique(alz_long$year_seq))
resid_var_by_year <- base_var * exp(2 * delta * anni)

# Calcola la varianza marginale prevista (intercetta + slope + residuo)
var_marginale <- sapply(anni, function(t) {
  var_intercept + 2 * t * cov_int_slope + t^2 * var_slope + resid_var_by_year[which(anni == t)]
})

df_var <- data.frame(
  year_seq = anni,
  var_residua = resid_var_by_year,
  var_marginale = var_marginale
)

# Grafico
library(ggplot2)
ggplot(df_var, aes(x = year_seq)) +
  geom_line(aes(y = var_residua, color = "Varianza residua (dal modello)"), size = 1) +
  geom_line(aes(y = var_marginale, color = "Varianza marginale predetta"), size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Varianze di BPRS nel tempo (modello con varExp)",
    x = "Anno (year_seq)",
    y = "Varianza di BPRS",
    color = "Tipo di varianza"
  ) +
  scale_color_manual(values = c("Varianza residua (dal modello)" = "steelblue",
                                "Varianza marginale predetta" = "darkred"))





















vc <- VarCorr(lme_terzo_modello)

var_intercept <- as.numeric(vc[1, "Variance"])
var_slope <- as.numeric(vc[2, "Variance"])
cov_int_slope <- as.numeric(vc[2, "Corr"]) * sqrt(var_intercept * var_slope)

# Estrai varianze residue specifiche per anno
base_var <- as.numeric(vc[nrow(vc), "Variance"])
scales <- coef(lme_terzo_modello$modelStruct$varStruct, unconstrained = FALSE)
names(scales) <- gsub(".*=", "", names(scales))
resid_var_by_year <- base_var * c(1, scales[order(names(scales))]^2)  # include baseline

# Calcolo varianza marginale corretta
anni <- sort(unique(alz_long$year_seq))
var_marginale <- sapply(seq_along(anni), function(i) {
  t <- anni[i]
  var_intercept + 2 * t * cov_int_slope + t^2 * var_slope + resid_var_by_year[i]
})

df_var <- data.frame(year_seq = anni, varianza = var_marginale)

ggplot(df_var, aes(x = year_seq, y = varianza)) +
  geom_line(color = "steelblue", size = 1.2) +
  theme_minimal() +
  labs(title = "Varianza marginale predetta (corretta) di BPRS nel tempo",
       x = "Anno (year_seq)",
       y = "Varianza marginale")


## Stima varianza predetta

pred <- predict(lme_terzo_modello, level = 0)  # predizione marginale
res <- residuals(lme_terzo_modello, level = 0)

alz_long$pred <- pred
alz_long$res <- res

# Varianza totale prevista (pred + residui)
var_pred_by_year <- alz_long %>%
  group_by(year_seq) %>%
  summarise(var_pred = var(pred + res, na.rm = TRUE))

ggplot(var_pred_by_year, aes(x = year_seq, y = var_pred)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Varianza totale predetta (empirica dal modello)",
       x = "Anno", y = "Varianza")


var_by_year <- alz_long %>%
  group_by(year_seq) %>%
  summarise(var_empirica = var(bprs, na.rm = TRUE))

df_var <- data.frame(year_seq = anni, var_predetta = var_marginale)

df_confronto <- left_join(var_by_year, df_var, by = "year_seq")

df_confronto <- df_confronto %>%
  mutate(var_empirica_z = scale(var_empirica),
         var_predetta_z = scale(var_predetta))

alz_long$resid <- residuals(lme_terzo_modello, level = 0)

var_resid_by_year <- alz_long %>%
  group_by(year_seq) %>%
  summarise(var_resid = var(resid, na.rm = TRUE))

ggplot(df_confronto, aes(x = year_seq)) +
  geom_line(aes(y = var_empirica, color = "Empirica"), size = 1.2) +
  geom_point(aes(y = var_empirica, color = "Empirica"), size = 2) +
  geom_line(aes(y = var_predetta, color = "Predetta"), size = 1.2, linetype = "dashed") +
  geom_point(aes(y = var_predetta, color = "Predetta"), size = 2, shape = 17) +
  theme_minimal() +
  labs(
    title = "Varianza di BPRS nel tempo: dati vs modello",
    x = "Anno (year_seq)",
    y = "Varianza di BPRS",
    color = "Tipo di varianza"
  ) +
  scale_color_manual(values = c("Empirica" = "black", "Predetta" = "steelblue")) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )








## 1️⃣ Varianza empirica totale nei dati
var_by_year <- alz_long %>%
  group_by(year_seq) %>%
  summarise(var_empirica = var(bprs, na.rm = TRUE))

## 2️⃣ Varianza empirica dei residui del modello
alz_long$resid <- residuals(lme_terzo_modello, level = 0)

var_resid_by_year <- alz_long %>%
  group_by(year_seq) %>%
  summarise(var_resid = var(resid, na.rm = TRUE))

## 3️⃣ Varianza marginale predetta dal modello
vc <- VarCorr(lme_terzo_modello)

var_intercept <- as.numeric(vc[1, "Variance"])
var_slope <- as.numeric(vc[2, "Variance"])
cov_int_slope <- as.numeric(vc[2, "Corr"]) * sqrt(var_intercept * var_slope)

# Varianze residue specifiche per anno (se hai varIdent)
base_var <- as.numeric(vc[nrow(vc), "Variance"])
scales <- coef(lme_terzo_modello$modelStruct$varStruct, unconstrained = FALSE)
names(scales) <- gsub(".*=", "", names(scales))
resid_var_by_year <- base_var * c(1, scales[order(names(scales))]^2)

anni <- sort(unique(alz_long$year_seq))
var_marginale <- sapply(seq_along(anni), function(i) {
  t <- anni[i]
  var_intercept + 2 * t * cov_int_slope + t^2 * var_slope + resid_var_by_year[i]
})

df_var <- data.frame(year_seq = anni, var_predetta = var_marginale)

## 4️⃣ Unisci tutto
df_confronto <- full_join(var_by_year, var_resid_by_year, by = "year_seq") %>%
  full_join(df_var, by = "year_seq") %>%
  tidyr::pivot_longer(cols = -year_seq, names_to = "tipo", values_to = "varianza")

## 5️⃣ Grafico comparativo
ggplot(df_confronto, aes(x = year_seq, y = varianza, color = tipo)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Confronto tra varianze di BPRS nel tempo",
    subtitle = "Empirica (dati), Residui (dal modello) e Predetta (teorica)",
    x = "Anno (year_seq)",
    y = "Varianza di BPRS",
    color = "Tipo di varianza"
  ) +
  scale_color_manual(
    values = c(
      "var_empirica" = "black",
      "var_resid" = "darkred",
      "var_predetta" = "steelblue"
    ),
    labels = c(
      "Varianza empirica (dati)",
      "Varianza residui (dal modello)",
      "Varianza marginale predetta"
    )
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom"
  )
