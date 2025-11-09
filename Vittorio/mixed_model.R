library(nlme)
library(ggplot2)
library(dplyr)
library(ggrepel)  
library(broom)
library(nlme)

# Modello OLS 

#Preliminary Mean Structure
mod_ols <- lm(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  na.action = na.exclude
)

# Modello lme 

#Preliminary Random-eﬀects Structure
# Residui OLS
df_ols <- alz_long %>%
  mutate(
    resid = resid(mod_ols),        # residui dal modello OLS
    resid2 = resid(mod_ols)^2,              # residui quadratici
    sample = as.factor(sample)     # ID paziente come fattore
  ) %>%
  select(sample, year, resid, resid2) %>%
  filter(is.finite(resid), is.finite(year))

# Seleziona 30 soggetti con più osservazioni
set.seed(123)  # per riproducibilità
top30 <- sample(unique(df_ols$sample), 30)

df_ols_30 <- df_ols %>% filter(sample %in% top30)

# Grafico profili residui
ggplot(df_ols_30, aes(x = year, y = resid, group = sample, color = sample)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  labs(title = "Profili residui OLS (30 pazienti)",
       x = "Anno", y = "Residui OLS") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",        # mostra legenda con i colori dei pazienti
    strip.text = element_text(face = "bold")
  )

# ---> RANDOM INTERCEPT

# Grafico con lisciamento unico
ggplot(df_ols, aes(x = year, y = resid2)) +
  geom_point(alpha = 0.15, color = "grey40", size = 1) +
  geom_smooth(se = FALSE, method = "loess",
              color = "#1f77b4", linewidth = 1.3) +
  labs(title = "Smooth residui OLS",
       x = "Anno", y = "Residui OLS") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

#leggermente inclinata quindi ---> RANDOM SLOPES

#procediamo quindi con la stima del modello mixed effects

mod_lme <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  random = ~ year | sample,
  na.action = na.exclude
)

# Estrazione D e sigma^2 dall'LME
# Prendiamo la matrice dei random effects (D) e l'errore residuo
vc <- VarCorr(mod_lme)

# Varianze
d00 <- as.numeric(vc["(Intercept)", "Variance"])
d11 <- as.numeric(vc["year", "Variance"])
sigma2 <- as.numeric(vc["Residual", "Variance"])

# Correlazione intercetta-slope
rho <- as.numeric(vc["year", "Corr"])

# Covarianza = correlazione * sqrt(var_intercetta * var_slope)
d01 <- rho * sqrt(d00 * d11)

# Matrice D
D <- matrix(c(d00, d01, d01, d11), nrow = 2)

# Componenti di D
d00 <- D[1, 1]                 # Var(random intercept)
d11 <- D[2, 2]                 # Var(random slope year)
d01 <- D[1, 2]                 # Cov(intercept, year)

# Griglia di anni su cui valutare la funzione di varianza
t_seq <- sort(unique(model.frame(mod_ols)$year))

var_fun <- d00 + 2 * d01 * t_seq + d11 * (t_seq^2) + sigma2

df_var <- data.frame(
  year = t_seq,
  resid2 = var_fun,
  modello = "Fitted variance (LME)"
)

# Grafico: linea continua per OLS smussata, linea tratteggiata per varianza fittata LME
ggplot() +
  # curva smussata dei residui OLS (empirica)
  geom_smooth(data = df_ols, aes(x = year, y = resid2),
              method = "loess", se = FALSE,
              linewidth = 1.3, color = "#1f77b4") +
  # curva teorica della varianza dal modello LME
  geom_line(data = df_var, aes(x = year, y = resid2),
            linewidth = 1.2, linetype = "dashed", color = "#d62728") +
  # etichette finali direttamente sulle curve
  geom_text_repel(
    data = df_ols %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "OLS (smussato)"),
    color = "#1f77b4", hjust = 0, nudge_x = 1, size = 4
  ) +
  geom_text_repel(
    data = df_var %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "LME (varianza fittata)"),
    color = "#d62728", hjust = 0, nudge_x = 1, size = 4
  ) +
  labs(
    x = "Anno",
    y = "Varianza / residui al quadrato"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
#curve molto vicine, sembra fare molto bene



#PROVA
#modello con solo random intercept e non random slopes
#modello con solo solo intercetta casuale

# 1. Modello OLS
mod_ols <- lm(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long
)

# 2. Modello LME con solo intercetta casuale
mod_lme1 <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  random = ~ 1 | sample,
  na.action = na.exclude
)


# 3. Varianza fittata dal modello LME (solo intercetta casuale)
vc <- VarCorr(mod_lme1)

d00 <- as.numeric(vc["(Intercept)", "Variance"])   # Var(random intercept)
sigma2 <- as.numeric(vc["Residual", "Variance"])   # Var(residuo)

var_fun <- d00 + sigma2   # costante nel tempo

df_var <- data.frame(
  year = sort(unique(df_ols$year)),
  resid2 = var_fun,
  modello = "Fitted variance (LME)"
)

# 5. Grafico: curva smussata OLS vs linea costante LME
ggplot() +
  # curva smussata dei residui OLS (empirica)
  geom_smooth(data = df_ols, aes(x = year, y = resid2),
              method = "loess", se = FALSE,
              linewidth = 1.3, color = "#1f77b4") +
  # linea costante della varianza fittata LME
  geom_hline(yintercept = var_fun,
             linetype = "dashed", color = "#d62728", linewidth = 1.2) +
  labs(x = "Anno", y = "Varianza / residui al quadrato") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

#linea orizzontale come ci si aspettava, ma non va bene --> RANDOM SLOPES NECESSARIA

#se ne può discutere, incremeneto minimo se si guarda al grafico con scala più grande (forse non ha s)