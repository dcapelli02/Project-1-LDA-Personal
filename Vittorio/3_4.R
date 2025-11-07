install.packages("nlme")
library("nlme")

#multivariate model
alz_long$year_fac <- as.factor(alz_long$year)

mult_model_1 <- gls(
  bprs ~ 
    year + sex + trial + job + age  + inkomen + adl + wzc + cdrsb_base + ab_base + tau_base +(sex+job+wzc)*year ,
  correlation = corAR1(form = ~ year | sample),  # type=ARH(1)
  weights = varIdent(form = ~ 1 | year),  # eterogeneità delle varianze
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)


alz_long <- alz_long %>%
  group_by(sample) %>%
  arrange(year) %>%
  mutate(time_idx = row_number()) %>%
  ungroup()

mult_model_2 <- gls(
  bprs ~ 
    sex +trial+ age + adl + wzc + cdrsb_base + year+(sex+wzc+age)*year ,
  correlation = corSymm(form = ~ time_idx | sample),  # type=TOEPH
  weights = varIdent(form = ~ 1 | year),  # eterogeneità delle varianze
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

prova <- gls(
  bprs ~ 
    sex +trial+ age + adl + wzc + cdrsb_base +(sex+wzc+age)*year ,
  correlation = corSymm(form = ~ time_idx | sample),  # type=TOEPH
  weights = varIdent(form = ~ 1 | year),  # eterogeneità delle varianze
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

mult_model_3 <- gls(
  bprs ~ 
    sex +trial+ age + adl + wzc + cdrsb_base +(wzc+age)*year ,
  correlation = corSymm(form = ~ time_idx | sample),  # type=TOEPH
  weights = varIdent(form = ~ 1 | year),  # eterogeneità delle varianze
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

summary(mult_model_2)

anova( mult_model_2, prova)

#STEPAIC PROCEDURE IS BETTER, DANIELE'S CODE



#cov and cor analysis
var_by_year <- alz_long %>%
  group_by(year) %>%
  summarise(var_bprs = var(bprs, na.rm = TRUE))

ggplot(var_by_year, aes(x = year, y = var_bprs)) +
  geom_line() + geom_point(size=2) +
  theme_minimal() +
  labs(title = "Variance of BPRS over time")


cov_matrix_bprs <- cov(alz[, c(18:24)], use = "pairwise.complete.obs")
round(cov_matrix_bprs, 2)

heatmap(cov_matrix_bprs, main = "Covariance matrix of BPRS")


cor_matrix_bprs <- cor(alz[, c(18:24)], use = "pairwise.complete.obs")
round(cor_matrix_bprs, 2)

heatmap(cor_matrix_bprs, main = "Correlation matrix of BPRS")



#INFORMATIVE DROPOUT

#I find the latest reading for each patient
ultima_rilevazione <- alz_long %>%
  group_by(patid) %>%
  summarise(last_time = max(year))

#  indidcator 1 = there is untile the last one observation
dati <- alz_long %>%
  left_join(ultima_rilevazione, by = "patid") %>%
  mutate(last_visit = ifelse(year == last_time, 1, 0))

ultimo_anno <- max(alz_long$year, na.rm = TRUE)

bprs_baseline <- alz_long %>%
  group_by(patid) %>%
  summarise(
    # first value of BPRS
    BPRS_start = bprs[!is.na(bprs) & year == min(year[!is.na(bprs)])][1],
    # latest reading
    last_year = max(year[!is.na(bprs)], na.rm = TRUE)
  ) %>%
  mutate(
    completed = ifelse(last_year == ultimo_anno, 1, 0)
  )

bprs_baseline %>%
  group_by(completed) %>%
  summarise(
    mean_BPRS = mean(BPRS_start, na.rm = TRUE),
    sd_BPRS = sd(BPRS_start, na.rm = TRUE),
    n = n()
  )

t_test <- t.test(BPRS_start ~ completed, data = bprs_baseline)

#conclusion:
#there is an informative drop out, who start the study with high levels of bprs drop out earlier than those who start with low levels of bprs


#two stage model
#Transformation of the time scale to linearize the profiles?
ggplot(alz_rist_long, aes(x = year, y = bprs, group = patid)) +
  geom_line(alpha = 0.3) +
  stat_smooth(aes(group = 1), method = "lm", se = FALSE) +
  labs(title = "Profili Individuali e Trend Medio")
#we don't need it

#outcome scale trasformation?
ggplot(alz_long, aes(x = bprs)) +
  geom_histogram(bins = 30) +
  labs(title = "Distribution of Outcome Variable")
  
#no, its normal

  library(tidyverse)
# stage 1
# Per ogni paziente, fitto un modello separato
stage1_models <- alz_long %>%
  group_by(patid) %>%
  do(model = lm(bprs ~ year, data = .))

## Calcola prima il numero di osservazioni per soggetto
n_obs_data <- alz_long %>%
  group_by(patid) %>%
  summarise(n_obs = sum(!is.na(bprs)), .groups = "drop")


# Estrai intercept, slope e R² dai modelli individuali
individual_fits <- stage1_models %>%
  rowwise() %>%  # serve perché ogni riga ha un modello diverso
  mutate(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    r2_linear = summary(model)$r.squared
  ) %>%
  ungroup() %>%
  select(patid, intercept, slope, r2_linear) %>%
  left_join(n_obs_data, by = "patid") %>%  # Ora la join funziona
  filter(!is.na(r2_linear))

r2_overall <- mean(individual_fits$r2_linear, na.rm = TRUE)

ggplot(individual_fits, aes(x = n_obs, y = r2_linear)) +
  geom_point(alpha = 0.6, color = "blue", size = 2) +
  geom_hline(yintercept = r2_overall, linetype = "dashed", color = "red", size = 1) +
  labs(
    title = "Scatterplot of R² under Linear Model",
    x = "Number n_i of measurements",
    y = "Coefficient Ri²"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

#linear model not bad, all the subject have a good fit 

# stage 2, SEE DANIELE'S CODE, HE DID THE STEPAIC!!!!
a=left_join(alz_long, individual_fits, by = "patid")
#Modello per slope
model_slope <- lm(slope ~ sex+wzc+age,
                  data = a)

summary(model_slope)


# Modello per l'intercetta
model_intercept <- lm(intercept ~ sex +trial+ age + adl + wzc + cdrsb_base,
                      data = a)

summary(model_intercept)


