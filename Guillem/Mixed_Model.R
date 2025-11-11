# MIXED MODEL ANALYSIS

library(nlme)
library(dplyr)
library(ggplot2)
library(tibble)


# Fit the mixed_effects model
# Random intercept and slope per patient
lme_model <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job + adl + wzc + cdrsb_base + abpet_base + taupet_base) * year,
  random = ~ year | patid,                       # random intercept + slope per patient
  correlation = corAR1(form = ~ year | patid),   
  data = alz_long,
  method = "REML",
  na.action = na.exclude
)

summary(lme_model)

#Extract the random effects from each patient
ranef_df <- ranef(lme_model) %>% 
  rownames_to_column("patid") %>%
  rename(random_intercept = `(Intercept)`, random_slope = year) %>%
  mutate(patid = as.numeric(patid))

head(ranef_df)


# Add predicted values from the mixed-effects model to the dataset
alz_pred <- alz_long %>%
  mutate(
    # Level = 1 includes random effects. This means that the predictions are subject-specific
    #Each patient gets their own intercept and slope values.
    bprs_pred = predict(lme_model, level = 1),   

    # Level = 0 includes only fixed effects. This means that the predictions are population-level, ignoring indivudual differences
    bprs_pred_fixed = predict(lme_model, level = 0)
  )


# Plot the predictions
ggplot(alz_pred, aes(x = year, y = bprs, group = patid)) +
  geom_point(alpha = 0.5) +
  
  # Line for level = 1, includes random effects
  geom_line(aes(y = bprs_pred), color = "red", alpha = 0.7) +  
  
  # Line for level = 0, includes only fixed effects
  geom_line(aes(y = bprs_pred_fixed), color = "blue", linetype = "dashed", alpha = 0.7) +
  
  labs(
    title = "Mixed-effects model: observed vs predicted",
    subtitle = "Red = subject-specific (with random effects), Blue = population-level (fixed effects)",
    x = "Year",
    y = "BPRS"
  ) +
  theme_minimal(base_size = 13)

# Sample 10 so the graphic can be better understood
set.seed(123)
sample_patients <- sample(unique(alz_pred$patid), 10)

ggplot(alz_pred %>% filter(patid %in% sample_patients),
       aes(x = year, y = bprs, group = patid)) +
  
  # Plots the observed BPRS points
  geom_point(aes(color = as.factor(patid)), alpha = 0.7, size = 2) +   
  
  # Line for level = 1, includes random effects
  geom_line(aes(y = bprs_pred), color = "red", size = 1) +             
  
  # Line for level = 0, includes only fixed effects
  geom_line(aes(y = bprs_pred_fixed), color = "blue", linetype = "dashed", size = 1) +  # population-level
  
  labs(
    title = "Mixed-effects model: observed vs predicted",
    subtitle = "Red = subject-specific (random effects), Blue dashed = population-level (fixed effects)",
    x = "Year",
    y = "BPRS",
    color = "Patient ID"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")



# =========================
# 5. Inspect random effects distributions
# =========================

# Random intercepts
# Difference between a patient’s baseline value and the overall population baseline
ggplot(ranef_df, aes(x = random_intercept)) +
  geom_histogram(bins = 20, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Random Intercepts Distribution", x = "Intercept", y = "Count") +
  theme_minimal()

# How much a patient’s trajectory over time deviates from the average population trend
ggplot(ranef_df, aes(x = random_slope)) +
  geom_histogram(bins = 20, fill = "red", color = "black", alpha = 0.7) +
  labs(title = "Random Slopes Distribution", x = "Slope (Year)", y = "Count") +
  theme_minimal()

# Intercept vs Slope 
# How patient starting points relate to their change over time
ggplot(ranef_df, aes(x = random_intercept, y = random_slope)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  labs(title = "Random Intercept vs Slope", x = "Intercept", y = "Slope") +
  theme_minimal()

