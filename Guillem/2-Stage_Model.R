# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(nlme)  

# Stage 1: Linear regressions for every subject
stage1_models <- alz_long %>%
  group_by(patid) %>%
  group_map(~ lm(bprs ~ year, data = .x), .keep = TRUE)

#bprs_it = b_0i + b_1i * year_it + e_it

# Extract coefficients for each subject
beta_coefs <- map_dfr(stage1_models, coef) %>%
  mutate(patid = alz_long %>% group_by(patid) %>% group_keys() %>% pull(patid))

head(beta_coefs)

# Stage 2: Explain variability in intercepts and slopes using baseline covariates

# Merge baseline covariates
subject_data <- alz_long %>%
  select(patid, trial, sex, age, edu, bmi, inkomen, job, adl, wzc, abpet_base, taupet_base, cdrsb_base) %>%
  distinct() %>%
  left_join(beta_coefs, by = "patid")

# subject_data has all the covariates + intercept + slope for each pacient

# Reshape data for GLS
stage2_data <- subject_data %>%
  # `(Intercept)` -> Stage 1 intercept; year -> Stage 1 slope
  select(patid, `(Intercept)`, year, trial, sex, age, edu, bmi, inkomen, job, adl, wzc, 
         abpet_base, taupet_base, cdrsb_base) %>%
  pivot_longer(cols = c(`(Intercept)`, year), 
               names_to = "param",     #new colum with intercept and year values
               values_to = "beta")     #contains the value of the coefficient

#Convertteix param en un factor, així el model sap que és una varaible categòrica
stage2_data$param <- factor(stage2_data$param) 


# Stage 2 GLS model with correlation between intercept and slope
stage2_gls_model <- gls(
  beta ~ param + trial + sex + age + edu + bmi + inkomen + job + adl + wzc + 
    abpet_base + taupet_base + cdrsb_base,
  data = stage2_data,
  correlation = corSymm(form = ~ 1 | patid),   # correlation structure per patient
  weights = varIdent(form = ~ 1 | param),     # allow different variances for intercept/slope
  method = "ML",
  na.action = na.exclude
)

summary(stage2_gls_model)


cat("paramyear = -68.3 -> slope much smaller than the intercept")
cat("Correlation intercept & slope = -0.977 -> strong negative correlation: higher intercept -> lower slope")
cat("Trial effects -> mainly affect intercept")
cat("Age = 0.786 -> older patients have higher beta")
cat("JobJob = -2.44 -> employment reduces beta")
cat("WZCResidence = 0.928 -> living in a care center increases beta")
cat("cdrsb_base = -0.0199 -> higher baseline cdrsb -> faster bprs decline")





# Plot Stage the Predictions!
ggplot(alz_long, aes(x = year, y = bprs, group = patid)) +
  geom_point(alpha = 0.5) +
  geom_line(data = subject_data %>%
              rowwise() %>%
              mutate(years = list(0:6)) %>%
              unnest(cols = c(years)) %>%
              mutate(bprs_pred = `(Intercept)` + year * years), # Calculates the Prediction!!!
            aes(x = years, y = bprs_pred, group = patid),
            color = "red", alpha = 0.5) +
  labs(title = "Stage 1: Subject-specific linear regressions for BPRS") +
  theme_minimal()




#LINEAR MODEL (Worse because it doesn't take into account Sigma and D)


# Stage 1

# Fas un model lineal de bprs que depengui de l'any.
stage1_models <- alz_long %>%
  group_by(patid) %>%
  group_map(~ lm(bprs ~ year, data = .x), .keep = TRUE)

# Extreus els coeficients per a cada pacient
beta_coefs <- map_dfr(stage1_models, coef) %>%
  mutate(patid = alz_long %>% group_by(patid) %>% group_keys() %>% pull(patid))
  

head(beta_coefs)

# Stage 2: Explain variability in intercepts and slopes using baseline covariates

# Merge baseline covariates
subject_data <- alz_long %>%
  select(patid, trial, sex, age, edu, bmi, inkomen, job, adl, wzc, abpet_base, taupet_base, cdrsb_base) %>%
  distinct() %>%
  left_join(beta_coefs, by = "patid")

# Intercept model
#Fas un model lineal per explicar l'intercept dels models lineals dels coeficients anteriors
stage2_intercept_model <- lm(`(Intercept)` ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + 
                        abpet_base + taupet_base + cdrsb_base,
                      data = subject_data)

summary(stage2_intercept_model)

cat("Most significant variables for the intercept: intercept, trial(most of them), age, bmi, job, wzc and cdrsb_base(**).")

# Slope model
#Fas un model lineal per explicar la slope dels models lineals dels coeficients anteriors
stage2_slope_model <- lm(year ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + 
                    abpet_base + taupet_base + cdrsb_base,
                  data = subject_data %>%
                    rename(year = `year`))  

summary(stage2_slope_model)

cat("Most significant variables for the slope: cdrsb_base (and  with (*): age, job, adl and wzc")


