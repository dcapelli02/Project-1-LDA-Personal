library(ggplot2)
library(gridExtra)

# Define the outcome variable
outcome <- "bprs"

# Compute summary statistics: AUC, endpoint, increment
summary_stats <- alz_long %>%
  group_by(patid) %>%
  summarise(
    trial = first(trial),
    sex = first(sex),
    age = first(age),
    edu = first(edu),
    bmi = first(bmi),
    inkomen = first(inkomen),
    job = first(job),
    adl = first(adl),
    wzc = first(wzc),
    abpet_base = first(abpet_base),
    taupet_base = first(taupet_base),
    cdrsb_base = first(cdrsb_base),
    
    # Baseline Value
    y0 = first(.data[[outcome]]),
    
    # Last Available Measurement
    yini = last(na.omit(.data[[outcome]])),
    
    # Last Time Point Observed
    tmax = last(na.omit(year)),
    
    # Increment
    increment = yini - y0,
    
    # Area Under the Curve (AUC) using trapezoidal rule
    AUC = sum(diff(year) * (head(.data[[outcome]], -1) + tail(.data[[outcome]], -1)) / 2, na.rm = TRUE)
  )

summary(summary_stats)


#Build a linear model to see how baseline factors explain each person's overall summary statistics:

#AUC: Area Under the Curve
model_auc <- lm(AUC ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                data = summary_stats)

summary(model_auc)

cat("The most significant ones are: sex, edu, inkomen, job, adl, wzc, abpet_base and taupet_base (***)")

#Endpoints
model_endpoint <- lm(yini ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                    data = summary_stats)

summary(model_endpoint)

cat("The most significant ones are: age, HigherEdu, inkomen, job, adl, wzc, abpet_base and cdrsb_base(***)")

#Endpoints (Depending also on the bprs baseline) Linear Regression including bprs as a predictor
model_ancova <- lm(yini ~ y0 + trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                   data = summary_stats)

summary(model_ancova)

cat("The most significant ones are: inkomen, job, adl, wzc, abpet_base and cdrsb_base(***) [y0 highly significant (**)]")

#Increment
model_increment <- lm(increment ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                      data = summary_stats)

summary(model_increment)

cat("The most significant ones are: inkomen, job, adl, wzc, abpet_base and cdrsb_base(***)")


AIC(model_auc, model_endpoint, model_increment, model_ancova)


# Create predicted vs observed plots for all models
p1 <- ggplot(summary_stats, aes(x = AUC, y = fitted(model_auc))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed: AUC",
       x = "Observed AUC", 
       y = "Predicted AUC") +
  theme_minimal()

p2 <- ggplot(summary_stats, aes(x = yini, y = fitted(model_endpoint))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed: Endpoint",
       x = "Observed Endpoint", 
       y = "Predicted Endpoint") +
  theme_minimal()

p3 <- ggplot(summary_stats, aes(x = yini, y = fitted(model_ancova))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed: Endpoint (ANCOVA)",
       x = "Observed Endpoint", 
       y = "Predicted Endpoint") +
  theme_minimal()

p4 <- ggplot(summary_stats, aes(x = increment, y = fitted(model_increment))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed: Increment",
       x = "Observed Increment", 
       y = "Predicted Increment") +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, ncol = 2)

# Residuals vs Fitted plots for all models
par(mfrow = c(2, 2))

plot(fitted(model_auc), residuals(model_auc),
     main = "Residuals vs Fitted: AUC",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

plot(fitted(model_endpoint), residuals(model_endpoint),
     main = "Residuals vs Fitted: Endpoint",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

plot(fitted(model_ancova), residuals(model_ancova),
     main = "Residuals vs Fitted: Endpoint (ANCOVA)",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

plot(fitted(model_increment), residuals(model_increment),
     main = "Residuals vs Fitted: Increment",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")




