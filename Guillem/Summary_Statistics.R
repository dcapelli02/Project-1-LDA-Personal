
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


#Endpoints
model_endpoint <- lm(yini ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                    data = summary_stats)

summary(model_endpoint)


#Endpoints (Depending also on the bprs baseline) Linear Regression including bprs as a predictor
model_ancova <- lm(yini ~ y0 + trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                   data = summary_stats)

summary(model_ancova)


#Increment
model_increment <- lm(increment ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                      data = summary_stats)

summary(model_increment)



AIC(model_auc, model_endpoint, model_increment, model_ancova)
BIC(model_auc, model_endpoint, model_increment, model_ancova)


plot(summary_stats$AUC, fitted(model_auc), 
     main = "Predicted vs Observed AUC", 
     xlab = "Observed", 
     ylab = "Predicted")
abline(0, 1, col = "red")
