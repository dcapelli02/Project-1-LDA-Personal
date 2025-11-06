# Carregar paquets
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggcorrplot)
library(nlme)

# Llegir dades
alz <- read_sas("C:/Users/win11/Documents/Guillem/Erasmus/Assignatures/Longitudinal Data Analysis/Project-1-LDA-Personal/alzheimer25.sas7bdat")

# Convertir variables qualitatives a factors
alz$trial <- factor(alz$trial)
alz$sex <- factor(alz$sex, levels = c(0, 1), labels = c("male", "female"))
alz$edu <- factor(alz$edu, levels = c(1,2,3,4),
                  labels = c("Primary Ed.", "Lower Secondary Ed.", "Upper Secondary Ed.", "Higher Ed."))
alz$job <- factor(alz$job, levels = c(0, 1), labels = c("No job", "Job"))
alz$wzc <- factor(alz$wzc, levels = c(0, 1), labels = c("Home", "Residence"))

alz$abpet_base <- alz$abpet0
alz$taupet_base <- alz$taupet0
alz$cdrsb_base <- alz$cdrsb0

# Transformar a long format
alz_long <- pivot_longer(
  alz,
  cols = matches("^(bprs|cdrsb|abpet|taupet)\\d+$"),
  names_to = c(".value", "year"),
  names_pattern = "(bprs|cdrsb|abpet|taupet)(\\d+)"
)

# Convertir year a numèric
alz_long$year <- as.numeric(as.character(alz_long$year))

vars_long <- c("cdrsb", "bprs", "abpet", "taupet")

for (var in vars_long) {
  
  cat("\n==============================\n")
  cat("Variable:", var, "\n")
  
  # Mean & SD
  summary_by_year <- alz_long %>%
    group_by(year) %>%
    summarise(
      mean_val = mean(.data[[var]], na.rm = TRUE),
      sd_val   = sd(.data[[var]], na.rm = TRUE)
    )
  
  # Gràfic de mitjana amb ribbon
  print(
    ggplot(summary_by_year, aes(x = year, y = mean_val)) +
      geom_line(color = "red", size = 1.2) +
      geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                  alpha = 0.2, fill = "blue") +
      labs(title = paste(var, "evolution"), x = "Year", y = var) +
      theme_minimal()
  )
  
  # Gràfic de SD
  print(
    ggplot(summary_by_year, aes(x = year, y = sd_val)) +
      geom_line(color = "darkgreen", size = 1.2) +
      geom_point(color = "darkgreen", size = 2) +
      labs(title = paste(var, "SD evolution (variance structure)"), 
           x = "Year", y = paste("SD of", var)) +
      theme_minimal()
  )
  
  # Correlation
  wide_df <- alz_long %>%
    select(patid, year, .data[[var]]) %>%
    pivot_wider(names_from = year, values_from = var)
  
  cor_mat <- cor(wide_df[,-1], use = "pairwise.complete.obs")
  print(cor_mat)
  
  print(
    ggcorrplot(cor_mat, 
               method = "circle", 
               lab = TRUE, 
               title = paste("Correlation structure of", var, "over years"))
  )
}


#Multivariate Model


model_gls_1 <- gls(
  bprs ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + cdrsb_base + 
    abpet_base + taupet_base + (trial + sex + age + edu + bmi + inkomen + 
        job + adl + wzc + cdrsb_base + abpet_base + taupet_base) * year,
  
  correlation = corAR1(form = ~ year | patid),
  weights = varIdent(form = ~ 1 | year),
  data = alz_long,
  method = "ML",
  na.action = na.exclude
)


model_gls_2 <- gls(
  bprs ~ trial + sex + age + edu + job + adl + wzc + cdrsb_base + 
    (sex + age + job + adl + wzc + cdrsb_base) * year,
  
  correlation = corAR1(form = ~ year | patid),
  weights = varIdent(form = ~ 1 | year),
  data = alz_long,
  method = "ML",
  na.action = na.exclude
)

model_gls_3 <- gls(
  bprs ~ trial + sex + age + adl + job + (wzc + cdrsb_base) * year,
  
  correlation = corAR1(form = ~ year | patid),
  weights = varIdent(form = ~ 1 | year),
  data = alz_long,
  method = "ML",
  na.action = na.exclude
)

#Without Trial - It has larger BIC 
model_gls_4 <- gls(
  bprs ~ sex + age + adl + job + (wzc + cdrsb_base) * year,
  
  correlation = corAR1(form = ~ year | patid),
  weights = varIdent(form = ~ 1 | year),
  data = alz_long,
  method = "ML",
  na.action = na.exclude
)

anova(model_gls_1, model_gls_2, model_gls_3, model_gls_4)


#alz_long <- alz_long %>%
#  arrange(patid, year) %>%
#  group_by(patid) %>%
#  mutate(year_seq = row_number()) %>%
#  ungroup()



