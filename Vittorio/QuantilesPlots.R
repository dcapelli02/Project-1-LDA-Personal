#Variabile discretized accoridng to quantiles, maybe better

#1
alz_long<- alz_long %>%
  mutate(age_disc = ntile(age, 4))

ggplot(alz_long, aes(x = year, y = bprs, group = age_disc, color = age_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age_disc without variance")


#2
alz_long<- alz_long %>%
  mutate(bmi_disc = ntile(bmi, 4))

ggplot(alz_long, aes(x = year, y = bprs, group = bmi_disc, color = bmi_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt bmi_disc without variance")


#3
alz_long<- alz_long %>%
  mutate(inkomen_disc = ntile(inkomen, 4))

ggplot(alz_long, aes(x = year, y = bprs, group = inkomen_disc, color = inkomen_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt inkomen_disc without variance")


#4
alz_long<- alz_long %>%
  mutate(adl_disc = ntile(adl, 4))
ggplot(alz_long, aes(x = year, y = bprs, group = adl_disc, color = adl_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt adl_disc without variance")

#5
alz_long<- alz_long %>%
  mutate(cdrsb_disc_base = ntile(cdrsb_base, 4))
ggplot(alz_long, aes(x = year, y = bprs, group = cdrsb_disc_base, color = cdrsb_disc_base)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt cdrsb_disc_base without variance")

#6
alz_long<- alz_long %>%
  mutate(ab_disc_base = ntile(ab_base, 4))
ggplot(alz_long, aes(x = year, y = bprs, group = ab_disc_base, color = ab_disc_base)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt ab_disc_base without variance")

#7
alz_long<- alz_long %>%
  mutate(tau_disc_base = ntile(tau_base, 4))
ggplot(alz_long, aes(x = year, y = bprs, group = tau_disc_base, color = tau_disc_base)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt tau_disc_base without variance")
