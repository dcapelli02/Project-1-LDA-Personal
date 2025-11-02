install.packages("nlme")
library("nlme")

alz_long$year_fac <- as.factor(alz_long$year)

mult_model_1 <- gls(
  bprs ~ 
    sex + trial + job + age  + inkomen + adl + wzc + cdrsb_base + ab_base + tau_base +(sex+job+wzc+ cdrsb_base+ tau_base)*year ,
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
    sex + trial + job + age  + inkomen + adl + wzc + cdrsb_base + ab_base + tau_base +(sex+job+wzc+ cdrsb_base+ tau_base)*year ,
  correlation = corSymm(form = ~ time_idx | sample),  # type=TOEPH
  weights = varIdent(form = ~ 1 | year),  # eterogeneità delle varianze
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_1, mult_model_2)