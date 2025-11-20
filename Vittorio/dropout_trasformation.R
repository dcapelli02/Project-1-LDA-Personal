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
t_test
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
