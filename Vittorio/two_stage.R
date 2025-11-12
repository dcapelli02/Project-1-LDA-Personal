library(metafor)

res.list <- lmList(bprs ~ year | sample, data = alz_long)
b <- lapply(res.list, coef)
V <- lapply(res.list, vcov)

estm <- rep(c("intercept","slope"), length(b))
subj <- rep(names(b), each=2)


b <- unlist(b)
V <- bldiag(V)

# Prendi solo la prima osservazione di ogni soggetto
covariate_data <- alz_long[!duplicated(alz_long$sample), 
                           c("sample", "trial", "sex", "age", "edu", 
                             "bmi", "inkomen", "job", "adl", "wzc", 
                             "ab_base", "tau_base", "cdrsb_base")]
# Allinea all'ordine dei soggetti di lmList
subj_names_ordered <- names(res.list)
covariate_data_aligned <- covariate_data[match(subj_names_ordered, covariate_data$sample), ]

# Crea vettori ripetuti per ogni stima (intercept/slope)
trial_rep   <- rep(covariate_data_aligned$trial, each = 2)
sex_rep     <- rep(covariate_data_aligned$sex, each = 2)
age_rep     <- rep(covariate_data_aligned$age, each = 2)
edu_rep     <- rep(covariate_data_aligned$edu, each = 2)
bmi_rep     <- rep(covariate_data_aligned$bmi, each = 2)
inkomen_rep <- rep(covariate_data_aligned$inkomen, each = 2)
job_rep     <- rep(covariate_data_aligned$job, each = 2)
adl_rep     <- rep(covariate_data_aligned$adl, each = 2)
wzc_rep     <- rep(covariate_data_aligned$wzc, each = 2)
ab_base_rep   <- rep(covariate_data_aligned$ab_base, each = 2)
tau_base_rep  <- rep(covariate_data_aligned$tau_base, each = 2)
cdrsb_base_rep   <- rep(covariate_data_aligned$cdrsb_base,  each = 2)
length(coef(res.list[[1]]))

subj <- rep(subj_names_ordered, each = 2)

length(b)           # 2 * n
length(trial_rep)   # deve essere uguale a length(b)
length(subj) 

k <- length(b)
V_identity <- diag(1, k)

res2_full <- rma.mv(b ~ estm + 
                      estm:trial_rep + estm:sex_rep + estm:age_rep + estm:edu_rep + 
                      estm:bmi_rep + estm:inkomen_rep + estm:job_rep + estm:adl_rep + 
                      estm:wzc_rep + estm:ab_base_rep + estm:tau_base_rep + estm:cdrsb_base_rep - 1, 
                    random = ~ estm | subj, 
                    V = V_identity,
                    struct = "UN")
res2_full
