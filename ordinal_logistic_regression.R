# ## Package load
# library(tidyverse)
# library(gtsummary)
# library(survival)
# library(survminer)
# library(gridExtra)
# library(ggplot2)
# library(lubridate)
# library(broom)
# library(purrr)
# library(car)
# library(pROC)
# library(Cairo)
# library(officer)
# library(rvg)
# library(cowplot)
# library(FSA)  # For Dunn test
# library(ggsignif)
# library(clipr)
# # library(cardx)

# # read data file
# chestCT <- read.csv("Chest_CT_0719.csv")

# # check variables starts with "Date" or ends with "date"
# date_columns <- grep("^Date|date$", names(chestCT))
# # as.Date
# if (length(date_columns) > 0) {
#   chestCT[, date_columns] <- lapply(chestCT[, date_columns], as.Date)
# } else {
#   print("No columns starting with 'Date' or ending with 'date' found.")
# }

# ordinal logistic regression ####
# LVI##
fig2 <- chestCT %>%
  filter(Diagnosis_check == 'ALS') %>%
  filter(Hospital_ID != 00000000) %>% # unknown onset
  filter(Time_TIV_death < 365*20) %>% # remove extreme long survival to satisfy the proportional hazard assumption of RMI
  filter(Final_Diagnosis_JSK != "Others") %>% # FAS. There is no evidence for the other segments
  filter(Chest_CT_indication_category != "Underlying_resp_ds") %>% # underlying respiratory disease
  mutate(FVC_filter = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), FVC, NA_real_)) %>%
  mutate(across(everything(), ~replace(., . == "", NA))) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>%
  mutate(Ht = ifelse(is.na(Ht), mean(Ht, na.rm = TRUE), Ht)) %>%
  mutate(CT_age = round((Date_CT - Date_birth)/365,1),
         onset_age = round((Date_onset - Date_birth)/365,1),
         onset_to_CT = round((Date_CT - Date_onset)/365*12,1),
         fu_duration = round((Date_censor_fu - First_visit_date)/365*12,1),
         diagnostic_delay = round((Date_dx_google - Date_onset)/365*12,1)) %>%
  mutate(Gastrostomy = if_else(!is.na(Date_Gastrostomy), 1, 0),
         NIV = if_else(!is.na(Date_NIV), 1, 0),
         TIV = if_else(!is.na(Date_TIV), 1, 0),
         Death = if_else(!is.na(Date_death), 1, 0)) %>%
  mutate(
    Cr_filter = if_else(between(as.numeric(Cr_date - Date_CT, units = "days"), -30, 30) & (is.na(Cr_hydration) | Cr_hydration == "" | Cr_hydration == " "), Cr, NA_real_),
    UA_filter = if_else(between(as.numeric(UA_date - Date_CT, units = "days"), -30, 30) & (is.na(UA_hydration) | UA_hydration == "" | UA_hydration == " "), UA, NA_real_),
    CK_filter = if_else(between(as.numeric(CK_date - Date_CT, units = "days"), -30, 30) & (is.na(CK_hydration) | CK_hydration == "" | CK_hydration == " "), CK, NA_real_),
    CO2_filter = if_else(between(as.numeric(CO2_date - Date_CT, units = "days"), -30, 30) & (is.na(ABGA_on_MV) | ABGA_on_MV == "" | ABGA_on_MV == " "), CO2, NA_real_),
    O2_filter = if_else(between(as.numeric(O2_date - Date_CT, units = "days"), -30, 30) & (is.na(ABGA_on_MV) | ABGA_on_MV == "" | ABGA_on_MV == " "), O2, NA_real_),
    Ca_filter = if_else(between(as.numeric(Ca_date - Date_CT, units = "days"), -30, 30), Ca, NA_real_),
    P_filter = if_else(between(as.numeric(P_date - Date_CT, units = "days"), -30, 30), P, NA_real_),
    Wt_filter = if_else(between(as.numeric(Wt_CT_date - Date_CT, units = "days"), -90, 90), Wt_CT, NA_real_),
    ALSFRS_filter = if_else(between(as.numeric(ALSFRS_closest_CT_date - Date_CT, units = "days"), -90, 90), ALSFRS_closest_CT, NA_real_),
    FVC_filter = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), FVC, NA_real_),
    FVC_measured = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), FVC_measured, NA_real_),
    FVC_supine = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), FVC_supine, NA_real_),
    MIP_filter = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), MIP, NA_real_),
    MEP_filter = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), MEP, NA_real_),
    SNIP_filter = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), SNIP, NA_real_),
    PEF_filter = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), PEF, NA_real_)
  ) %>%
  mutate(Bulbar_filter = if_else(between(as.numeric(ALSFRS_closest_CT_date - Date_CT, units = "days"), -30, 30), Bulbar, NA_real_),
         Fine_filter = if_else(between(as.numeric(ALSFRS_closest_CT_date - Date_CT, units = "days"), -30, 30), Fine, NA_real_),
         Gross_filter = if_else(between(as.numeric(ALSFRS_closest_CT_date - Date_CT, units = "days"), -30, 30), Gross, NA_real_),
         Resp_filter = if_else(between(as.numeric(ALSFRS_closest_CT_date - Date_CT, units = "days"), -30, 30), Resp, NA_real_)) %>%
  mutate(Ca_P = Ca_filter/P_filter) %>%
  mutate(BMI_filter = round(Wt_filter/(Ht/100)^2,1),
         BMI_pre = round(Wt_pre/(Ht/100)^2,1)) %>%
  mutate(Wt_change = Wt_pre - Wt_filter,
         BMI_change = BMI_pre - BMI_filter) %>%
  mutate(Onset_region = case_when(
    Onset_region %in% c("C", "LS") ~ "Limb",
    TRUE ~ Onset_region)) %>%
  mutate(Final_Diagnosis_JSK = case_when(
    Final_Diagnosis_JSK %in% c("ALS", "fALS") ~ "ALS",
    TRUE ~ Final_Diagnosis_JSK)) %>%
  mutate(Time = as.numeric(ALSFRS_closest_CT_date - Date_onset)/365*12) %>%
  mutate(ALSFRS_rate = if_else(between(as.numeric(ALSFRS_closest_CT_date - Date_CT, units = "days"), -90, 90), ALSFRS_closest_CT, NA_real_)) %>%
  mutate(Progression_rate = round((48-ALSFRS_rate)/Time, 2)) %>%
  mutate(across(where(~ inherits(., "difftime")), as.numeric)) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) 

fig2$Stage_final <- factor(fig2$Stage_final, ordered = TRUE)

library(MASS)
# LVI##
model <- polr(Stage_final ~ LVI + Sex + CT_age, data = fig2, Hess = TRUE)

summary(model)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 

model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

# RMI##
model <- polr(Stage_final ~ SMI + Sex + CT_age, data = fig2, Hess = TRUE)

summary(model)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 

model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)


# FVC##
model <- polr(Stage_final ~ FVC_filter + Sex + CT_age, data = fig2, Hess = TRUE)

summary(model)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 

model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)


# Muscle_HU##
model <- polr(Stage_final ~ Muscle_HU + Sex + CT_age, data = fig2, Hess = TRUE)

summary(model)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 

model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)


# SFI##
model <- polr(Stage_final ~ SFI + Sex + CT_age, data = fig2, Hess = TRUE)

summary(model)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 

model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)


# SF_HU##
model <- polr(Stage_final ~ SF_HU + Sex + CT_age, data = fig2, Hess = TRUE)

summary(model)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 

model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)




# Univariable check ####
model <- polr(Stage_final ~ LVI, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- polr(Stage_final ~ SMI, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- polr(Stage_final ~ FVC_filter, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- polr(Stage_final ~ Muscle_HU, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- polr(Stage_final ~ SFI, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- polr(Stage_final ~ SF_HU, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)



model <- polr(Stage_final ~ CT_age, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)


model <- polr(Stage_final ~ Sex, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- polr(Stage_final ~ Onset_region, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- polr(Stage_final ~ ALSFRS_filter, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- polr(Stage_final ~ Progression_rate, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- polr(Stage_final ~ BMI_filter, data = fig2, Hess = TRUE)
coefs <- coef(summary(model))
p_values <- pnorm(abs(coefs[, "t value"]), lower.tail = FALSE) * 2
p_values_rounded <- round(p_values, 4)
coefs <- cbind(coefs, "p value" = p_values_rounded)
print(coefs) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)


