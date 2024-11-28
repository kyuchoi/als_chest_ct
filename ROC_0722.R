library(tidyverse)
library(gtsummary)
library(survival)
library(survminer)
library(gridExtra)
library(ggplot2)
library(lubridate)
library(broom)
library(purrr)
library(car)
library(pROC)
library(Cairo)
library(officer)
library(rvg)
library(cowplot)

library(SurvMetrics)



# read data file
chestCT <- read.csv("Chest_CT_0719.csv")

# check variables starts with "Date" or ends with "date"
date_columns <- grep("^Date|date$", names(chestCT))
# as.Date
if (length(date_columns) > 0) {
  chestCT[, date_columns] <- lapply(chestCT[, date_columns], as.Date)
} else {
  print("No columns starting with 'Date' or ending with 'date' found.")
}

## inclusion criteria and refine variable names
included <- chestCT %>%
  filter(Stage_final != 0) %>%
  filter(Diagnosis_check == 'ALS') %>%
  filter(Hospital_ID != 00000000) %>% # unknown onset
  filter(Time_TIV_death < 365*20) %>% # remove extreme long survival to satisfy the proportional hazard assumption of RMI
  filter(Final_Diagnosis_JSK != "Others") %>% # FAS. There is no evidence for the other segments
  filter(Chest_CT_indication_category != "Underlying_resp_ds") %>% # underlying respiratory disease
  mutate(across(everything(), ~replace(., . == "", NA)))


# redefine variables
included <- included %>%
  mutate(Ht = ifelse(is.na(Ht), mean(Ht, na.rm = TRUE), Ht)) %>%
  mutate(CT_age = round((Date_CT - Date_birth)/365,1),
         onset_age = round((Date_onset - Date_birth)/365,1),
         onset_to_CT = round((Date_CT - Date_onset)/365*12,1),
         fu_duration = round((Date_censor_fu - First_visit_date)/365*12,1)) %>%
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



## ROC ####
tdROC <- included %>% 
  filter(Stage_final %in% c('1','2','3')) %>%
  filter(onset_to_CT <= 12) %>%
  rename(RMI = SMI)%>%
  mutate(onset_CT_day = Date_CT - Date_onset,
         CT_TIV_death = Time_TIV_death - onset_CT_day)

tdROC <- tdROC[, c("CT_TIV_death", "Event_TIV_death", "LVI", "RMI", "FVC")]

tdROC_LVI <- tdROC[!is.na(tdROC$LVI),]
tdROC_SMI <- tdROC[!is.na(tdROC$RMI),]
tdROC_FVC <- tdROC[!is.na(tdROC$FVC),]

# variable lists
variables <- c("LVI", "RMI", "FVC")
data_list <- list(tdROC_LVI, tdROC_SMI, tdROC_FVC)
colors <- c("#D54C4C", "#228B22", "#4169E1")


#Time dependent ROC: ############################################################################################################
roc_results <- list()
predict_times <- seq(365, 365*4, by = 365) # (1-year, 2-year, 3-year, 4-year)
# calculate ROC for each variable
for (i in 1:length(variables)) {
  var <- variables[i]
  data <- data_list[[i]]
  cox.lp <- predict(coxph(Surv(CT_TIV_death, Event_TIV_death) ~ get(var), data = data), type = "lp")
  roc_results[[var]] <- list()
  # calculate ROC for each time point
  for (time in predict_times) {
    response <- ifelse(data$Event_TIV_death == 1 & data$CT_TIV_death <= time, 1, 0)
    roc_results[[var]][[paste0("ROC_", time)]] <- roc(response, cox.lp, directions="<")
  }
}

par(mfrow = c(2, 2)) 
for (time in predict_times) {
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "1-Specificity", ylab = "Sensitivity", 
       main = paste(time / 365, "years"))
  abline(0, 1, lty = 3)  
  auc_values <- c() 
  roc_curves <- list() 
  for (i in 1:length(variables)) {
    var <- variables[i]
    roc_curve <- roc_results[[var]][[paste0("ROC_", time)]]
    lines(1-roc_curve$specificities, roc_curve$sensitivities, type = "l", col = colors[i], lwd = 2)
    auc_values <- c(auc_values, round(roc_curve$auc, 3))
  }
  legend_text <- paste0(variables, " (AUC=", auc_values, ")")
  # legend(x = 0.4, y = 0.1, , legend = legend_text, col = colors, lty = 1, lwd = 2, bty = "o")
  legend(x = 0.5, y = 0.3, legend = legend_text, col = colors, lty = 1, lwd = 2, bty = "o")

  # p-value 
  p_values <- c()
  roc_test <- roc.test(roc_results$LVI[[paste0("ROC_", time)]], roc_results$FVC[[paste0("ROC_", time)]], method = "delong", va1=1, var2=1)
  p_values <- c(p_values, paste("LVI vs. FVC: ", round(roc_test$p.value, 3)))
  roc_test <- roc.test(roc_results$RMI[[paste0("ROC_", time)]], roc_results$FVC[[paste0("ROC_", time)]], method = "delong", va1=1, var2=1)
  p_values <- c(p_values, paste("RMI vs. FVC: ", round(roc_test$p.value, 3)))
  roc_test <- roc.test(roc_results$LVI[[paste0("ROC_", time)]], roc_results$RMI[[paste0("ROC_", time)]], method = "delong", va1=1, var2=1)
  p_values <- c(p_values, paste("LVI vs. RMI: ", round(roc_test$p.value, 3)))
  legend("right", legend = p_values, bty = "n", title="DeLong's test (p-value)")
}

#Time dependent AUC, iAUC: ########################################################################################################################

par(mfrow = c(1, 1)) 
tmax <- 365 * 4 
plot(NULL, xlim = c(0.35, tmax / 365), ylim = c(0.4, 1), xlab = "Follow-up Time (years)", ylab = "AUC", main = "Time dependent AUC graph")
abline(0.5, 0, lty = 3)
# calculate ROC for each variable
iauc_values <- c()
for (i in 1:length(variables)) {
  var <- variables[i]
  data <- data_list[[i]]
  cox.lp <- predict(coxph(Surv(CT_TIV_death, Event_TIV_death) ~ get(var), data = data), type = "lp")
  roc_results[[var]] <- list(AUC = numeric(), S = list())
  # calculate ROC for each time point
  utimes <- unique(data$CT_TIV_death[data$Event_TIV_death == 1 & data$CT_TIV_death <= tmax]) # predict when the event occurs
  predict_times <- sort(utimes)
  predict_times <- predict_times[predict_times >= 200 & predict_times <= 365*4] # from 1-year to 4-year
  for (time in predict_times) {
    response <- ifelse(data$Event_TIV_death == 1 & data$CT_TIV_death <= time, 1, 0)
    roc_curve <- roc(response, cox.lp, directions="<")
    roc_results[[var]]$AUC <- c(roc_results[[var]]$AUC, roc_curve$auc)
  }
  # ref: https://github.com/bhklab/survcomp/blob/master/R/iauc.comp.R
  diffs <- c(predict_times[1], predict_times[2:length(predict_times)] - predict_times[1:(length(predict_times) - 1)])
  iauc_values <- c(iauc_values, sum(diffs * roc_results[[var]]$AUC) / max(predict_times))
  lines(predict_times / 365, roc_results[[var]]$AUC, type = "l", col = colors[i], lwd = 2)
}

legend_text <- paste0(variables, " (iAUC=", round(iauc_values, 3), ")")
legend("bottomright", legend = legend_text, col = colors, lty = 1, lwd = 2, bty = "o")


# Brier Score, IBS: ########################################################################################################################
tmax <- 365 * 4 
brier_results <- list()
ibs_values <- c()
# Brier Score calculation function
calculate_brier_score <- function(time_point, surv_object, pred_matrix, time_points) {
  # Prediction probability at a specific point in time
  preds_at_time <- pred_matrix[, which(time_points == time_point)]
  Brier(surv_object, preds_at_time, time_point)
}
for (i in 1:length(variables)) {
  var <- variables[i]
  data <- data_list[[i]]
  
  # cox model
  object <- Surv(data$CT_TIV_death, data$Event_TIV_death)
  fit <- coxph(as.formula(paste("object ~", var)), data = data)
  # predict when the event occurs
  predict_times <- sort(unique(data$CT_TIV_death[data$Event_TIV_death == 1]))
  predict_times <- predict_times[predict_times <= tmax]  # Limit to 4 years
  sp_matrix <- matrix(NA, nrow = nrow(data), ncol = length(predict_times))
  for (j in 1:length(predict_times)) {
    surv_fit <- survfit(fit, newdata = data, se.fit = FALSE)
    surv_probs <- summary(surv_fit, times = predict_times[j])$surv
    sp_matrix[, j] <- surv_probs
  }
  brier_score <- sapply(predict_times, calculate_brier_score, surv_object = object, pred_matrix = sp_matrix, time_points = predict_times)
  brier_results[[var]] <- data.frame(time = predict_times, brier_score = brier_score)
  ibs_values <- c(ibs_values, round(IBS(object, sp_matrix, predict_times), 3))
}

# plotting
plot(NULL, xlim = c(0.35, tmax / 365), ylim = c(0, 1), xlab = "Follow-up Time (years)", ylab = "Brier Score", main = "Brier Scores")
abline(h = 0.25, lty = 3)
for (i in 1:length(variables)) {
  var <- variables[i]
  brier <- brier_results[[var]]
  lines(brier$time / 365, brier$brier_score, type = "l", col = colors[i], lwd = 2)
}

legend_text <- paste0(variables, " (IBS=", ibs_values, ")")
legend("topright", legend = legend_text, col = colors, lty = 1, lwd = 2, bty = "o")
