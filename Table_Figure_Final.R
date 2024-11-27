## Package load
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
library(FSA)  # For Dunn test
library(ggsignif)
library(clipr)
# library(cardx)

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
  filter(Final_Diagnosis_JSK != "Others") %>% 
  filter(Chest_CT_indication_category != "Underlying_resp_ds") %>% # underlying respiratory disease
  mutate(across(everything(), ~replace(., . == "", NA)))


# redefine variables
included <- included %>%
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

#### Table1: Characteristics of study population ####
## check normal distribution
# Ensure all selected continuous variables are numeric and remove rows with NA.
continuous_vars <- c("CT_age", "onset_to_CT", "diagnostic_delay","ALSFRS_filter","Progression_rate",
                     "BMI_filter", "FVC_filter", "CO2_filter", "O2_filter", "Cr_filter", "UA_filter", "CK_filter", "Ca_P",
                     "LVI", "SMI")

# Divide variables into chunks of 4 for manageable plotting
chunks <- split(continuous_vars, ceiling(seq_along(continuous_vars)/4))

# Process each chunk separately
plot_grids <- lapply(chunks, function(vars) {
  # Generate plots for each variable in the chunk
  plot_list <- lapply(vars, function(var) {
    data <- included[[var]]
    
    # Histogram
    p1 <- ggplot(included, aes_string(x = var)) +
      geom_histogram(bins = 30, fill = "blue", color = "black") +
      ggtitle(paste("Histogram of", var))
    
    # QQ Plot
    p2 <- ggplot() +
      geom_qq(aes(sample = data)) +
      geom_qq_line(aes(sample = data)) +
      ggtitle(paste("QQ Plot of", var))
    
    # Shapiro-Wilk test if applicable
    shapiro_test <- if(length(data) > 5000) {
      list(p.value = NA, message = "Sample size too large for Shapiro-Wilk test")
    } else {
      shapiro.test(data)
    }
    
    # Print Shapiro-Wilk test result
    print(paste(var, "Shapiro-Wilk test p-value:", shapiro_test$p.value))
    
    # Return plots
    list(p1, p2)
  })
  
  # Flatten the list of lists into a single list of plots
  plot_list <- do.call(c, plot_list)
  
  # Combine all plots of this chunk into a grid
  do.call(grid.arrange, c(plot_list, ncol = 2))
})


# Table
table1 <- included %>%
  select(Sex, CT_age, onset_to_CT, Final_Diagnosis_JSK, Level_certainty_JSK, Onset_region,
         ALSFRS_filter, Progression_rate,Stage_final,
         BMI_filter, BMI_change,
         CO2_filter, O2_filter, Cr_filter, UA_filter, CK_filter, Ca_P,
         FVC_filter, MIP_filter, MEP_filter, SNIP_filter, PEF_filter,
         NIV, TIV, Gastrostomy, Riluzole, 
         NIV_at_CT, TIV_at_CT, Gastrostomy_at_CT,
         Death, fu_duration,diagnostic_delay,
         Chest_CT_indication_category) %>%
  tbl_summary(
    by = NULL,  # no grouping
    statistic = all_continuous() ~ "{median} ({p25}, {p75})",  # apply statistics value for every continuous variable
    type = all_continuous() ~ "continuous2",  # set type as continuous variable
    digits = all_continuous() ~ 1  
  ) %>%
  add_n()


table1 <- included %>%
  select(Sex, CT_age, onset_to_CT, Final_Diagnosis_JSK, Level_certainty_JSK, Onset_region,
         ALSFRS_filter, Progression_rate,Stage_final,
         BMI_filter, BMI_change,
         CO2_filter, O2_filter, Cr_filter, UA_filter, CK_filter, Ca_P,
         FVC_filter, MIP_filter, MEP_filter, SNIP_filter, PEF_filter,
         NIV, TIV, Gastrostomy, Riluzole, 
         NIV_at_CT, TIV_at_CT, Gastrostomy_at_CT,
         Death, fu_duration,diagnostic_delay,
         Chest_CT_indication_category) %>%
  tbl_summary(
    by = NULL,  # no grouping
    statistic = all_continuous() ~ "{mean} ({sd})",  # apply mean and variance for every continuous variable
    type = all_continuous() ~ "continuous2",  # set type as continuous variable
    digits = all_continuous() ~ 1  
  ) %>%
  add_n()

table1  



table1_tibble <- as_tibble(table1)
# library(clipr)
write_clip(table1_tibble)


#### Figure 2: Boxplot by stage ####
fig2 <- chestCT %>%
  filter(Diagnosis_check == 'ALS') %>%
  filter(Hospital_ID != 00000000) %>% # unknown onset
  filter(Time_TIV_death < 365*20) %>% # remove extreme long survival to satisfy the proportional hazard assumption of RMI
  filter(Final_Diagnosis_JSK != "Others") %>% # FAS. There is no evidence for the other segments
  filter(Chest_CT_indication_category != "Underlying_resp_ds") %>% # underlying respiratory disease
  mutate(FVC_filter = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), FVC, NA_real_)) %>%
  mutate(across(everything(), ~replace(., . == "", NA))) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>%
  select(LVI, SMI, FVC_filter, Stage_final) %>%
  rename(FVC = FVC_filter,
         RMI = SMI) 

# *all-in-one ####
variable_order <- c("LVI", "RMI", "FVC")
fig2 <- fig2 %>%
  pivot_longer(cols = -Stage_final, names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = variable_order))


p <- ggplot(fig2, aes(x = factor(Stage_final), y = value)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~variable, scales = "free_y", nrow = 2, ncol = 3, as.table = TRUE) +
  labs(x = "", y = "") +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 8), 
    panel.grid.major = element_blank(),    
    panel.grid.minor = element_blank(),   
    axis.ticks = element_line(color = "black", size = 0.5),  
    axis.line = element_line(color = "black"),  
    panel.spacing = unit(0.8, "cm"),  # Adjust spacing between panels
    axis.text.x = element_text(size = 7),  # Adjust x-axis text size
    axis.text.y = element_text(size = 7)   # Adjust y-axis text size
  )

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 6.5, height = 3))
print(doc, target = "plot.pptx")



#### Table S1 (A) ####
by_stage <-  chestCT %>%
  filter(Diagnosis_check == 'ALS') %>%
  filter(Hospital_ID != 00000000) %>% # unknown onset
  filter(Time_TIV_death < 365*20) %>% # remove extreme long survival to satisfy the proportional hazard assumption of RMI
  filter(Final_Diagnosis_JSK != "Others") %>% # FAS. There is no evidence for the other segments
  filter(Chest_CT_indication_category != "Underlying_resp_ds") %>% # underlying respiratory disease
  mutate(FVC_filter = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), FVC, NA_real_)) %>%
  mutate(across(everything(), ~replace(., . == "", NA))) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>%
  select(LVI, SMI, FVC_filter, Muscle_HU, SFI, SF_HU, Stage_final) %>%
  rename(FVC = FVC_filter,
         RMI = SMI) %>%
  tbl_summary(
    by = Stage_final,  # no grouping
    statistic = all_continuous() ~ "{median} ({p25}, {p75})",  # apply statistics value for every continuous variable
    type = all_continuous() ~ "continuous2",  # set type as continuous variable
    digits = all_continuous() ~ 0  
  ) %>%
  add_n()%>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 3))

by_stage  

table_stage_tibble <- as_tibble(by_stage)
# library(clipr)
write_clip(table_stage_tibble)


#### Fig S1: Changes in Respiratory and Body Composition Parameters According to King’s Stage ####

figS1 <- chestCT %>%
  filter(Diagnosis_check == 'ALS') %>%
  filter(Hospital_ID != 00000000) %>% # unknown onset
  filter(Time_TIV_death < 365*20) %>% # remove extreme long survival to satisfy the proportional hazard assumption of RMI
  filter(Final_Diagnosis_JSK != "Others") %>% # FAS. There is no evidence for the other segments
  filter(Chest_CT_indication_category != "Underlying_resp_ds") %>% # underlying respiratory disease
  mutate(FVC_filter = if_else(between(as.numeric(PFT_date - Date_CT, units = "days"), -90, 90), FVC, NA_real_)) %>%
  mutate(across(everything(), ~replace(., . == "", NA))) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>%
  select(Muscle_HU, SFI, SF_HU, Stage_final) %>%
  mutate(SF_HU = if_else(Stage_final == '4b' & SF_HU < -250 | SF_HU >0, NA_real_, SF_HU)) %>%
  mutate(SFI = if_else(Stage_final == '1' & SFI > 3000, NA_real_, SFI)) 


# *all-in-on ####
variable_order <- c("Muscle_HU", "SFI", "SF_HU")
figS1 <- figS1 %>%
  pivot_longer(cols = -Stage_final, names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = variable_order))


p <- ggplot(figS1, aes(x = factor(Stage_final), y = value)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~variable, scales = "free_y", nrow = 2, ncol = 3, as.table = TRUE) +
  labs(x = "", y = "") +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 8), 
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.line = element_line(color = "black"), 
    panel.spacing = unit(0.8, "cm"),  # Adjust spacing between panels
    axis.text.x = element_text(size = 7),  # Adjust x-axis text size
    axis.text.y = element_text(size = 7)   # Adjust y-axis text size
  )

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 6.5, height = 3))
print(doc, target = "plot.pptx")


#### Figure 3: Correlation Plots Between Clinical Variables and LVI or RMI  ####

## * LVI####
fig_LVI <- included %>%
  select(LVI, FVC_filter, MIP_filter, MEP_filter, SNIP_filter,
         ALSFRS_filter, Bulbar_filter, Fine_filter, Gross_filter, Resp_filter, Progression_rate,
         Wt_pre, Wt_filter, Wt_change, BMI_filter, BMI_change,
         CO2_filter, O2_filter, Cr_filter, UA_filter, CK_filter, Ca_filter, P_filter, Ca_P,
         severe_bulbar_dysfunction) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA)))  # Convert Inf to NA (Inf occurs when Height is missing)

## LVI - ALSFRS
model <- lm(ALSFRS_filter ~ LVI, data = fig_LVI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]


fig_LVI_clean <- fig_LVI %>%
  drop_na(LVI, ALSFRS_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, ALSFRS_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) + 
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max(fig_LVI_clean$ALSFRS_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")


p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## LVI - FVC 
model <- lm(FVC_filter ~ LVI, data = fig_LVI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]


fig_LVI_clean <- fig_LVI %>%
  drop_na(LVI, FVC_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, FVC_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "FVC") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) + 
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, by = 25)) +
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max(fig_LVI_clean$FVC_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.2e", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p

doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")


## LVI - O2 
model <- lm(O2_filter ~ LVI, data = fig_LVI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_LVI_clean <- fig_LVI %>%
  drop_na(LVI, O2_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, O2_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "PaO2") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) + 
  scale_y_continuous(limits = c(min(fig_LVI_clean$O2_filter), max(fig_LVI_clean$O2_filter, na.rm = TRUE)), breaks = seq(0, 120, by = 20)) +
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max(fig_LVI_clean$O2_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")


## LVI - CO2 
model <- lm(CO2_filter ~ LVI, data = fig_LVI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_LVI_clean <- fig_LVI %>%
  drop_na(LVI, CO2_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, CO2_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "PaCO2") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) + 
  scale_y_continuous(limits = c(30, max(fig_LVI_clean$CO2_filter, na.rm = TRUE)), breaks = seq(0, 75, by = 15)) +
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max(fig_LVI_clean$CO2_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.4f", -sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")






## * RMI####
fig_SMI <- included %>%
  select(SMI, FVC_filter, MIP_filter, MEP_filter, SNIP_filter,
         ALSFRS_filter, Bulbar_filter, Fine_filter, Gross_filter, Resp_filter, Progression_rate,
         Wt_pre, Wt_filter, Wt_change, BMI_filter, BMI_change,
         CO2_filter, O2_filter, Cr_filter, UA_filter, CK_filter, Ca_filter, P_filter, Ca_P,
         severe_bulbar_dysfunction) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA)))  # Convert Inf to NA (Inf occurs when Height is missing)


## RMI - ALSFRS
model <- lm(ALSFRS_filter ~ SMI, data = fig_SMI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI %>%
  drop_na(SMI, ALSFRS_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, ALSFRS_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1000, by = 250)) + 
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max(fig_SMI_clean$ALSFRS_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")


## RMI - FVC
model <- lm(FVC_filter ~ SMI, data = fig_SMI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI %>%
  drop_na(SMI, FVC_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, FVC_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "FVC") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(250, 1000, by = 250)) + 
  scale_y_continuous(limits = c(0, max(fig_SMI_clean$FVC_filter, na.rm = TRUE)), breaks = seq(0, 100, by = 25)) +
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max(fig_SMI_clean$FVC_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.2e", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")




## RMI - O2 
model <- lm(O2_filter ~ SMI, data = fig_SMI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI %>%
  drop_na(SMI, O2_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, O2_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "PaO2") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1000, by = 250)) + 
  scale_y_continuous(limits = c(min(fig_SMI_clean$O2_filter, na.rm = TRUE), max(fig_SMI_clean$O2_filter, na.rm = TRUE)), breaks = seq(0, 120, by = 20)) +
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max(fig_SMI_clean$O2_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## RMI - CO2 
model <- lm(CO2_filter ~ SMI, data = fig_SMI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI %>%
  drop_na(SMI, CO2_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, CO2_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "PaCO2") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1000, by = 250)) + 
  scale_y_continuous(limits = c(30, max(fig_SMI_clean$CO2_filter, na.rm = TRUE)), breaks = seq(0, 75, by = 15)) +
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max(fig_SMI_clean$CO2_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.4f", -sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")





#### Fig S2: Correlation Plots Between Blood Biomarkers (Cr, CK, UA, Ca/P) and LVI or RMI  ####
## LVI - Cr 
model <- lm(Cr_filter ~ LVI, data = fig_LVI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_LVI_clean <- fig_LVI %>%
  drop_na(LVI, Cr_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, Cr_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "Cr") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) + 
  scale_y_continuous(limits = c(min(fig_LVI_clean$Cr_filter), max(fig_LVI_clean$Cr_filter, na.rm = TRUE)), breaks = seq(0, 1.5, by = 0.5)) +
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max(fig_LVI_clean$Cr_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.2e", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## LVI - CK 
model <- lm(CK_filter ~ LVI, data = fig_LVI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_LVI_clean <- fig_LVI %>%
  drop_na(LVI, CK_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, CK_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "CK") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) +
  coord_cartesian(ylim = c(0, 800)) +
  scale_y_continuous(limits = c(-100, max(fig_LVI_clean$CK_filter, na.rm = TRUE)), breaks = seq(0, 600, by = 300)) +
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = 800, 
           label = sprintf("R = %.2f, p = %.2e", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## LVI - UA 
model <- lm(UA_filter ~ LVI, data = fig_LVI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_LVI_clean <- fig_LVI %>%
  drop_na(LVI, UA_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, UA_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "UA") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by =500)) + 
  scale_y_continuous(limits = c(min(fig_LVI_clean$UA_filter), max(fig_LVI_clean$UA_filter, na.rm = TRUE))) +
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max(fig_LVI_clean$UA_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.3f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## LVI - Ca/P 
model <- lm(Ca_P ~ LVI, data = fig_LVI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_LVI_clean <- fig_LVI %>%
  drop_na(LVI, Ca_P)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, Ca_P)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "Ca/P") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) + 
  scale_y_continuous(limits = c(min(fig_LVI_clean$Ca_P), max(fig_LVI_clean$Ca_P, na.rm = TRUE))) +
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max(fig_LVI_clean$Ca_P, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.3f", -sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")





## RMI - Cr 
model <- lm(Cr_filter ~ SMI, data = fig_SMI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI %>%
  drop_na(SMI, Cr_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, Cr_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "Cr") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1250, by = 250)) + 
  scale_y_continuous(limits = c(min(fig_SMI_clean$Cr_filter), max(fig_SMI_clean$Cr_filter, na.rm = TRUE)), breaks = seq(0, 1.5, by = 0.5)) +
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max(fig_SMI_clean$Cr_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.2e", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## RMI - CK 
model <- lm(CK_filter ~ SMI, data = fig_SMI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI %>%
  drop_na(SMI, CK_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, CK_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "CK") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1250, by = 250)) + 
  coord_cartesian(ylim = c(0, 800)) +
  scale_y_continuous(limits = c(-100, max(fig_SMI_clean$CK_filter, na.rm = TRUE)), breaks = seq(0, 600, by = 300)) +
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = 800, 
           label = sprintf("R = %.2f, p = %.2e", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## RMI - UA 
model <- lm(UA_filter ~ SMI, data = fig_SMI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI %>%
  drop_na(SMI, UA_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, UA_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "UA") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1250, by = 250)) + 
  scale_y_continuous(limits = c(min(fig_SMI_clean$UA_filter), max(fig_SMI_clean$UA_filter, na.rm = TRUE))) +
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max(fig_SMI_clean$UA_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.2e", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## RMI - Ca/P 
model <- lm(Ca_P ~ SMI, data = fig_SMI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI %>%
  drop_na(SMI, Ca_P)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, Ca_P)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "Ca/P") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1250, by = 250)) + 
  scale_y_continuous(limits = c(min(fig_SMI_clean$Ca_P), max(fig_SMI_clean$Ca_P, na.rm = TRUE))) +
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max(fig_SMI_clean$Ca_P, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.3f", -sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")




#### Table 2: Cox regression ####
cox <- included %>%
  filter(Stage_final %in% c('1','2','3')) 

cox_TIV_death <- cox %>%
  mutate(FVC_supine_pred = round(FVC_supine/(FVC_measured/(FVC_filter/100))*100,0)) %>%
  select(Time_TIV_death, Event_TIV_death, 
         LVI, SMI, SFI, Muscle_HU, SF_HU, volume,
         onset_age, Sex, Stage_final, BMI_filter, BMI_change, Cr_filter, CK_filter, UA_filter, Ca_P, Slope2023, Preslope,
         FVC_filter, MIP_filter, MEP_filter, SNIP_filter, PEF_filter, Resp_filter,
         CT_age,onset_to_CT,
         Onset_region, Final_Diagnosis_JSK, Level_certainty_JSK, Wt_pre, Wt_filter, CO2_filter, O2_filter,
         Gastrostomy, NIV, Gastrostomy_at_CT, NIV_at_CT,
         ALSFRS_filter,
         Riluzole, Edaravone, Progression_rate,
         diagnostic_delay,
         FVC_supine_pred) %>%
  rename(RMI = SMI)

# Univariable####
model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ SFI, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Muscle_HU, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)


model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ SF_HU, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ CT_age, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)


model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ onset_to_CT, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ diagnostic_delay, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Final_Diagnosis_JSK, data = cox_TIV_death)

cox_TIV_death$Onset_region <- factor(cox_TIV_death$Onset_region)
cox_TIV_death$Onset_region <- relevel(cox_TIV_death$Onset_region, ref = "Limb")
model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Onset_region, data = cox_TIV_death) # R 3.96e-05
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ ALSFRS_filter, data = cox_TIV_death)
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Progression_rate, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Stage_final, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ BMI_filter, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ BMI_change, data = cox_TIV_death)
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ CO2_filter, data = cox_TIV_death)
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ O2_filter, data = cox_TIV_death)
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Cr_filter, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model) # Violation of the proportional hazards assumption

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ CK_filter, data = cox_TIV_death)
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)


model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ UA_filter, data = cox_TIV_death)
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Ca_P, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ FVC_filter, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ FVC_supine_pred, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ MIP_filter, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ MEP_filter, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ SNIP_filter, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ PEF_filter, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ LVI, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ RMI, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model)


coxph(Surv(Time_TIV_death, Event_TIV_death) ~ NIV_at_CT, data = cox_TIV_death)
coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Gastrostomy_at_CT, data = cox_TIV_death)
coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Riluzole, data = cox_TIV_death)




# Multivariable####
model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ CT_age + Sex + SFI, data = cox_TIV_death)


model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ CT_age + Sex + Onset_region + Stage_final + SFI, data = cox_TIV_death)

vif(model)
cox.zph(model)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ CT_age + Sex + Onset_region + BMI_filter + Stage_final + Muscle_HU, data = cox_TIV_death)
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)


model1 <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex + CT_age + Progression_rate + BMI_filter + LVI, data = cox_TIV_death)
model_summary <- tidy(model1, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model1)
vif(model1)

model2 <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex + CT_age + Progression_rate + BMI_filter + RMI, data = cox_TIV_death)
model_summary <- tidy(model2, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model2)
vif(model2)

model3 <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex + CT_age + BMI_filter + FVC_filter, data = cox_TIV_death)
model_summary <- tidy(model3, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model3)
vif(model3)



#### Table S2: in subgroup with chest CT within 1 year of onset ####
cox <- included %>%
  filter(Stage_final %in% c('1','2','3')) %>%
  filter(Chest_CT_indication_category == 'Etiology_wu')
cox_TIV_death <- cox %>%
  mutate(FVC_supine_pred = round(FVC_supine/(FVC_measured/(FVC_filter/100))*100,0)) %>%
  select(Time_TIV_death, Event_TIV_death, 
         LVI, SMI, SFI, Muscle_HU, SF_HU, volume,
         onset_age, Sex, Stage_final, BMI_filter, BMI_change, Cr_filter, CK_filter, UA_filter, Ca_P, Slope2023, Preslope,
         FVC_filter, MIP_filter, MEP_filter, SNIP_filter, PEF_filter, Resp_filter,
         CT_age,onset_to_CT,
         Onset_region, Final_Diagnosis_JSK, Level_certainty_JSK, Wt_pre, Wt_filter, CO2_filter, O2_filter,
         Gastrostomy, NIV, Gastrostomy_at_CT, NIV_at_CT,
         ALSFRS_filter,
         Riluzole, Edaravone, Progression_rate,
         diagnostic_delay,
         FVC_supine_pred) %>%
  rename(RMI = SMI)
# use cox_TIV_death data from Table 2 code
# Multivariable
model1 <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex + CT_age + Stage_final + BMI_filter + LVI, data = cox_TIV_death)
model_summary <- tidy(model1, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model1)
vif(model1)

model2 <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex + CT_age +  Stage_final + BMI_filter + RMI, data = cox_TIV_death)
model_summary <- tidy(model2, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)
cox.zph(model2)
vif(model2)

model3 <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex + CT_age + Stage_final +  BMI_filter + FVC_filter, data = cox_TIV_death)
model_summary <- tidy(model3, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)




#### Table S3: in subgroup with bulbar symptom ####
cox <- included %>%
  filter(Stage_final %in% c('1','2','3')) %>%
  filter(Onset_region == 'B')

cox_TIV_death <- cox %>%
  select(Time_TIV_death, Event_TIV_death, 
         LVI, SMI, SFI, Muscle_HU, SF_HU, volume,
         onset_age, Sex, Stage_final, BMI_filter, BMI_change, Cr_filter, CK_filter, UA_filter, Ca_P, Slope2023, Preslope,
         FVC_filter,
         CT_age,onset_to_CT,
         Onset_region, Final_Diagnosis_JSK, Level_certainty_JSK, Wt_pre, Wt_filter, CO2_filter, O2_filter,
         Gastrostomy, NIV,
         ALSFRS_filter,
         Riluzole, Edaravone, Progression_rate, diagnostic_delay) %>%
  rename(RMI = SMI)


# Univariable
model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ CT_age, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ diagnostic_delay, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Stage_final, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ BMI_filter, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ FVC_filter, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ LVI, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ RMI, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ SFI, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Muscle_HU, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ SF_HU, data = cox_TIV_death) 
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)


# Multivariable
coxph(Surv(Time_TIV_death, Event_TIV_death) ~ CT_age + Stage_final + RMI, data = cox_TIV_death)


model1 <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex + CT_age + Stage_final + BMI_filter+ LVI, data = cox_TIV_death)
model_summary <- tidy(model1, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model2 <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex + CT_age + Stage_final + BMI_filter+ RMI, data = cox_TIV_death)
model_summary <- tidy(model2, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)

model3 <- coxph(Surv(Time_TIV_death, Event_TIV_death) ~ Sex + CT_age + Stage_final + BMI_filter + FVC_filter, data = cox_TIV_death)
model_summary <- tidy(model3, conf.int = TRUE, exponentiate = TRUE) 
model_summary <- model_summary %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))
print(model_summary)



cox.zph(model1)
vif(model1)





#### Fig 4: Survival curve #### 
cox <- included %>%
  filter(Stage_final %in% c('1','2','3'))

# Calculate the medians of LVI and SMI
median_LVI <- median(cox$LVI, na.rm = TRUE)
median_SMI <- median(cox$SMI, na.rm = TRUE)
median_FVC_filter <- median(cox$FVC_filter, na.rm = TRUE)

# Create binary groups based on medians
cox_curve <- cox %>%
  mutate(LVI_group = ifelse(LVI >= median_LVI, "High LVI", "Low LVI"),
         SMI_group = ifelse(SMI >= median_SMI, "High RMI", "Low RMI"),
         FVC_group = ifelse(FVC_filter >= median_FVC_filter, "High FVC", "Low FVC")) %>%
  mutate(Time_TIV_death = Time_TIV_death/365*12)




# LVI
km_LVI <- survfit(Surv(Time_TIV_death, Event_TIV_death) ~ LVI_group, data= cox_curve)
survdiff(Surv(Time_TIV_death, Event_TIV_death) ~ LVI_group, data= cox_curve)

p <- ggsurvplot(km_LVI, 
                conf.int = FALSE,  
                legend.labs = c("High LVI", "Low LVI"),
                legend = c(0.8, 0.9),
                font.legend = 8,
                font.x = 6, 
                font.y = 6,
                font.tickslab = c(13, 'black'),
                ylab = "Survival probability", 
                xlab = "Time (Month)",
                xlim = c(0,125),
                break.time.by = 25,
                legend.title = "",
                palette = c("salmon", "royalblue"),
                size = 0.5,
                pval = TRUE,
                pval.size = 3,
                pval.coord = c(10, 0.1),
                tables.theme = theme_cleantable(),
                ggtheme = theme_minimal() +
                  theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(color = "black")), 
                risk.table = TRUE,
                risk.table.height = 0.2,
                risk.table.y.text = FALSE, 
                risk.table.title = "Number at risk",
                risk.table.fontsize = 2.1)

p$plot <- p$plot + 
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6))  

p$ table <- p$table + theme(plot.title = element_text(size=8))

p
surv_plot <- p$plot
risk_table <- p$table + theme_cleantable()

# Align the plots based on their x-axis
aligned_plots <- align_plots(surv_plot, risk_table, align = 'v', axis = 'lr')

# Combine the aligned plots
combined_plot <- plot_grid(aligned_plots[[1]], aligned_plots[[2]], ncol = 1, rel_heights = c(3, 1))

# Create PowerPoint and add the combined plot
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = combined_plot), location = ph_location(width = 3.15, height = 2.8))
print(doc, target = "plot.pptx")


# RMI
km_SMI <- survfit(Surv(Time_TIV_death, Event_TIV_death) ~ SMI_group, data= cox_curve)
survdiff(Surv(Time_TIV_death, Event_TIV_death) ~ SMI_group, data= cox_curve)

p <- ggsurvplot(km_SMI, 
                conf.int = FALSE,  
                legend.labs = c("High RMI", "Low RMI"),
                legend = c(0.8, 0.9),
                font.legend = 8,
                font.x = 6, 
                font.y = 6,
                font.tickslab = c(13, 'black'),
                ylab = "Survival probability", 
                xlab = "Time (Month)",
                xlim = c(0,125),
                break.time.by = 25,
                legend.title = "",
                palette = c("salmon", "royalblue"),
                size = 0.5,
                pval = TRUE,
                pval.size = 3,
                pval.coord = c(10, 0.1),
                tables.theme = theme_cleantable(),
                ggtheme = theme_minimal() +
                  theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(color = "black")), 
                risk.table = TRUE,
                risk.table.height = 0.2,
                risk.table.y.text = FALSE, 
                risk.table.title = "Number at risk",
                risk.table.fontsize = 2.1)

p$plot <- p$plot + 
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6))  

p$ table <- p$table + theme(plot.title = element_text(size=8))

p
surv_plot <- p$plot
risk_table <- p$table + theme_cleantable()

# Align the plots based on their x-axis
aligned_plots <- align_plots(surv_plot, risk_table, align = 'v', axis = 'lr')

# Combine the aligned plots
combined_plot <- plot_grid(aligned_plots[[1]], aligned_plots[[2]], ncol = 1, rel_heights = c(3, 1))

# Create PowerPoint and add the combined plot
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = combined_plot), location = ph_location(width = 3.15, height = 2.8))
print(doc, target = "plot.pptx")



#### Fig S4: Survival curve in subgroup with chest CT within 1 year of onset#### 
cox <- included %>%
  filter(Stage_final %in% c('1','2','3')) %>%
  filter(Chest_CT_indication_category == 'Etiology_wu')

# Calculate the medians of LVI and SMI
median_LVI <- median(cox$LVI, na.rm = TRUE)
median_SMI <- median(cox$SMI, na.rm = TRUE)
median_FVC_filter <- median(cox$FVC_filter, na.rm = TRUE)

# Create binary groups based on medians
cox_curve <- cox %>%
  mutate(LVI_group = ifelse(LVI >= median_LVI, "High LVI", "Low LVI"),
         SMI_group = ifelse(SMI >= median_SMI, "High RMI", "Low RMI"),
         FVC_group = ifelse(FVC_filter >= median_FVC_filter, "High FVC", "Low FVC")) %>%
  mutate(Time_TIV_death = Time_TIV_death/365*12)


# LVI
km_LVI <- survfit(Surv(Time_TIV_death, Event_TIV_death) ~ LVI_group, data= cox_curve)
survdiff(Surv(Time_TIV_death, Event_TIV_death) ~ LVI_group, data= cox_curve)

p <- ggsurvplot(km_LVI, 
                conf.int = FALSE,  
                legend.labs = c("High LVI", "Low LVI"),
                legend = c(0.8, 0.9),
                font.legend = 8,
                font.x = 6, 
                font.y = 6,
                font.tickslab = c(13, 'black'),
                ylab = "Survival probability", 
                xlab = "Time (Month)",
                xlim = c(0,125),
                break.time.by = 25,
                legend.title = "",
                palette = c("salmon", "royalblue"),
                size = 0.5,
                pval = TRUE,
                pval.size = 3,
                pval.coord = c(10, 0.1),
                tables.theme = theme_cleantable(),
                ggtheme = theme_minimal() +
                  theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(color = "black")), 
                risk.table = TRUE,
                risk.table.height = 0.2,
                risk.table.y.text = FALSE, 
                risk.table.title = "Number at risk",
                risk.table.fontsize = 2.1)

p$plot <- p$plot + 
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6))  

p$ table <- p$table + theme(plot.title = element_text(size=8))

p
surv_plot <- p$plot
risk_table <- p$table + theme_cleantable()

# Align the plots based on their x-axis
aligned_plots <- align_plots(surv_plot, risk_table, align = 'v', axis = 'lr')

# Combine the aligned plots
combined_plot <- plot_grid(aligned_plots[[1]], aligned_plots[[2]], ncol = 1, rel_heights = c(3, 1))

# Create PowerPoint and add the combined plot
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = combined_plot), location = ph_location(width = 3.15, height = 2.8))
print(doc, target = "plot.pptx")


# RMI
km_SMI <- survfit(Surv(Time_TIV_death, Event_TIV_death) ~ SMI_group, data= cox_curve)
survdiff(Surv(Time_TIV_death, Event_TIV_death) ~ SMI_group, data= cox_curve)

p <- ggsurvplot(km_SMI, 
                conf.int = FALSE,  
                legend.labs = c("High RMI", "Low RMI"),
                legend = c(0.8, 0.9),
                font.legend = 8,
                font.x = 6, 
                font.y = 6,
                font.tickslab = c(13, 'black'),
                ylab = "Survival probability", 
                xlab = "Time (Month)",
                xlim = c(0,125),
                break.time.by = 25,
                legend.title = "",
                palette = c("salmon", "royalblue"),
                size = 0.5,
                pval = TRUE,
                pval.size = 3,
                pval.coord = c(10, 0.1),
                tables.theme = theme_cleantable(),
                ggtheme = theme_minimal() +
                  theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(color = "black")), 
                risk.table = TRUE,
                risk.table.height = 0.2,
                risk.table.y.text = FALSE, 
                risk.table.title = "Number at risk",
                risk.table.fontsize = 2.1)

p$plot <- p$plot + 
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6))  

p$ table <- p$table + theme(plot.title = element_text(size=8))

p
surv_plot <- p$plot
risk_table <- p$table + theme_cleantable()

# Align the plots based on their x-axis
aligned_plots <- align_plots(surv_plot, risk_table, align = 'v', axis = 'lr')

# Combine the aligned plots
combined_plot <- plot_grid(aligned_plots[[1]], aligned_plots[[2]], ncol = 1, rel_heights = c(3, 1))

# Create PowerPoint and add the combined plot
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = combined_plot), location = ph_location(width = 3.15, height = 2.8))
print(doc, target = "plot.pptx")

#### Fig S7: correlation with respiratory ALSFRS-R depending on severe bulbar sx ####
fig_LVI <- included %>%
  select(LVI, FVC_filter, MIP_filter, MEP_filter, SNIP_filter,
         ALSFRS_filter, Bulbar_filter, Fine_filter, Gross_filter, Resp_filter, Progression_rate,
         Wt_pre, Wt_filter, Wt_change, BMI_filter, BMI_change,
         CO2_filter, O2_filter, Cr_filter, UA_filter, CK_filter, Ca_filter, P_filter, Ca_P,
         severe_bulbar_dysfunction) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA)))  # Convert Inf to NA (Inf occurs when Height is missing)

fig_SMI <- included %>%
  select(SMI, FVC_filter, MIP_filter, MEP_filter, SNIP_filter,
         ALSFRS_filter, Bulbar_filter, Fine_filter, Gross_filter, Resp_filter, Progression_rate,
         Wt_pre, Wt_filter, Wt_change, BMI_filter, BMI_change,
         CO2_filter, O2_filter, Cr_filter, UA_filter, CK_filter, Ca_filter, P_filter, Ca_P,
         severe_bulbar_dysfunction) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA)))  # Convert Inf to NA (Inf occurs when Height is missing)

fig_FVC <- included %>%
  select(FVC_filter, MIP_filter, MEP_filter, SNIP_filter,
         ALSFRS_filter, Bulbar_filter, Fine_filter, Gross_filter, Resp_filter, Progression_rate,
         Wt_pre, Wt_filter, Wt_change, BMI_filter, BMI_change,
         CO2_filter, O2_filter, Cr_filter, UA_filter, CK_filter, Ca_filter, P_filter, Ca_P,
         severe_bulbar_dysfunction,
         Hospital_ID) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA)))  # Convert Inf to NA (Inf occurs when Height is missing)


## with severe_bulbar dysfunction 
fig_LVI_bulbar <- fig_LVI %>%
  filter(severe_bulbar_dysfunction == 1)
fig_SMI_bulbar <- fig_SMI %>%
  filter(severe_bulbar_dysfunction == 1)
fig_FVC_bulbar <- fig_FVC %>%
  filter(severe_bulbar_dysfunction == 1)
## without severe bulbar dysfunction 
fig_LVI_no_bulbar <- fig_LVI %>%
  filter(severe_bulbar_dysfunction == 0)
fig_SMI_no_bulbar <- fig_SMI %>%
  filter(severe_bulbar_dysfunction == 0)
fig_FVC_no_bulbar <- fig_FVC %>%
  filter(severe_bulbar_dysfunction == 0)


## * LVI, total ####
model <- lm(Resp_filter ~ LVI, data = fig_LVI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_LVI_clean <- fig_LVI %>%
  drop_na(LVI, Resp_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, Resp_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "Respiratory ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) + 
  scale_y_continuous(limits = c(min(fig_LVI_clean$Resp_filter), max(fig_LVI_clean$Resp_filter, na.rm = TRUE)*1.2), breaks = seq(0, 12, by = 3)) +
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max(fig_LVI_clean$Resp_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = -1, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")







## LVI, bulbar dysfunction 
model <- lm(Resp_filter ~ LVI, data = fig_LVI_bulbar)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_LVI_clean <- fig_LVI_bulbar %>%
  drop_na(LVI, Resp_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, Resp_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "Respiratory ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) + 
  coord_cartesian(ylim = c(0, max(fig_LVI_clean$Resp_filter, na.rm = TRUE) * 1.2)) +
  scale_y_continuous(breaks = seq(0, 12, by = 3)) +
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max(fig_LVI_clean$Resp_filter, na.rm = TRUE) * 1.15, 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## LVI, no bulbar dysfunction
model <- lm(Resp_filter ~ LVI, data = fig_LVI_no_bulbar)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_LVI_clean <- fig_LVI_no_bulbar %>%
  drop_na(LVI, Resp_filter)

min_LVI <- min(fig_LVI_clean$LVI, na.rm = TRUE)
max_LVI <- max(fig_LVI_clean$LVI, na.rm = TRUE)
max_Resp <- max(fig_LVI_clean$Resp_filter, na.rm = TRUE)

p <- ggplot(fig_LVI_clean, aes(LVI, Resp_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "LVI", y = "Respiratory ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_LVI, max_LVI), breaks = seq(0, 2000, by = 500)) + 
  coord_cartesian(ylim = c(0, max_Resp * 1.2)) + 
  scale_y_continuous(breaks = seq(0, 12, by = 3)) + 
  annotate("text", x = min(fig_LVI_clean$LVI, na.rm = TRUE), 
           y = max_Resp * 1.15, 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")





## * RMI, total####
model <- lm(Resp_filter ~ SMI, data = fig_SMI)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI %>%
  drop_na(SMI, Resp_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, Resp_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "Respiratory ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1000, by = 250)) + 
  scale_y_continuous(limits = c(min(fig_SMI_clean$Resp_filter), max(fig_SMI_clean$Resp_filter, na.rm = TRUE)*1.2), breaks = seq(0, 12, by = 3)) +
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max(fig_SMI_clean$Resp_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = -1, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## RMI, bulbar dysfunction
model <- lm(Resp_filter ~ SMI, data = fig_SMI_bulbar)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI_bulbar %>%
  drop_na(SMI, Resp_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, Resp_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "Respiratory ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1000, by = 250)) + 
  coord_cartesian(ylim = c(0, max(fig_SMI_clean$Resp_filter, na.rm = TRUE) * 1.2)) +
  scale_y_continuous(breaks = seq(0, 12, by = 3)) +
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max(fig_SMI_clean$Resp_filter, na.rm = TRUE) * 1.15, 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## RMI, no bulbar dysfunction 
model <- lm(Resp_filter ~ SMI, data = fig_SMI_no_bulbar)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_SMI_clean <- fig_SMI_no_bulbar %>%
  drop_na(SMI, Resp_filter)

min_SMI <- min(fig_SMI_clean$SMI, na.rm = TRUE)
max_SMI <- max(fig_SMI_clean$SMI, na.rm = TRUE)
max_Resp <- max(fig_SMI_clean$Resp_filter, na.rm = TRUE)

p <- ggplot(fig_SMI_clean, aes(SMI, Resp_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "RMI", y = "Respiratory ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_SMI, max_SMI), breaks = seq(0, 1000, by = 250)) + 
  coord_cartesian(ylim = c(0, max_Resp * 1.2)) + 
  scale_y_continuous(breaks = seq(0, 12, by = 3)) + 
  annotate("text", x = min(fig_SMI_clean$SMI, na.rm = TRUE), 
           y = max_Resp * 1.15, 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## * FVC, total ####
model <- lm(Resp_filter ~ FVC_filter, data = fig_FVC)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_FVC_clean <- fig_FVC %>%
  drop_na(FVC_filter, Resp_filter)

min_FVC <- min(fig_FVC_clean$FVC_filter, na.rm = TRUE)
max_FVC <- max(fig_FVC_clean$FVC_filter, na.rm = TRUE)

p <- ggplot(fig_FVC_clean, aes(FVC_filter, Resp_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "FVC", y = "Respiratory ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_FVC, max_FVC), breaks = seq(0, max_FVC, by = 25)) + 
  scale_y_continuous(limits = c(min(fig_FVC_clean$Resp_filter), max(fig_FVC_clean$Resp_filter, na.rm = TRUE)*1.2), breaks = seq(0, 12, by = 3)) +
  annotate("text", x = min(fig_FVC_clean$FVC_filter, na.rm = TRUE), 
           y = max(fig_FVC_clean$Resp_filter, na.rm = TRUE), 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = -1, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")



## FVC, bulbar dysfunction
model <- lm(Resp_filter ~ FVC_filter, data = fig_FVC_bulbar)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_FVC_clean <- fig_FVC_bulbar %>%
  drop_na(FVC_filter, Resp_filter)

min_FVC <- min(fig_FVC_clean$FVC_filter, na.rm = TRUE)
max_FVC <- max(fig_FVC_clean$FVC_filter, na.rm = TRUE)

p <- ggplot(fig_FVC_clean, aes(FVC_filter, Resp_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "FVC", y = "Respiratory ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_FVC, max_FVC), breaks = seq(0, max_FVC, by = 25)) + 
  coord_cartesian(ylim = c(0, max(fig_FVC_clean$Resp_filter, na.rm = TRUE) * 1.2)) +
  scale_y_continuous(breaks = seq(0, 12, by = 3)) +
  annotate("text", x = min(fig_FVC_clean$FVC_filter, na.rm = TRUE), 
           y = max(fig_FVC_clean$Resp_filter, na.rm = TRUE) * 1.15, 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")

## FVC, no bulbar dysfunction 
model <- lm(Resp_filter ~ FVC_filter, data = fig_FVC_no_bulbar)
summary_model <- summary(model)
r_value <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

fig_FVC_clean <- fig_FVC_no_bulbar %>%
  drop_na(FVC_filter, Resp_filter)

min_FVC <- min(fig_FVC_clean$FVC_filter, na.rm = TRUE)
max_FVC <- max(fig_FVC_clean$FVC_filter, na.rm = TRUE)
max_Resp <- max(fig_FVC_clean$Resp_filter, na.rm = TRUE)

p <- ggplot(fig_FVC_clean, aes(FVC_filter, Resp_filter)) + 
  geom_point(size = 0.2) + 
  geom_smooth(method = 'lm', se = TRUE, color = 'black', fill = 'grey', linewidth = 0.5) + 
  labs(x = "FVC", y = "Respiratory ALSFRS-R") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) +
  scale_x_continuous(limits = c(min_FVC, max_FVC), breaks = seq(0, max_FVC, by = 25)) + 
  coord_cartesian(ylim = c(0, max_Resp * 1.2)) + 
  scale_y_continuous(breaks = seq(0, 12, by = 3)) + 
  annotate("text", x = min(fig_FVC_clean$FVC_filter, na.rm = TRUE), 
           y = max_Resp * 1.15, 
           label = sprintf("R = %.2f, p = %.4f", sqrt(r_value), p_value), 
           hjust = 0, vjust = 0, size = 2.7, color = "black")

p
doc <- read_pptx()
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- ph_with(doc, dml(ggobj = p), location = ph_location(width = 3.15, height = 2.1))
print(doc, target = "plot.pptx")


