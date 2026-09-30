############################################################
## Author: Amir J. Lueth
## Purpose: Published analysis code sample
## Topic: OPE exposures and preeclampsia
##
## NOTE:
## This script is adapted from an original published analysis.
## Cohort-specific file paths, identifiers, record-level
## decisions, and confidential data details have been removed.
############################################################

library(tidyverse)
library(survival)
library(broom)

## ---------------------------------------------------------
## Expected analysis variables
## ---------------------------------------------------------
## participant_id        unique participant identifier
## pe_status             1 = preeclampsia, 0 = no preeclampsia
## pe_ga                 gestational age at preeclampsia diagnosis
## delivery_ga           gestational age at delivery
## sampling_weight       study sampling weight
## age                   maternal age
## bmi                   prepregnancy BMI
## education             maternal education
## race_ethnicity        race/ethnicity
## insurance             insurance category
## parity                parity
##
## ope_1_iqr ... ope_5_iqr
##     continuous OPE biomarkers scaled to one IQR
##
## ope_6_detect ... ope_8_detect
##     OPE biomarkers modeled as detected/not detected
##
## Replace the object below with an approved analysis dataset.
# pe_data <- readRDS("data/approved_analysis_data.rds")


## ---------------------------------------------------------
## 1. Time-to-event outcome construction
## ---------------------------------------------------------

pe_analysis <- pe_data %>%
  mutate(
    event_time = if_else(
      pe_status == 1 & !is.na(pe_ga),
      as.numeric(pe_ga),
      as.numeric(delivery_ga)
    ),
    event = as.integer(pe_status == 1)
  )


## ---------------------------------------------------------
## 2. IQR-standardized exposure variables
## ---------------------------------------------------------
## The original analysis used weighted IQRs for continuous OPE
## biomarkers so estimates corresponded to an IQR increase.

continuous_exposures_raw <- paste0("ope_", 1:5)

pe_design <- svydesign(
  id = ~0,
  weights = ~sampling_weight,
  data = pe_analysis
)

for (x in continuous_exposures_raw) {

  q <- svyquantile(
    as.formula(paste0("~", x)),
    design = pe_design,
    quantiles = c(0.25, 0.75),
    ci = FALSE,
    se = FALSE,
    na.rm = TRUE
  )

  q_values <- as.numeric(q[[1]])
  exposure_iqr <- q_values[2] - q_values[1]

  pe_analysis[[paste0(x, "_iqr")]] <-
    as.numeric(pe_analysis[[x]]) / exposure_iqr
}


## ---------------------------------------------------------
## 3. Adjusted weighted Cox proportional hazards models
## ---------------------------------------------------------

continuous_exposures <- paste0("ope_", 1:5, "_iqr")
binary_exposures <- paste0("ope_", 6:8, "_detect")
exposures <- c(continuous_exposures, binary_exposures)

run_cox_model <- function(exposure, data) {

  model_formula <- as.formula(
    paste0(
      "Surv(event_time, event) ~ ", exposure,
      " + age + bmi + education + race_ethnicity + insurance + parity"
    )
  )

  fit <- coxph(
    model_formula,
    weights = sampling_weight,
    robust = TRUE,
    data = data
  )

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure) %>%
    mutate(exposure = exposure)
}

pe_results <- map_dfr(
  exposures,
  ~run_cox_model(.x, pe_analysis)
) %>%
  select(exposure, estimate, conf.low, conf.high, p.value)

print(pe_results)

# write_csv(pe_results, "output/preeclampsia_cox_results.csv")


## ---------------------------------------------------------
## 4. Sensitivity-analysis structure
## ---------------------------------------------------------
## Examples evaluated in the original workflow included:
## - visit-specific exposure models;
## - analyses excluding smokers;
## - adjustment for prior preeclampsia history;
## - analyses of preeclampsia phenotype/timing.

# Example:
# pe_no_smoking <- pe_analysis %>%
#   filter(smoking != "Yes")
#
# sensitivity_results <- map_dfr(
#   exposures,
#   ~run_cox_model(.x, pe_no_smoking)
# )
