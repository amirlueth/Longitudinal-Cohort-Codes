############################################################
## Author: Amir J. Lueth
## Purpose: Published analysis code sample
## Topic: OPE exposures and longitudinal blood pressure
##
## NOTE:
## This script is adapted from an original published analysis.
## Cohort-specific file paths, identifiers, record-level
## decisions, and confidential data details have been removed.
############################################################

library(tidyverse)
library(survey)
library(broom)

## ---------------------------------------------------------
## Expected analysis variables
## ---------------------------------------------------------
## participant_id        unique participant identifier
## gest_age              gestational age at BP measurement
## delivery_ga           gestational age at delivery
## sbp                   systolic blood pressure
## dbp                   diastolic blood pressure
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
# bp_data <- readRDS("data/approved_analysis_data.rds")


## ---------------------------------------------------------
## 1. Blood pressure quality control
## ---------------------------------------------------------

bp_qc <- bp_data %>%
  mutate(
    ga_week = floor(as.numeric(gest_age)),
    delivery_week = floor(as.numeric(delivery_ga)),
    ga_difference = if_else(
      ga_week == delivery_week,
      as.numeric(gest_age),
      as.numeric(delivery_week - ga_week)
    )
  ) %>%
  # Remove measurements recorded essentially at delivery.
  filter(!ga_difference %in% c(-1, -2, -3, -4, -5, -6, -7, 1))


## ---------------------------------------------------------
## 2. Gestational-age-specific BP outlier screening
## ---------------------------------------------------------

bp_reference <- bp_qc %>%
  group_by(ga_week) %>%
  summarise(
    mean_sbp = mean(sbp, na.rm = TRUE),
    sd_sbp   = sd(sbp, na.rm = TRUE),
    mean_dbp = mean(dbp, na.rm = TRUE),
    sd_dbp   = sd(dbp, na.rm = TRUE),
    .groups = "drop"
  )

bp_qc <- bp_qc %>%
  left_join(bp_reference, by = "ga_week") %>%
  mutate(
    sbp_z = (sbp - mean_sbp) / sd_sbp,
    dbp_z = (dbp - mean_dbp) / sd_dbp,
    sbp_outlier = abs(sbp_z) > 4,
    dbp_outlier = abs(dbp_z) > 4
  )

# In the original workflow, flagged trajectories were reviewed before
# exclusions were finalized. The code sample applies the pre-specified
# +/- 4 SD rule after that review step.
bp_analysis <- bp_qc %>%
  filter(
    (is.na(sbp_outlier) | !sbp_outlier) &
    (is.na(dbp_outlier) | !dbp_outlier)
  )


## ---------------------------------------------------------
## 3. Repeated-measure weighting
## ---------------------------------------------------------

visit_count <- bp_analysis %>%
  filter(!is.na(sbp) | !is.na(dbp)) %>%
  count(participant_id, name = "n_bp_measurements")

bp_analysis <- bp_analysis %>%
  left_join(visit_count, by = "participant_id") %>%
  mutate(
    bp_weight = sampling_weight * (1 / n_bp_measurements)
  ) %>%
  filter(!is.na(bp_weight))

bp_design <- svydesign(
  id = ~participant_id,
  weights = ~bp_weight,
  data = bp_analysis
)


## ---------------------------------------------------------
## 4. Adjusted weighted longitudinal BP models
## ---------------------------------------------------------
## Gestational age is modeled using linear and quadratic terms.
## Models are adjusted for maternal age, prepregnancy BMI,
## education, race/ethnicity, insurance, and parity.

continuous_exposures <- paste0("ope_", 1:5, "_iqr")
binary_exposures <- paste0("ope_", 6:8, "_detect")
exposures <- c(continuous_exposures, binary_exposures)

run_bp_model <- function(outcome, exposure, design) {

  model_formula <- as.formula(
    paste0(
      outcome, " ~ ", exposure,
      " + gest_age + I(gest_age^2)",
      " + age + bmi + education + race_ethnicity + insurance + parity"
    )
  )

  fit <- svyglm(
    model_formula,
    design = design,
    na.action = na.omit
  )

  tidy(fit, conf.int = TRUE) %>%
    filter(term == exposure) %>%
    mutate(
      outcome = outcome,
      exposure = exposure
    )
}

sbp_results <- map_dfr(
  exposures,
  ~run_bp_model("sbp", .x, bp_design)
)

dbp_results <- map_dfr(
  exposures,
  ~run_bp_model("dbp", .x, bp_design)
)

bp_results <- bind_rows(sbp_results, dbp_results) %>%
  select(outcome, exposure, estimate, conf.low, conf.high, p.value)

print(bp_results)

# write_csv(bp_results, "output/longitudinal_bp_results.csv")


## ---------------------------------------------------------
## 5. Example sensitivity analysis
## ---------------------------------------------------------
## In the original analysis, visit-specific exposure models were
## evaluated to address temporality. The same modeling structure can
## be applied after replacing pregnancy-average exposure terms with
## visit-specific exposure terms.

# visit1_design <- update(bp_design, ...)
# visit1_results <- map_dfr(
#   visit1_exposures,
#   ~run_bp_model("sbp", .x, visit1_design)
# )
