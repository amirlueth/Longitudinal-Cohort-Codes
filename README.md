# Published Epidemiologic Analysis Code Sample

**Author:** Amir J. Lueth, PhD, MPH

## Overview
This repository is a code sample adapted from analyses I conducted for published epidemiologic research on organophosphate ester (OPE) exposures, longitudinal blood pressure, and preeclampsia.

The purpose of this repository is to demonstrate my analytic workflow and coding style while protecting cohort-specific and institutional information. The original analyses were conducted in the LIFECODES Pregnancy Cohort. This public-facing version removes participant identifiers, internal file paths, record-level decisions, collaborator comments, and other cohort-specific implementation details.

**No participant-level data are included.**

## What this sample demonstrates
- longitudinal blood pressure data preparation and quality control;
- construction and application of sampling weights;
- repeated-measures modeling;
- time-to-event outcome construction;
- weighted Cox proportional hazards regression;
- exposure scaling using interquartile ranges;
- sensitivity-analysis structure;
- reproducible extraction and organization of model results.

## Repository structure

```text
R/
  01_longitudinal_bp_analysis.R
  02_preeclampsia_survival_analysis.R
data/
  README.md
output/
  README.md
```

## Data
The original research data cannot be distributed. The scripts therefore assume an analysis-ready data object with generalized variable names. Variable names are documented in the scripts so the analytic logic can be reviewed without access to the underlying cohort.

## Notes
This repository is intended as a professional code sample rather than a fully executable replication package. The analytic logic is adapted from my original research code, but confidential and cohort-specific details have been removed or generalized.

## Methods represented
The original workflow included:
- inverse-probability/sampling weights;
- gestational-age-specific blood pressure cleaning and outlier review;
- weighted repeated-measures analysis of systolic and diastolic blood pressure;
- pregnancy-averaged and visit-specific exposure analyses;
- weighted Cox proportional hazards models for preeclampsia;
- sensitivity analyses.

## Contact
Amir J. Lueth, PhD, MPH
