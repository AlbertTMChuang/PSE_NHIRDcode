# Sodium-Glucose Co-Transporter 2 Inhibitors and Risk of Post-Stroke Epilepsy
This repository contains SAS code supporting the study "Sodium-Glucose Co-Transporter 2 Inhibitors and Risk of Post-Stroke Epilepsy."
The study evaluated whether initiation of sodium-glucose co-transporter 2 (SGLT2) inhibitors was associated with the risk of post-stroke epilepsy among adults with type 2 diabetes and a recent cerebrovascular event, 
using dipeptidyl peptidase 4 (DPP4) inhibitors as the active comparator.

The program performs the following steps:
1. Identifies baseline medication variables beginning with `med_` and covariate variables beginning with `cov_` in `final.cohort`.
2. Applies propensity-score fine stratification using 50 strata and the average treatment effect estimand.
3. Estimates propensity scores using logistic regression.
4. Calculates unstabilized and stabilized inverse probability of treatment weights.
5. Combines the cohort, fine-stratification weights, and IPTW variables into `final.weight`.
6. Generates a baseline-characteristics table using the `%table1` macro.
7. Calculates weighted event counts, person-years, and incidence rates per 100 person-years.
8. Fits weighted Cox proportional hazards models for convulsions, epilepsy, and their composite outcome.
9. Combines outcome-specific incidence-rate and hazard-ratio results into summary datasets.
