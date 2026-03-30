# Roche ADS Programmer Coding Assessment

## Overview
This repository contains my solutions to the Roche Analytical Data Science Programmer coding assessment. The assessment covers SDTM dataset creation, ADaM dataset creation, and clinical reporting using open-source R packages from the Pharmaverse ecosystem.

## Repository Structure
```
├── question_1_sdtm/
│   ├── question_1.R          # SDTM DS domain creation script
│   ├── ds_domain.csv         # Output DS dataset
│   └── ds_domain_log.txt     # Evidence of error-free execution
│
├── question_2_adam/
│   ├── question_2.R          # ADaM ADSL creation script
│   ├── adsl.csv              # Output ADSL dataset
│   └── adsl_log.txt          # Evidence of error-free execution
│
├── question_3_tlg/
│   ├── question_3_ae_summary.R       # AE summary table script
│   ├── question_3_visualizations.R   # AE visualizations script
│   ├── ae_summary_table.html         # TEAE summary table output
│   ├── ae_severity_plot.png          # AE severity distribution plot
│   ├── ae_top10_plot.png             # Top 10 AEs with 95% CI plot
│   └── ae_tlg_log.txt                # Evidence of error-free execution
```

## Questions Summary

### Question 1: SDTM DS Domain Creation
- **Package**: `{sdtm.oak}`, `{dplyr}`
- **Input**: `pharmaverseraw::ds_raw`
- **Output**: DS domain with 12 required variables (STUDYID, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD, DSCAT, VISITNUM, VISIT, DSDTC, DSSTDTC, DSSTDY)
- **Key decisions**: Controlled terminology mapping using study_ct with extended manual mapping for terms not directly covered; DSSTDY derived by joining DM for RFSTDTC; USUBJID aligned to DM format

### Question 2: ADaM ADSL Dataset Creation
- **Package**: `{admiral}`, `{dplyr}`, `{tidyr}`, `{lubridate}`
- **Input**: `pharmaversesdtm::dm`, `vs`, `ex`, `ds`, `ae`
- **Output**: ADSL with 26 variables including custom derivations
- **Key derivations**:
  - AGEGR9/AGEGR9N: Age grouping (<18, 18-50, >50)
  - TRTSDTM/TRTSTMF: Treatment start datetime with time imputation from EX
  - ITTFL: Randomization flag (Y/N) based on ARM population
  - LSTAVLDT: Last known alive date from VS, AE, DS, and EX sources

### Question 3: TLG - Adverse Events Reporting
- **Packages**: `{gtsummary}`, `{gt}`, `{ggplot2}`, `{dplyr}`
- **Input**: `pharmaverseadam::adae`, `pharmaverseadam::adsl`
- **Outputs**:
  - TEAE summary table by SOC and Preferred Term with treatment columns and total
  - AE severity distribution stacked bar chart by treatment arm
  - Top 10 most frequent AEs forest plot with 95% Clopper-Pearson confidence intervals

## How to Run

### Prerequisites
```r
install.packages(c(
  "admiral", "sdtm.oak", "pharmaverseraw", "pharmaversesdtm",
  "pharmaverseadam", "gt", "ggplot2", "gtsummary",
  "dplyr", "tidyr", "lubridate"
))
```

### Execution
Run each script independently in order:
1. `question_1_sdtm/question_1.R`
2. `question_2_adam/question_2.R`
3. `question_3_tlg/question_3_ae_summary.R`
4. `question_3_tlg/question_3_visualizations.R`

## Environment
- R version: 4.5
- Platform: Posit Cloud
- AI assistance: Claude (Anthropic) used as coding assistant per assessment guidelines
