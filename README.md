# Salivary Biomarkers for Classifying Acute Physical Fatigue

Reproducible analysis code and data accompanying:

> Lindsey B, Bowden K, Shaul Y, Petricoin E, Caswell SV, Alhammad R, Elayadi AN,
> Roberts BM, Martin JR. **Salivary Biomarkers for Classifying Acute Physical
> Fatigue: A Comparative, Exploratory Study.** *bioRxiv* preprint (2025).
> https://doi.org/10.1101/2025.06.04.657971. Under revision at *Translational
> Journal of the American College of Sports Medicine*.

See [CITATION.cff](CITATION.cff) for full citation details.

## Overview

This repository reproduces every table, figure, and statistic reported in the
manuscript: comparison of targeted salivary stress biomarkers, salivary
proteins identified via untargeted proteomics, and a combined multi-platform
panel for classifying acute physical fatigue, using nested leave-one-subject-
out cross-validation (LOSOCV).

## Repository Structure

```
├── Data/
│   ├── Biomarker_Metadata_Saliva.csv            targeted biomarker + performance data
│   └── Fatigue_Study_10_Samples(Proteins).csv   untargeted proteomics data
├── Scripts/
│   ├── PROSPER_Fatigue_Biomarker_Analysis.R     main analysis (Tables 1-2, Figures 2-3, all panels)
│   └── PROSPER_Imputation_Sensitivity.R         imputation sensitivity analysis (companion, run independently)
├── Results/                                     outputs from running the scripts above (populated on run)
├── PROSPER_Fatigue_Biomarker_Analysis.Rproj     RStudio project file — open this first
├── LICENSE                                      code license (MIT)
├── DATA_LICENSE                                 data license (CC BY 4.0)
└── CITATION.cff
```

`Biomarker_Metadata_Saliva.csv` is a filtered extract of the full study
database, containing only the participant ID, timepoint, targeted salivary
biomarker, and physical performance columns used in this analysis. The full
database also includes item-level responses to several psychological/sleep
survey instruments and pilot/protocol-development data not part of the
reported study sample; those are not included here.

## How to Run

1. Open `PROSPER_Fatigue_Biomarker_Analysis.Rproj` in RStudio (double-click it,
   or File → Open Project). This sets the working directory correctly for the
   `here`-based relative paths used throughout — running the scripts without
   opening this project first will cause `here::i_am()` path-resolution
   errors.
2. Install the required packages (see **Reproducibility** below).
3. Run `Scripts/PROSPER_Fatigue_Biomarker_Analysis.R` to reproduce Table 1,
   Table 2, Figures 2-3, permutation importance, and the SVM-vs-logistic
   regression comparison. Expected run time: ~30-90 minutes, dominated by the
   protein and combined-panel nested screens (~1,350+ mixed-effects models
   fit per fold).
4. Run `Scripts/PROSPER_Imputation_Sensitivity.R` independently to reproduce
   the imputation sensitivity analysis. Expected run time: ~15-20 minutes.
5. Outputs are written to `Results/` (CSVs are version-controlled in this
   repo for review; TIFF figures are `.gitignore`d since they're large and
   fully reproducible from the CSVs/code — see `.gitignore`).

## Reproducibility

Analyses were run under:

- **R version 4.5.2**

Package versions confirmed directly from console output at analysis time:

| Package | Version |
|---|---|
| glmnet | 4.1-10 |
| dplyr | 1.1.4 |
| readr | 2.1.6 |
| forcats | 1.0.1 |
| stringr | 1.6.0 |
| ggplot2 | 4.0.3 |
| tibble | 3.3.1 |
| lubridate | 1.9.4 |
| tidyr | 1.3.2 |
| purrr | 1.2.1 |

Additional packages used (`here`, `e1071`, `pROC`, `emmeans`, `broom`,
`future`, `furrr`, `lme4`, `lmerTest`, `readxl`, `ggsignif`, `ggnewscale`)
were loaded successfully at analysis time but their exact version numbers
were not captured in the retained console output.

**For exact reproducibility**, run the following once inside this project and
commit the resulting `renv.lock` file:

```r
install.packages("renv")
renv::init()      # scans installed packages and creates renv.lock
renv::snapshot()  # pins exact versions
```

Anyone re-running this analysis can then restore the exact package versions
used with `renv::restore()`. Alternatively, run `sessionInfo()` after
sourcing the main script for a complete, exact version list.

Given the small sample size (n=10) and the iterative nature of some model
fits (e.g., `lme4` mixed-effects models, which can occasionally report
`boundary (singular) fit` warnings — expected and non-fatal at this sample
size), minor numeric differences at the level of floating-point precision may
occur across R/package versions or operating systems. These do not affect
any value as reported (rounded) in the manuscript.

## Data Availability Statement

The data and analysis code that support the findings of this study are
openly available. The complete analysis pipeline, targeted biomarker and
proteomic datasets, and all reported results are publicly available at
https://github.com/blindse3/salivary-biomarkers-fatigue.

## License

- **Code** (`Scripts/`): [MIT License](LICENSE)
- **Data** (`Data/`, `Results/`): [CC BY 4.0](DATA_LICENSE)

## Citation

See [CITATION.cff](CITATION.cff), or use GitHub's "Cite this repository"
button once this repo is published.

## Contact

Corresponding author: Dr. Joel Martin, George Mason University —
jmarti38@gmu.edu
