# What You Miss Won't Move You
## Awareness Connects Respiratory Change to Subjective Arousal

This repository contains analysis code for a five-study empirical paper
investigating whether conscious detection of breathing changes is required
for those changes to influence subjective arousal.

**Preregistration (Study 5):** https://osf.io/r6zja/overview?view_only=4edce0bc913947d3a8491bfbdeb0deb3

**OSF archive (data, RDS objects, trial-by-trial QC files):** https://osf.io/g7rdb/overview?view_only=00d7f593dad94c3b9c4c34b994b8c162

**Task code (Studies 1–5):** https://anonymous.4open.science/r/BCAT2026-Tasks-9B85/

---

## Replication instructions

### Step 1 — Download the code and data

**1a. Download this repository**

Click the green **Code** button at the top of this GitHub page and choose
**Download ZIP**. Unzip it somewhere convenient (e.g. your Documents folder).
You should now have a folder — call it your *root folder* — that contains an
`Analysis/` subfolder.

**1b. Download the Data folder from OSF**

1. Go to: https://osf.io/g7rdb/files/osfstorage?view_only=00d7f593dad94c3b9c4c34b994b8c162
2. Find the **`Data`** folder in the file browser.
3. Click the **three dots (⋯)** to the right of the `Data` folder and choose
   **Download as zip**.
4. Once downloaded, unzip it. You should get a folder called `Data` containing
   CSV and RDS files.
5. Move or copy the `Data` folder into your root folder (the same folder that
   contains `Analysis/`).

When set up correctly your folder should look like this:

```
BCAT2026/          ← this is your root folder (it can be named anything)
├── Analysis/
└── Data/
```

### Step 2 — Tell the script where to find your files

Open `Analysis/MainAnalysis.R` in RStudio (or any text editor). Near the top
you will find this line:

```r
BASE_DIR <- "."
```

Replace `"."` with the full path to your root folder. For example:

- **Windows:** `BASE_DIR <- "C:/Users/YourName/Documents/BCAT2026"`
- **Mac/Linux:** `BASE_DIR <- "/Users/YourName/Documents/BCAT2026"`

**How to find your path:**
- *Windows:* Open the root folder in File Explorer, click in the address bar,
  and copy the path shown. Replace any backslashes (`\`) with forward slashes (`/`).
- *Mac:* Right-click the root folder in Finder, hold the Option key, and choose
  "Copy … as Pathname".
- *RStudio shortcut:* Open any file inside the root folder in RStudio, then run
  `dirname(dirname(rstudioapi::getActiveDocumentContext()$path))` in the Console.

**Confirm everything is in the right place:**

```r
file.exists(file.path(BASE_DIR, "Analysis", "MainAnalysis.R"))
file.exists(file.path(BASE_DIR, "Data"))
```

Both should print `TRUE`.

### Step 3 — Run the analysis

```r
source("Analysis/MainAnalysis.R")
```

Or open `Analysis/MainAnalysis.R` in RStudio and click **Source**.

`MainAnalysis.R` runs all five studies in order, writes result CSVs to
`Results/`, generates all manuscript figures to `Figures/`, and builds
the formatted table documents to `Tables/`. Total runtime is approximately
45–90 minutes depending on hardware; Bayesian models in Studies 4 and 5
are the bottleneck.

---

## Repository structure

```
/
├── Analysis/
│   ├── MainAnalysis.R               Entry point — source this to run everything
│   ├── utils.R                      Shared helper functions and d' formula
│   ├── theme_bcat.R                 ggplot2 theme and colour palette
│   ├── meta_analysis.R              Random-effects meta-analysis across studies
│   │
│   ├── analysis_arousal.R           H4A–H4C: arousal gating tests (hits vs misses)
│   ├── analysis_belt.R              Study 5 belt compliance, salience, regime comparison
│   ├── analysis_hbd.R               Heartbeat detection: cross-modal dissociation (S6)
│   ├── analysis_individual_differences.R  Exploratory individual-differences (S5)
│   ├── analysis_maia.R              H3A–H3B: MAIA dissociation; MAIA × gating moderation
│   ├── analysis_miss_baseline.R     Bayesian null tests: missed trials vs no-change baseline
│   ├── analysis_s4_entrainment.R    Belt entrainment and Breath-over-Visual advantage (S4)
│   ├── analysis_s7_maia_selfesteem.R  MAIA specificity controls: self-esteem, trait self-doubt
│   ├── analysis_study3_attraction_mediation.R  Study 3 misattribution mediation
│   ├── analysis_study5_exploratory.R  Study 5 exploratory analyses
│   ├── analysis_tce.R               Sensitivity analyses (regime, matched magnitude, prior);
│   │                                Study 5 completer check
│   ├── analysis_val_detection.R     Change² dose-response; detection accuracy models
│   ├── analysis_val_pilot_studies.R Studies 1A/1B and 2 procedure comparison
│   ├── analysis_val_thresholds.R    Threshold descriptives, test-block H5 validity, ICC
│   ├── belt_salience_followup.R     Belt-salience independence and direction compliance
│   ├── test_block_accuracy.R        Test block: d' (3AFC), detection accuracy,
│   │                                direction asymmetry by study/group/session
│   ├── test_block_arousal.R         Test block: Bayesian mediation of salience → arousal
│   │
│   ├── fig_accuracy.R               Figure 1: psychometric functions
│   ├── fig_arousal.R                Figure 2: arousal gating across studies
│   ├── fig_regime_comparison.R      Figure S3: TCE hits/misses regime sensitivity
│   ├── fig_staircase.R              Figure S1: staircase convergence
│   │
│   ├── Build_Main_Tables.R          Main manuscript tables
│   ├── Build_Reliability_Tables.R   Reliability tables
│   ├── Build_Supplementary_Tables.R Supplementary tables (ST1–ST8)
│   │
│   ├── DataCleaning/                Study-level cleaning scripts
│   │   ├── run_all_cleaning.R
│   │   ├── study1_clean.R
│   │   ├── study2_clean.R
│   │   ├── study3_clean.R
│   │   ├── study4_clean.R
│   │   └── study5_clean.R
│   │
│   ├── ScaleReliability/            Scale reliability preparation
│   │   ├── Study1_PrepScales.R
│   │   ├── Study2_PrepScales.R
│   │   ├── Study3_PrepScales.R
│   │   ├── Study4_PrepScales.R
│   │   └── Study5_PrepScales.R
│   │
│   └── Study5/                      Study 5 physiological processing pipeline
│       ├── study5_processing_README.md   ← read before running physio scripts
│       ├── breath_pipeline.R
│       ├── analysis_study5.R
│       └── Intero2025_*.R           Individual processing steps
│
├── Data/                            Not in git — download from OSF archive
├── Results/                         Not in git — generated by MainAnalysis.R
├── Figures/                         Not in git — generated by MainAnalysis.R
└── Tables/                          Not in git — generated by MainAnalysis.R
```

---

## Requirements

- R ≥ 4.2.0
- Key packages: `brms`, `lme4`, `lmerTest`, `broom.mixed`, `mediation`,
  `metafor`, `tidyverse`, `purrr`, `patchwork`, `flextable`, `officer`,
  `BayesFactor`, `MuMIn`, `signal`, `psycho`

Install all dependencies:

```r
install.packages(c(
  "brms", "lme4", "lmerTest", "broom.mixed", "mediation", "metafor",
  "tidyverse", "purrr", "patchwork", "readr", "tibble",
  "flextable", "officer", "BayesFactor", "MuMIn", "signal", "psycho"
))
```

---

## Data

Processed summary CSVs (one row per participant per study) and RDS data objects
are archived on OSF and are read directly by `MainAnalysis.R`. Raw PsychoPy
output, Qualtrics exports, and physiological recordings are also on OSF.

Participant IDs have been anonymized. No identifying information is present
in the data files.

---

## Deviations from preregistration

Deviations from the Study 5 preregistration are documented in the
Supplementary Materials, available on the OSF archive:
https://osf.io/g7rdb/overview?view_only=00d7f593dad94c3b9c4c34b994b8c162

---

## License

Code: MIT License  
Data: CC-BY 4.0
