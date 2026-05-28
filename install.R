install.packages(c(
  "tidyverse", "readxl",
  "lme4", "lmerTest",
  "emmeans", "broom.mixed", "MuMIn", "irr",
  "BayesFactor", "ppcor", "car",
  "brms", "posterior",
  "metafor", "mediation",
  "patchwork", "ggeffects", "scales",
  "flextable", "officer"
), Ncpus = parallel::detectCores())
