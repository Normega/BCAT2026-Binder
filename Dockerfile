FROM rocker/binder:4.4.2

USER root

# Install additional R packages not already in rocker/binder
# Using Posit Package Manager binaries for Ubuntu 24.04 (noble)
RUN Rscript -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/all/__linux__/noble/latest')); \
    install.packages(c( \
      'lmerTest', 'emmeans', 'broom.mixed', 'MuMIn', 'irr', \
      'BayesFactor', 'ppcor', 'car', \
      'brms', 'posterior', \
      'metafor', 'mediation', \
      'patchwork', 'ggeffects', \
      'flextable', 'officer' \
    ), Ncpus = parallel::detectCores())"

USER ${NB_USER}
