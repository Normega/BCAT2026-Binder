FROM rocker/binder:4.4.2

USER root

# Install additional R packages not already in rocker/binder
# pkgType=binary avoids source compilation
RUN Rscript -e " \
  options(repos = c(CRAN = 'https://packagemanager.posit.co/all/__linux__/noble/latest')); \
  install.packages(c( \
    'lmerTest', 'emmeans', 'broom.mixed', 'MuMIn', 'irr', \
    'BayesFactor', 'ppcor', 'car', \
    'brms', 'posterior', \
    'metafor', 'mediation', \
    'patchwork', 'ggeffects', \
    'flextable', 'officer' \
  ), Ncpus = parallel::detectCores())"

# Copy the repository into the image
COPY --chown=rstudio:rstudio . /home/rstudio/

# Run postBuild as rstudio so Results/ is owned by the right user
USER rstudio
WORKDIR /home/rstudio
RUN bash postBuild
