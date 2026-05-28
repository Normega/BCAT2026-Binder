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

# Download pre-computed brms model objects (baked into the image layer)
RUN cd /home/rstudio && bash postBuild

USER rstudio
WORKDIR /home/rstudio
