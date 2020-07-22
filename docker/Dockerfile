# Note: this must be built in the GillespieSSA2 root directory, not in the docker directory

FROM rocker/tidyverse:4.0.2

ARG GITHUB_PAT

ARG GITHUB_BRANCH
ENV GITHUB_BRANCH="${GITHUB_BRANCH:-master}"

RUN Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")' && \
  Rscript -e "remotes::install_github('rcannood/GillespieSSA2@$GITHUB_BRANCH', dependencies = TRUE, repos = 'https://cran.rstudio.com')"
  
ADD vignettes/*.Rmd /home/rstudio/vignettes/
