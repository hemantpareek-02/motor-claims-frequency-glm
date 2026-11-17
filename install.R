# Installs the packages needed to run analysis.Rmd
pkgs <- c("insuranceData", "dplyr", "tidyr", "ggplot2", "scales", "rmarkdown", "knitr")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
# MASS ships with R (recommended package), so no install is needed.
