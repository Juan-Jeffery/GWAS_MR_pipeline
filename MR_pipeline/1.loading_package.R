# ─────────────── load package ───────────────
load_packages <- function() {
  packages <- c(
    "TwoSampleMR",
    "gwasglue",
    "gwasvcf",
    "VariantAnnotation",
    "MendelianRandomization",
    "dplyr",
    "ggplot2",
    "CMplot",
    "ieugwasr",
    "data.table",
    "R.utils"
  )

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Package '%s' not found. Installing...", pkg))
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }

}
