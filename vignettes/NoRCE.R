## ----style, echo=FALSE, results="asis", message=FALSE, warnings = FALSE----
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE,fig.width=6, fig.height=6 )

## ----set size------------------------------------------------------------
knitr::opts_knit$set(width = 80)

## ----Install, eval=FALSE, echo=TRUE, include=TRUE------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("NoRCE")

## ----Load, eval=TRUE, echo=TRUE, include=TRUE----------------------------
library(NoRCE)

## ----Load package, eval=TRUE, echo=TRUE, include=TRUE--------------------
library(NoRCE)

