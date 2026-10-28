# irtDFI

Dynamic fit index cutoffs for polytomous item response models fitted with
`mirt` (graded, GPCM, Rasch/PCM). Derives model-specific cutoffs for
RMSEA, SRMSR, CFI and TLI (C2-based) by simulating from the fitted model
(Level 0) and from versions with increasing local dependence (Levels 1-L).

```r
# install.packages("remotes")
remotes::install_github("USERNAME/irtDFI")
library(irtDFI); library(mirt); data(sim_grm)   # computer-generated data
m <- mirt(sim_grm[, 1:10], 1, itemtype = "graded")
r <- dfi_irt(m, nrep = 200, L = 3, cores = 1)
r; dfi_plot(r, "RMSEA")
```

Status: version 0.1.0, under development. `sim_grm` is computer-generated;
no human respondents. Simulation scripts are in `inst/scripts/`.
