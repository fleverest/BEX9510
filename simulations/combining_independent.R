#!/usr/bin/env Rscript

# This file replicates the first simulation study in Section 7 of the paper
# "E-values:Calibration, combination, and applications".

library(reshape2)

#### Calibrators
source("simulations/lib/calibrate.R")

#### merging functions for e-values and p-values
source("simulations/lib/merge.R")


#### e-values

# Calucates e-values using the product of likelihood ratios for the simple
# hypotheses H0 and H1.
e_product <- function(x, δ) cumprod(dnorm(x, δ, 1) / dnorm(x, 0, 1))

# Calcualtes e-values using the average of likelihood rations for the simple
# hypotheses H0 and H1.
e_average <- function(x, δ) {
  cumsum(dnorm(x, δ, 1) / dnorm(x, 0 ,1)) / seq_along(x)
}

# Calculates the e-values using the universal test martingale.
e_universal <- function(x) {
  sapply(seq_along(x), \(k) {
    1 / sqrt(k + 1) * exp(1 / (2 * k + 2) * sum(x[1:k]) ^ 2)
  })
}


#### p-values

# Calculates p-values by applying the std normal cdf to obtain a p-value for
# each observation, and then combining them using fishers method
p_fisher <- function(x) {
  sapply(seq_along(x), \(k) {
    fisher_merge_p(pnorm(x[1:k]))
  })
}

# Calculates p-values by applying the std normal cdf to obtain a p-value per
# observation, then combines them using simes' method
p_simes <- function(x) {
  sapply(seq_along(x), \(k) {
    simes_merge_p(pnorm(x[1:k]))
  })
}

# Calculate p-values by applying the std-normal cdf to obtain a p-value per
# observation, then combines them by applying the bonferroni correction.
p_bonferroni <- function(x) {
  p_values <- pnorm(x)
  sapply(seq_along(x), \(k) {
    bonferroni_merge_p(p_values[1:k])
  })
}


#### Conduct simulations based on e-values

# In this simulation study, we generate observations from the Gaussian model
# rnorm(μ, 1) with unknown μ. The null hypothesis is H0: μ=0, and the
# alternative is H1: μ=δ.

# "Global" simulation parameters
seed <- 123456789
μ <- 0
δ <- -0.1
n_observations <- 10000
n_seeds <- 100

set.seed(seed)
product <- melt(
  replicate(
    n_seeds,
    e_product(rnorm(n_observations, δ, 1), δ)
  )
)
names(product) <- c("observations", "simulation_idx", "e_value")
product[, "method"] <- "product"
product[, "p_value"] <- calibrate_e_to_p(product[, "e_value"])

set.seed(seed)
universal <- melt(
  replicate(
    n_seeds,
    e_universal(rnorm(n_observations, δ, 1))
  )
)
names(universal) <- c("observations", "simulation_idx", "e_value")
universal[, "method"] <- "universal"
universal[, "p_value"] <- calibrate_e_to_p(universal[, "e_value"])

#### Conduct simulations based on p-values.

# The two plots in figure 1 use the same underlying p-values, calculated by
# combining per-observation p-values using Fisher's method
set.seed(seed)
p_values <- melt(
  replicate(
    n_seeds,
    p_fisher(rnorm(n_observations, δ, 1))
  )
)

# The line labelled "Fisher" employs the approximation of an e-value
# given by the reciprocal of the p-values.
fisher <- p_values
names(fisher) <- c("observations", "simulation_idx", "p_value")
fisher[, "method"] <- "fisher"
fisher[, "e_value"] <- reciprocal(fisher[, "p_value"])

# Each value combines the p-values to obtain either an e-value (or a
# family-wise bound in the case of fisher_VS)
fisher_vs <- p_values
names(fisher_vs) <- c("observations", "simulation_idx", "p_value")
fisher_vs[, "method"] <- "fisher_vs"
fisher_vs[, "e_value"] <- vs_p(p_values[, "value"])

data <- do.call(
  rbind,
  list(
    product,
    universal,
    fisher,
    fisher_vs
  )
)

write.csv2(data, "results/data/combining_independent_p_and_e_fig1.csv", row.names = FALSE)

# "Global" simulation parameters
seed <- 987654321
μ <- 0
δ <- -0.1
n_observations <- 1000
n_seeds <- 1000

# Using averaging
set.seed(seed)
average <- melt(
  replicate(
    n_seeds,
    e_average(rnorm(n_observations, δ, 1), δ)
  )
)
names(average) <- c("observations", "simulation_idx", "e_value")
average[, "method"] <- "average"
average[, "p_value"] <- calibrate_e_to_p(average[, "e_value"])

# Simes method
set.seed(seed)
simes <- melt(
  replicate(
    n_seeds,
    p_simes(rnorm(n_observations, δ, 1))
  )
)
names(simes) <- c("observations", "simulation_idx", "p_value")
simes[, "method"] <- "simes"
simes[, "e_value"] <- reciprocal(simes[, "p_value"])

simes_vs <- simes
simes_vs[, "method"] <- "simes_vs"
simes_vs[, "e_value"] <- vs_p(simes_vs[, "p_value"])

set.seed(seed)
bonferroni <- melt(
  replicate(
    n_seeds,
    p_bonferroni(rnorm(n_observations, δ, 1))
  )
)
names(bonferroni) <- c("observations", "simulation_idx", "p_value")
bonferroni[, "method"] <- "bonferroni"
bonferroni[, "e_value"] <- reciprocal(bonferroni[, "p_value"])

bonferroni_vs <- bonferroni
bonferroni_vs[, "method"] <- "bonferroni_vs"
bonferroni_vs[, "e_value"] <- vs_p(bonferroni_vs[, "p_value"])

data <- do.call(
  rbind,
  list(
    average,
    simes,
    simes_vs,
    bonferroni,
    bonferroni_vs
  )
)

write.csv2(
  data,
  "results/data/combining_independent_p_and_e_fig2.csv",
  row.names = FALSE
)