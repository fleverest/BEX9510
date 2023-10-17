#!/usr/bin/env Rscript

# This file replicates the second simulation study in Section 7 of the paper
# "E-values:Calibration, combination, and applications".

library(reshape2)

source("simulations/lib/merge.R")
source("simulations/lib/calibrate.R")

# In this study, we conduct 20 hypothesis tests.
# Each hypothesis is identical. H0: μ=0, H1: μ=δ=-4.
# In the first 10 hypothesis tests, H1 is true.
# In the last 10 hypothesis tests, H0 is true.
# To test each hypothesis, we take a single observation from the underlying
# distribution (N(-4,1) for the first 10 and N(0, 1) for the last 10). Then, we
# use the p-value given by the std-Normal CDF.

δ <- -4

#### e-values 
e <- function(x) 0.5 * exp(x * δ - δ ^ 2 / 2) + 0.5
e2 <- function(x) exp(x * δ - δ ^ 2 / 2)

#### p-values
p <- function(x) pnorm(x, 0, 1)

# Adjust e and p values
source("simulations/lib/adjust.R")

#### Simulations for figure 3
n_seeds <- 1000
n_hypotheses <- 20
n_true <- 10
n_false <- 10
seed <- 543219876

set.seed(seed)
average <- melt(
  replicate(
    n_seeds,
    algorithm1(e(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1)))
  )
)
names(average) <- c("hypothesis", "simulation_idx", "e_value")
average[, "method"] <- "Average"
average[, "p_value"] <- calibrate_e_to_p(average[, "e_value"])
average[, "bound"] <- "Bona fide"

set.seed(seed)
average2 <- melt(
  replicate(
    n_seeds,
    algorithm1(e2(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1)))
  )
)
names(average2) <- c("hypothesis", "simulation_idx", "e_value")
average2[, "method"] <- "Average2"
average2[, "p_value"] <- calibrate_e_to_p(average2[, "e_value"])
average2[, "bound"] <- "Bona fide"


set.seed(seed)
product <- melt(
  replicate(
    n_seeds,
    algorithm2(e(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1)))
  )
)
names(product) <- c("hypothesis", "simulation_idx", "e_value")
product[, "method"] <- "Product"
product[, "p_value"] <- calibrate_e_to_p(product[, "e_value"])
product[, "bound"] <- "Bona fide"

set.seed(seed)
simes <- melt(
  replicate(
    n_seeds,
    simes_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    )
  )
)
names(simes) <- c("hypothesis", "simulation_idx", "p_value")
simes[, "method"] <- "Simes"
simes[, "e_value"] <- reciprocal(simes[, "p_value"])
simes[, "bound"] <- "Global UB"

simes_vs <- simes
simes_vs[, "method"] <- "Simes VS"
simes_vs[, "e_value"] <- vs_p(simes_vs[, "p_value"])
simes_vs[, "bound"] <- "VS UB"

set.seed(seed)
bonferroni <- melt(
  replicate(
    n_seeds,
    bonferroni_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    )
  )
)
names(bonferroni) <- c("hypothesis", "simulation_idx", "p_value")
bonferroni[, "method"] <- "Bonferroni"
bonferroni[, "e_value"] <- reciprocal(bonferroni[, "p_value"])
bonferroni[, "bound"] <- "Global UB"

bonferroni_vs <- bonferroni
bonferroni_vs[, "method"] <- "Bonferroni VS"
bonferroni_vs[, "e_value"] <- vs_p(bonferroni_vs[, "p_value"])
bonferroni_vs[, "bound"] <- "VS UB"

set.seed(seed)
fisher <- melt(
  replicate(
    n_seeds,
    fisher_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    )
  )
)
names(fisher) <- c("hypothesis", "simulation_idx", "p_value")
fisher[, "method"] <- "Fisher"
fisher[, "e_value"] <- reciprocal(fisher[, "p_value"])
fisher[, "bound"] <- "Global UB"

fisher_vs <- fisher
fisher_vs[, "method"] <- "Fisher VS"
fisher_vs[, "e_value"] <- vs_p(fisher_vs[, "p_value"])
fisher_vs[, "bound"] <- "VS UB"

data <- do.call(
  rbind,
  list(
    average,
    average2,
    product,
    simes,
    simes_vs,
    bonferroni,
    bonferroni_vs,
    fisher,
    fisher_vs
  )
)
write.csv2(
  data,
  "results/data/multiple_hypothesis_testing_fig3.csv",
  row.names = FALSE
)

#### Simulations for figure 4
n_seeds <- 1000
n_hypotheses <- 200
n_true <- 100
n_false <- 100
seed <- 242219226

set.seed(seed)
average <- melt(
  replicate(
    n_seeds,
    algorithm1(e(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1)))
  )
)
names(average) <- c("hypothesis", "simulation_idx", "e_value")
average[, "method"] <- "Average"
average[, "p_value"] <- calibrate_e_to_p(average[, "e_value"])
average[, "bound"] <- "Bona fide"

set.seed(seed)
average2 <- melt(
  replicate(
    n_seeds,
    algorithm1(e2(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1)))
  )
)
names(average2) <- c("hypothesis", "simulation_idx", "e_value")
average2[, "method"] <- "Average2"
average2[, "p_value"] <- calibrate_e_to_p(average2[, "e_value"])
average2[, "bound"] <- "Bona fide"


set.seed(seed)
product <- melt(
  replicate(
    n_seeds,
    algorithm2(e(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1)))
  )
)
names(product) <- c("hypothesis", "simulation_idx", "e_value")
product[, "method"] <- "Product"
product[, "p_value"] <- calibrate_e_to_p(product[, "e_value"])
product[, "bound"] <- "Bona fide"

set.seed(seed)
simes <- melt(
  replicate(
    n_seeds,
    simes_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    )
  )
)
names(simes) <- c("hypothesis", "simulation_idx", "p_value")
simes[, "method"] <- "Simes"
simes[, "e_value"] <- reciprocal(simes[, "p_value"])
simes[, "bound"] <- "Global UB"

simes_vs <- simes
simes_vs[, "method"] <- "Simes VS"
simes_vs[, "e_value"] <- vs_p(simes_vs[, "p_value"])
simes_vs[, "bound"] <- "VS UB"

set.seed(seed)
bonferroni <- melt(
  replicate(
    n_seeds,
    bonferroni_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    )
  )
)
names(bonferroni) <- c("hypothesis", "simulation_idx", "p_value")
bonferroni[, "method"] <- "Bonferroni"
bonferroni[, "e_value"] <- reciprocal(bonferroni[, "p_value"])
bonferroni[, "bound"] <- "Global UB"

bonferroni_vs <- bonferroni
bonferroni_vs[, "method"] <- "Bonferroni VS"
bonferroni_vs[, "e_value"] <- vs_p(bonferroni_vs[, "p_value"])
bonferroni_vs[, "bound"] <- "VS UB"

set.seed(seed)
fisher <- melt(
  replicate(
    n_seeds,
    fisher_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    )
  )
)
names(fisher) <- c("hypothesis", "simulation_idx", "p_value")
fisher[, "method"] <- "Fisher"
fisher[, "e_value"] <- reciprocal(fisher[, "p_value"])
fisher[, "bound"] <- "Global UB"

fisher_vs <- bonferroni
fisher_vs[, "method"] <- "Fisher VS"
fisher_vs[, "e_value"] <- vs_p(fisher_vs[, "p_value"])
fisher_vs[, "bound"] <- "VS UB"

data <- do.call(
  rbind,
  list(
    average,
    average2,
    product,
    simes,
    simes_vs,
    bonferroni,
    bonferroni_vs,
    fisher,
    fisher_vs
  )
)
write.csv2(
  data,
  "results/data/multiple_hypothesis_testing_fig4.csv",
  row.names = FALSE
)


#### SORTING THE PVALUES BEFORE PLOTTING
sort_fn <- function(p, decreasing) {
  c(
    sort(
      p[1:n_true],
      decreasing = decreasing
    ),
    sort(
      p[(n_true+1):(n_true + n_false)],
      decreasing = decreasing
    )
  )
}
n_seeds <- 1000
n_hypotheses <- 20
n_true <- 10
n_false <- 10
seed <- 34328165

set.seed(seed)
average_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(algorithm1(e(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))), decreasing=TRUE)
  )
)
names(average_sorted) <- c("hypothesis", "simulation_idx", "e_value")
average_sorted[, "method"] <- "Average"
average_sorted[, "p_value"] <- calibrate_e_to_p(average_sorted[, "e_value"])
average_sorted[, "bound"] <- "Bona fide"

set.seed(seed)
average2_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(algorithm1(e2(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))),decreasing = TRUE)
  )
)
names(average2_sorted) <- c("hypothesis", "simulation_idx", "e_value")
average2_sorted[, "method"] <- "Average2"
average2_sorted[, "p_value"] <- calibrate_e_to_p(average2_sorted[, "e_value"])
average2_sorted[, "bound"] <- "Bona fide"


set.seed(seed)
product_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(algorithm2(e(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))), decreasing=TRUE)
  )
)
names(product_sorted) <- c("hypothesis", "simulation_idx", "e_value")
product_sorted[, "method"] <- "Product"
product_sorted[, "p_value"] <- calibrate_e_to_p(product_sorted[, "e_value"])
product_sorted[, "bound"] <- "Bona fide"

set.seed(seed)
simes_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(simes_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    ), decreasing=FALSE)
  )
)
names(simes_sorted) <- c("hypothesis", "simulation_idx", "p_value")
simes_sorted[, "method"] <- "Simes"
simes_sorted[, "e_value"] <- reciprocal(simes_sorted[, "p_value"])
simes_sorted[, "bound"] <- "Global UB"

simes_vs_sorted <- simes_sorted
simes_vs_sorted[, "method"] <- "Simes VS"
simes_vs_sorted[, "e_value"] <- vs_p(simes_vs_sorted[, "p_value"])
simes_vs_sorted[, "bound"] <- "VS UB"

set.seed(seed)
bonferroni_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(bonferroni_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    ), decreasing=FALSE)
  )
)
names(bonferroni_sorted) <- c("hypothesis", "simulation_idx", "p_value")
bonferroni_sorted[, "method"] <- "Bonferroni"
bonferroni_sorted[, "e_value"] <- reciprocal(bonferroni_sorted[, "p_value"])
bonferroni_sorted[, "bound"] <- "Global UB"

bonferroni_vs_sorted <- bonferroni_sorted
bonferroni_vs_sorted[, "method"] <- "Bonferroni VS"
bonferroni_vs_sorted[, "e_value"] <- vs_p(bonferroni_vs_sorted[, "p_value"])
bonferroni_vs_sorted[, "bound"] <- "VS UB"

set.seed(seed)
fisher_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(fisher_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    ), decreasing=FALSE)
  )
)
names(fisher_sorted) <- c("hypothesis", "simulation_idx", "p_value")
fisher_sorted[, "method"] <- "Fisher"
fisher_sorted[, "e_value"] <- reciprocal(fisher_sorted[, "p_value"])
fisher_sorted[, "bound"] <- "Global UB"

fisher_vs_sorted <- fisher_sorted
fisher_vs_sorted[, "method"] <- "Fisher VS"
fisher_vs_sorted[, "e_value"] <- vs_p(fisher_vs_sorted[, "p_value"])
fisher_vs_sorted[, "bound"] <- "VS UB"

data <- do.call(
  rbind,
  list(
    average_sorted,
    average2_sorted,
    product_sorted,
    simes_sorted,
    simes_vs_sorted,
    bonferroni_sorted,
    bonferroni_vs_sorted,
    fisher_sorted,
    fisher_vs_sorted
  )
)
write.csv2(
  data,
  "results/data/multiple_hypothesis_testing_fig3_sorted.csv",
  row.names = FALSE
)

#### SORTING THE PVALUES BEFORE PLOTTING
n_seeds <- 1000
n_hypotheses <- 200
n_true <- 10
n_false <- 10
seed <- 54723765

set.seed(seed)
average_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(algorithm1(e(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))), decreasing=TRUE)
  )
)
names(average_sorted) <- c("hypothesis", "simulation_idx", "e_value")
average_sorted[, "method"] <- "Average"
average_sorted[, "p_value"] <- calibrate_e_to_p(average_sorted[, "e_value"])
average_sorted[, "bound"] <- "Bona fide"

set.seed(seed)
average2_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(algorithm1(e2(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))),decreasing = TRUE)
  )
)
names(average2_sorted) <- c("hypothesis", "simulation_idx", "e_value")
average2_sorted[, "method"] <- "Average2"
average2_sorted[, "p_value"] <- calibrate_e_to_p(average2_sorted[, "e_value"])
average2_sorted[, "bound"] <- "Bona fide"


set.seed(seed)
product_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(algorithm2(e(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))), decreasing=TRUE)
  )
)
names(product_sorted) <- c("hypothesis", "simulation_idx", "e_value")
product_sorted[, "method"] <- "Product"
product_sorted[, "p_value"] <- calibrate_e_to_p(product_sorted[, "e_value"])
product_sorted[, "bound"] <- "Bona fide"

set.seed(seed)
simes_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(simes_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    ), decreasing=FALSE)
  )
)
names(simes_sorted) <- c("hypothesis", "simulation_idx", "p_value")
simes_sorted[, "method"] <- "Simes"
simes_sorted[, "e_value"] <- reciprocal(simes_sorted[, "p_value"])
simes_sorted[, "bound"] <- "Global UB"

simes_vs_sorted <- simes_sorted
simes_vs_sorted[, "method"] <- "Simes VS"
simes_vs_sorted[, "e_value"] <- vs_p(simes_vs_sorted[, "p_value"])
simes_vs_sorted[, "bound"] <- "VS UB"

set.seed(seed)
bonferroni_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(bonferroni_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    ), decreasing=FALSE)
  )
)
names(bonferroni_sorted) <- c("hypothesis", "simulation_idx", "p_value")
bonferroni_sorted[, "method"] <- "Bonferroni"
bonferroni_sorted[, "e_value"] <- reciprocal(bonferroni_sorted[, "p_value"])
bonferroni_sorted[, "bound"] <- "Global UB"

bonferroni_vs_sorted <- bonferroni_sorted
bonferroni_vs_sorted[, "method"] <- "Bonferroni VS"
bonferroni_vs_sorted[, "e_value"] <- vs_p(bonferroni_vs_sorted[, "p_value"])
bonferroni_vs_sorted[, "bound"] <- "VS UB"

set.seed(seed)
fisher_sorted <- melt(
  replicate(
    n_seeds,
    sort_fn(fisher_closure(
      p(rnorm(n_hypotheses, c(rep(-4, n_false), rep(0, n_true)), 1))
    ), decreasing=FALSE)
  )
)
names(fisher_sorted) <- c("hypothesis", "simulation_idx", "p_value")
fisher_sorted[, "method"] <- "Fisher"
fisher_sorted[, "e_value"] <- reciprocal(fisher_sorted[, "p_value"])
fisher_sorted[, "bound"] <- "Global UB"

fisher_vs_sorted <- fisher_sorted
fisher_vs_sorted[, "method"] <- "Fisher VS"
fisher_vs_sorted[, "e_value"] <- vs_p(fisher_vs_sorted[, "p_value"])
fisher_vs_sorted[, "bound"] <- "VS UB"

data <- do.call(
  rbind,
  list(
    average_sorted,
    average2_sorted,
    product_sorted,
    simes_sorted,
    simes_vs_sorted,
    bonferroni_sorted,
    bonferroni_vs_sorted,
    fisher_sorted,
    fisher_vs_sorted
  )
)
write.csv2(
  data,
  "results/data/multiple_hypothesis_testing_fig4_sorted.csv",
  row.names = FALSE
)