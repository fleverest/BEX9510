# This file contains a selection of (i)e-merging and (i)p-mergine functions from
# the 2020 paper "E-values: Calibration, combination, and applications", by
# Vladimir Vovk and Ruodu Wang.

#### p-merging

# The bonferroni correction for dependent p-values
bonferroni_merge_p <- function(p) min(length(p) * min(p), 1)

# Simes' method for combining positively-dependent p-values
simes_merge_p <- function(p) min(sort(p) * length(p) / seq_along(p))

#### ip-merging

# The fisher
fisher_merge_p <- function(p, log.p = FALSE) {
  pchisq(
    -2 * sum(log(p)),
    df = 2 * length(p),
    lower.tail = FALSE,
    log.p = log.p
  )
}

#### e-merging excluded because one can just use mean(), prod()...