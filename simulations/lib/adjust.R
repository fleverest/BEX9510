# e-values

# Algorithm 1 for combining (not necessarily independent) e-values
algorithm1 <- function(e_values) {
  pi <- order(e_values)
  s <- cumsum(e_values[pi])
  sapply(
    seq_along(e_values), \(k) {
      min(
        sapply(
          seq_len(k),
          \(i) mean(c(e_values[pi[k]], s[seq_len(i)]))
        )
      )
    }
  )[order(pi)]
}

# Algorithm 2 for combining sequential e-values
algorithm2 <- function(e_values) {
  e_values * prod(e_values[e_values < 1])
}


# p-values

# This function implements the p-value adjustments described in the paper
# "Adjusted P-Values for Simultaneous Inference" by S. Paul Wright.
simes_closure <- function(p_values) {
  n <- length(p_values)
  # Sort p_values and keep permutation to unsort later
  pi <- order(order(p_values))
  p_values <- sort(p_values)
  a <- p_values

  for (m in seq(n, 2, by = -1)) {
    c_values <- numeric(n)
    # "Very least significant subset of size m"
    vlssm <- seq(n - m + 1, n, by = 1)
    cmin <- min(m * p_values[vlssm] / (vlssm + m - n))
    c_values[vlssm] <- cmin

    if (m < n) {
      c_values[seq(1, n - m)] <- sapply(
        seq(1, n - m),
        \(i) min(cmin, m * p_values[i])
      )
    }
    a[a < c_values] <- c_values[a < c_values]
  }
  return(a[pi])
}

bonferroni_closure <- function(p) {
  n <- length(p)
  pi <- order(order(p))
  p <- sort(p)
  pnew <- numeric(n)
  for (i in seq_len(n)) {
    pnew[i] <- min(1, max(sapply(seq_len(i), \(j) {(n - j + 1) * p[j]})))
  }
  pnew[pi]
}

fisher_closure <- function(p) {
  n <- length(p)
  pi <- order(order(p))
  p <- sort(p)
  pnew <- numeric(n)
  for (k in seq_len(n)) {
    pnew[k] <- max(sapply(seq(k, n), \(j) {fisher_merge_p(p[c(k, seq(j, n))])}))
  }
  pnew[pi]
}