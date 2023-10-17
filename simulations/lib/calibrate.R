# This file contains a selection of e-to-p and p-to-e calibrators from the 2020
# paper "E-values: Calibration, combination, and applications", by Vladimir Vovk
# and Ruodu Wang.

# The admissible e-to-p calibrator.
calibrate_e_to_p <- Vectorize(\(e) min(1, 1 / e))

# The p-to-e calibrators in the family f_κ
f_kappa <- function(p, kappa) kappa * p ^ (kappa - 1)

# The p-to-e calibrator obtained by integrating the simple family over κ.
calibrate_p_to_e <- function(p) (1 - p + p * log(p)) / (p * log(p)^2)

# The Vovke-Sellke bound on the simple family of p-to-e calibrators.
vs_p <- Vectorize(\(p) if (p <= exp(-1)) -exp(-1) / (p * log(p)) else 1)

# The upper bound on all p-to-e callibrators: the reciprocal of a p-value.
reciprocal <- \(p) 1 / p