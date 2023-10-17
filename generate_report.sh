#!/bin/bash

# Run all simulations
./simulations/combining_independent.R
./simulations/multiple_hypothesis_testing.R

# Generate plots
./plots/figure_1.R
./plots/figure_1_pvals.R
./plots/figure_2.R
./plots/figure_2_pvals.R
./plots/figure_3.R
./plots/figure_3_box.R
./plots/figure_3_sorted.R
./plots/figure_4.R
./plots/figure_4_box.R
./plots/figure_4_sorted.R

# TODO: Render quarto presentation