# Methodological considerations for semialgebraic hypothesis testing with incomplete U-statistics

V.1.0 Update: July/2025

This repository contains code implementing methods described in the paper:
[1] _Methodological considerations for semialgebraic hypothesis testing with incomplete U-statistics_, by D Barnhill, M Garrote-López, E Gross, M Hill,
B Kagy, J A Rhodes, AND J Z Zhang, available at ...

A portion of the R and c++ code in this repository is adapted from the implementation in the TestGGM R package, which is available at [`github.com/NilsSturma/TestGGM`](https://github.com/NilsSturma/TestGGM). We gratefully acknowledge the authors of TestGGM for their work and allowing us to use their code.

This repository includes code implemented in `R`, `C++` and `Macaulay2`.


## Code for hypothesis testing with incomplete U-statistics:
The code is based on the following files:

- **[`src/SDL-test.R`](src/SDL-test.R)**  
  Contains the main high-level function `test_U_stat`, which runs the SDL test.  
  └ **[`src/SDL-test_functions.cpp`](src/SDL-test_functions.cpp)**: Additional C++ functions for implementing the SDL test.

- **[`src/main_incompleteU-statistics_trinomials.R`](src/main_incompleteU-statistics_trinomials.R)**  
  Additional scripts for incomplete U-statistics for the trinomial models (see Section 3 in [1]).

- **[`src/main-incompleteU-statistics-CFN.R`](src/run-CFN-experiment.R)**
  Code for running experiments using the CFN 4-leaf tree model (see Section 4 in [1]).  
  └ **[`src/CFN-model_utils.R`](src/CFN-model_utils.R)**: Functions for creating and modifying CFN-type models.

- **[`src/plot_examles.R`](src/run-CFN-experiment.R)**
  Code for ploting ....
  └ **[`src/plot_utils.R`](src/CFN-model_utils.R)**: Additional functions for plotting.
  └ CFN data in **[`data/big_simulation_CFN.R`](data/big_simulation_CFN)**: Data...


## Additional Macaulay2 code:
* [`singular_locus_CFN.m2`](singular_locus_CFN.m2) This file can be used to compute the singular locus of the CFN model on a 4-leaf binary tree.
* [`docs/generating_sets_CFN.m2`](docs/generating_sets_CFN.m2) This file contains the code presented in Appendix C.1, used to compute the 5 different constrain sets (CDD, CDM, PDM, CDR and PDR) used thoughout the paper. It also contains the list of invariants for the CFN model for the 3 tree topologies of 4-leaf binary trees.


To discuss: Add docs/TestModels.pdf? or similar R code?


