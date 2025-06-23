# Semialgebraic Hypothesis Testing with Incomplete U-Statistics: Practical Issues

V.1.0 Update: Apr/2025

This repository contains code implementing methods described in the paper "Methodological considerations for semialgebraic hypothesis testing with incomplete U-statistics" (available at ...).

A portion of the R and c++ code in this repository is adapted from the implementation in the TestGGM R package, which is available at https://github.com/NilsSturma/TestGGM. We gratefully acknowledge the authors of TestGGM for their valuable work.

This repository includes code implemented in `R`, `C++` and `Macaulay2`.


## Code for hypothesis testing with incomplete U-statistics:
The code is based on the following files:

- **[`src/SDL-test.R`](src/SDL-test.R)**  
  Contains the main high-level function `test_U_stat`, which runs the SDL test.  
  └ **[`src/SDL-test_functions.cpp`](src/SDL-test_functions.cpp)**: Additional C++ functions for implementing the SDL test.

- **[`src/main_incompleteU-statistics_trinomials.R`](src/main_incompleteU-statistics_trinomials.R)**  
  Additional scripts for incomplete U-statistics in coalescent models (details TBD).

- **[`src/main-incompleteU-statistics-CFN.R`](src/run-CFN-experiment.R)**
  Code for running experiments using the CFN 4-leaf tree model.  
  └ **[`src/CFN-model_utils.R`](src/CFN-model_utils.R)**: Functions for creating and modifying CFN-type models.

- **[`src/plot_examles.R`](src/run-CFN-experiment.R)**
  Code for ploting ....
  └ **[`src/plot_utils.R`](src/CFN-model_utils.R)**: Additional functions for plotting.



## Additional Macaulay2 code:
* [`singular_locus_CFN.m2`](singular_locus_CFN.m2) This file can be used to compute the singular locus of the CFN model on a 4-leaf binary tree.
* [`generating_sets_CFN.m2`](generating_sets_CFN.m2) This file contains the code presented in Appendix C.1, used to compute the 5 different constrain sets (CDD, CDM, PDM, CDR and PDR) used thoughout the paper.
* [`list_constraints_3topologies_CFN.m2`](list_constraints_3topologies.m2) List of invariants for the CFN model on the 3 tree topologies of 4-leaf binary trees.

* <span style="color:red"> To discuss: Add docs/TestModels.pdf? or similar R code? </span>.
