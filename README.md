# Semialgebraic Hypothesis Testing with Incomplete U-Statistics: Practical Issues

V.1.0 Update: Apr/2025

This repository contains code implementing methods described in the paper "Methodological considerations for semialgebraic hypothesis testing with incomplete U-statistics" (available at ...).

A portion of the R and c++ code in this repository is adapted from the implementation in the TestGGM R package, which is available at https://github.com/NilsSturma/TestGGM. We gratefully acknowledge the authors of TestGGM for their valuable work.

This repository includes code implemented in R, C++ and Macaulay2.


## Code for hypothesis testing with incomplete U-statistics:
The code is based on the following files:

- **[`src/SDL-test.R`](src/SDL-test.R)**  
  Contains the main high-level function `test_U_stat`, which runs the SDL test.  
  └ **[`SDL-test_functions.cpp`](SDL-test_functions.cpp)**: Additional C++ functions for implementing the SDL test.

* _src/SDL-test.R_
* _src/SDL-test_functions.cpp_
* _src/data_plots_functions.R_
* 



## Additional Macaulay2 code:
* _singular_locus_CFN.m2_ This file can be used to compute the singular locus of the CFN model on a 4-leaf binary tree.
* _generating_sets_CFN.m2_ This file contains the code presented in Appendix C.1, used to compute the 5 different constrain sets (CDD, CDM, PDM, CDR and PDR) used thoughout the paper.
