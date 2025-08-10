# Methodological considerations for semialgebraic hypothesis testing with incomplete U-statistics

V.1.0 Update: July/2025

This repository contains code implementing methods described in the paper:

[1] *Methodological considerations for semialgebraic hypothesis testing with incomplete U-statistics*,  
D. Barnhill, M. Garrote-López, E. Gross, M. Hill, B. Kagy, J. A. Rhodes, and J. Z. Zhang.  
Available at [arxiv.org/abs/2507.13531](https://arxiv.org/abs/2507.13531).


A portion of the `R` and `C++` code in this repository is adapted from the implementation in the TestGGM R package, which is available at [`github.com/NilsSturma/TestGGM`](https://github.com/NilsSturma/TestGGM). We gratefully acknowledge the authors of TestGGM for their work and allowing us to use their code.

This repository includes code implemented in `R`, `C++` and `Macaulay2`.

## Code for hypothesis testing with incomplete U-statistics:
The code is based on the following files:

- **[`src/SDL-test.R`](src/SDL-test.R)**  
  Contains the main high-level function `test_U_stat`, which runs the SDL test.  
  └ **[`src/SDL-test_functions.cpp`](src/SDL-test_functions.cpp)**: Additional C++ functions for implementing the SDL test.

- **[`src/trinomial-models_main_SDL-test.R`](src/trinomial-models_main_SDL-test.R)**
  Code for running experiments for the trinomial models (see Section 3 in [1]).  
  └ **[`src/trinomial-models_utils.R`](src/trinomial-models_utils.R)**: Additional functions for generating data and plotting trinomial models.  

- **[`src/CFN-model_main_SDL-test.R`](src/CFN-model_main_SDL-test.R)**
  This file contains code for how to use the SDL test with the CFN 4-leaf tree model (see Section 4 in [1]).  
  └ **[`src/CFN-model_utils.R`](src/CFN-model_utils.R)**: Additional functions for working with the CFN model.

- **[`src/plot-CFN-results.R`](src/plot-CFN-results.R)** 
  This file contains code with examples illustrating how to generate the p-value histograms and treespace plots in Section 4 of [1].  
  └ The dataset used for this is Collection 1, see **[`the data readme`](data/README.md)**.


## Additional Macaulay2 code:
- **[`docs/singular_locus_CFN.m2`](docs/singular_locus_CFN.m2)** This file can be used to compute the singular locus of the CFN model on a 4-leaf binary tree.

- **[`docs/generating_sets_CFN.m2`](docs/generating_sets_CFN.m2)** This file contains the code presented in Appendix C.1, used to compute the 5 different constrain sets (CDD, CDM, PDM, CDR and PDR) used thoughout the paper. It also contains the list of invariants for the CFN model for the 3 tree topologies of 4-leaf binary trees.


To discuss: Add docs/TestModels.pdf? or similar R code?


