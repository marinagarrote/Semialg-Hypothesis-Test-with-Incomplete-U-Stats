#-------------------------------------------------------------------------------
# CFN-model_main_SDL-test.R -- This file contains example code illustrating
# how we implement the SDL method, using data generated according to the CFN
# 4-leaf tree model, to test hypotheses about the topology of the tree
# parameter (Section 4 in the paper)
# -------------------------------------------------------------------------------
#
#
# HOW TO RUN: copy-paste the code into R session with working directory
# 'src/'.

# ------------------------------------------------- #
# ------ Load packages & auxiliary functions ------ #
# ------------------------------------------------- #
path = ""
source(paste0(path, "SDL-test.R"))
source(paste0(path, "CFN-model_utils.R"))



# ------------------------------------------------- #
# ---------- Choose model and parameters ---------- #
# ------------------------------------------------- #

# Here we'll use the same parameter names as we use in the paper (see also
# SDL-test.R for a glossary).
n = 1000 
m = 12
N = 1000
n1 = 80
A = 5000 
s = 100
r = 0
original_constraints = T

## Choose generating set (i.e. the specific choice of polynomial constraints
## defining the submodel)
generating_set = "CDD" # alternatives: CDR, CDM, PDM, PDR
use_internal_branch_inequality = TRUE # equation 4.3 in the paper

# Save the information about the generating set in a format is readable by the
# functions in SDL-test_functions.cpp
model_type = make_CFN_model_string(generating_set,use_internal_branch_inequality)

# Set the number of cores
numCores = c(6,6)
# Commentary: there are two steps in the computation of the SDL test which
# are paralellized in our implmentation: (see steps 5 and 6 in the function
# test_U_stat). The variable numCores specifies how many cores to use in
# those two steps. On our hardware, setting both values to 6 was optimal.

# Next, we'll define some null hypotheses to be tested. In this case, we will
# simultaneously test three null hypotheses named H_12, H_13, and H_14. Here,
# H_12 is the null hypothesis that the data comes from the CFN model on a tree
# with topology 12|34. H_13 and H_14 correspond to tests of topologies 13|24
# and 14|23 respectively.
H_12 = paste0(model_type,"-12|34")
H_13 = paste0(model_type,"-13|24")
H_14 = paste0(model_type,"-14|23")
# Note that the choice of generating polynomials used to represent the above
# hypotheses is encoded in the string `model_type`.

# Next, we set the true model parameters which will be used to generate the
# data.
true_topology = "12|34" # (alternatives: 13|24, 14|23)
branch_lengths = c(.2, .3, .24, .1, .17) # (can be any nonnegative numbers)
# Commentary: The true tree parameter consists of an unrooted binary quartet
# topology together with a vector of five branch length parameters, which
# specify the expected number of mutations occuring leaf 1, leaf 2, leaf 3,
# leaf 4, and the internal branch (in that order). In this example, the true
# tree topology is "12|34" so the null hypothesis H_12 is a true null
# hypothesis, whereas H_13 and H_13 are not.



# ------------------------------------------------- #
# ------ Generate data and compute p-values ------- #
# ------------------------------------------------- #

# Next, we generate data according to the CFN process with tree parameter
# specified above.
CFN_distribution = compute_CFN_distribution(branch_lengths, true_topology) 
X = t(rmultinom(n, size=1, prob = CFN_distribution))
# Commentary: Here, the data X consists of n iid multinomial trials with
# distribution specified by the vector CFN_distribution. In our code, we
# represent the outcome of each multnomial trial with one of the standard
# basis vectors e₁, e₂, ..., e₈. This is why X is a matrix with n rows
# (=samples) and eight columns. The biological interpretation of the data is
# that each time the CFN process is run, there are 8 possible site pattern
# probabilities that may be observed at the four leaves of the tree (namely,
# xxxx, xxxy, xxyx, xxyy, xyxx, xyxy, xyyx, xyyy), and the standard basis
# vector tells us which one was observed. 

# Finally, we now compute a p-value (and t-statistic) for each null
# hypothesis:
pvals_12 = test_U_stat(H_12, X, m, N, A, n1, s, original_constraints, r,
                       H_parallel = T, G_parallel = T, numCores = numCores)

pvals_13 = test_U_stat(H_13, X, m, N, A, n1, s, original_constraints, r,
                       H_parallel = T, G_parallel = T, numCores = numCores)

pvals_14 = test_U_stat(H_14, X, m, N, A, n1, s, original_constraints, r,
                       H_parallel = T, G_parallel = T, numCores = numCores)

pvals_12 # this p-value is probably big
pvals_13 # this p-value one is probably small
pvals_14 # this p-value one is probably small

