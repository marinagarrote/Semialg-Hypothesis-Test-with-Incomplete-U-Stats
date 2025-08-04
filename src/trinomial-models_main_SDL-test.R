# Main file for the SDL test on trinomial models used for the simulations in
# Section 3 in `Methodological considerations for semialgebraic  
#               hypothesis testing with incomplete U-statistics`

# ------------------------------------------------- #
# ------ Load packages & auxiliary functions ------ #
# ------------------------------------------------- #
## change to your local path
path = ""
source(paste0(path, "SDL-test.R"))
source(paste0(path, "trinomial-models_utils.R"))


# ------------------------------------------------- #
# ---------- Choose model and parameters ---------- #
# ------------------------------------------------- #

## Possible models: model-1, model-2, model-3, model-4, HW
model <- "model-1"

## Fix the SDL test parameters
n = 300   #  n: Sample size
m = 6     #  m: Number of arguments taken by the kernel function h
N = 1000  #  N: Average number of terms in the sum defining the incomplete U-statistic.
n1 = 80   #  n1: Number of divide-and-conquer estimators to construct.
E = 500   #  E: Number of bootstrap iterations.

num_permutations = 0      #  Total number of permutations: permutations_symetrization_h(model, m, 1)
#                            - if num_permutations == 0: No permutationon the subsets Xm when computing h.
#                            - if num_permutations == 1: It takes one random permutation of the rows of Xm when computing h.
original_constraints = T  #  original_constraints: (T/F) Whether or not to add the original contrains to the combex combinations
num_cc_constraints = 10   #  num_cc_constraints: Number of extra contrains obtained by combex combination of the original ones

numCores <- 4  #  numCores: Number of cores used if any function is executed using parallelization.
#              #   Use detectCores(logical = FALSE) to get the number of physical cores


# ------------------------------------------------- #
# --- Algorithm Complexity & Computational time --- #
# ------------------------------------------------- #

# ## Uncomment the following lines to check the algorithm Complexity & Computational time
# ## Complexity
# iterations_H_Computation(model, m, N, num_permutations) ## Number of times we compute the kernel function h() for subsets of the data in calculate_H
# iterations_G_Computation(model, n, m, n1, num_permutations) ## Number of times we compute the kernel function h() for subsets of the data in calculate_H

# # times_ustat: Returns the min/mean/max time to compute the functions compute_H and compute_G using or not parallelization for rep_times times.
# #              This can be useful to decide whether use or not paralellization in test_U_stat()
# library(microbenchmark)
# rep_times = 10  # Number of times to run the functions
# times_ustat(model, n, m, N, n1, numCores, rep_times, num_permutations, original_constraints, num_cc_constraints)


# ------------------------------------------------- #
# ------- Compute p-value using the SDL test ------ #
# ------------------------------------------------- #

# Coordinates (x,y,z) with x + y + z = 1 for trinomial models. 
# For convinience we take x + y + z = n
x = 10
y = 45
z = n - x - y

X = generate_dataset_trinomial(x, y, z) # Randomly generate dataset X from coordinates (x,y,z).

pstat = test_U_stat(model, X, m, N, E, n1, num_permutations, original_constraints, num_cc_constraints, 
                    H_parallel = F, G_parallel = F, numCores = numCores)
p_value = pstat$PVAL


# ------------------------------------------------- #
# ------- Compute p-value using the SDL test ------ #
# ------------------------------------------------- #

## Uncomment to produce a simplex plot

# initializeSimplexPlot(model, n, m, N, n1, num_permutations, original_constraints, num_cc_constraints)
# step_size = 5
# for (i in seq(10, n-2, step_size)) { 
#   for (j in seq(10, n-i-1, step_size)) {
#     # Data
#     X = generate_dataset_trinomial(i, j, n - i - j) # Generate data
# 
#     # Compute p-value
#     pstat = test_U_stat(model, X, m, N, E, n1, num_permutations, original_constraints, num_cc_constraints,
#                         H_parallel = F, G_parallel = F, numCores = numCores)
#     p_value = pstat$PVAL
# 
#     # Plot the p-values in the simplex
#     pcolor = setcolors[1 + sum(p_value < alphas)] # determine color for plotting
#     simplexPoint(c(x,y,z), pch = "*", col = pcolor) # plot point
#   }
# }



