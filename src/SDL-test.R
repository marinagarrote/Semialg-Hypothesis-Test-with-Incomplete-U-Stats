# SDL-test.R -- contains the main high-level function 'test_U_stat' for
# implementing the SDL test.

library(Rcpp)
library(parallel)

# load the C++ functions
sourceCpp(paste0(path, "SDL-test_functions.cpp"))


# COMMENTARY: In this code, we have attempted to mirror the notation of Sturma et
# al as closely as possible. Here are some important definitions:

# ------------------------------------------------------------------------------ 
# -------------------- Arguments of function "test_U_stat" ---------------------
# ------------------------------------------------------------------------------

# model = a string identifying which model to test:
#         - Accepted trinomial models are 'cut1', 'T1', 'cut', 'T3', and 'HW'. 
#         - For testing CFN models, see documentation in the file
#           'CFN-model_main_SDL-test.R'

# X = the data, taking the form of n row vectors. Each row vector is a standard
#     basis vector, indicating what was observed for that data point. 
#     - For trinomial models: rows are of the form (0,0,1), (0,1,0) and (1,0,0), 
#       and indicate which topology was observed for a particular gene.
#     - For CFN models, they are (0,0,0,0,0,0,0,1), ..., (1,0,0,0,0,0,0,0),
#       and indicate the observed site pattern for a given DNA site.
#       For more information, see commentary in 'CFN-model_main_SDL-test.R'

# m = number of arguments taken by the kernel function h. Determines the number
#     of columns in the matrix 'indices'. Needs to be divisble by the maximum
#     degree of the semi-algebraic constraints defining the submodel being tested.
 
# N = computational budget parameter. This is the average number of terms in 
#     the sum defining the incomplete U-statistic. (In particular, we randomly 
#     choose on average N indices from I_nm, where I_nm is set of m-tuples of 
#     {1,...,n}). User specified (up to a maximum of 0.7*(n choose m).)
 
# A = the number of bootstrap replicates to use when computing the p-values.
#     It just needs to be a large number.
 
# n1 = number of divide-and-conquer estimators to construct. Determines the
#      cardinality of S1. SDL suggest using n1 = n.

# num_permutations = the number of (randomly-chosen) permutations to use to
#                    symmetrize the kernel function. In the paper, this
#                    quantity is denoted s. If m is small, the kernel can be
#                    fully symmetrized by taking num_permutations = m!. If
#                    possible, we recommend doing this.

# original_constraints: A boolean flag, either T or F. If set to T (the
#                       default), then the test of the null hypothesis will
#                       use the polynomial constraints in the user-specified
#                       generating set. On the other hand, if
#                       original_constraints=F, then only convex combinations
#                       of those polynomial constraints are used.

# Ncc = the number of additional inequalities to be used which are obtained by
#       taking random convex combinations of the user-supplied polynomial
#       inequalities. This is a nonnegative number. It is denoted by r in the
#       paper.

# H_parallel: A boolean flag to enable parallel computation of the U-statistic sum.

# G_parallel: A boolean flag to enable parallel computation of the
#             divide-and-conquer estimates.

# numCores = a pair of positive integers, indiciating the number of cores to use
#            when computing H and G, respectively. This is an optional parameter
#            to be used if you want to parallelize the computation. We found
#            numCores = c(6,6) worked good, along with taking H_parallel = T and
#            G_parallel = T


# ------------------------------------------------------------------------------  
# -------------------- Other important notation/definitions --------------------
# ------------------------------------------------------------------------------
# n = number of data samples

# S1 = an index set used to construct divide-and-conquer estimators. (Note: we 
#      can choose this to be an arbitrary size n1 subset of {1,...,n}. As such,
#      in this code, we take S1 = {1,...,n1}.)

# p = total number of semi-algebraic constraints (i.e. inequality constraints).
#     In particular, the output of h is a vector length p.

# H = An intermediate matrix used to store the terms of the incomplete
#     U-statistic. This is an N_hat x p matrix. There are N_hat nonzero terms in
#     the sum defining the incomplete U-statistic. For each tuple (i₁,...,iₘ)
#     corresponding to one of these nonzero terms, we have a row vector
#     h(X_i₁,...,X_iₘ) which has length p. These are the N_hat rows of H. The
#     incompelete U statistic is then obtained by taking the column sums of H.

# G = Similar to H, but for the divide-and-conquer estimators.

# cov_diag = a (p x p)-diagonal matrix whose entries are of the form
#            \hat{\sigma}^{2}_{j}, which is defined on p.12 of Sturma et al.

# s : see entry for num_permutations above

# r : see entry for Ncc above

# ------------------------------------------------------------------------------ 
# -------------------- Main function to apply the SDL test ---------------------
# ------------------------------------------------------------------------------ 

test_U_stat <- function(model, X, m, N, A, n1, num_permutations, original_constraints,
                        Ncc, H_parallel = F, G_parallel = T, numCores=1){
  ## This is the main high-level function which implements the SDL test.

  ## Step 0. Determine number of cores to use for parallel processing
  if(length(numCores) == 1){
    numCoresH = numCores; # number of cores to use when computing H
    numCoresG = numCores  # number of cores to use when computing G
  }else{
    numCoresH = numCores[1]; numCoresG = numCores[2]
  }
  if(numCoresH == 1){H_parallel=F}
  if(numCoresG == 1){G_parallel=F}
  
  ## Step 1. Initial definitions
  n = nrow(X) # number of samples
  N = min(0.7*choose(n,m), N) # ensure N isn't too big
  
  ## Step 2. Determine N_hat by Bernoulli sampling 
  N_hat = stats::rbinom(1, choose(n,m), (N / choose(n,m)))
  
  ## Step 3. Create the list of tuples to be used in computing the incomplete U-statistic
  indices = do.call(rbind, random_combs(n, m, N_hat)[[1]])
  permuted_indices <- indices[sample(nrow(indices)), , drop = FALSE]
  
  ## Step 4. Compute a matrix of convex combinations by calling the C++ function
  ## convex_combination_coefficients)
  cc_coefficients = convex_combination_coefficients(model, original_constraints, Ncc)

  ## Step 5. Compute H.
  if(H_parallel){
    H = calculate_H_parallel(numCoresH, model, X, permuted_indices, num_permutations, cc_coefficients)
  }else{
    ## call the C++ function calculate_H
    H = calculate_H(model, X, permuted_indices, num_permutations, cc_coefficients)
  }
  
  ## Step 6. Compute matrix G on a size n1 subset of samples.
  if(G_parallel){
    G = calculate_G_parallel_outer(numCoresG, model, X, n1, m, num_permutations, cc_coefficients)
  }else{
    G = calculate_G(model, X, n1, m, num_permutations, cc_coefficients)
  }
  
  ## Step 7. Compute some simple statistics of H and G.
  H_mean = colMeans(H) # incomplete U-statistic based on Bernoulli sampling (Eq 2.1 in the SDL paper)
  H_centered = t(t(H) - H_mean)
  G_mean = colMeans(G)
  G_centered = t(t(G) - G_mean)
  
  ## Step 8. The diagonal of the approximate variance of H. 
  cov_H_diag = colSums(H_centered**2) / N_hat
  cov_G_diag = colSums(G_centered**2) / n1
  cov_diag = m**2 * cov_G_diag + (n/N) * cov_H_diag
  standardizer = cov_diag**(-1/2) # Vector for standardizing
  
  ## Step 9. Compute test statistic
  marginal_stats = sqrt(n) * H_mean
  test_stat =  max(standardizer * marginal_stats) # This is \mathcal{T} in SDL (p 13)
  
  ## Step 10. Do bootstrapping
  bootstrap_res = bootstrap_U(A, H_centered, G_centered)
  #print("Boostrap done!")
  U_A = bootstrap_res[[1]]
  U_B = bootstrap_res[[2]]
  U = m*U_A + sqrt(n/N) * U_B # I think each row(?) of U is U^{\#}_{n,N}, that
                              # the rows of U_A = U^{\#}_{n_1,g} and the rows of
                              # U_B = U^{\#}_{n_1,h}. I should changes these to
                              # Usharp, Usharp_g, and Usharp_h?
  U_standardized = t(t(U) * standardizer) # I think this is W
  results = matrixStats::rowMaxs(U_standardized) 
  
  ## Step 11. Compute p-value using bootstrap results
  pval = (1 + sum(results >= test_stat)) / (1+A)
  return(list("PVAL"=pval, "TSTAT"=test_stat))
}


# ------------------------------------------------------------------------------ 
# ---------------------------- Auxiliary functions  ----------------------------
# ------------------------------------------------------------------------------

call_calculate_H <- function(permuted_indices, f_model, f_data, f_num_permutations, f_cc_coefficients) {
  ## Call the C++ function calculuate_H
  calculate_H(f_model, f_data, permuted_indices, f_num_permutations, f_cc_coefficients)
}

call_calculate_G <- function(S1, f_model, f_data, f_m, f_num_permutations, f_cc_coefficients) {
  ## Call the C++ function calculuate_G
  calculate_G_S1(f_model, f_data, S1, f_m, f_num_permutations, f_cc_coefficients)
}
                                                                                
calculate_H_parallel <- function(numCores, model, X, permuted_indices, num_permutations, cc_coefficients){
  ## compute the U-statistic sum in parallel
  permuted_indices_list <- partition_matrix(permuted_indices, numCores)
  results <- mclapply(permuted_indices_list, call_calculate_H, 
                      f_model = model, f_data= X, f_num_permutations=num_permutations,
                      f_cc_coefficients=cc_coefficients, mc.cores = numCores);
  H <- do.call(rbind, results)
  return(H)
}

calculate_G_parallel_outer <- function(numCores, model, X, n1, m, num_permutations, cc_coefficients){
  ## compute the divide-and-conquer estimates in parallel
  n = nrow(X)
  S1 = sample(n, n1, replace = F)
  S1_list <- partition_vector(S1, numCores)
  results <- mclapply(S1_list, call_calculate_G, f_model = model, f_data = X, f_m = m,
                      f_num_permutations = num_permutations, f_cc_coefficients = cc_coefficients,
                      mc.cores = numCores);
  G <- do.call(rbind, results)
  return(G)
}




