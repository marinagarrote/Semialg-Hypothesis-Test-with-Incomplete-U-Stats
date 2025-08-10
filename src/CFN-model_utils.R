#' Function to compute the CFN Probability Model for a 4-Leaf Tree
#'
#' Computes the CFN probability vector corresponding to a given 4-leaf tree
#' topology and edge parameters, based on the formula from Semple and Steel's
#' book 'Phylogenetics' (Corollary 8.6.6, p. 201).
#'
#' @param branch_lengths A numeric vector of length 5, corresponding to the
#'   edge parameters for [leaf1, leaf2, leaf3, leaf4, internalBranch]. The
#'   values are interpreted as evolutionary distances, i.e., the expected
#'   number of mutations along the edge.
#' 
#' @param quartet_topology A string specifying the unrooted quartet topology.
#'   Must be one of `"12|34"`, `"13|24"`, or `"14|23"`.
#'
#' @return A numeric vector of length 8 representing the computed site pattern
#'   probabilities under the CFN process (on a tree with branch lenghths
#'   `branch_lengths' and unrooted topology `quartet_topology`), with order:
#'   pxxxx, pxxxy, pxxyx, pxxyy, pxyxx, pxyxy, pxyyx, pxyyy
#'
#' @examples
#' branch_lengths <- c(0.2, 0.3, 0.4, 0.5, 0.6)
#' compute_probability_vector(branch_lengths, "12|34")
#'
#' @export
compute_CFN_distribution <- function(branch_lengths, quartet_topology) {
  # First, convert branch lengths to Hadamard parameters
  θ <- exp(-2*branch_lengths)

  # Assign assign values 1,2,3,4 to the variables j1, j2, k1, and k2 according
  # to the given quartet topology.
  if (quartet_topology == "12|34") {
    j1 <- 1; j2 <- 2;
    k1 <- 3; k2 <- 4
  } else if (quartet_topology == "13|24") {
    j1 <- 1; j2 <- 3;
    k1 <- 2; k2 <- 4
  } else if (quartet_topology == "14|23") {
    j1 <- 1; j2 <- 4;
    k1 <- 2; k2 <- 3
  } else {
    stop("quartet_topology must be one of '12|34', '13|24', or '14|23'")
  }
  # Define site pattern list
  S <- list(
    c(1, 1, 1, 1), c(1, 1, 1,-1), c(1, 1,-1, 1), c(1, 1,-1,-1),
    c(1,-1, 1, 1), c(1,-1, 1,-1), c(1,-1,-1, 1), c(1,-1,-1,-1)
  )
  # Apply formula from Corollary 8.6.6.
  prob_vector <- sapply(1:8, function(i) {
    (1 +
       S[[i]][j1] * S[[i]][j2] * θ[j1] * θ[j2] +
       S[[i]][k1] * S[[i]][k2] * θ[k1] * θ[k2] +
       S[[i]][j1] * S[[i]][k1] * θ[j1] * θ[k1] * θ[5] +
       S[[i]][j1] * S[[i]][k2] * θ[j1] * θ[k2] * θ[5] +
       S[[i]][j2] * S[[i]][k1] * θ[j2] * θ[k1] * θ[5] +
       S[[i]][j2] * S[[i]][k2] * θ[j2] * θ[k2] * θ[5] +
       S[[i]][1] * S[[i]][2] * S[[i]][3] * S[[i]][4] * θ[1] * θ[2] * θ[3] * θ[4]) / 8
  })
  return(prob_vector)
}



# Function to save the information about the generating set in a string format
# that is readable by the functions in SDL-test_functions.cpp
make_CFN_model_string <- function(generating_set, use_internal_branch_inequality){
  if (use_internal_branch_inequality==FALSE) {
    if (generating_set == "CDD") { return("CFN-BD-NI----NL---")
    } else if (generating_set == "CDR") { return("CFN-BR-NI----NL---")
    } else if (generating_set == "CDM") { return("CFN-BM-NI----NL---")
    } else if (generating_set == "PDR") { return("CFN-UR-NI----NL---")
    } else if (generating_set == "PDM") { return("CFN-UM-NI----NL---")}
  }
  else if (use_internal_branch_inequality==TRUE) {
    if (generating_set == "CDD") { return("CFN-BD-BI1H--NL---")
    } else if (generating_set == "CDR") { return("CFN-BR-BI1H--NL---")
    } else if (generating_set == "CDM") { return("CFN-BM-BI1H--NL---")
    } else if (generating_set == "PDR") { return("CFN-UR-BI1H--NL---")
    } else if (generating_set == "PDM") { return("CFN-UM-BI1H--NL---")}
  }
}

