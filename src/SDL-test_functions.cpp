// SDL-test_functions.cpp -- auxilliary functions for implementing the SDL test.

#include <Rcpp.h>
#include <random> 
#include <iostream>
#include <vector>
#include <algorithm>

using namespace Rcpp;

/* Document map (Functions marked with * are exported)

1. AUXILLIARY FUNCTIONS
   
   A) FUNCTIONS FOR DOING STUFF WITH VECTORS AND MATRICES
   - generate_sequence
   - extract_rows
   - get_entries
   - col_mean
   - partition_matrix*
   - partition_vector*
   - printIntegerVector
   - printNumericMatrix

   B) FUNCTIONS FOR READING MODEL STRINGS
   - is_a_CFN_model*: test if a model string is for CFN
   - getInvariantType*: extract the part of the string that specifies which algebraic invariants to use
   - getInternalType*: extract the part of the string that specifies which internal branch inequality to use
   - getLeafType*: extract the part of the string that specifies which which leaf inequalities to use
   - getQuartetTopology*: extract the part of the string that specifies which topology is to be tested
   
   C) PERMUTATIONS
   - random_combs*: randomly samples a specified number of unique subsets of size 'k' from the {1,...,n}

   D) FUNCTIONS FOR PARTITIONING {1,...,n}
   - log_factorial*: compute log(n!)
   - total_number_partitions*
   - concatenate_partitions*
   - get_random_k_partitions_efficient*

2. MODEL-SPECIFIC FUNCTIONS

   A) PARAMETER FUNCTIONS
   - value_s*: returns the degree of the model
   - value_p*: returns total number of inequalities defining the model

   B) KERNEL FUNCTIONS
   - h_tilde_cut1
   - h_tilde_T1
   - h_tilde_cut
   - h_tilde_T3
   - h_tilde_HW
   - h_tilde_CFN

3. FUNCTIONS TO COMPUTE INCOMPLETE U-STATISTIC

   A) FUNCTIONS TO IMPLMENT CONVEX COMBINATIONS OF CONSTRAINTS
      - convex_combination*
      - dirichlet_sample*
      - convex_combination_coefficients*

   B) MAIN KERNEL AND KERNEL SUMMATION FUNCTIONS
      - h*: the symmetrized kernel
      - calculate_H*

4. FUNCTIONS TO COMPUTE DIVIDE-AND-CONQUER ESTIMATORS
   - indices_ĝ_i1*
   - ĝ_i1
   - calculate_G_S1*
   - calculate_G*

5. BOOTSTRAP FUNCTION
   - bootstrap_U*


GRAPH OF SOME KEY DEPENDENCIES

h  (the symmetrized kernel function)
└─ convex_combination
    ├── value_p
    └── h_tilde
        └── h_tilde_MODEL  (i.e., any function in section 2B)
            └── value_s



TODO (as of 2024-12-04):
- [ ] check the correctness of function dirichlet_sample_cpp
- [ ] improve commentary for section 1D
- [ ] add commentary for section 4
- [ ] add commentary for section 5          
*/


////////////////////////////////////////////// 
/////////// 1. AUXILLIARY FUNCTIONS //////////
////////////////////////////////////////////// 


/////////////////////////////////////////////////////////////////////
////// 1A. FUNCTIONS FOR DOING STUFF WITH VECTORS AND MATRICES //////
/////////////////////////////////////////////////////////////////////

IntegerVector generate_sequence(int length, int initial_value){
  // Generate an integer vector of the specified length:
  // [initial_value, initial_value+1,...,initial_value + length - 1]
  IntegerVector v(length);
  for (int i=0; i < length ; i++){
    v[i] = initial_value + i;
  }
  return v;
}

NumericMatrix extract_rows(NumericMatrix A, IntegerVector indices) {
  int n = indices.size();
  int m = A.cols();
  NumericMatrix B(n, m);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      B(i, j) = A(indices(i)-1, j);
    }
  }
  return B;
}

IntegerVector get_entries(IntegerVector v, IntegerVector indices) {
  int n = indices.size();
  IntegerVector result(n);
  
  for (int i = 0; i < n; ++i) {
    result[i] = v[indices[i] - 1]; 
  }
  
  return result;
}

NumericVector col_mean(NumericMatrix A){
  int n = A.rows();
  int m = A.cols();
  
  NumericVector mean(m);
  for(int j = 0; j < m; ++j){
    double col_sum = 0.0;
    for(int i = 0; i < n; ++i){
      col_sum += A(i,j);
    }
    mean(j) = col_sum/n;
  }
  return mean;
}

// [[Rcpp::export]]
List partition_matrix(NumericMatrix M, int parts) {
  if(parts < 1) parts = 1;
  int nrows = M.nrow();
  int k = ceil((double)nrows / parts); 
  List M_list(parts);
  
  for (int i = 0; i < parts - 1; i++) {
    M_list[i] = M(Range((i * k), (i * k + k - 1)), _); 
  }
  M_list[parts - 1] = M(Range((parts - 1) * k, nrows - 1), _);
  
  return M_list;
}

// [[Rcpp::export]]
List partition_vector(NumericVector v, int parts) {
  int n = v.size();
  int k = floor(n / parts);
  List v_list(parts);
  
  for (int i = 0; i < parts - 1; i++) {
    v_list[i] = v[Range((i * k), (i * k + k - 1))]; 
  }
  v_list[parts - 1] = v[Range((parts - 1) * k, n - 1)]; 
  
  return v_list;
}

// Printing functions
void printIntegerVector(const IntegerVector& v) {
  for (int i = 0; i < v.size(); ++i) {
    Rcpp::Rcout << v(i) << " ";
  }
  Rcpp::Rcout << std::endl;
}

void printNumericMatrix(const NumericMatrix& M) {
  for (int i = 0; i < M.rows(); ++i) {
    for (int j = 0; j < M.cols(); ++j) {
      Rcpp::Rcout << M(i, j) << " ";
    }
    Rcpp::Rcout << std::endl;
  }
  Rcpp::Rcout << std::endl;
}


/////////////////////////////////////////////////////
////// 1B. FUNCTIONS FOR READING MODEL STRINGS //////
/////////////////////////////////////////////////////

// [[Rcpp::export]]
bool is_a_CFN_model(const std::string& model) {
  // Check if the first 3 characters are "CFN"
  return model.size() >= 3 && model.substr(0, 3) == "CFN";
}

// [[Rcpp::export]]
std::string getInvariantType(const std::string& model) { return model.substr(4, 2); }

// [[Rcpp::export]]
std::string getInternalType(const std::string& model) { return model.substr(7, 5); }

// [[Rcpp::export]]
std::string getLeafType(const std::string& model) { return model.substr(13, 5); }

// [[Rcpp::export]]
std::string getQuartetTopology(const std::string& model) {
  // Extract the last 5 letters from a string
  if (model.length() >= 5) {
    return model.substr(model.length() - 5);
  } else {
    return "";  // Return an empty string if the input is too short
  }
}


////////////////////////////////////
///////// 1C. PERMTUATIONS /////////
////////////////////////////////////

// [[Rcpp::export]]
List random_combs(int n, int k, int nr){
  // randomly sample 'nr' unique subsets of size 'k' from {1,...,n}
  std::set<std::set<int>> sub_sets;

  // initialize the vecto {1,...,n}
  IntegerVector v(n);
  for(int i=0; i < n; i++){
    v[i] = i+1;
  }
  
  std::set<int> s;
  IntegerVector w(k);

  // generate unique subsets
  while(sub_sets.size() < (unsigned)nr){
    w = sample(v,k);
    for (int i = 0; i < k; i++) {
      s.insert(w[i]);
    }

    // add the subset to the set of subsets
    sub_sets.insert(s);
    s.clear(); 
  }

  // convert output to a list and return
  List res = List::create(sub_sets);
  return(res);
}


/////////////////////////////////////////////////////////////////
///////// 1D. FUNCTIONS FOR PARTITIONING {1,...,n} //// ///////// 
/////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double log_factorial(int n) {
  // Auxilliary function for 'total_number_partitions'. Computes log(n!)
  double log_fact = 0;
  for (int i = 1; i <= n; ++i) {
    log_fact += std::log(i);
  }
  return log_fact;
}

// [[Rcpp::export]]
long total_number_partitions(int n, int k) {
  // Auxilliary function for 'get_random_k_partitions_efficient'
  if (k == 0 || n < k) return 0;  // edge cases
  int partition_size = n / k;
  // Use logarithmic scale to avoid overflow
  double log_total = log_factorial(n) - k * log_factorial(partition_size);
  // Convert back from log scale, handle large numbers carefully
  return static_cast<long>(std::round(std::exp(log_total)));
}

// [[Rcpp::export]]
std::vector<int> concatenate_partitions(std::vector<std::vector<int>> partition){
  // Auxilliary function for 'get_random_k_partitions_efficient'
  std::vector<int> result;
  for (const auto& part : partition) {
    result.insert(result.end(), part.begin(), part.end());
  }
  return result;
}

// [[Rcpp::export]]
IntegerMatrix get_random_k_partitions_efficient(int n, int k, int num_partitions) {
  
  if(num_partitions == 0){
    IntegerMatrix noPermutation(1,n);
    noPermutation(0,_) = generate_sequence(n, 1);
    return(noPermutation);
  }
  
  // Ensure that the partition is possible
  if (n % k != 0) {
    std::cerr << "Partition is not possible, n is not divisible by k" << std::endl;
    return IntegerMatrix(0, 0); // Return an empty matrix
  }
  
  int partition_size = n / k;
  long total_partitions = total_number_partitions(n, k);
  //int num_partitions = std::floor(total_partitions*prop);
  //num_partitions = std::min(std::max(num_partitions, 1), total_partitions);
  num_partitions = static_cast<int>(std::min(std::max(static_cast<long>(num_partitions), 1L), total_partitions));
  
  std::random_device rd;
  std::mt19937 gen(rd());
  
  std::set<std::vector<std::vector<int>>> unique_partitions; // To store unique partitions
  IntegerMatrix random_combinations(num_partitions, n); // Initialize with zero rows
  int row = 0;
  
  while (row < num_partitions) {
    std::vector<int> elements(n);
    std::iota(elements.begin(), elements.end(), 1); // Fill with 1 to n
    std::shuffle(elements.begin(), elements.end(), gen); // Shuffle randomly

    std::vector<std::vector<int>> partition;
    for (int j = 0; j < k; ++j) {
      std::vector<int> part(partition_size);
      for (int m = 0; m < partition_size; ++m) {
        part[m] = elements[j * partition_size + m];
      }
      std::sort(part.begin(), part.end());
      partition.push_back(part);
    }

    // Insert into set to check for uniqueness
    if (unique_partitions.insert(partition).second) {
      // If the partition is unique, add it to the matrix
      std::vector<int> part_conc =  concatenate_partitions(partition);
      random_combinations(row, _) = Rcpp::IntegerVector(part_conc.begin(), part_conc.end());
      ++row;
    }
  }

  return random_combinations;
}


//////////////////////////////////////////////////
/////////// 2. MODEL-SPECFIC FUNCTIONS ///////////
//////////////////////////////////////////////////

//////////////////////////////////////////////////////////
///////// 2A. MODEL-SPECIFIC PARAMETER FUNCTIONS ///////// 
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
int value_s(std::string model){
  // Input: a model string
  //
  // Output: a number s, which equals the maximum degree of the polynomials
  // inequalities defining the submodel.
  if(model == "model-1" || model == "model-2"){ return(1);
  }else if(model == "HW" || is_a_CFN_model(model)){ return(2);
  }else if(model == "model-3"){return(3);
  }else if(model == "model-4"){return(5);
  }else{
    Rcerr << "Unknown model: " << model << std::endl;
    Rcerr << " -> Accepted models: cut1, T1, cut, T3, CFN, HW." << model << std::endl;
    return(0);
  }
}

// [[Rcpp::export]]
int value_p(std::string model){
  // Input: a model string
  //
  // Output: a number p, which equals the number of polynomial inequalities
  // defining the submodel, i.e., p = 2*equalities + inequalities.
  if(model == "model-1" || model == "model-3" || model == "HW") {return(2);
  }else if(model == "model-2"){ return(3);
  }else if(model == "model-4"){ return(5);
  }else if(is_a_CFN_model(model)){ 
    int num_internal_constraints; int num_leaf_constraints;
    std::string internal_type=getInternalType(model);
    std::string leaf_type=getLeafType(model);

    if(internal_type == "UI4--"){num_internal_constraints = 4;
    }else if(internal_type == "UI2H-"){num_internal_constraints = 2;
    }else if(internal_type == "UI2NH"){num_internal_constraints = 2;
    }else if(internal_type == "BI1H-"){num_internal_constraints = 1;
    }else if(internal_type == "BI1NH"){num_internal_constraints = 1;
    }else if(internal_type == "BI1--"){num_internal_constraints = 1;
    }else if(internal_type == "NI---"){num_internal_constraints = 0;}

    if(leaf_type== "UL8--"){num_leaf_constraints = 8;
    }else if (leaf_type == "BL4--"){num_leaf_constraints = 4;
    }else if (leaf_type == "NL---"){num_leaf_constraints = 0;}
    
    return 4+num_internal_constraints+num_leaf_constraints;
    
  }else{
    Rcerr << "Unknown model: " << model << std::endl;
    Rcerr << " -> Accepted models: cut1, T1, cut, T3, CFN, HW." << model << std::endl;
    return(0);
  }
}


///////////////////////////////////////////////////////
///////// 2B. MODEL-SPECIFIC KERNEL FUNCTIONS ///////// 
///////////////////////////////////////////////////////

////// h_tilde functions for TRINOMIAL MODELS ////// 
// Model 1
NumericVector h_tilde_cut1(NumericMatrix Xm){
  NumericVector mean_Xm = col_mean(Xm);
  double eq1 = mean_Xm(1) - mean_Xm(2);
  return NumericVector::create(eq1, -eq1);
}

// Model 2
NumericVector h_tilde_T1(NumericMatrix Xm){
  NumericVector mean_Xm = col_mean(Xm);
  double eq1 = mean_Xm(1) - mean_Xm(2);
  double ineq1 = (double)1/3 - mean_Xm(0);
  return NumericVector::create(ineq1, eq1, -eq1);
}

// Model 3
NumericVector h_tilde_cut(std::string model, NumericMatrix Xm, IntegerVector perm){
  int eta = perm.size()/value_s(model);
  NumericVector mean_Xm_1 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1))));
  NumericVector mean_Xm_2 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1+eta))));
  NumericVector mean_Xm_3 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1+2*eta))));
  double factor1 = mean_Xm_1(0) - mean_Xm_1(1);
  double factor2 = mean_Xm_2(1) - mean_Xm_2(2);
  double factor3 = mean_Xm_3(0) - mean_Xm_3(2);
  double eq1 = factor1*factor2*factor3;
  return NumericVector::create(eq1, -eq1);
}

// Model 4
NumericVector h_tilde_T3(std::string model, NumericMatrix Xm, IntegerVector perm){ 
  int eta = perm.size()/value_s(model);
  NumericVector mean_Xm_1 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1))));
  NumericVector mean_Xm_2 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1+eta))));
  NumericVector mean_Xm_3 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1+2*eta))));
  NumericVector mean_Xm_4 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1+3*eta))));
  NumericVector mean_Xm_5 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1+4*eta))));
  
  double eq1 = (mean_Xm_1(0) - mean_Xm_1(1))*
               (mean_Xm_2(0) - mean_Xm_2(2))*
               (mean_Xm_3(1) - mean_Xm_3(2));
  
  double ineq1 = (mean_Xm_1(0) - mean_Xm_1(1))*(mean_Xm_2(0) - mean_Xm_2(1))*
                 (mean_Xm_3(0) - mean_Xm_3(2))*(mean_Xm_4(0) - mean_Xm_4(2))*
                 ((double)1/3 - mean_Xm_5(0));
  
  double ineq2 = (mean_Xm_1(1) - mean_Xm_1(0))*(mean_Xm_2(1) - mean_Xm_2(0))*
                 (mean_Xm_3(1) - mean_Xm_3(2))*(mean_Xm_4(1) - mean_Xm_4(2))*
                 ((double)1/3 - mean_Xm_5(1));
  
  double ineq3 = (mean_Xm_1(2) - mean_Xm_1(0))*(mean_Xm_2(2) - mean_Xm_2(0))*
                 (mean_Xm_3(2) - mean_Xm_3(1))*(mean_Xm_4(2) - mean_Xm_4(1))*
                 ((double)1/3 - mean_Xm_5(2));
  
  return NumericVector::create(eq1, -eq1, ineq1, ineq2, ineq3);
}

// Hardy-Weinberg
NumericVector h_tilde_HW(std::string model, NumericMatrix Xm, IntegerVector perm){
  int eta = perm.size()/value_s(model);
  NumericVector mean_Xm_1 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1))));
  NumericVector mean_Xm_2 = col_mean(extract_rows(Xm, get_entries(perm, generate_sequence(eta,1+eta))));
   
  double monomial1 = mean_Xm_1(1)*mean_Xm_2(1);
  double monomial2 = mean_Xm_1(0)*mean_Xm_2(2);
  double eq1 = monomial1 - 4*monomial2;
  
  return NumericVector::create(eq1, -eq1);
  
}


/////// h_tilde functions for the CFN MODEL /////// 
// [[Rcpp::export]]
NumericVector h_tilde_CFN(std::string model, NumericMatrix Xm, IntegerVector perm){
  // This is the nonsymmetric estimator with which the symmetric kernel h is defined for use when testing
  // CFN submodels. (We sum this over many, possibly randomly chosen permutations to get a symmetrized
  // kernel function).
  // 
  // Inputs:
  //
  // Xm -- this is the data, consist of m rows. Each row is a standard basis vector e₁,...,e₈. Each
  // standard basis vector corresponds to one of the 8 site patterns (aaaa aaab aaba aabb abaa abab abba
  // abbb), where a,b∈{0,1} are distinct nucleotides.
  //
  // perm -- an element of the symmetric group on {1,...,m} which is applied to the data. We use it to
  // permute the rows of X.
  //
  // model -- a string which indicates which polynomial inequalities to use. This must conform to the
  // following tagging scheme:
  //
  //    model = "CFN-XX-YYYYY-ZZZZZ-T"
  //
  // where:
  //
  //    XX ∈ {BD, UR, BR, UM, BM}
  //    YYYYY ∈ {UI4--, UI2H-, UI2NH, BI1H-, BI1NH, BI1--, NI---}
  //    ZZZZZ ∈ {UL8--, BL4--, NL---}
  //    T ∈ {12|34, 13|24, 14|23}
  //
  // For example, "CFN-BD-BI1H--NL----12|34" is a valid input for 'model'. Yes, the extra dashes are
  // necessary, as YYYYY and ZZZZZ are assumed to consist of 5 characters. There is no error-checking in
  // this function, so make sure your input is correct!
  //
  // COMMENTARY: XX specifies the algebraic invariant type, YYYYY specifies which internal-branch
  // inequalities to use, ZZZZZ specifies which leaf inequalities to use, and T specifies the tree
  // topology of the null hypothesis being tested. See the file compute_invariants.m2 for detailed
  // definitions of these tags.

  // FIRST, INITIALIZE SOME VARIABLES
  std::string topology = getQuartetTopology(model); // possible values: "12|34", "13|24", or "14|23"
  std::string invariant_type = getInvariantType(model); // possible values: "BD", "UR", "BR", "UM", and
                                                        // "BM"
  std::string leaf_type = getLeafType(model); // possible values: "NL---", "UL8--", and "BL4--"
  std::string internal_type = getInternalType(model); // possible values: "UI4--", "UI2H-", "UI2NH",
                                                      // "BI1H-" "BI1NH", "BI1--" and "NI---"
  
  // SECOND, ESTIMATE THE VECTOR OF SITE PATTERN PROBABILITIES
  // Specifically, we'll estimate the vector (paaaa,paaab,paaba,paabb,pabaa,pabab,pabba,pabbb). Here,
  // paaaa is the probability that a single sample has pattern 0000 or 1111; paaab is the probability of
  // 0001 or 1110, and so forth.
  int m = Xm.nrow();  // number of rows of Xm
  int d = Xm.ncol();  // number of columns of Xm. This is always 8 for CFN models.
  int s = 2; // s is maximum degree of the polynomials defining the submodel, which for the CFN model on
             // 4-leaf trees is always 2. 
  int eta = m/s; 
  // Next, construct two (=s) independent estimates p1 and p2 of the vector (paaaa,...,pabbb) using two
  // disjoint subsets of rows of Xm, each of size η:=m/s. For example, p1[0] and p2[0] are both estimates
  // of paaaa, p1[1] and p2[1] are both estimates of paaab, etc. Since p1 and p2 are constructed using
  // disjoint sets of samples, they are independent, which is important for the kernel function to be
  // unbiased.
  NumericVector p1(d); // 8 is the number of columns of Xm
  NumericVector p2(d);
  for(int i = 0; i < d; ++i){
    for(int j = 0; j < eta; ++j){
      p1(i) = p1(i) + Xm(perm(j) - 1, i);
      p2(i) = p2(i) + Xm(perm(eta + j) - 1, i);
    };
    p1(i) = p1(i)/eta;
    p2(i) = p2(i)/eta;
  };
  // COMMENTARY: Note that in the above, we have permuted the rows of Xm using the input permutation
  // 'perm'; this will allow us to achieve symmetrization (either partial or complete) of h, which is
  // constructed by summing over many calls of h_tilde_CFN, using a different permutation each time.
  
  // THIRD, COMPUTE THE ALGEBRAIC INVARIANTS
  double eq1; double eq2; double tmp1; double tmp2;
  if(invariant_type == "BD"){ // case 1: bounded determinental invariants
    if(topology == "12|34"){
      // q1111 - q0011*q1100 = 0
      // q0101*q1010 - q0110*q1001 = 0
      eq1 = p1[0]-p1[1]-p1[2]+p1[3]-p1[4]+p1[5]+p1[6]-p1[7]-
        (p1[0]-p1[1]-p1[2]+p1[3]+p1[4]-p1[5]-p1[6]+p1[7])*(p2[0]+p2[1]+p2[2]+p2[3]-p2[4]-p2[5]-p2[6]-p2[7]); 
      eq2 = (p1[0]+p1[1]-p1[2]-p1[3]+p1[4]+p1[5]-p1[6]-p1[7])*(p2[0]-p2[1]+p2[2]-p2[3]-p2[4]+p2[5]-p2[6]+p2[7])-
        (p1[0]+p1[1]-p1[2]-p1[3]-p1[4]-p1[5]+p1[6]+p1[7])*(p2[0]-p2[1]+p2[2]-p2[3]+p2[4]-p2[5]+p2[6]-p2[7]); 
    }else if(topology == "13|24"){
      // q1111 - q0101*q1010 = 0
      // q0011*q1100 - q0110*q1001 = 0
      eq1 = p1[0]-p1[1] -p1[2]+p1[3]-p1[4]+p1[5]+p1[6]-p1[7]-
        (p1[0]+p1[1]-p1[2]-p1[3]+p1[4]+p1[5]-p1[6]-p1[7])*(p1[0]-p1[1]+p1[2]-p1[3]-p1[4]+p1[5]-p1[6]+p1[7]); 
      eq2 = (p1[0]-p1[1]-p1[2]+p1[3]+p1[4]-p1[5]-p1[6]+p1[7])*(p2[0]+p2[1]+p2[2]+p2[3]-p2[4]-p2[5]-p2[6]-p2[7])-
        (p1[0]+p1[1]-p1[2]-p1[3]-p1[4]-p1[5]+p1[6]+p1[7])*(p2[0]-p2[1]+p2[2]-p2[3]+p2[4]-p2[5]+p2[6]-p2[7]); 
    }else if(topology == "14|23"){
      // q1111 - q0110*q1001 = 0
      // q0011*q1100 - q0101*q1010 = 0
      eq1 = p1[0]-p1[1]-p1[2]+p1[3]-p1[4]+p1[5]+p1[6]-p1[7]-
        (p1[0]+p1[1]-p1[2]-p1[3]-p1[4]-p1[5]+p1[6]+p1[7])*(p2[0]-p2[1]+p2[2]-p2[3]+p2[4]-p2[5]+p2[6]-p2[7]); 
      eq2 = (p1[0]-p1[1]-p1[2]+p1[3]+p1[4]-p1[5]-p1[6]+p1[7])*(p2[0]+p2[1]+p2[2]+p2[3]-p2[4]-p2[5]-p2[6]-p2[7])-
        (p1[0]+p1[1]-p1[2]-p1[3]+p1[4]+p1[5]-p1[6]-p1[7])*(p2[0]-p2[1]+p2[2]-p2[3]-p2[4]+p2[5]-p2[6]+p2[7]); 
    }
  }else if(invariant_type == "UR" || invariant_type == "BR"){ // case 2: rank constraints
    if(     topology == "12|34"){
      tmp1 = p1[2]*p2[4]-p1[3]*p2[5]-p1[0]*p2[6]+p1[1]*p2[7];
      tmp2 = p1[1]*p2[4]-p1[0]*p2[5]-p1[3]*p2[6]+p1[2]*p2[7];
    }else if(topology == "13|24"){
      tmp1 = p1[2]*p2[4]-p1[3]*p2[5]-p1[0]*p2[6]+p1[1]*p2[7];
      tmp2 = p1[1]*p2[2]-p1[0]*p2[3]-p1[5]*p2[6]+p1[4]*p2[7];
    }else if(topology == "14|23"){
      tmp1 = p1[1]*p2[4]-p1[0]*p2[5]-p1[3]*p2[6]+p1[2]*p2[7];
      tmp2 = p1[1]*p2[2]-p1[0]*p2[3]-p1[5]*p2[6]+p1[4]*p2[7];}
    if(invariant_type == "UR"){eq1 = tmp1; eq2 = tmp2;}
    if(invariant_type == "BR"){eq1 = tmp1 + tmp2; eq2 = tmp1 - tmp2;}
  }else if(invariant_type == "UM" || invariant_type == "BM"){ // case 3: mingens 
    if(topology == "12|34"){
      // p0010*p0100-p0011*p0101+p0001*p0110+p0010*p0110+p0011*p0110+p0100*p0110+p0101*p0110+p0110^2+p0001*p0111+p0110*p0111-p0110 = 0
      // p0001*p0100+p0001*p0101+p0010*p0101+p0011*p0101+p0100*p0101+p0101^2-p0011*p0110+p0101*p0110+p0010*p0111+p0101*p0111-p0101 = 0 
      tmp1 = p1[2]*p2[4]-p1[3]*p2[5]+p1[1]*p2[6]+p1[2]*p2[6]+p1[3]*p2[6]+p1[4]*p2[6]+p1[5]*p2[6]+p1[6]*p2[6]+p1[1]*p2[7]+p1[6]*p2[7]-p1[6];
      tmp2 = p1[1]*p2[4]+p1[1]*p2[5]+p1[2]*p2[5]+p1[3]*p2[5]+p1[4]*p2[5]+p1[5]*p2[5]-p1[3]*p2[6]+p1[5]*p2[6]+p1[2]*p2[7]+p1[5]*p2[7]-p1[5];
    }else if(topology == "13|24"){
      // p0010*p0100-p0011*p0101+p0001*p0110+p0010*p0110+p0011*p0110+p0100*p0110+p0101*p0110+p0110^2+p0001*p0111+p0110*p0111-p0110 = 0
      // p0001*p0010+p0001*p0011+p0010*p0011+p0011^2+p0011*p0100+p0011*p0101+p0011*p0110-p0101*p0110+p0011*p0111+p0100*p0111-p0011 = 0
      tmp1 = p1[2]*p2[4]-p1[3]*p2[5]+p1[1]*p2[6]+p1[2]*p2[6]+p1[3]*p2[6]+p1[4]*p2[6]+p1[5]*p2[6]+p1[6]*p2[6]+p1[1]*p2[7]+p1[6]*p2[7]-p1[6];
      tmp2 = p1[1]*p2[2]+p1[1]*p2[3]+p1[2]*p2[3]+p1[3]*p2[3]+p1[3]*p2[4]+p1[3]*p2[5]+p1[3]*p2[6]-p1[5]*p2[6]+p1[3]*p2[7]+p1[4]*p2[7]-p1[3];
    }else if(topology == "14|23"){
      // p0001*p0100+p0001*p0101+p0010*p0101+p0011*p0101+p0100*p0101+p0101^2-p0011*p0110+p0101*p0110+p0010*p0111+p0101*p0111-p0101 = 0
      // p0001*p0010+p0001*p0011+p0010*p0011+p0011^2+p0011*p0100+p0011*p0101+p0011*p0110-p0101*p0110+p0011*p0111+p0100*p0111-p0011 = 0
      tmp1 = p1[1]*p2[4]+p1[1]*p2[5]+p1[2]*p2[5]+p1[3]*p2[5]+p1[4]*p2[5]+p1[5]*p2[5]-p1[3]*p2[6]+p1[5]*p2[6]+p1[2]*p2[7]+p1[5]*p2[7]-p1[5];
      tmp2 = p1[1]*p2[2]+p1[1]*p2[3]+p1[2]*p2[3]+p1[3]*p2[3]+p1[3]*p2[4]+p1[3]*p2[5]+p1[3]*p2[6]-p1[5]*p2[6]+p1[3]*p2[7]+p1[4]*p2[7]-p1[3];
    }
    if(invariant_type == "UM"){eq1 = tmp1; eq2 = tmp2;}
    else if(invariant_type == "BM"){eq1 = tmp1 + tmp2; eq2 = tmp1 - tmp2;}
  }
  
  // FOURTH, COMPUTE THE INTERNAL BRANCH INEQUALITIES (IF ANY)  [this part could be implemented more efficiently]
  double tmp_internal_ineq1; double tmp_internal_ineq2; double tmp_internal_ineq3; double tmp_internal_ineq4;
  double internal_ineq1; double internal_ineq2; double internal_ineq3; double internal_ineq4;
  if(internal_type != "NI---"){
    if(topology == "12|34"){
      // q1010*q0101 - q1100*q0011 ≤ 0
      // q0110*q1001 - q1100*q0011 ≤ 0
      // q1010*q0101 - q1111 ≤ 0
      // q0110*q1001 - q1111 ≤ 0
      tmp_internal_ineq1 = p1[1]*p2[2]-p1[0]*p2[3]-p1[1]*p2[4]+p1[0]*p2[5]+p1[3]*p2[6]-p1[5]*p2[6]-p1[2]*p2[7]+p1[4]*p2[7];
      tmp_internal_ineq2 = p1[1]*p2[2]-p1[0]*p2[3]-p1[2]*p2[4]+p1[3]*p2[5]+p1[0]*p2[6]-p1[5]*p2[6]-p1[1]*p2[7]+p1[4]*p2[7]; 
      tmp_internal_ineq3 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[1]*p2[2]-p1[2]*p2[2]-2*p1[0]*p2[3]+p1[3]*p2[3]-2*p1[1]*p2[4]+       
        2*p1[2]*p2[4]-p1[4]*p2[4]+2*p1[0]*p2[5]-2*p1[3]*p2[5]+ p1[5]*p2[5]-2*p1[0]*p2[6]+2*p1[3]*p2[6]-2*p1[5]*p2[6]+
        p1[6]*p2[6]+2*p1[1]*p2[7]-2*p1[2]*p2[7]+2*p1[4]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]+p1[4]-p1[5]-p1[6]+p1[7];   
      tmp_internal_ineq4 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[1]*p2[2]-p1[2]*p2[2]-2*p1[0]*p2[3]+p1[3]*p2[3]+2*p1[1]*p2[4]-
        2*p1[2]*p2[4]-p1[4]*p2[4]-2*p1[0]*p2[5]+2*p1[3]*p2[5]+p1[5]*p2[5]+2*p1[0]*p2[6]-2*p1[3]*p2[6]-2*p1[5]*p2[6]+
        p1[6]*p2[6]-2*p1[1]*p2[7]+2*p1[2]*p2[7]+2*p1[4]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]+p1[4]-p1[5]-p1[6]+p1[7];   
    }else if(topology == "13|24"){
      // q0011*q1100 - q0101*q1010 ≤ 0
      // q0110*q1001 - q0101*q1010 ≤ 0
      // q0011*q1100 - q1111 ≤ 0
      // q0110*q1001 - q1111 ≤ 0
      tmp_internal_ineq1 = -p1[1]*p2[2]+p1[0]*p2[3]+p1[1]*p2[4]-p1[0]*p2[5]-p1[3]*p2[6]+p1[5]*p2[6]+p1[2]*p2[7]-p1[4]*p2[7];
      tmp_internal_ineq2 =  p1[1]*p2[4]-p1[2]*p2[4]-p1[0]*p2[5]+p1[3]*p2[5]+p1[0]*p2[6]-p1[3]*p2[6]-p1[1]*p2[7]+p1[2]*p2[7];
      tmp_internal_ineq3 = p1[0]*p2[0]-p1[1]*p2[1]-2*p1[1]*p2[2]-p1[2]*p2[2]+2*p1[0]*p2[3]+p1[3]*p2[3]+2*p1[1]*p2[4]+
        2*p1[2]*p2[4]-p1[4]*p2[4]-2*p1[0]*p2[5]-2*p1[3]*p2[5]+p1[5]*p2[5]-2*p1[0]*p2[6]-2*p1[3]*p2[6]+2*p1[5]*p2[6]+
        p1[6]*p2[6]+2*p1[1]*p2[7]+2*p1[2]*p2[7]-2*p1[4]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]+p1[4]-p1[5]-p1[6]+p1[7];   
      tmp_internal_ineq4 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[1]*p2[2]-p1[2]*p2[2]-2*p1[0]*p2[3]+p1[3]*p2[3]+2*p1[1]*p2[4]-
        2*p1[2]*p2[4]-p1[4]*p2[4]-2*p1[0]*p2[5]+2*p1[3]*p2[5]+p1[5]*p2[5]+2*p1[0]*p2[6]-2*p1[3]*p2[6]-2*p1[5]*p2[6]+
        p1[6]*p2[6]-2*p1[1]*p2[7]+2*p1[2]*p2[7]+2*p1[4]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]+p1[4]-p1[5]-p1[6]+p1[7];   
    }else if(topology == "14|23"){
      // q0011*q1100 - q0110*q1001 ≤ 0
      // q0101*q1010 - q0110*q1001 ≤ 0
      // q0011*q1100 - q1111 ≤ 0
      // q0101*q1010 - q1111 ≤ 0
      tmp_internal_ineq1 = -p1[1]*p2[2]+p1[0]*p2[3]+p1[2]*p2[4]-p1[3]*p2[5]-p1[0]*p2[6]+p1[5]*p2[6]+p1[1]*p2[7]-p1[4]*p2[7]; 
      tmp_internal_ineq2 = -p1[1]*p2[4]+p1[2]*p2[4]+p1[0]*p2[5]-p1[3]*p2[5]-p1[0]*p2[6]+p1[3]*p2[6]+p1[1]*p2[7]-p1[2]*p2[7]; 
      tmp_internal_ineq3 = p1[0]*p2[0]-p1[1]*p2[1]-2*p1[1]*p2[2]-p1[2]*p2[2]+2*p1[0]*p2[3]+p1[3]*p2[3]+2*p1[1]*p2[4]+
        2*p1[2]*p2[4]-p1[4]*p2[4]-2*p1[0]*p2[5]-2*p1[3]*p2[5]+p1[5]*p2[5]-2*p1[0]*p2[6]-2*p1[3]*p2[6]+2*p1[5]*p2[6]+
        p1[6]*p2[6]+2*p1[1]*p2[7]+2*p1[2]*p2[7]-2*p1[4]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]+p1[4]-p1[5]-p1[6]+p1[7];   
      tmp_internal_ineq4 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[1]*p2[2]-p1[2]*p2[2]-2*p1[0]*p2[3]+p1[3]*p2[3]-2*p1[1]*p2[4]+
        2*p1[2]*p2[4]-p1[4]*p2[4]+2*p1[0]*p2[5]-2*p1[3]*p2[5]+p1[5]*p2[5]-2*p1[0]*p2[6]+2*p1[3]*p2[6]-2*p1[5]*p2[6]+
        p1[6]*p2[6]+2*p1[1]*p2[7]-2*p1[2]*p2[7]+2*p1[4]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]+p1[4]-p1[5]-p1[6]+p1[7];   
    }
    if(internal_type == "UI4--"){ // all 4 inequalities
      internal_ineq1 = tmp_internal_ineq1;
      internal_ineq2 = tmp_internal_ineq2;
      internal_ineq3 = tmp_internal_ineq3;
      internal_ineq4 = tmp_internal_ineq4;
    }else if(internal_type == "UI2H-"){ // unbalance homogeneous
      internal_ineq1 = tmp_internal_ineq1;
      internal_ineq2 = tmp_internal_ineq2;
    }else if(internal_type == "UI2NH"){ // unbalanced non-homogeneous
      internal_ineq1 = tmp_internal_ineq3;
      internal_ineq2 = tmp_internal_ineq4;
    }else if(internal_type == "BI1H-"){ // balanced single homogeneous
      internal_ineq1 = tmp_internal_ineq1 + tmp_internal_ineq2;
    }else if(internal_type == "BI1NH"){ // balanced single non-homogeneous
      internal_ineq1 = tmp_internal_ineq3 + tmp_internal_ineq4;
    }else if(internal_type == "BI1--"){ // balanced, using all 4
      internal_ineq1 =  tmp_internal_ineq1 + tmp_internal_ineq2 + tmp_internal_ineq3 + tmp_internal_ineq4;
    }
  }
  
  // FIFTH, COMPUTE THE LEAF INEQUALITIES (IF ANY)
  double leaf_ineq1; double leaf_ineq2; double leaf_ineq3; double leaf_ineq4; double leaf_ineq5; double leaf_ineq6; double leaf_ineq7; double leaf_ineq8;
  if(leaf_type == "UL8--" || leaf_type == "BL4--"){
    if(topology == "12|34"){
      leaf_ineq1 = p1[0]*p2[0]+2*p1[0]*p2[1]+p1[1]*p2[1]-p1[2]*p2[2]-2*p1[2]*p2[3]-p1[3]*p2[3]+2*p1[2]*p2[4]+2*p1[3]*p2[4]-p1[4]*p2[4]+2*p1[2]*p2[5]+2*p1[3]*p2[5]-
        2*p1[4]*p2[5]-p1[5]*p2[5]-2*p1[0]*p2[6]-2*p1[1]*p2[6]+p1[6]*p2[6]-2*p1[0]*p2[7]-2*p1[1]*p2[7]+2*p1[6]*p2[7]+p1[7]*p2[7]-p1[0]-p1[1]+p1[2]+p1[3]+p1[4]+p1[5]-p1[6]-p1[7]; // -q0110 + q1010*q1100 ≤ 0
      leaf_ineq2 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[0]*p2[2]+p1[2]*p2[2]-2*p1[1]*p2[3]-p1[3]*p2[3]+2*p1[1]*p2[4]+2*p1[3]*p2[4]-p1[4]*p2[4]-2*p1[0]*p2[5]-2*p1[2]*p2[5]+
        p1[5]*p2[5]+2*p1[1]*p2[6]+2*p1[3]*p2[6]-2*p1[4]*p2[6]-p1[6]*p2[6]-2*p1[0]*p2[7]-2*p1[2]*p2[7]+2*p1[5]*p2[7]+p1[7]*p2[7]-p1[0]+p1[1]-p1[2]+p1[3]+p1[4]-p1[5]+p1[6]-p1[7]; // -q0101 + q1001*q1100 ≤ 0
      leaf_ineq3 = p1[0]*p2[0]+2*p1[0]*p2[1]+p1[1]*p2[1]-p1[2]*p2[2]-2*p1[2]*p2[3]-p1[3]*p2[3]-2*p1[0]*p2[4]-2*p1[1]*p2[4]+p1[4]*p2[4]-2*p1[0]*p2[5]-2*p1[1]*p2[5]+
        2*p1[4]*p2[5]+p1[5]*p2[5]+2*p1[2]*p2[6]+2*p1[3]*p2[6]-p1[6]*p2[6]+2*p1[2]*p2[7]+2*p1[3]*p2[7]-2*p1[6]*p2[7]-p1[7]*p2[7]-p1[0]-p1[1]+p1[2]+p1[3]-p1[4]-p1[5]+p1[6]+p1[7]; // -q1010 + q0110*q1100 ≤ 0
      leaf_ineq4 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[0]*p2[2]+p1[2]*p2[2]-2*p1[1]*p2[3]-p1[3]*p2[3]-2*p1[0]*p2[4]-2*p1[2]*p2[4]+p1[4]*p2[4]+2*p1[1]*p2[5]+2*p1[3]*p2[5]-
        p1[5]*p2[5]-2*p1[0]*p2[6]-2*p1[2]*p2[6]+2*p1[4]*p2[6]+p1[6]*p2[6]+2*p1[1]*p2[7]+2*p1[3]*p2[7]-2*p1[5]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]-p1[2]+p1[3]-p1[4]+p1[5]-p1[6]+p1[7]; // -q1001 + q0101*q1100 ≤ 0
      leaf_ineq5 = p1[0]*p2[0]-p1[1]*p2[1]-2*p1[0]*p2[2]+p1[2]*p2[2]+2*p1[1]*p2[3]-p1[3]*p2[3]+2*p1[0]*p2[4]-2*p1[2]*p2[4]+p1[4]*p2[4]-2*p1[1]*p2[5]+2*p1[3]*p2[5]-
        p1[5]*p2[5]-2*p1[0]*p2[6]+2*p1[2]*p2[6]-2*p1[4]*p2[6]+p1[6]*p2[6]+2*p1[1]*p2[7]-2*p1[3]*p2[7]+2*p1[5]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]-p1[2]+p1[3]-p1[4]+p1[5]-p1[6]+p1[7]; // -q1001 + q0011*q1010 ≤ 0
      leaf_ineq6 = p1[0]*p2[0]-p1[1]*p2[1]-2*p1[0]*p2[2]+p1[2]*p2[2]+2*p1[1]*p2[3]-p1[3]*p2[3]+2*p1[1]*p2[4]-2*p1[3]*p2[4]-p1[4]*p2[4]-2*p1[0]*p2[5]+2*p1[2]*p2[5]+
        p1[5]*p2[5]-2*p1[1]*p2[6]+2*p1[3]*p2[6]+2*p1[4]*p2[6]-p1[6]*p2[6]+2*p1[0]*p2[7]-2*p1[2]*p2[7]-2*p1[5]*p2[7]+p1[7]*p2[7]-p1[0]+p1[1]-p1[2]+p1[3]+p1[4]-p1[5]+p1[6]-p1[7]; // -q0101 + q0110*q0011 ≤ 0
      leaf_ineq7 = p1[0]*p2[0]-2*p1[0]*p2[1]+p1[1]*p2[1]-p1[2]*p2[2]+2*p1[2]*p2[3]-p1[3]*p2[3]+2*p1[0]*p2[4]-2*p1[1]*p2[4]+p1[4]*p2[4]-2*p1[0]*p2[5]+2*p1[1]*p2[5]-
        2*p1[4]*p2[5]+p1[5]*p2[5]-2*p1[2]*p2[6]+2*p1[3]*p2[6]-p1[6]*p2[6]+2*p1[2]*p2[7]-2*p1[3]*p2[7]+2*p1[6]*p2[7]-p1[7]*p2[7]-p1[0]-p1[1]+p1[2]+p1[3]-p1[4]-p1[5]+p1[6]+p1[7]; // -q1010 + q0011*q1001 ≤ 0
      leaf_ineq8 = p1[0]*p2[0]-2*p1[0]*p2[1]+p1[1]*p2[1]-p1[2]*p2[2]+2*p1[2]*p2[3]-p1[3]*p2[3]+2*p1[2]*p2[4]-2*p1[3]*p2[4]-p1[4]*p2[4]-2*p1[2]*p2[5]+2*p1[3]*p2[5]+
        2*p1[4]*p2[5]-p1[5]*p2[5]-2*p1[0]*p2[6]+2*p1[1]*p2[6]+p1[6]*p2[6]+2*p1[0]*p2[7]-2*p1[1]*p2[7]-2*p1[6]*p2[7]+p1[7]*p2[7]-p1[0]-p1[1]+p1[2]+p1[3]+p1[4]+p1[5]-p1[6]-p1[7]; // -q0110 + q0101*q0011 ≤ 0
    }else if(topology == "13|24"){
      leaf_ineq1 = p1[0]*p2[0]+2*p1[0]*p2[1]+p1[1]*p2[1]-p1[2]*p2[2]-2*p1[2]*p2[3]-p1[3]*p2[3]+2*p1[2]*p2[4]+2*p1[3]*p2[4]-p1[4]*p2[4]+2*p1[2]*p2[5]+2*p1[3]*p2[5]-
        2*p1[4]*p2[5]-p1[5]*p2[5]-2*p1[0]*p2[6]-2*p1[1]*p2[6]+p1[6]*p2[6]-2*p1[0]*p2[7]-2*p1[1]*p2[7]+2*p1[6]*p2[7]+p1[7]*p2[7]-p1[0]-p1[1]+p1[2]+p1[3]+p1[4]+p1[5]-p1[6]-p1[7]; // -q0110 + q1010*q1100 ≤ 0
      leaf_ineq2 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[1]*p2[2]-p1[2]*p2[2]-2*p1[0]*p2[3]+p1[3]*p2[3]+2*p1[0]*p2[4]-2*p1[3]*p2[4]+p1[4]*p2[4]-2*p1[1]*p2[5]+2*p1[2]*p2[5]-
        p1[5]*p2[5]+2*p1[1]*p2[6]-2*p1[2]*p2[6]+2*p1[5]*p2[6]-p1[6]*p2[6]-2*p1[0]*p2[7]+2*p1[3]*p2[7]-2*p1[4]*p2[7]+p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]-p1[4]+p1[5]+p1[6]-p1[7]; // -q0011 + q1001*q1010 ≤ 0
      leaf_ineq3 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[0]*p2[2]+p1[2]*p2[2]-2*p1[1]*p2[3]-p1[3]*p2[3]-2*p1[0]*p2[4]-2*p1[2]*p2[4]+p1[4]*p2[4]+2*p1[1]*p2[5]+2*p1[3]*p2[5]-
        p1[5]*p2[5]-2*p1[0]*p2[6]-2*p1[2]*p2[6]+2*p1[4]*p2[6]+p1[6]*p2[6]+2*p1[1]*p2[7]+2*p1[3]*p2[7]-2*p1[5]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]-p1[2]+p1[3]-p1[4]+p1[5]-p1[6]+p1[7]; // -q1001 + q0101*q1100 ≤ 0
      leaf_ineq4 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[1]*p2[2]-p1[2]*p2[2]-2*p1[0]*p2[3]+p1[3]*p2[3]-2*p1[0]*p2[4]+2*p1[3]*p2[4]+p1[4]*p2[4]+2*p1[1]*p2[5]-2*p1[2]*p2[5]-
        p1[5]*p2[5]-2*p1[1]*p2[6]+2*p1[2]*p2[6]+2*p1[5]*p2[6]-p1[6]*p2[6]+2*p1[0]*p2[7]-2*p1[3]*p2[7]-2*p1[4]*p2[7]+p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]-p1[4]+p1[5]+p1[6]-p1[7]; // -q0011 + q0101*q0110 ≤ 0
      leaf_ineq5 = p1[0]*p2[0]+2*p1[0]*p2[1]+p1[1]*p2[1]-2*p1[0]*p2[2]-2*p1[1]*p2[2]+p1[2]*p2[2]-2*p1[0]*p2[3]-2*p1[1]*p2[3]+2*p1[2]*p2[3]+p1[3]*p2[3]-p1[4]*p2[4]-
        2*p1[4]*p2[5]-p1[5]*p2[5]+2*p1[4]*p2[6]+2*p1[5]*p2[6]-p1[6]*p2[6]+2*p1[4]*p2[7]+2*p1[5]*p2[7]-2*p1[6]*p2[7]-p1[7]*p2[7]-p1[0]-p1[1]-p1[2]-p1[3]+p1[4]+p1[5]+p1[6]+p1[7]; // -q1100 + q0110*q1010 ≤ 0
      leaf_ineq6 = p1[0]*p2[0]-p1[1]*p2[1]-2*p1[0]*p2[2]+p1[2]*p2[2]+2*p1[1]*p2[3]-p1[3]*p2[3]+2*p1[0]*p2[4]-2*p1[2]*p2[4]+p1[4]*p2[4]-2*p1[1]*p2[5]+2*p1[3]*p2[5]-
        p1[5]*p2[5]-2*p1[0]*p2[6]+2*p1[2]*p2[6]-2*p1[4]*p2[6]+p1[6]*p2[6]+2*p1[1]*p2[7]-2*p1[3]*p2[7]+2*p1[5]*p2[7]-p1[7]*p2[7]-p1[0]+p1[1]-p1[2]+p1[3]-p1[4]+p1[5]-p1[6]+p1[7]; // -q1001 + q0011*q1010 ≤ 0
      leaf_ineq7 = p1[0]*p2[0]-2*p1[0]*p2[1]+p1[1]*p2[1]+2*p1[0]*p2[2]-2*p1[1]*p2[2]+p1[2]*p2[2]-2*p1[0]*p2[3]+2*p1[1]*p2[3]-2*p1[2]*p2[3]+p1[3]*p2[3]-p1[4]*p2[4]+
        2*p1[4]*p2[5]-p1[5]*p2[5]-2*p1[4]*p2[6]+2*p1[5]*p2[6]-p1[6]*p2[6]+2*p1[4]*p2[7]-2*p1[5]*p2[7]+2*p1[6]*p2[7]-p1[7]*p2[7]-p1[0]-p1[1]-p1[2]-p1[3]+p1[4]+p1[5]+p1[6]+p1[7]; // -q1100 + q0101*q1001 ≤ 0
      leaf_ineq8 = p1[0]*p2[0]-2*p1[0]*p2[1]+p1[1]*p2[1]-p1[2]*p2[2]+2*p1[2]*p2[3]-p1[3]*p2[3]+2*p1[2]*p2[4]-2*p1[3]*p2[4]-p1[4]*p2[4]-2*p1[2]*p2[5]+2*p1[3]*p2[5]+
        2*p1[4]*p2[5]-p1[5]*p2[5]-2*p1[0]*p2[6]+2*p1[1]*p2[6]+p1[6]*p2[6]+2*p1[0]*p2[7]-2*p1[1]*p2[7]-2*p1[6]*p2[7]+p1[7]*p2[7]-p1[0]-p1[1]+p1[2]+p1[3]+p1[4]+p1[5]-p1[6]-p1[7]; // -q0110 + q0101*q0011 ≤ 0
    }else if(topology == "14|23"){
      leaf_ineq1 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[0]*p2[2]+p1[2]*p2[2]-2*p1[1]*p2[3]-p1[3]*p2[3]+2*p1[1]*p2[4]+2*p1[3]*p2[4]-p1[4]*p2[4]-2*p1[0]*p2[5]-2*p1[2]*p2[5]+
        p1[5]*p2[5]+2*p1[1]*p2[6]+2*p1[3]*p2[6]-2*p1[4]*p2[6]-p1[6]*p2[6]-2*p1[0]*p2[7]-2*p1[2]*p2[7]+2*p1[5]*p2[7]+p1[7]*p2[7]-p1[0]+p1[1]-p1[2]+p1[3]+p1[4]-p1[5]+p1[6]-p1[7]; // -q0101 + q1001*q1100 ≤ 0
      leaf_ineq2 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[1]*p2[2]-p1[2]*p2[2]-2*p1[0]*p2[3]+p1[3]*p2[3]+2*p1[0]*p2[4]-2*p1[3]*p2[4]+p1[4]*p2[4]-2*p1[1]*p2[5]+2*p1[2]*p2[5]-
        p1[5]*p2[5]+2*p1[1]*p2[6]-2*p1[2]*p2[6]+2*p1[5]*p2[6]-p1[6]*p2[6]-2*p1[0]*p2[7]+2*p1[3]*p2[7]-2*p1[4]*p2[7]+p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]-p1[4]+p1[5]+p1[6]-p1[7]; // -q0011 + q1001*q1010 ≤ 0
      leaf_ineq3 = p1[0]*p2[0]+2*p1[0]*p2[1]+p1[1]*p2[1]-p1[2]*p2[2]-2*p1[2]*p2[3]-p1[3]*p2[3]-2*p1[0]*p2[4]-2*p1[1]*p2[4]+p1[4]*p2[4]-2*p1[0]*p2[5]-2*p1[1]*p2[5]+
        2*p1[4]*p2[5]+p1[5]*p2[5]+2*p1[2]*p2[6]+2*p1[3]*p2[6]-p1[6]*p2[6]+2*p1[2]*p2[7]+2*p1[3]*p2[7]-2*p1[6]*p2[7]-p1[7]*p2[7]-p1[0]-p1[1]+p1[2]+p1[3]-p1[4]-p1[5]+p1[6]+p1[7]; // -q1010 + q0110*q1100 ≤ 0
      leaf_ineq4 = p1[0]*p2[0]-p1[1]*p2[1]+2*p1[1]*p2[2]-p1[2]*p2[2]-2*p1[0]*p2[3]+p1[3]*p2[3]-2*p1[0]*p2[4]+2*p1[3]*p2[4]+p1[4]*p2[4]+2*p1[1]*p2[5]-2*p1[2]*p2[5]-
        p1[5]*p2[5]-2*p1[1]*p2[6]+2*p1[2]*p2[6]+2*p1[5]*p2[6]-p1[6]*p2[6]+2*p1[0]*p2[7]-2*p1[3]*p2[7]-2*p1[4]*p2[7]+p1[7]*p2[7]-p1[0]+p1[1]+p1[2]-p1[3]-p1[4]+p1[5]+p1[6]-p1[7]; // -q0011 + q0101*q0110 ≤ 0
      leaf_ineq5 = p1[0]*p2[0]+2*p1[0]*p2[1]+p1[1]*p2[1]-2*p1[0]*p2[2]-2*p1[1]*p2[2]+p1[2]*p2[2]-2*p1[0]*p2[3]-2*p1[1]*p2[3]+2*p1[2]*p2[3]+p1[3]*p2[3]-p1[4]*p2[4]-
        2*p1[4]*p2[5]-p1[5]*p2[5]+2*p1[4]*p2[6]+2*p1[5]*p2[6]-p1[6]*p2[6]+2*p1[4]*p2[7]+2*p1[5]*p2[7]-2*p1[6]*p2[7]-p1[7]*p2[7]-p1[0]-p1[1]-p1[2]-p1[3]+p1[4]+p1[5]+p1[6]+p1[7]; // -q1100 + q0110*q1010 ≤ 0
      leaf_ineq6 = p1[0]*p2[0]-p1[1]*p2[1]-2*p1[0]*p2[2]+p1[2]*p2[2]+2*p1[1]*p2[3]-p1[3]*p2[3]+2*p1[1]*p2[4]-2*p1[3]*p2[4]-p1[4]*p2[4]-2*p1[0]*p2[5]+2*p1[2]*p2[5]+
        p1[5]*p2[5]-2*p1[1]*p2[6]+2*p1[3]*p2[6]+2*p1[4]*p2[6]-p1[6]*p2[6]+2*p1[0]*p2[7]-2*p1[2]*p2[7]-2*p1[5]*p2[7]+p1[7]*p2[7]-p1[0]+p1[1]-p1[2]+p1[3]+p1[4]-p1[5]+p1[6]-p1[7]; // -q0101 + q0110*q0011 ≤ 0
      leaf_ineq7 = p1[0]*p2[0]-2*p1[0]*p2[1]+p1[1]*p2[1]+2*p1[0]*p2[2]-2*p1[1]*p2[2]+p1[2]*p2[2]-2*p1[0]*p2[3]+2*p1[1]*p2[3]-2*p1[2]*p2[3]+p1[3]*p2[3]-p1[4]*p2[4]+
        2*p1[4]*p2[5]-p1[5]*p2[5]-2*p1[4]*p2[6]+2*p1[5]*p2[6]-p1[6]*p2[6]+2*p1[4]*p2[7]-2*p1[5]*p2[7]+2*p1[6]*p2[7]-p1[7]*p2[7]-p1[0]-p1[1]-p1[2]-p1[3]+p1[4]+p1[5]+p1[6]+p1[7]; // -q1100 + q0101*q1001 ≤ 0
      leaf_ineq8 = p1[0]*p2[0]-2*p1[0]*p2[1]+p1[1]*p2[1]-p1[2]*p2[2]+2*p1[2]*p2[3]-p1[3]*p2[3]+2*p1[0]*p2[4]-2*p1[1]*p2[4]+p1[4]*p2[4]-2*p1[0]*p2[5]+2*p1[1]*p2[5]-
        2*p1[4]*p2[5]+p1[5]*p2[5]-2*p1[2]*p2[6]+2*p1[3]*p2[6]-p1[6]*p2[6]+2*p1[2]*p2[7]-2*p1[3]*p2[7]+2*p1[6]*p2[7]-p1[7]*p2[7]-p1[0]-p1[1]+p1[2]+p1[3]-p1[4]-p1[5]+p1[6]+p1[7]; // -q1010 + q0011*q1001 ≤ 0
    }
  }

  // FINALLY, DEFINE THE OUTPUT VECTOR, WHICH DEPENDS ON WHICH CASE WE ARE IN
  if(leaf_type == "NL---"){
    if(internal_type == "UI4--"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,internal_ineq2,internal_ineq3,internal_ineq4);}
    if(internal_type == "UI2H-"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,internal_ineq2);}
    if(internal_type == "UI2NH"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,internal_ineq2);}
    if(internal_type == "BI1H-"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1);}
    if(internal_type == "BI1NH"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1);}
    if(internal_type == "BI1--"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1);}
    if(internal_type == "NI---"){ return NumericVector::create(eq1,-eq1,eq2,-eq2);}
  }
  if(leaf_type == "UL8--"){
    if(internal_type == "UI4--"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,internal_ineq2,internal_ineq3,internal_ineq4,leaf_ineq1,leaf_ineq2,leaf_ineq3,leaf_ineq4,leaf_ineq5,leaf_ineq6,leaf_ineq7,leaf_ineq8);}
    if(internal_type == "UI2H-"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,internal_ineq2,leaf_ineq1,leaf_ineq2,leaf_ineq3,leaf_ineq4,leaf_ineq5,leaf_ineq6,leaf_ineq7,leaf_ineq8);}
    if(internal_type == "UI2NH"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,internal_ineq2,leaf_ineq1,leaf_ineq2,leaf_ineq3,leaf_ineq4,leaf_ineq5,leaf_ineq6,leaf_ineq7,leaf_ineq8);}
    if(internal_type == "BI1H-"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,leaf_ineq1,leaf_ineq2,leaf_ineq3,leaf_ineq4,leaf_ineq5,leaf_ineq6,leaf_ineq7,leaf_ineq8);}
    if(internal_type == "BI1NH"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,leaf_ineq1,leaf_ineq2,leaf_ineq3,leaf_ineq4,leaf_ineq5,leaf_ineq6,leaf_ineq7,leaf_ineq8);}
    if(internal_type == "BI1--"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,leaf_ineq1,leaf_ineq2,leaf_ineq3,leaf_ineq4,leaf_ineq5,leaf_ineq6,leaf_ineq7,leaf_ineq8);}
    if(internal_type == "NI---"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,leaf_ineq1,leaf_ineq2,leaf_ineq3,leaf_ineq4,leaf_ineq5,leaf_ineq6,leaf_ineq7,leaf_ineq8);}
  }
  if(leaf_type == "BL4--"){
    if(internal_type == "UI4--"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,internal_ineq2,internal_ineq3,internal_ineq4,leaf_ineq1+leaf_ineq2,leaf_ineq3+leaf_ineq4,leaf_ineq5+leaf_ineq6,leaf_ineq7+leaf_ineq8);}
    if(internal_type == "UI2H-"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,internal_ineq2,leaf_ineq1+leaf_ineq2,leaf_ineq3+leaf_ineq4,leaf_ineq5+leaf_ineq6,leaf_ineq7+leaf_ineq8);}
    if(internal_type == "UI2NH"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,internal_ineq2,leaf_ineq1+leaf_ineq2,leaf_ineq3+leaf_ineq4,leaf_ineq5+leaf_ineq6,leaf_ineq7+leaf_ineq8);}
    if(internal_type == "BI1H-"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,leaf_ineq1+leaf_ineq2,leaf_ineq3+leaf_ineq4,leaf_ineq5+leaf_ineq6,leaf_ineq7+leaf_ineq8);}
    if(internal_type == "BI1NH"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,leaf_ineq1+leaf_ineq2,leaf_ineq3+leaf_ineq4,leaf_ineq5+leaf_ineq6,leaf_ineq7+leaf_ineq8);}
    if(internal_type == "BI1--"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,internal_ineq1,leaf_ineq1+leaf_ineq2,leaf_ineq3+leaf_ineq4,leaf_ineq5+leaf_ineq6,leaf_ineq7+leaf_ineq8);}
    if(internal_type == "NI---"){ return NumericVector::create(eq1,-eq1,eq2,-eq2,leaf_ineq1+leaf_ineq2,leaf_ineq3+leaf_ineq4,leaf_ineq5+leaf_ineq6,leaf_ineq7+leaf_ineq8);}
  }
  
  return NULL;
}


//////////// Main h_tilde function //////////// 
// [[Rcpp::export]]
NumericVector h_tilde(std::string model, NumericMatrix Xm, IntegerVector perm){
  // this function determines which h_tilde_MODEL function to use
  if(model == "model-1"){
    return h_tilde_cut1(Xm);
  }else if(model == "model-2"){
    return h_tilde_T1(Xm);
  }else if(model == "model-3"){
    return h_tilde_cut(model, Xm, perm);
  }else if(model == "model-4"){
    return h_tilde_T3(model, Xm, perm);
  }else if(model == "HW"){
    return h_tilde_HW(model, Xm, perm);
  }else if(is_a_CFN_model(model)){
    return h_tilde_CFN(model, Xm, perm);
  }else{
    Rcout << "Unknown model: " << model << std::endl;
    Rcout << " -> Accepted models: cut1, T1, cur, T3, CFN, HW." << std::endl;
    return(NULL);
  }
}


//////////////////////////////////////////////////////////////////////// 
/////////// 3. FUNCTIONS TO COMPUTE INCOMPLETE U-STATISTIC /////////////
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
/////////// 3A. FUNCTIONS TO IMPLEMENT CONVEX COMBINATIONS OF CONSTAINTS /////////// 
////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
NumericVector convex_combination(NumericVector h, NumericMatrix cc_coefficients) {
  int n = h.size();
  int n_cc = cc_coefficients.nrow();
  NumericVector convex_combinations(n_cc);
  // Generate random convex combinations and calculate convex combinations of constraints
  for (int j = 0; j < n_cc; j++) {
    double conv = 0;
    for (int i = 0; i < n; i++) conv += cc_coefficients(j,i)*h(i);
    convex_combinations(j) = conv;
  }
  return(convex_combinations);
}

// [[Rcpp::export]]
NumericVector dirichlet_sample(int n, NumericVector alpha) {
  // Auxilliary function to 'convex_combination_coefficients'. This function draws samples from a
  // Dirichlet distribution (this is used to draw convex combinations uniformly at random)
  //
  // Input:
  //   - n: Number of samples to generate.
  //   - alpha: a length K vector of shape parameters for the Dirichlet distribution.
  //
  // The output is a NumericMatrix of size n x K, where each row is a Dirichlet sample.
  int K = alpha.size();               // Number of parameters (dimension)
  NumericMatrix gamma_samples(n, K);     // Vector to store Gamma samples
  
  // Step 1: Sample from Gamma distribution
  for(int j = 0; j< n; ++j){
    double sum_gamma = 0.0;
    for (int i = 0; i < K; i++) {
      gamma_samples[i] = R::rgamma(alpha[i], 1.0);  // Shape = alpha[i], Scale = 1.0
      sum_gamma += gamma_samples[i];                // Accumulate the sum
    }
    
    // Step 2: Normalize the Gamma samples
    for (int i = 0; i < K; i++) {
      gamma_samples[i] /= sum_gamma;   // Normalize by the total sum
    }
  }
  return gamma_samples;  // This is thes Dirichlet samples
}

// [[Rcpp::export]]
NumericMatrix convex_combination_coefficients(std::string model, bool original_constraints, int num_cc_constraints) {
  int model_constraints = value_p(model);
  // this function gets called by the R function test_U_stat
  
  // if(num_cc_constraints == 0) return();
  int ini;
  NumericMatrix coefficients;
  
  if(original_constraints || num_cc_constraints == 0){
    coefficients = NumericMatrix(model_constraints + num_cc_constraints, model_constraints);
    for (int i = 0; i < model_constraints; i++){
      coefficients(i,_) = rep(0,model_constraints);
      coefficients(i,i) = 1;
    }
    ini = model_constraints;
  }else{
    coefficients = NumericMatrix(num_cc_constraints, model_constraints);
    ini = 0;
  }
  
  // Generate random convex combinations and calculate convex combinations of constraints
  for (int j = 0; j < num_cc_constraints; j++) {
    // pick convex combinations uniformly at random and store them in a matrix:
    NumericVector c = dirichlet_sample(1, rep(1.0, model_constraints)); 
    coefficients(ini + j, _) = c;
  }
  return(coefficients);
}


//////////////////////////////////////////////////////////////////////// 
/////////// 3B. MAIN KERNEL AND KERNEL SUMMATION FUNCTIONS /////////////
////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
NumericVector h(std::string model, NumericMatrix Xm, int num_permutations, NumericMatrix cc_coefficients){
  int m = Xm.nrow();
  int s = value_s(model);
  //Rcout << m << " " << s << " " << num_permutations << std::endl;
  IntegerMatrix random_permutations = get_random_k_partitions_efficient(m, s, num_permutations);
  num_permutations = random_permutations.nrow();
  //num_permutations = std::min(std::max(num_permutations, 1), total_number_partitions(m, s));
  int ncols_h = cc_coefficients.nrow(); 
  //Rcout << ncols_h << std::endl;
  NumericMatrix total_h_tilde(num_permutations, ncols_h);
  for(int i = 0; i < num_permutations; ++i){
    IntegerVector perm = random_permutations(i,_);
    //printIntegerVector(perm);
    total_h_tilde(i,_) = convex_combination(h_tilde(model, Xm, perm), cc_coefficients);
  }
  //printNumericMatrix(total_h_tilde);
  NumericVector H = col_mean(total_h_tilde);
  return H;
}

// [[Rcpp::export]]
NumericMatrix calculate_H(std::string model, NumericMatrix X, IntegerMatrix permuted_indices, int num_permutations, NumericMatrix cc_coefficients){
  // Compute the incomplete U statistic by summing over many randomly-chosen intstances of the kernel
  // function evaluated at randoms subsets of the data.
  int N_hat = permuted_indices.nrow();
  
  int ncolsH = cc_coefficients.nrow();
  NumericMatrix H(N_hat, ncolsH);
  
  for(int i = 0; i < N_hat; ++i){
    NumericMatrix Xm = extract_rows(X, permuted_indices(i,_)); 
    H(i,_) = h(model, Xm, num_permutations, cc_coefficients);
  }
  return H;
}


//////////////////////////////////////////////////////////////////////////////
/////////// 4. FUNCTIONS TO COMPUTE DIVIDE-AND-CONQUER ESTIMATORS /////////// 
//////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
IntegerMatrix indices_ĝ_i1(int n, int m, int i1){
  int K = floor((n-1)/(m-1));
  
  IntegerVector comp_i1 = generate_sequence(n, 1);
  comp_i1.erase(i1-1);
  
  IntegerMatrix indices(K, m);
  
  int aux = 0;
  for(int k = 0; k < K; ++k){
    indices(k, 0) = i1;
    for(int j = 1; j < m; ++j) indices(k, j) = comp_i1(aux+j-1);
    aux += m-1;
  }
  
  return indices;
}

NumericVector ĝ_i1(std::string model, NumericMatrix X, int i1, int m, int num_permutations, NumericMatrix cc_coefficients){
   int n = X.nrow();
   int K = floor((n-1)/(std::max(1,m-1)));
   
   IntegerVector comp_i1 = generate_sequence(n, 1);
   //printIntegerVector(comp_i1);
   comp_i1.erase(i1-1);
   
   int ncolsg_i1 = cc_coefficients.nrow();
   NumericMatrix total_g_i1(K, ncolsg_i1);
   
   int aux = 0;
   for(int k = 0; k < K; ++k){
     IntegerVector indices(m);
     indices(0) = i1;
     for(int j = 1; j < m; ++j) indices(j) = comp_i1(aux+j-1);
     aux += m-1;
     //printIntegerVector(indices);
     NumericMatrix Xm = extract_rows(X, indices);
     total_g_i1(k,_) = h(model, Xm, num_permutations, cc_coefficients);
   }
   
   NumericVector g_i1= col_mean(total_g_i1);
   return g_i1;
 }

// [[Rcpp::export]]
NumericMatrix calculate_G_S1(std::string model, NumericMatrix X, IntegerVector S1, int m, int num_permutations, NumericMatrix cc_coefficients){
  
  int n1 = S1.size();
  
  int ncolsG = cc_coefficients.nrow();
  NumericMatrix G(n1, ncolsG);
  
  for(int i=0; i < n1; ++i){
    int i1 = S1(i);
    G(i, _) = ĝ_i1(model, X, i1, m, num_permutations, cc_coefficients);
  }
  return(G);
}

// [[Rcpp::export]]
NumericMatrix calculate_G(std::string model, NumericMatrix X, int n1, int m, int num_permutations, NumericMatrix cc_coefficients){
  int n = X.nrow();
  IntegerVector S1 = sample(n, n1, false);

  return(calculate_G_S1(model, X, S1, m, num_permutations, cc_coefficients));
}


/////////////////////////////////////////////
////////// 5. BOOTSTRAP FUNCTION ////////////
/////////////////////////////////////////////

// [[Rcpp::export]]
List bootstrap_U(int num_replicates, NumericMatrix H_centered, NumericMatrix G_centered){
  // num_replicates = number of bootstrap replicates to do. This is A in the SDL paper.
  int N_hat = H_centered.nrow();
  int n = G_centered.nrow();  // rename to the more correct n1?
  int p = H_centered.ncol();
  
  NumericVector epsilons((N_hat + n));
  NumericVector colsums_B(p);
  NumericMatrix U_B(num_replicates,p);
  NumericVector colsums_A(p);
  NumericMatrix U_A(num_replicates,p);
  
  for (int i = 0; i < num_replicates; i++){
    epsilons = rnorm((N_hat+n),0,1);
    // 
    for (int j = 0; j < p; j++){
      colsums_B[j] = sum(H_centered(_,j) * epsilons[Range(0, (N_hat-1))]);
    }
    U_B(i,_) = (1/sqrt(N_hat)) * colsums_B;
    //
    for (int j = 0; j < p; j++){
      colsums_A[j] = sum(G_centered(_,j) * epsilons[Range(N_hat, ((N_hat+n)-1))]);
    }
    U_A(i,_) = (1/sqrt(n)) * colsums_A;
  }
  List res = List::create(U_A, U_B);
  return res;
}



