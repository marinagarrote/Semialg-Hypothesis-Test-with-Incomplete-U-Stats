# ------------------------------------------------- #
# ----------------- Generate DATA ----------------- #
# ------------------------------------------------- #

generate_dataset_trinomial <- function(x, y, z){
    ## Input: a triple of counts (x, y, z) of integers with at least one non-zero entry.
    ##
    ## Output: A data matrix X of the same form as that output by the function
    ##         generate_random_CF_data. X has 3 columns and x+y+z rows. Each row
    ##         is a 0/1 vector which encodes a tree topology, s.t. the column sums
    ##         are (x,y,z)
  
    if (x > 0) 
    data1 = t(matrix(c(1, 0, 0), 3, x))
    else
    data1 = NULL
    if (y > 0)
    data2 = t(matrix(c(0, 1, 0), 3, y))
    else
    data2 = NULL
    if (z > 0)
    data3 = t(matrix(c(0, 0, 1), 3, z))
    else
    data3 = NULL
    treedata = rbind(data1, data2, data3)
    reorder = sample(x+y+z) #randomize order
    treedata = treedata[reorder,]
    return(treedata)
}



# ------------------------------------------------- #
# ----------------- Simplex Plots ----------------- #
# ------------------------------------------------- #
library(MSCquartets)

### Simplex
alphas = c(.10, .05, .01) # p-value cutoffs between regions
setcolors = c("purple", "blue", "green", "red")

initializeSimplexPlot <- function(model, n, m, N, n1, num_permutations, original_constraints, num_cc_constraints){

    maintitle = paste0("SDL Test for ", model , "\n\n")
    titletext= paste0("n=", n, "; m=", m, "; N=", N, "; n1=", n1, "; Num permutations=", num_permutations,  "\n",
                    "# RCC constraints=", num_cc_constraints, "; original constraints included=", original_constraints)

    if (model ==  "model-1") {
        simplexPrepare("T1", main = maintitle, titletext = titletext)
        simplexSegment(c(1, 0, 0), c(0,1,1), lwd = 2)
    }else if(model ==  "model-2") {
        simplexPrepare("T1", main = maintitle, titletext = titletext)
    }else if(model ==  "model-3") {
        simplexPrepare("cut", main = maintitle, titletext = titletext)
    }else if(model ==  "model-4") {
        simplexPrepare("T3", main = maintitle, titletext = titletext)
    }else if(model ==  "point") {
        simplexPrepare("T3", main = maintitle, titletext = titletext)
    }
} 
