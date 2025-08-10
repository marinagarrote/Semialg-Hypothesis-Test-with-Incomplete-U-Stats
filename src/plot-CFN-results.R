## main_plot_CFN.R -- This file contains code to draw treespace plots and
##                    p-value histograms for tests of the CFN model.
##
##
## How to run: first extract the file ../data/collection-1.zip. Then copy-paste
## code form this file line-by-line into R session with working directiory
## 'src/'
##
## Assumed structure of the working directory from which to run this code:
## ├── gitrepo/
## │   ├── src/
## │   │   ├── ...
## │   │   ├── plot-CFN-results.R
## │   │   ├── plot-CFN-results_utils.R
## │   │   └── ...
## │   ├── ...
## │   └── data/
## │       ├── ...
## │       └── collection-1/
## │           ├── ...
## │           ├── XYZ/
## │           │   ├── XYZ__PARAMETERS1
## │           │   │   ├── parameters.txt
## │           │   │   ├── results_pvals_1
## │           │   │   ├── results_pvals_2
## │           │   │   ├── ...
## │           │   │   ├── results_tstats_1
## │           │   │   ├── results_tstats_2
## │           │   │   └── ...
## │           │   ├── XYZ__PARAMETERS2
## │           |   └── ...
## │           └── ...
##
## Here, XYZ ranges over the set {CDD, PDR, PDM, CDR, CDM}. See the
## ../data/readme.md for details.


#_______________________________________________________________________________
#
# Load necessary packages and define helper functions for plotting
# _______________________________________________________________________________

library(gridExtra)  # Load the gridExtra library for arranging plots
library(data.table) # Load the data.table library for efficient data manipulation

# Function to read data from multiple files and optionally merge them
read_data <- function(data_path, merge = T){
  dataList = list()  # Initialize an empty list to store data from each file
  files = list.files(data_path, pattern = "results_pvals")  # Get a list of files matching the pattern
  # Loop through each file
  for(i in 1:length(files)){
    load(file=paste0(data_path, "/", files[i]))  # Load data from file (assumes data is named 'data_i_pvals')
    data = data_i_pvals 
    data = do.call(rbind, data)  # Combine data if it's in a list format
    dataList[[i]] = data  # Add the data to the list
  }
  # Merge the dataframes if merge = TRUE, otherwise return the list of dataframes
  if(merge){ 
    return(merge_dataframes(dataList))  # Call the merge_dataframes function to combine the data
  }else{ 
    return(dataList) 
  }  
}

# Function to merge a list of dataframes into a single data.table. 
merge_dataframes <- function(dataList) {
  library(data.table)
  combined_df <- rbindlist(lapply(dataList, function(df) df[, 1:6, drop = FALSE]))
  # COMMENTARY: we only save the first six columns. Some datasets have
  # additional columns (svd scores) that we don't need.
  return(combined_df)
}

# Function to plot histograms of p-values for a specific topology
plot_pvalue_histogram <- function(data, topology, title=paste0("P-values when testing topology ", topology)){
  # Create a copy of the data to avoid modifying the original
  data_copy <- copy(data)  
  # Define a named vector to map topology names to column names
  topologies <- c("12|34" = "pval_12", "13|24" = "pval_13", "14|23" = "pval_14")  
  # Rename the relevant p-value column to "pval"
  setnames(data_copy, topologies[topology], "pval")
  # Plot a single histogram for all experiments
  hist(data_copy[, pval], main = title, xlab = "p-values")  
}

# Function to find the maximum pvalue in each row, and add this information to
# the dataframe as a new column called 'max_top'. The value of 'max_top'
# equals 1,2, or 3 depending on whether the largest pvalue corresponded to
# topology 12|34, 13|24, or 14|23, respectively.
find_max_pval <- function(data){
  # Define a named vector to map column names to numeric values
  top_names <- c(pval_12 = 1, pval_13 = 2, pval_14 = 3)
  # Find the column name with the maximum value for each row
  data[, max_top := names(.SD)[apply(.SD, 1, which.max)], 
       .SDcols = c("pval_12", "pval_13", "pval_14")]
  # Convert column names to numeric values, and return output
  data[, max_top := top_names[max_top]]
  return(data)
}

# Function to plot a treespace visualization showing the success rate of
# inferring a specified toplogy.
plot_treespace <- function(data, topology, plot_title){
  # Create a copy of the data
  data_copy <- copy(data) 
  # Define a named vector to map topology names to numeric values
  topologies <- c("12|34" = 1, "13|24" = 2, "14|23" = 3)
  # Calculate the number of successes for the specified topology, and rename
  # the relevant column to "success_count"
  result_df <- aggregate(max_top ~ a + b, data = data_copy,
                         FUN = function(x) sum(as.numeric(x) == topologies[topology])) 
  setnames(result_df, "max_top", "success_count")
  # Calculate the success percentage
  total_num_pvals = nrow(data_copy)
  succ <- paste0("Success rate: ", 100*sum(result_df$success_count)/total_num_pvals, "%")  
  # Construct the image plot
  # Get unique values for the x- and y-axis
  x <- unique(result_df$b)
  y <- unique(result_df$a)
  # Create a matrix from the "top" values
  z <- array(result_df$success_count, dim=c(length(x),length(y)))  
  z_transpose <- t(z)
  # Paint the plot
  p <- image(x, y, z_transpose, col=gray.colors(100, start=1, end=0),
             main=paste0(plot_title," - ", succ), xlab='', ylab='') 
}

#_______________________________________________________________________________
#
# Examples of how to make histograms and treespace plots
#_______________________________________________________________________________

# First, specify the path to the data, and then read the data into an R
# dataframe. For these examples, we'll use Collection 1 data obtained using the
# CDD generating set, found in the following folder:
data_folder="CDD/CDD__n=1000,m=12,N=1000,n1=80,A=5000,s=100,r=20"
data <- read_data(file.path("../data/collection-1",data_folder), merge = T)

# EXAMPLE 1. Create a histogram of p-values aggregated over all choices of tree
# parameters (a,b) in Collection 1.
plot_pvalue_histogram(data,"12|34")
plot_pvalue_histogram(data,"13|24")
plot_pvalue_histogram(data,"14|23")

# EXAMPLE 2. Create treespace plots from data in Collection 1.
find_max_pval(data) # adds a column to 'data' indicating which topology has the
                    # largest p-value. This is needed to run plot_treespace().
plot_treespace(data, "12|34", "Example treespace plot for topology 12|34")
plot_treespace(data, "12|34", "Example treespace plot for topology 13|24")
plot_treespace(data, "12|34", "Example treespace plot for topology 14|23")



