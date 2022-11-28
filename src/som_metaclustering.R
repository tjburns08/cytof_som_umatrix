# Date: June 20, 2018
# Procedure: Meta-cluster a SOM and then color the map by the meta-cluster
#   identities
# Purpose: to get a visual of how the meta-clusters fit on to the U-matrix

##################### SETUP ##################### 
library(tidyverse)
library(Sconify)

##################### FUNCTIONS ##################### 

heatmap_wrapper <- function(matrix) {
    # Quick way to utilize pheatmap for SOM-related visualizations
    # Args:
    #   matrix: the matrix to be used in the visualization. Elements need to
    #       be the expression values of interest (eg. density or marker level)
    # Returns:
    #   pheatmap object without any clustering for visualization
    pheatmap(matrix, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE)
}


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
    # Credit: https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}


clust_to_matrix <- function(id, nrow, ncol) {
    # Given a vector of cluster IDs, wrap up into the same 2-D structure of
    #   the SOM, such that its plot can be compared to the U-Matrix
    # Args:
    #   id: vector of cluster id
    #   nrow: number of rows of the SOM
    #   ncol: number of columns of the SOM
    result <- matrix(id, nrow, ncol) %>% t()
    return(result)
}


randomize_cluster_id <- function(id) {
    # Create a randomized vector between 1 and the max number in the vector
    tmp <- sample(1:max(id))
    
    # Switch the input vector to be the elements of the ranomdized vector
    result <- lapply(id, function(i) {
        return(tmp[i])
    }) %>% unlist
    
    # return it
    return(result)
}

##################### PIPELINE ##################### 

setwd("/Volumes/Samsung_T5/new.mac/drfz.external/large_som/results/2018-06-18/data_objects")

# Get out the list
som_list <- list.files() %>% 
    .[grep("rds", .)] %>% 
    .[grep("80.to.100", .)] %>% 
    readRDS()

# Get out the markers
input <- readRDS("nd_cells_input_markers.rds")

# Get out the 100 x 100 SOM for further work
som <- som_list$`100`
som_exp <- som$Weights %>% as.tibble()
names(som_exp) <- input

# The clustering and ID generation
n <- nrow(som_exp)
hc <- hclust(dist(som_exp))
hc_id <- cutree(hc, k = 10) %>% randomize_cluster_id()

# The heatmap
id_som <- clust_to_matrix(hc_id, 100, 100)
heatmap_wrapper(id_som)

# The loop for this
clust_titration <- c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000, 2000, 5000)
n <- nrow(som_exp)
nrow <- ncol <- sqrt(nrow(som_exp))

clust_results <- lapply(clust_titration, function(i) {
    curr_id <- cutree(hc, k = i) %>% randomize_cluster_id()
    curr_id_som <- clust_to_matrix(curr_id, nrow, ncol)
    h <- heatmap_wrapper(curr_id_som)
    save_pheatmap_pdf(h, paste("heatmap", "nclust", i, "pdf", sep = "."))
})

# Extra credit: Density matrix rather than U matrix
k_titration <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
knn_de_results <- lapply(k_titration, function(i) {
    som_nn <- Sconify::Fnn(som_exp, names(som_exp), k = i)[[2]] %>% as.tibble()
    som_dist <- apply(som_nn, 1, mean) %>% clust_to_matrix(., nrow, ncol)
    h <- heatmap_wrapper(som_dist)
    save_pheatmap_pdf(h, paste("heatmap", "knn_de", i, "pdf", sep = "."))
})


