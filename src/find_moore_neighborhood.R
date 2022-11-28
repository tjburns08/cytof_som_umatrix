# Date: June 22, 2018
# Procedure: to find the identity of the moore neighborhood of the SOM
# Purpose: to look at SOM fidelity

################## SETUP ################## 
library(tidyverse)
library(Umatrix)
library(pheatmap)

################## FUNCTIONS ##################

GetMoore <- function(element, input_matrix) {
    # Returns the rest of the moore neighborhood outside of the center element
    # Args:
    #   element: the center element of the input matrix
    #   input_matrix: the input matrix
    # Returns:
    #   a vector of the 8 remaining members of the moore neighborhood
    locale <- which(input_matrix == element, arr.ind = TRUE)
    if(locale[1] == 100 || locale[2] == 100) {
        print("100 found")
        return(0)
    }
    
    if(locale[1] == 1 || locale[2] == 1) {
        print("1 found")
        return(0)
    }
    
    neighbor_coords <- list(c(locale[1] + 1, locale[2]),
                      c(locale[1] - 1, locale[2]),
                      c(locale[1], locale[2] + 1),
                      c(locale[1], locale[2] - 1),
                      c(locale[1] + 1, locale[2] + 1),
                      c(locale[1] + 1, locale[2] - 1),
                      c(locale[1] - 1, locale[2] + 1),
                      c(locale[1] - 1, locale[2] - 1)) 
    
    neighbors <- lapply(neighbor_coords, function(i) {
        input_matrix[i[1], i[2]]
    }) %>% unlist %>%
        unname
    
    return(neighbors)
}

heatmap_wrapper <- function(matrix) {
    # Quick way to utilize pheatmap for SOM-related visualizations
    # Args:
    #   matrix: the matrix to be used in the visualization. Elements need to
    #       be the expression values of interest (eg. density or marker level)
    # Returns:
    #   pheatmap object without any clustering for visualization
    pheatmap(matrix, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE,
             show_rownames = FALSE,
             show_colnames = FALSE)
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

################## PIPELINE ##################

# Get out the 100 x 100 SOM
setwd("/Volumes/Samsung_T5/new.mac/drfz.external/large_som/results/2018-06-18/data_objects")
som <- readRDS("som.results.n.100k.dim.80.to.100.rds")[[3]]

# Get out the respective cells
cells <- readRDS("nd_cells.100k.rds")

# Get out the identities
side <- 100

# The format by which the U matrix is laid out, as determined emperically
id_matrix <- 1:side**2 %>%
    matrix(data = ., nrow = side, ncol =  side) %>% 
    t %>%
    as.tibble

# For each part of the matrix, find the moore neighborhood
# Lay the matrix into a vector
id_vec <- 1:max(id_matrix)

nn_moore <- lapply(id_vec, function(i) {
    neighbor <- GetMoore(i, id_matrix)
}) %>% do.call(rbind, .) %>%
    as.tibble

# Find K-nearest neighbors
som_exp <- som$Weights %>% as.tibble()
nn <- Sconify::Fnn(som_exp, names(som_exp), k = 8)$nn.index %>% as.tibble

# Calculate percentages of the union for each. We expect the outsides to be
#   zero
nn_sim <- lapply(seq(nrow(nn)), function(i) {
    curr <- nn[i,] %in% nn_moore[i,] %>% 
        as.numeric() %>% 
        sum()
    curr <- curr/8
    return(curr)
}) %>% unlist() 

# Now convert the vector into the umatrix format and heatmap it
nn_sim_matrix <- matrix(nn_sim, nrow = side, ncol = side) %>% 
    t %>% 
    as.tibble

moore_hm <- heatmap_wrapper(nn_sim_matrix)
save_pheatmap_pdf(moore_hm, "knn_vs_nn_moore.pdf")

# Get the moore neighborhood out of it (minus the center)

