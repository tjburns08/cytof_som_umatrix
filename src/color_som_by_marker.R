# Date: June 19, 2018
# Procedure: color the SOM results by marker expression
# Purpose: to be able to make direct comparisons between itself and t-SNE

################### SETUP ################### 
library(tidyverse)
library(pheatmap)

################### FUNCTIONS ################### 

heatmap_wrapper <- function(matrix) {
    # Quick way to utilize pheatmap for SOM-related visualizations
    # Args:
    #   matrix: the matrix to be used in the visualization. Elements need to
    #       be the expression values of interest (eg. density or marker level)
    # Returns:
    #   pheatmap object without any clustering for visualization
    p <- pheatmap(matrix, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE
             #color = colorRampPalette(c("white","black"))(100)
    )
    return(p)
}


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
# Credit: https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    png(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

exp_to_matrix <- function(som_exp, marker, nrow, ncol) {
    # Places marker expression format into the proper data for SOM visualization
    # Args:
    #   som_exp: expression tibble of each SOM node 
    #   marker: The marker of interest, in vector format of expression across
    #       SOM nodes, to compare with the U-Matrix visualizations
    #   ncol: the number of columns in the SOM
    #   nrow: the number of rows in the SOM
    result <- som_exp[[marker]] %>% matrix(., nrow, ncol) %>% t()
    return(result)
}


################### DATA PRE-PROCESSING ################### 

# Working directory containing the SOM results
# For first run
#setwd("/Volumes/Samsung_T5/new.mac/drfz.external/large_som/results/2018-06-18/data_objects")

# Get out the list, for first run
# som_list <- list.files() %>% 
#    .[grep("rds", .)] %>% 
#    .[grep("80.to.100", .)] %>% 
#    readRDS()

# For ESOM 
som_list <- list.files() %>% 
    .[grep("rds", .)] %>% 
    .[grep("dimrange", .)] %>% 
    readRDS()

# Get out whatever dimension you need
dim <- 1000 %>% as.character()
som <- som_list[[dim]]
som_exp <- som$Weights %>% as.tibble()

# Get out the markers
input <- readRDS("nd_cells_input_markers.rds")
names(som_exp) <- input

# Get the cell selection abundance for each node
node_to_cell <- som$BestMatches %>% as.tibble()
names(node_to_cell) <- c("cell_id", "som_row", "som_col")

################### COLOR SOM BY ABUNDANCE ################### 

# Cell abundance for each node of the SOM
sparse <- TRUE
if(sparse == TRUE) {
    som_table <- matrix(0, nrow = 1000, ncol = 1000)
    for(i in seq(nrow(node_to_cell))) {
        curr <- node_to_cell[i,]
        som_table[curr[[2]], curr[[3]]] <- som_table[curr[[2]], curr[[3]]] + 1
    }
} else {
    som_table <- table(node_to_cell$som_row, node_to_cell$som_col) 
}

# Populate the matrix of zeros with the values of the table
som_abundance <- lapply(seq(ncol(som$Umatrix)), function(i) {
    result <- som_table[,i] %>% unname()
}) %>% do.call(cbind, .)

dim(som_abundance) # Should be what was originally specified

################### HEATMAPS ###################

# Plot it
# The Umatrix
hu <- heatmap_wrapper(som$Umatrix)
save_pheatmap_pdf(hu, "heatmap_umatrix.pdf")

# Abundance of cells
ha <- heatmap_wrapper(som_abundance)
save_pheatmap_pdf(ha, "heatmap_abundance.pdf")

# Thresholded abundance of cells
som_abundance_thr <- ifelse(som_abundance > 0, 1, 0)
ht <- heatmap_wrapper(som_abundance_thr)
save_pheatmap_pdf(ht, "heatmap_abundance_0_1_thr.pdf")

################### MODIFIED UMATRIX CONSTRUCTION ###################

# Modifying the U-matrix generation function (see mod_umatrix_maker.R)
source('/Volumes/Samsung_T5/new.mac/drfz.external/large_som/src/mod_umatrix_maker.R')

umatrix_euclidean_mean <- umatrixForEsom_adapted(som$Weights, 
                                       som$Lines, 
                                       som$Columns, 
                                       som$Toroid,
                                       dist_measure = "euclidean",
                                       neighbor_agg = "mean")

umatrix_euclidean_max <- umatrixForEsom_adapted(som$Weights, 
                                                som$Lines, 
                                                som$Columns, 
                                                som$Toroid,
                                                dist_measure = "euclidean",
                                                neighbor_agg = "max")

umatrix_cheb_mean <- umatrixForEsom_adapted(som$Weights, 
                                            som$Lines, 
                                            som$Columns, 
                                            som$Toroid,
                                            dist_measure = "chebyshev",
                                            neighbor_agg = "mean")

umatrix_cheb_max <- umatrixForEsom_adapted(som$Weights, 
                                           som$Lines, 
                                           som$Columns, 
                                           som$Toroid,
                                           dist_measure = "chebyshev",
                                           neighbor_agg = "max")

emean <- heatmap_wrapper(umatrix_euclidean_mean)
save_pheatmap_pdf(emean, "euclidean_mean_umatrix.pdf")

emax <- heatmap_wrapper(umatrix_euclidean_max)
save_pheatmap_pdf(emax, "euclidean_max_umatrix.pdf")

cmean <- heatmap_wrapper(umatrix_cheb_mean)
save_pheatmap_pdf(cmean, "chebyshev_mean_umatrix.pdf")

cmax <- heatmap_wrapper(umatrix_cheb_max)
save_pheatmap_pdf(cmax, "chebyshev_max_umatrix.pdf")

################### COLOR SOM BY MARKER ################### 

# Get out the marker expression per node
som_markers <- som$Weights %>% as.tibble()
names(som_markers) <- input

# Wrap a given marker in this matrix into the 100 x 100 format and test
for(i in input) {
    curr <- exp_to_matrix(som_markers, i, 1000, 1000)
    p <- heatmap_wrapper(curr)
    save_pheatmap_pdf(p, paste("heatmap_exp", i, "png", sep = "."), width = 1000, height = 1000)
}




