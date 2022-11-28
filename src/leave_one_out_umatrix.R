# Date: June 20, 2018
# Procedure: leave-one-out U matrix calculation
# Purpose: to determine if leaving out markers reveals additional information

################### SETUP ################### 
library(tidyverse)

################### FUNCTIONS ###################

heatmap_wrapper <- function(matrix, title) {
    # Quick way to utilize pheatmap for SOM-related visualizations
    # Args:
    #   matrix: the matrix to be used in the visualization. Elements need to
    #       be the expression values of interest (eg. density or marker level)
    # Returns:
    #   pheatmap object without any clustering for visualization
    pheatmap(matrix, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE,
             main = title)
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
################### PIPELINE ###################

source('/Volumes/Samsung_T5/new.mac/drfz.external/large_som/src/mod_umatrix_maker.R')
setwd("/Volumes/Samsung_T5/new.mac/drfz.external/large_som/results/2018-06-18/data_objects")

# Get out the list
som_list <- list.files() %>% 
    .[grep("rds", .)] %>% 
    .[grep("80.to.100", .)] %>% 
    readRDS()

# Get out the 100 x 100 SOM for further work
som <- som_list$`100`
som_exp <- som$Weights %>% as.tibble()

# Get out the markers
input <- readRDS("nd_cells_input_markers.rds")
names(som_exp) <- input

# The loop and the list
umatrix_loo_list <- lapply(seq(length((input))), function(i) {
    umatrix <- umatrixForEsom_adapted(Weights = som$Weights[,-i], 
                                      Lines = som$Lines, 
                                      Columns = som$Columns, 
                                      Toroid = FALSE)
    print(i)
    return(umatrix)
})

names(umatrix_loo_list) <- paste("without", input, sep = ".")
saveRDS(umatrix_loo_list, "umatrix_list_loo_input_markers.rds")

# Visualization
for(i in seq(length(umatrix_loo_list))) {
    curr <- umatrix_loo_list[[i]]
    curr_hm <- heatmap_wrapper(curr, title = paste("umatrix", names(umatrix_loo_list[i]), sep = " "))
    save_pheatmap_pdf(curr_hm, filename = paste("heatmap_umatrix", names(umatrix_loo_list[i]), "pdf", sep = "."))
}
