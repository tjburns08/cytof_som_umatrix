# Date: June 18, 2018
# Procedure: make a large self organizing map
# Purpose: to extend the size of Flow-SOM experimentally, to expand on some
#   of the work pioneered by Quentin Rouchon in the laboratory of Yvan Saeys. 

#################### SETUP #################### 
library(flowCore)
library(tidyverse)
library(Umatrix)
library(Sconify)
library(reshape)

#################### PIPELINE ####################

# Test on the iris dataset
#test <- Umatrix::esomTrain(Data = as.matrix(iris[,2:4]), Columns = 20, Lines = 20, Toroid = FALSE)

# Test on the TSORA healthy control
setwd("/Volumes/Samsung_T5/new.mac/drfz.external/large_som/data/tsora")
nd_markers <- ParseMarkers("markers.csv")
surface <- nd_markers[[1]]

# Get out the cells input
ncells <- 10000
nd_fcs <- list.files() %>% .[grep("fcs", .)]
set.seed(384302)
nd_cells <- Sconify::ProcessMultipleFiles(nd_fcs, 
                                          numcells = ncells, 
                                          input = surface)

# Prepare the data for the SOM mapping (requires a matrix)
nd_matrix <- as.matrix(nd_cells[,surface])

# The x and y size of the SOM in the following loop
# dim_titration <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 1000) 
dim_titration <- c(10, 100, 1000)

# Build the SOM of the size of the aforementioned dimensions
set.seed(10284)
som_results <- lapply(dim_titration, function(i) {
    nrow <- i
    ncol <- i
    
    # The SOM built from the Umatrix package. NOT a torous. 
    nd_som <- Umatrix::esomTrain(Data = nd_matrix, 
                                 Columns = ncol, 
                                 Lines = nrow, 
                                 Toroid = FALSE)
    
    # Make a plot of the SOM U-matrix and save it using pre-made function
    # p <- pheatmap(res_esomtrain = nd_som)
    # ggsave(filename = paste("matrix", nrow, "x", ncol, "ncells", ncells, "png", sep = "."))
    
    return(nd_som)
})

# Save the Rtsne results
nd_tsne <- Rtsne::Rtsne(nd_cells[,surface], verbose = TRUE)$Y %>% as.tibble()
names(nd_tsne) <- c("bh-SNE1", "bh-SNE2")
nd_cells <- bind_cols(nd_cells, nd_tsne)
saveRDS(nd_cells, paste("nd_cells", ncells, "rds", sep = "."))
write.csv(nd_cells, paste("nd_cells", ncells, "csv", sep = "."))
saveRDS(surface, "nd_cells_input_markers.rds")

# Save the SOM results
names(som_results) <- dim_titration
saveRDS(som_results, paste("som.results.n", 
        ncells, "dimrange", 
        dim_titration[1], 
        dim_titration[length(dim_titration)], 
        "rds", 
        sep = "."))


