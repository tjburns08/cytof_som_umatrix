## Introduction

[FlowSOM](https://bioconductor.org/packages/release/bioc/html/FlowSOM.html), developed by Sofie VanGassen in the laboratory of Yvan Saeys, is a widely used algorithm in CyTOF data, as it has been found to be [both fast and competent at detecting rare cell subsets](https://pubmed.ncbi.nlm.nih.gov/27992111/). 

FlowSOM uses a self-organizing map (SOM) to quickly subset cells into clusters. The SOM, which is shaped like a 2-D grid, can be unraveled and visualized to give the users better intuition around the clustering tool and the dataset at large. This visualization is called the U-Matrix, and as of 2022, this data structure is not widely used in the CyTOF community. 

This project, done in the laboratory of Yvan Saeys in 2018, was an attempt to explore the U-Matrix using large SOMs. FlowSOM's default is 10x10, but here I looked at up to 100x100 SOMs for datasets of 100k cells. 

The data used in this project was CyTOF PBMCs from the laboratory of Henrik Mei at the German Rheumatism Research Center (DRFZ), a laboratory known for producing top quality CyTOF data. 

## How to navigate this project

To look at the punchline, look at the pdf in doc/. To look at the code that generated the analysis and visualizations that ended up in the aforementioned pdf, go to src/. We note that we don't use the FlowSOM package here. We use the [Umatrix package from CRAN](https://cran.r-project.org/web/packages/Umatrix/index.html). 

Specific scripts of interest, in src/:

make_som_dim_titration.R: this is where you'll see how the SOMs are made. 
mod_umatrix_maker.R: this is a helper function I didn't make, but I modified, to make the U-matrix. 
som_mataclustering.R: I do the meta-clustering step to see how this looks on the U-Matrix, and I do a KNN density estimation as well. 
find_moore_neighborhood.R: I compare the moore neighborhood of the SOM to the KNN of each node of the trained SOM in feature space. 
leave_one_out_umatrix.R: I remove one marker from the original input markers and determine if this has an impact on how the U-matrix looks. 
color_som_by_marker.R: I produce the U-matrix and color it by marker as one would a dimension reduction embedding. 
external_code/: code I didn't write, but either used or modified throughout the project. 
