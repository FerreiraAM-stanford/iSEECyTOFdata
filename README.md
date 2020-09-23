# Visualizing CyTOF data with iSEE
Anne-Maud Ferreira - Statistics department, Stanford University

## Version information

R version: R version 4.0.2 (2020-06-22) 
Bioconductor version: 3.11

library(FlowRepositoryR)
library(flowCore)
library(iSEE)
library(CyTOFtoolbox)
library(CytoGLMM)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(tidyverse)
library(CATALYST)

## Data

### Load the data

#### FCS files

Read the FCS files as a `flowSet` using the `flowCore` package.

```
# Read files
fcs_files <- list.files(path = data_path, pattern = "fcs")
# Read files as a flowSet
flowSet_data <- read.flowSet(files = fcs_files,
                             path = data_path)
```

#### Marker information

The marker file contains information about the markers used in the FCS files. It is composed of 3 columns corresponding to the fcs column names, the antigen and the marker class.

```
# Read file
markers_info <- read.table(paste0(data_path, "markers_info.txt"))
# Check
## Column names
colnames(markers_info)
## [1] "fcs_colname"  "antigen"      "marker_class"
## Marker class
unique(markers_info$marker_class)
## [1] "none"  "type"  "state"
```

#### Sample information

The samples information describes at least the file name, patient ID and condition. This file may contain any useful additional information for the analysis.

```
# Read file
samples_info <- read.table(paste0(data_path, "samples_info.txt"))
## NB: File name must be a character string and not a factor
# Check
## Column names
colnames(samples_info)
## [1] "filename"   "patient_id" "condition"
```

### Prepare the data

The last step is to create a `SingleCellExperiment` (SCE) containing the FCS file values, as well as the metadata. The input parameters are (1) the flowSet, (2) marker's information, (3) sample's information, (4) a list specifying the column names in the sample's information. 

Details about the data preparation function can be find by typing `?CATALYST::prepData`.

```
## Prepare the data
sce_data <- CATALYST::prepData(flowSet_data,
                               panel = markers_info,
                               md = samples_info,
                               md_cols = list(file = "filename",
                                              id = "filename",
                                              factors = c("patient_id", "condition")),
                               cofactor = 5) # automatically applied arcsinh transformation with cofactor = 5
# Class
class(aghaeepour_data)
## [1] "SingleCellExperiment"
## attr(,"package")
## [1] "SingleCellExperiment"
```

## Analysis

### Clustering analysis with `CATALYST`

To perform the clustering, we need to specify the type of the markers on which the clustering will be performed.

```
# Specify markers to use for clustering
type_markers(sce_data)
# Cluster
aghaeepour_data <- cluster(sce_data, 
                           features = type_markers(sce_data),
                           xdim = 10, ydim = 10, maxK = 20, 
                           verbose = FALSE, seed = 1234)
```

#### Delta area plot

Based on the delta area plot which represents the stability gained when using *k* clusters, we can determine the optimal number of cluster:

```
# Print the delta area plot
metadata(sce_data)$delta_area
```

#### Marker expression's heatmap

Once the optimal number of clusters is determined, the marker expressions in each cluster can be visualized in the following heatmap:

```
# For instance if we decide on 10 metaclusters
plotFreqHeatmap(sce_data, k = "meta10")
```

### Dimensionality reduction with UMAP

The dimensionality reduction (DR) is computed with UMAP, using the `uwot` package.

```
# Extract data from sce object
umap_data_input <- t(assays(sce_data)$exprs)
# UMAP dimensionality reduction
set.seed(1234)
dr_umap_data <- uwot::umap(umap_data_input, 
                               n_neighbors = 20, 
                               min_dist = 0.1)
```

The computed DR is added to the SCE object.

```
# Add the dimensionality reduction to the SCE
colnames(dr_umap_data) <- c("X1", "X2")
reducedDims(sce_data) <- list("UMAP" = dr_umap_data)
```

## Interactive visualization of the data

```
iSEE(sce_data)
```
