# Visualizing CyTOF data with iSEE

*Anne-Maud Ferreira - Statistics department, Stanford University*

## Version information

R version: R version 4.0.2 (2020-06-22) 

Bioconductor version: 3.11

```
library(flowCore)
library(iSEE)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(CATALYST)
library(uwot)
```

## Introduction

The goal of this README page is to show how to create a `SingleCellExperiment` (SCE) from the *.fcs* files and launch the `iSEE` shiny app to perform some data exploration.

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
class(sce_data)
## [1] "SingleCellExperiment"
## attr(,"package")
## [1] "SingleCellExperiment"
```

![SingleCellExperiment object structure](https://bioconductor.github.io/BiocWorkshops/202_Das_SingleCellRNASeq/SingleCellExperiment.png)

image from the **The Bioconductor 2018 Workshop Compilation** online book. 

## Analysis

### Clustering analysis with `CATALYST`

To perform the clustering with the `CATALYST` package, we need to specify the type of the markers on which the clustering will be performed.

```
# Specify markers to use for clustering
type_markers(sce_data)
# Cluster
sce_data <- cluster(sce_data, 
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

The `iSEE` package provides a shiny app with a combination of interactive panels. Just launch the app with the following command:

```
iSEE(sce_data)
```

**Disclosure**: All the information compiled below are from the `iSEE` package documentation available online at the following link: https://bioconductor.org/packages/release/bioc/html/iSEE.html under the *Documentation* section.

### Parameters

There are 3 types of parameters available in `iSEE` app:

- Data parameters: to control the parameters which are specific to each plot (e.g. dimensionality reduction input, marker selection);
- Visual parameters: to determine the visual of each plots (e.g. colors, facet, font);
- Selection parameters: to control the point selection that the plot will receive.

### Plots

The app have eight default panels. 

The five following panels were presented during the group meeting:

- Reduced dimension plot: to visualize the dimensionality reduction stored in the SCE object (e.g. UMAP). It is a great plot to visualize the different clusters, as well as explore any visual differences between conditions or patients;
- Row data table: to display the values of the `rowData` slot, i.e. marker metadata;
- Feature assay plot: to plot the assayed values, i.e. counts and exprs of each marker in the cells. This plot can be used to explore data in 2D plot (marker A versus marker B). It also allows to see the distribution of each marker with a violin plot;
- Column data table: to display the values of the `colData` slot, i.e. sample metadata;
- Heatmap: to overview the data for multiple features. This plot shows the assayed values, where the markers are the rows and the cells are the columns.

The 3 additional panels exist:

- Column data plots: sample metadata (colData slot);
- Row data plots: marker metadata  (rowData slot);
- Sample assay plots: assayed values (marker intensities) for a particular cell across the markers on the y-axis.

### Additional features

The app also allows:

- To remove and reformat plots;
- To download the plots;
- To extract R code to reproduce plots;
- To extract settings to reproduce panel visualization.

## References

Rue-Albrecht K, Marini F, Soneson C, Lun ATL (2018). “iSEE: Interactive SummarizedExperiment Explorer.” F1000Research, 7, 741. doi: 10.12688/f1000research.14966.1.

Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J, Jiang M, Finak G (2020). flowCore: flowCore: Basic structures for flow cytometry data. R package version 2.0.1.

Crowell H, Zanotelli V, Chevrier S, Robinson M (2020). CATALYST: Cytometry dATa anALYSis Tools. R package version 1.12.2, https://github.com/HelenaLC/CATALYST.

Lun A, Risso D (2020). SingleCellExperiment: S4 Classes for Single Cell Data. R package version 1.10.1.

https://bioconductor.github.io/BiocWorkshops/

Morgan M, Obenchain V, Hester J, Pagès H (2020). SummarizedExperiment: SummarizedExperiment container. R package version 1.18.2.

https://cran.r-project.org/web/packages/uwot/index.html
