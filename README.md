# PopDiff_MonocyteIAV
Repository associated with the following manuscript: Basal Expression of Interferon-Stimulated Genes Drives Population Differences in Monocyte Susceptibility to Influenza Infection

System requirements
-all scripts were run on linux CentOS 6.0 with R.3.0 or macOS Catalina with R.3.6.3

The following R libraries and their dependencies were installed from CRAN or Bioconductor:
library(scran)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(reshape2)
library(Matrix)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

Installation guide:
The provided scripts require no specific installation aside from R and required libraries. Libraries and dependencies can be installed in <2 hours.

Demo & Instructions for use:
Cell Metadata and miscellaneous data needed to reproduces single cell analyses from the manuscript can be found in the data Folder. 
The normalized count matrix is too large to be hosted normally on github and is being tracked using git LFS. In order to reproduce the results you will need to use the full normalized count matrix; if you simply wish to see the structure of the dataset you can instead look at a much smaller sample of the matrix (data/normalized_count_matrix_sample_2600monocytes_6845hostgenes.rds).

Scripts can be found in the scripts folder. Each script can be run independently of the others. These will generate results in the results folder.
Total expected run time (2-3 hours).

Expected outputs can be found in the results folder with 'compare' captured at the end of the file names prior to the file extension.
