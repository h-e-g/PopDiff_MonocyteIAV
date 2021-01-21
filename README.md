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
The provided scripts require no specific installation aside from R and required libraries.

Demo & Instructions for use:
Cell Metadata and miscellaneous data needed to reproduces single cell analyses from the manuscript can be found in the data Folder. 
The normalized count matrix is too large to be hosted on github and needs to be manually added into the data folder.

Scripts can be found in the scripts folder. Each script can be run independently of the others. These will generate results in the results folder.
Total expected run time (2-3 hours).

Expected outputs can be found in the results folder with 'compare' captured at the end of the file names prior to the file extension.
