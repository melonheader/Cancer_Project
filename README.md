# Cancer_Project
Partly rewritten R code for Master's thesis project

#GDC.Download.by.Project.R 
1.Downloads TCGA expression data using TCGA biolinks library

#GDC.Matr.Prepare.R
1.Prepares expression matrices from downloaded TCGA data and other auxiliary gene expression data 
which need to be specified
(in my case it is a gene expression for normal tissues from Gtex and Definitive Endoderm Cells 
from GEO)
2.Downloads and writes metadata for all TCGA projects and other samples.
3.Merges together all expression matrices

#GDC.Normalize.R
1.Normalizes (edgeR) merged gene expression matrix from previous script and writes it.

#GDC.Filter.R
1.Filter expression matrix by specified set of genes.
2.Merges filtered table with metadata.#

#GDC.Plot.R
1.Code to ggplot PCA results from generated expression tables using
