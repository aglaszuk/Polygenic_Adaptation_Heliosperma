
# Perform GO terms enrichment of differentially expressed genes.

We performed separate GO terms enrichment analysis for each ecotype pair using R. Fisher test statistics implemented in the Bioconductor package topGO v.2.34.0 (https://bioconductor.org/packages/release/bioc/html/topGO.html) were run with the algorithm “weight01” to test for over-representation of specific functions. 
The scripts make.R and functions.R were used for this step.
