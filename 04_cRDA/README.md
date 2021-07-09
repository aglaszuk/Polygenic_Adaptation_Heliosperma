# Expression divergence between ecotypes captured by conditioned redundancy analysis (cRDA).

As a table of response variables in the cRDA we used the matrix of read counts after filtering the reads using a mean count per million threshold of one as in the script for differential expression analysis. The cRDA includes a multiple regression step of gene expression on the explanatory variable(s). 
In our case, the RDA was conditioned to remove the effects of the geographic ecotype pair using the formula “~ ecotype + Condition(pair)”. 
In the second step, a principal component analysis (PCA) of the fitted values from the multiple regression is performed to produce canonical axes, based on which an ordination in the space of the explanatory variable is performed. 
Therefore, the first axis of the cRDA (Fig. 4 of the preprint) shows the variance explained by the constrained variable “ecotype”, while the second axis is the first component of the PCA nested into the RDA, representing the main axis of unconstrained variance.
The scripts cRDA.Rmd and cRDAfunctions.R show the single steps of this analysis.
