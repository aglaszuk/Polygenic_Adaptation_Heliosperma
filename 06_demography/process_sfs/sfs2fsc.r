#! /usr/bin/env Rscript

# Author: Aglaia Szukala, aglaia.szukala [at] univie.ac.at
# Title: sfs2fsc.r
# Script to convert the 2D-sfs as output by angsd into a 2D-sfs in fastsimcoal format + plotting as heatmap

# Usage: sfs2fsc.r inputfile.ml outdir n_0 n_1 id_0 id_1
# inputfile.ml is a derived allele 2D_sfs in ANGSD format
# n_0 number of entries for deme 0 (columns) in the 2D_sfs, i.e. the no. of individuals*2+1
# n_1 number of entries for deme 1 (rows) in the 2D_sfs, i.e. the no. of individuals*2+1
# id_0 and id_1 are the populations identifiers used by fastsimcoal (0,1,2,3)

# Load required packages
#library(RColorBrewer)
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape)

# Read in arguments
args <- commandArgs(trailingOnly=T)
sfs <- args[1]
outdir <- args[2]
n_0 <- as.numeric(args[3])
n_1 <-  as.numeric(args[4])
id_0 <- args[5]
id_1 <- args[6]

# Function to transform 2D site frequency spectrum output by ANGSD into fastsimcoal format
write.2D.fsc<-function(sfs,f.output)
{
  # Load the .ml file
  data <- matrix(as.numeric(read.table(sfs)), nrow = n_1, ncol = n_0)
  
  # Add row and column ID
  colnames(data) <- paste0(rep(paste0("d",id_0,"_"),n_0),c(0:(n_0-1)))
  rownames(data) <- paste0(rep(paste0("d",id_1,"_"),n_1),c(0:(n_1-1))) 
  # Create output
  cat("1 observations\n\t",file=paste0(outdir,f.output))
  write.table(data,file=paste0(outdir,f.output),
              append=T,sep="\t",col.names=T,row.names=T,quote=F)

}
write.2D.fsc(sfs,f.output = paste0("HELIO_jointMAFpop",id_1,"_",id_0,".obs"))

# Plot sfs as heatmap and save as pdf
# Load the .ml file
data <- matrix(as.numeric(read.table(sfs)), nrow = n_1, ncol = n_0)

# Add row and column ID
colnames(data) <- paste0(rep(paste0("d",id_0,"_"),n_0),c(0:(n_0-1)))
rownames(data) <- paste0(rep(paste0("d",id_1,"_"),n_1),c(0:(n_1-1))) 

# Prepare the data.frame for ggplot2
data_plot <- melt(data)
data_plot$pop0 <- factor(data_plot$X2, levels = colnames(data))
data_plot$pop1 <- factor(data_plot$X1, levels = rownames(data))
data_plot$freq <- data_plot$value

# Plot the heatmap
pdf(paste0(outdir,"HELIO_jointMAFpop",id_1,"_",id_0,".pdf"))
ggplot(data_plot, aes(x=pop0,y=pop1, fill=freq)) + 
  geom_tile()+
  theme_bw()+
  #scale_fill_distiller(palette = 'Spectral')+
  scale_fill_gradientn(colours = c("white", "cyan4", "darkred"),
                       values = scales::rescale(c(-0.05, -0.03, 0, 0.6, 0.8))) +
  xlab(paste0('pop',id_0)) +
  ylab(paste0('pop',id_1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

