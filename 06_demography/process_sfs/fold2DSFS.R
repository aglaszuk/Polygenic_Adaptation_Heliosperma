#! /usr/bin/env Rscript

# Original author: David Marques, davidalexander.marques [at] eawag.ch
# Modified by Aglaia Szukala, aglaia.szukala [at] univie.ac.at
# Title: fold2DSFS.R

# Usage: fold2DSFS.R filename_jointMAFpopX_Y.[obs/txt] outdir
# Folds 2-dimensional SFS in fastsimcoal2 format

library(ggplot2)
library(dplyr)
library(reshape)

# Parse file name
args<-commandArgs(trailingOnly=T)
base=sub("_joint.*","",args[1])
pops=strsplit(sub(".*pop","",args[1]),"[.]")[[1]][1]
suffix=substr(args[1], nchar(args[1])-3+1, nchar(args[1]))
outdir <- args[2]

# Read SFS
a=read.table(args[1],skip=1)

# Function to fold 2D-SFS
derived2maf=function(derived_sfs2d){
  n1=nrow(derived_sfs2d)
  n2=ncol(derived_sfs2d)
  maf_2dsfs=matrix(0,nrow=n1,ncol=n2)
  colnames(maf_2dsfs)=colnames(derived_sfs2d)
  rownames(maf_2dsfs)=rownames(derived_sfs2d)
  threshold_freq=0.5*(n1+n2-2)
  for(i in 0:(n1-1)){
    for(j in 0:(n2-1)){
      if(i+j < threshold_freq){
        maf_2dsfs[i+1,j+1]=maf_2dsfs[i+1,j+1]+derived_sfs2d[i+1,j+1]
      }
      else if(i+j == threshold_freq){
        maf_2dsfs[i+1,j+1]=maf_2dsfs[i+1,j+1]+(0.5*derived_sfs2d[i+1,j+1]+0.5*derived_sfs2d[n1-i,n2-j])
      }
      else{
        maf_2dsfs[n1-i,n2-j]=maf_2dsfs[n1-i,n2-j]+derived_sfs2d[i+1,j+1]
      }
    }
  }
  maf_2dsfs
}

# Apply to read 2D-SFS
b=derived2maf(a)
b[1,1] <- as.numeric(0)

# Round if observed SFS
# if(grepl("obs",suffix)){b=round(b)}

# Print new, folded 2D-SFS
cat(paste("1 observations","\n",
          paste(c("",colnames(b)),collapse="\t"),"\n",sep=""),file=paste(outdir,base,"_jointMAFpop",pops,".",suffix,sep=""))
suppressWarnings(write.table(b,paste(outdir,base,"_jointMAFpop",pops,".",suffix,sep=""),col.names=F,
                             quote=F,sep="\t",append=T))

# Prepare the data.frame for ggplot2
data_plot <- melt(b)
#data_plot$pop1 <- rep(rownames(a),length(colnames(a)))
data_plot$pop0 <- factor(data_plot$X2, levels = colnames(b))
data_plot$pop1 <- factor(data_plot$X1, levels = rownames(b))
data_plot$freq <- data_plot$value

# Plot sfs as heatmap
suppressWarnings(pdf(paste(outdir,base,"_jointMAFpop",pops,".pdf",sep="")))
ggplot(data_plot, aes(x=pop0,y=pop1, fill=freq)) + 
  geom_tile()+
  theme_bw()+
  #scale_fill_distiller(palette = 'Spectral')+
  scale_fill_gradientn(colours = c("white", "cyan4", "darkred"),
                       values = scales::rescale(c(-0.05, -0.03, 0, 0.6, 0.8))) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


