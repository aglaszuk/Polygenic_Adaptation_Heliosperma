#! /usr/bin/env Rscript

# Author: Aglaia Szukala, aglaia.szukala [at] univie.ac.at
# Title: fstByGene.r
# Usage: fstByGene.r angsd.fst.txt gene.choord.bed out_file
# angsd.fst.txt is the angsd output file with per site Fst 
# gene.choord.bed is a bed file with genes coordinates of all gene features in the genome 

args<-commandArgs(trailingOnly=T)
fst = read.table(args[1], header=F)
genes = read.table(args[2], header=F)
out_prefix <- args[3]
colnames(fst) <- c("chr","pos","alpha","alpha+beta")
colnames(genes) <- c("chr","start","end")
genes$fst <- NA

for (i in 1:nrow(genes))
{
  sites = fst[fst$chr %in% genes[i,1],] 
  sites = sites[sites$pos >= genes[i,2] & sites$pos <= genes[i,3],]
  genes$fst[i] = sum(sites$alpha)/sum(sites$`alpha+beta`)
}

write.table(genes,file=paste0(out_prefix,"_Genes_fst.txt"), sep="\t",
            quote = F)
