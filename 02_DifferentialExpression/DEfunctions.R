# Perform data trimming
trimDat = function (counts) {
  data = read.table(counts,header= T, row.names = 1)
  print(paste0("Number of samples: ", dim(data)[2]))
  print(paste0("Starting number of genes: ", dim(data)[1]))
  data=data[apply(cpm(data),1,function(x){!(mean(x)<1)}),]
  print(paste0("SNumber of genes retaned after trimming: ", dim(data)[1]))
  return(data)
}

# Produce edgeR data object
makeObj = function (data){
  #define group variables
  groupGLM <- factor(c(rep("1.H",3),rep("3.H",3),rep("4.H",3), # Populations
                       rep("5.H",3),rep("1.L",3),rep("3.L",3),
                       rep("4.L",3),rep("5.L",3)))
  ecotype <- c(rep("H",12),rep("L",12)) # Ecotypes
  locality <- c(rep("1",3),rep("3",3),rep("4",3),rep("5",3), # Geographic localities
                rep("1",3),rep("3",3),rep("4",3),rep("5",3))
  y <- DGEList(counts=data,group=groupGLM) # list-based edgeR data object
  y$samples$locality=locality # add locality and ecotype information
  y$samples$ecotype=ecotype
  samples <- y$samples[,c(4,5,1,2,3)] # reorder columns
  y$samples=samples
  return(y)
}

# Estimate dispersion and fit GLM
fitGLM = function(obj) {
  groupGLM <- factor(c(rep("1.H",3),rep("3.H",3),rep("4.H",3), # Populations
                       rep("5.H",3),rep("1.L",3),rep("3.L",3),
                       rep("4.L",3),rep("5.L",3)))
  ModelDesign=model.matrix(~0+groupGLM)
  DGE=estimateDisp(obj,design = ModelDesign,robust = T) 
  GLM=glmFit(DGE,design = ModelDesign)
}

# Find DEGs
# LRT test to identify DE genes between ecotypes and plot as MA plot
LRT_test = function (contr, glm, col, l){
  lrt = glmLRT(glm,contrast = contr)
  res = lrt$table
  res$padj=p.adjust(res$PValue,method = "BH")
  DEG <- sum(res$padj < 0.05) 
  DEG_Hup <- sum(res$padj < 0.05 & res$logFC < 0) 
  DEG_Lup <- sum(res$padj < 0.05 & res$logFC > 0) 
  print(paste0("Number of DEG: ", DEG))
  print(paste0("Number of genes underexpressed in L: ", DEG_Hup))
  print(paste0("Number of genes overexpressed in L: ", DEG_Lup))
  # save tables with DEG before and after multiple testing correction
  write.table(row.names(res)[res$PValue < 0.05 & res$logFC < 0],
              file=paste0("DEG_lists/Hup_",l,"_pval0.05.txt"), quote = F, 
              row.names = F, col.names = F)
  write.table(row.names(res)[res$PValue < 0.05 & res$logFC > 0],
              file=paste0("DEG_lists/Lup_",l,"_pval0.05.txt"), quote = F,
              row.names = F, col.names = F)
  write.table(row.names(res)[res$padj < 0.05 & res$logFC < 0],
              file=paste0("DEG_lists/Hup_",l,"_pval0.05.txt"), quote = F, 
              row.names = F, col.names = F)
  write.table(row.names(res)[res$padj < 0.05 & res$logFC > 0],
              file=paste0("DEG_lists/Lup_",l,"_pval0.05.txt"), quote = F,
              row.names = F, col.names = F)
  # plot MA plots and save as pdf
  pdf(paste0("plots/MA_plots/Loc",l,".pdf"))
  plot(res$logCPM,res$logFC, col="cornsilk2",
       xlab="average expression (log(cpm))",ylab="log(FC)",
       cex.main = 2,cex.lab=1.55,
       main=paste0(DEG, " DEG (FDR<0.05) in locality ", l))
  points(res$logCPM[res$padj<0.05 & res$logFC>0],
         res$logFC[res$padj<0.05 & res$logFC>0],
         col=col,pch=19, cex=1.5)
  points(res$logCPM[res$padj<0.05 & res$logFC<0],
         res$logFC[res$padj<0.05 & res$logFC<0],
         col=col,pch=17, cex=1.45)
  legend("topright",
         legend=c(paste(DEG_Hup," DEG underexp in L"),
                  paste(DEG_Lup," DEG overexp in L")),
         pch=c(17, 19), cex=1.5, bty = "n")
  dev.off()
}

# Find DEGs between pops
# LRT test to identify DE genes between populations and plot as MA plot
LRT_test_loc = function (contr, glm, l1,l2){
  lrt = glmLRT(glm,contrast = contr)
  res = lrt$table
  res$padj=p.adjust(res$PValue,method = "BH")
  DEG <- sum(res$padj < 0.05) 
  DEG_1up <- sum(res$padj < 0.05 & res$logFC < 0) 
  DEG_2up <- sum(res$padj < 0.05 & res$logFC > 0) 
  print(paste0("Number of DEG: ", DEG))
  print(paste0("Number of genes underexpressed in ",l1,": ", DEG_1up))
  print(paste0("Number of genes overexpressed in ",l2,": ", DEG_2up))
  pdf(paste0("plots/MA_plots/Loc",l1,"_",l2,".pdf"))
  plot(res$logCPM,res$logFC, col="cornsilk2", cex=0.5,
       xlab="average expression (log(cpm))",ylab="log(FC)",
       cex.main = 1.2,cex.lab=1,
       main=paste0(DEG, " DEG (FDR<0.05) between pop ", l1," and ",l2))
  points(res$logCPM[res$padj<0.05 & res$logFC>0],
         res$logFC[res$padj<0.05 & res$logFC>0],
         col="lightcyan4",pch=8, cex=0.5)
  points(res$logCPM[res$padj<0.05 & res$logFC<0],
         res$logFC[res$padj<0.05 & res$logFC<0],
         col="lightcyan4",pch=9, cex=0.5)
  legend("topright", 
         legend=c(paste0(DEG_1up," DEG up in ",l1), 
                  paste0(DEG_2up," DEG overexp in ",l2)),
         pch=c(8, 9), cex=1, bty = "n")
  dev.off()
}

# Produce a data table to be used as input in GOplot
makeGOplotTable <- function(contr,glm){
  lrt = glmLRT(glm,contrast = contr)
  res = lrt$table
  res$FDR = p.adjust(res$PValue,method = "BH")
  colnames(res) <- c("logFC","logCPM","F","PValue","FDR")
  resGOplot <- res[res$FDR<0.05,]
  # resGOplot <- cbind(rownames(resGOplot), 
  #                    data.frame(resGOplot, row.names=NULL))
  # colnames(resGOplot)[1] <- "ID"
  # colnames(resGOplot)[3] <- "AveExpr"
  # resGOplot = resGOplot[,c(1,2)]
  # resGOplot_Hup = resGOplot[resGOplot$logFC<0,]
  # resGOplot_Lup = resGOplot[resGOplot$logFC>0,]
  #resGOplot<-return(resGOplot)
  write.table(resGOplot, 
              file = paste0("GOplot_tables/DEG_",lrt$comparison,".txt"),
              #quote = F,
              row.names = T,
              col.names = T)
}

# Find overlpas between ecotype pairs
findOverlapsPair = function (contr1, contr2, col1, col2, l1,l2){
  lrt1 = glmLRT(GLM,contrast = contr1)
  lrt2 = glmLRT(GLM,contrast = contr2)
  res1 = lrt1$table
  res2 = lrt2$table
  res1$padj=p.adjust(res1$PValue,method = "BH")
  res2$padj=p.adjust(res2$PValue,method = "BH")
  Hup_1=row.names(res1)[res1$padj<0.05 & res1$logFC<0]
  Lup_1=row.names(res1)[res1$padj<0.05 & res1$logFC>0]
  Hup_2=row.names(res2)[res2$padj<0.05 & res2$logFC<0]
  Lup_2=row.names(res2)[res2$padj<0.05 & res2$logFC>0]
  L_overlap <- sum(res1$padj<0.05 & res1$logFC>0 & 
                     res2$padj<0.05 & res2$logFC>0)
  H_overlap <- sum(res1$padj<0.05 & res1$logFC<0 &
                     res2$padj<0.05 & res2$logFC<0)
  print(paste0("Number of overlapping DEG overexpressed in L: ",L_overlap))
  print(paste0("Number of overlapping DEG underexpressed in L: ",H_overlap))
  write.table(Hup_1,file = paste0("Hup_",lrt1$comparison,".txt"),
              quote = F, row.names = F, col.names = F)
  write.table(Lup_1,file = paste0("Lup_",lrt1$comparison,".txt"),
              quote = F, row.names = F, col.names = F)
  write.table(Hup_2,file = paste0("Hup_",lrt2$comparison,".txt"),
              quote = F, row.names = F, col.names = F)
  write.table(Lup_2,file = paste0("Lup_",lrt2$comparison,".txt"),
              quote = F, row.names = F, col.names = F)
  plot(res1$logFC,res2$logFC,
       xlim=c(-11,11),ylim=c(-11,11),
       xlab=paste0("H vs L in ", l1," (logFC)"),ylab=paste0("H vs L in ",l2," (logFC)"), 
       main=paste0("Ecotype pair ",l1," and ",l2),cex=0.7)
  abline(h=0,v=0,col="grey20",lty=2)
  abline(a=0,b=1,col="grey20",lty=2)
  abline(a=0,b=-1,col="grey20",lty=2)
  points(res1$logFC[res1$padj>0.05&res2$padj>0.05],
         res2$logFC[res1$padj>0.05&res2$padj>0.05],col="cornsilk2",pch=19,cex=0.7)
  points(res1$logFC[res1$padj<0.05&res2$padj>0.05],
         res2$logFC[res1$padj<0.05&res2$padj>0.05],col=col1,pch=19)
  points(res1$logFC[res1$padj>0.05&res2$padj<0.05],
         res2$logFC[res1$padj>0.05&res2$padj<0.05],col=col2,pch=19,cex=0.7)
  points(res1$logFC[res1$padj<0.05&res2$padj<0.05],
         res2$logFC[res1$padj<0.05&res2$padj<0.05],col="indianred3",pch=19,cex=0.7)
  legend("bottomright", legend=c(paste0("DEG in ",l1), paste0("DEG in ",l2),paste0("DEG in ",l1," ",l2)),
         pch=c(19,19,19), col=c(col1,col2,"indianred3"), cex=0.8)
}
