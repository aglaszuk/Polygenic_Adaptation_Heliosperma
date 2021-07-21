# Script to filter unlinked positions from a chromosome + position list
# The input will be a tab separated contig + position list 
# The output is a subset of this list including only SNPs that have a minimum distance (e.g. 10 Kb) between them

library(dplyr)

# Load list
data <- read.table("FourFoldDeg_ContigPos.txt")
data <- data.frame("col1" = data[,1], "col2" = as.numeric(data[,2])) 

val <- ave(data$col2, data$col1, FUN=function(x) c(0, diff(x)))
pos0 <- data.frame("col1" = data$col1, "pos_s" = as.character(apply(data, 1, function(x) paste0(x[1],"_",x[2]))))
pos0$pos_s <- as.character(pos0$pos_s)
pos0["pos_e"] <- NA
pos0$pos_e[2:nrow(pos0)] <- as.character(pos0$pos_s[1:(nrow(pos0)-1)])
equals <- apply(pos0, 1, function(x) strsplit(x[3],"_")[[1]][1] == strsplit(x[2],"_")[[1]][1])
equals[is.na(equals)] <- FALSE
val <- val[equals]
pos0 <- pos0[equals, 2:3]
pos <- as.character(apply(pos0, 1, function(x) paste0(x[1],"-",x[2])))


out <- as.data.frame(cbind(pos,val)) 
keep <- as.numeric(as.character(out[,2]))>=10000 #windowsize: 10 kb (change in case)
out<-out[keep,] #filter out all values that are not at least 10000 bp distant
out <- na.exclude(out)
write.table(out,file = "SNPs_windows_out_marta",row.names = F,col.names = F,quote = F, sep = "\t")
dev.off()
