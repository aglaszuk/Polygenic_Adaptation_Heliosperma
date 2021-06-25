# Build a table of explanatory variables
makeEnv = function(eco, loc, data){
  env <- data.frame(cbind(eco,loc))
  rownames(env)=rownames(data)
  colnames(env)=c("ecotype","locality")
  return(env)
}

# Perform conditioned RDA analysis
helio_rda = function(dat, env){
  helio.rda <- rda(dat ~ ecotype + Condition(locality), data=env, scale=T, na.action=na.omit)
  summary(helio.rda, display=NULL)
  print(RsquareAdj(helio.rda))
  summary(eigenvals(helio.rda, model = "constrained"))
  summary(eigenvals(helio.rda, model = "unconstrained"))
  return(helio.rda)
}  

# Plot RDA
plot_rda = function(rda){
  plot(rda, main="cRDA of expression data", 
     xlab="RDA1 - 1.8%", ylab="PC1 - 52.11%")
  points(scores(helio.rda)$sites[1:12,],pch=17, cex=1.2, bg="black")
  points(scores(helio.rda)$sites[1:3,],pch=17, cex=1, col="#A71E34")
  points(scores(helio.rda)$sites[4:6,],pch=17,cex=1, col="#00798C")
  points(scores(helio.rda)$sites[7:9,],pch=17, cex=1, col="#F59700")
  points(scores(helio.rda)$sites[10:12,],pch=17, cex=1, col="#6E5DAC")
  points(scores(helio.rda)$sites[13:24,],pch=19, cex=1.2, bg="black")
  points(scores(helio.rda)$sites[13:15,],pch=19, cex=1, col="#EBADB5")
  points(scores(helio.rda)$sites[16:18,],pch=19, cex=1, col="#B7E1D3")
  points(scores(helio.rda)$sites[19:21,],pch=19, cex=1, col="#FFD085")
  points(scores(helio.rda)$sites[22:24,],pch=19, cex=1, col="#CEACEC")
  col = c("grey","black")
  ordihull(helio.rda, env$ecotype, col=col, lwd=1)
  rp <- vector('expression',1)
  rp[2] <- substitute(expression(italic(P) == valueB), 
                    list(valueB = format(0.18, digits = 2)))[2]
  legend("bottomleft",legend = rp, bty = 'n')
}

# Transform RDA scores into z-scores
rda2zscore = function(scores){
  rda1_score_exp <- read.table(scores, header = TRUE)
  hist(rda1_score_exp$RDA1)
  dat_scaled <- scale(rda1_score_exp$RDA1, center = TRUE, scale = F) #standardization: center the values around zero
  zrda <- (dat_scaled - mean(dat_scaled)) / sd(dat_scaled)
  hist(zrda, main="transcripts z-score distribution", xlab = "z-score", xlim = c(-4,4))
  zrda= matrix(zrda, ncol = 1, nrow = length(zrda))
  row.names(zrda) <- row.names(rda1_score_exp)
  colnames(zrda) <- "z-score"
  return(zrda)
}

# Detect outlier loci
outliers = function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Make dataframe including gene name and loadings
rda_out = function(candidates){
  cand <- cbind.data.frame(rep(1,times=length(candidates)), names(candidates), unname(candidates))
  colnames(cand) <- c("axis","snp","loading")
  write.table(cand, "candidateLoci_cRDA.txt", quote = F, sep="\t", row.names = T)
  return(cand)
}