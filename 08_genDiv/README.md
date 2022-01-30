
# Compute pi and theta Watterson using ANGSD 

ANGSD genome-wide saf files were computed for alpine and montane ecotype populations from pair 1 and 3.
```
for i in list*; \
do \
angsd \
-b $i \
-anc path/to/genome.scf.fasta \
-out ${i/list/pop} \
-dosaf 1 \
-GL 2 \
-minQ 20 \
-minMapQ 30 \
-P 8 \
-skipTriallelic 1 \
-doMajorMinor 1; \
done
```
Obtain the maximum likelihood estimate of the folded SFS and calculate thetas for each site. For Tajima's D estimates we used only FFD sites by adding the -sites option.
```
for i in *.saf.idx; \
do \
realSFS \
$i \
-P 24 \
-fold 1 > ${i/.saf.idx/.sfs}; \
done

for i in *.sfs; \
do \
realSFS \
saf2theta \ 
${i/.sfs/.saf.idx} \
-outname ${i/.sfs/} \
-sfs $i \
-fold 1; \
done
```
Perform sliding window analysis
```
for i in *thetas.idx; \
do \
./misc/thetaStat \
do_stat $i \
-win 50000 \
-step 10000 \
-outnames $i.thetasWindow.gz; \
done
```
Each .pestPG file was modified to have the first and last position with data in each window in a own column.
```
for file in *pestPG; \
do \
sed 's/)(/       /g' $file | \
sed 's/)//g' | \
sed s'/(//g' | \
tr ',' '\t' > ${file/.thetas.idx.thetasWindow.gz.pestPG/.edit.pestPG}; \
done
```
The resulting files were processed in R to obtain the weighted mean estimates and standard deviations
```
analyze <- function(filename){
  # Read input pestPG file
  data <- read.table(file = filename, header=F) 
  colnames(data) <- c("indexStart","indexStop","firstPos_withData","lastPos_withData","WinStart","WinStop",
                      "Chr","WinCenter","tW","tP","tF","tH","tL","Tajima",
                      "fuf","fud","fayh","zeng","nSites")
  data$totalSites <- data$lastPos_withData-data$firstPos_withData
  # Exclude NA values
  data <- na.omit(data)
  # Divide each tP and tW value by the number of sites in the window
  data$tPWeight <- data$tP/data$totalSites
  data$tWWeight <- data$tW/data$totalSites
  # Compute mean values and standard deviation
  m_tP <- mean(data$tPWeight)
  sd_tP <- sd(data$tPWeight)
  m_tW <- mean(data$tWWeight)
  sd_tW <- sd(data$tWWeight)
  
  dat_out <- c(m_tP,m_tW,sd_tP,sd_tW)
  out <- as.data.frame(matrix(data = dat_out, nrow = 2, ncol = 2))
  colnames(out) <- c("estimate","sd")
  rownames(out) <- c("pi", "thetaW")
  return(out)
}

par(mfrow=c(4,4))
p1 <- analyze("popP1.edit.pestPG")
v1 <- analyze("popV1.edit.pestPG")
p3 <- analyze("popP3.edit.pestPG")
v3 <- analyze("popV3.edit.pestPG")
p4 <- analyze("popP4.edit.pestPG")
v4 <- analyze("popV4.edit.pestPG")
p5 <- analyze("popP5.edit.pestPG")
v5 <- analyze("popV5.edit.pestPG")
```
