
# Compute Fst by gene using ANGSD 

ANGSD genome-wide saf files were computed for alpine and montane ecotype populations from pair 1 and 3.
```
for counter in 1 3; \
do \
for i in list*$counter; \
do \
angsd \
-b $i \
-anc path/to/genome.scf.fasta \
-out path/to/outdir/${i/list/pop} \
-dosaf 1 \
-GL 2 \
-minQ 20 \
-minMapQ 30 \
-minInd 7 \
-P 8 \
-skipTriallelic 1 \
-doMajorMinor 1 \
-doCounts 1 \
-setMinDepthInd 4 \
-setMaxDepthInd 150; \
done; \
done
```
Compute folded site frequency spectrum
```
for counter in 1 3; \
do \
/path/to/realSFS \
popP$counter.saf.idx \ #P=H.pusillum (alpine ecotype)
popV$counter.saf.idx \ #V=H.veselskyi (montane ecotype)
-fold 1 \
-P 8 \
> popP$counter.popV$counter.ml; \
done
```
Prepare Fst for easy window analysis and output the Bhatia estimator (see [Bhatia et al (2013)](https://pubmed.ncbi.nlm.nih.gov/23861382/))
```
for counter in 1 3; \
do \
/path/to/realSFS \
fst index \
popP$counter.saf.idx \
popV$counter.saf.idx \
-sfs popP$counter.popV$counter.ml \
-fstout ./out/P$counterV$counter \
-whichFst 1; \
done
```
Print numerator and denominator of Bhatia Fst per site estimator 
```
/path/to/realSFS fst print P1V1.fst.idx > P1V1.fst.txt
/path/to/realSFS fst print P3V3.fst.idx > P3V3.fst.txt
```
Sum numerator and denominator across sites in a gene region and take the ratio of the sums to get the estimate. This can be done using the R script fstByGene.r 
```
Rscript fstByGene.r P1V1.fst.txt gene.choord.bed out_file_1
Rscript fstByGene.r P3V3.fst.txt gene.choord.bed out_file_3
```
The loci falling in the top 5% of the Fst distribution were defined as Fst outliers.
