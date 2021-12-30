
# Pipeline to compute by gene Fst using ANGSD 

Angsd Saf files were computed for both ecotype populations from pair 1 and 3.
```
for counter in 1 3; \
do \
for i in list*$counter; \
do \
angsd \
-b $i \
-anc ../../refGenome/genome.scf.fasta \
-out ../review_Fst/${i/list/pop} \
-dosaf 1 \
-GL 2 \
-minQ 20 \
-minMapQ 30 \
-minInd 8 \
-P 8 \
-skipTriallelic 1 \
-doMajorMinor 1 \
-doCounts 1 \
-setMinDepthInd 10 \
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
Prepare Fst for easy window analysis and output the Bhatia estimator
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
Print Fst files
```
/path/to/realSFS fst print P1V1.fst.idx > P1V1.fst.txt
/path/to/realSFS fst print P3V3.fst.idx > P3V3.fst.txt
```
