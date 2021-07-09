# Steps to perform differential expression analysis


First, we created a table of read counts using FeatureCounts v.1.6.3 (Liao et al. 2104). This step takes as input the individual bam files output by STAR.
```
 
 featureCounts \
 -a /path/to/genome_annotation_file.gff \
 -o counts.txt \
 -g ID \
 -t gene \
 -G /path/to/genome.fasta \
 --tmpDir /path/to/tmpDir \
 *.bam
 
  ```
The table of read counts output by FeatureCounts was further processed using the following command to keep only the first column (genes) and the columns with the individual samples, and remove leftovers of previous steps from the samples IDs.
```
grep -v "^#" counts.txt | \
cut -f 1,7- | \
sed 's/.map2.sortedByCoord.out.bam//g' > counts.edit.txt
```
This table was used as input matrix to perform differential expression analyses in R. 
The scripts DEanalysis.Rmd, DEfunctions.R and DEanalysis_dataInspect.Rmd were used for this task.
