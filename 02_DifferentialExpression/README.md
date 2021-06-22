# RNA-seq analysis pipeline

```
Create table of read counts using FeatureCounts v. 1.6.3. Takes as input the individual bam files mapped using STAR.
```
 
 featureCounts \
 -a /path/to/genome_annotation_file.gff \
 -o output.txt \
 -g ID \
 -t gene \
 -G /path/to/genome.fasta \
 --tmpDir /path/to/tmpDir \
 *.bam
 
  ```
