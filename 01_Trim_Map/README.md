# Trimming and mapping of RNA-Seq fastq files

To trim adapters from individual fastq files and perform quality trimming we used Trimmomatic v. 0.36 (Bolger et al. 2014).
```
for file in *.fq; \
do \
nice -n 19 java -Xmx32g \
-jar trimmomatic-0.36.jar \
SE \
-threads 8 \
-phred33 \
$file ${file/.fq/.trim.fq} \
ILLUMINACLIP:/path/to/adapters/TruSeq3-SE.fa:2:30:10 \
LEADING:13 \
TRAILING:13 \
SLIDINGWINDOW:4:15 \
MINLEN:36; \
done
```
Mapping of each individual reads to the reference genome was performed using STAR v. 2.6.0c (Dobin et al. 2013) in the first plus second pass mode. Previous to mapping, the reference genome was indexed using STAR. 
 ```
nice -n 19 \
STAR \
--runMode genomeGenerate \
--genomeDir /path/to/genome/dir  
--genomeFastaFiles /path/to/genome.scf.fasta \
--genomeSAindexNbases 14 \
--genomeChrBinNbits 14
```
Then first pass mapping was performed. De novo information on the intron-exon junctions is collected during this step and used to generate a new genome index for the 2nd pass mapping.
 ```
for file in *.fq; \
do nice -n 19 \
STAR \
--genomeDir /path/to/genome/dir \
--readFilesIn ./$file \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /path/to/1pass/${file/.fq/.map1} \
--outFilterMismatchNoverLmax 0.04 \
--limitBAMsortRAM 5617609531; \
done
```
A new genome index for the second pass mapping was generated and second pass mapping performed. 
 ```
nice -n 19 \
STAR \   \
--runMode genomeGenerate\ 
--genomeDir /path/to/genome/dir \
--genomeFastaFiles /path/to/genome.scf.fasta \
--genomeSAindexNbases 14 \
--genomeChrBinNbits 14 \
--sjdbOverhang 99 \
--sjdbFileChrStartEnd mapping/1pass/*.out.tab
 
for file in *.fq; \
do nice -n 19 \
STAR \
--genomeDir /path/to/genome/dir \
--readFilesIn ./$file \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /path/to/2pass/${file/.fq/.map2} \
--outFilterMismatchNoverLmax 0.04 \
--limitBAMsortRAM 12000000000 \
--sjdbFileChrStartEnd mapping/1pass/*.out.tab; \
done
 
```
