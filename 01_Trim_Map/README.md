# Trimming and mapping of RNA-Seq fastq files

Trim .fq files using Trimmomatic v. 0.36 (Bolger et al. 2014)
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
 Mapping to reference genome using STAR v. 2.6.0c (Dobin et al. 2013) first and second pass mode 
 ```
 for file in *.fq; \
 do \
 nice -n 19 \
 STAR \
 --genomeDir /path/to/genome/dir \
 --readFilesIn ./$file \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix /path/to/1pass/${file/.fq/.map1} \
 --outFilterMismatchNoverLmax 0.04 \
 --limitBAMsortRAM 5617609531 \
 --runThreadN 12; \
 done
 
 nice -n 19 \
 STAR \   \
 --runMode genomeGenerate\ 
 --runThreadN 4 \
 --genomeDir /path/to/genome/dir \
 --genomeFastaFiles /path/to/genome.scf.fasta \
 --genomeSAindexNbases 12 \
 --sjdbOverhang 99 \
 --limitGenomeGenerateRAM 53000000000
 
 for file in *.fq; \
 do nice -n 19 \
 STAR \
 --genomeDir /path/to/genome/dir \
 --readFilesIn ./$file \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix /path/to/2pass/${file/.fq/.map2} \
 --outFilterMismatchNoverLmax 0.04 \
 --limitBAMsortRAM 12000000000 \
 --runThreadN 12 \
 --sjdbFileChrStartEnd mapping/1pass/*.out.tab; \
 done
 
```
