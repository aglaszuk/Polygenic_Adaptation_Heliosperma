
# Estimate genotype likelihoods and analyze genetic diversity and structure

We used a the following pipeline to estimate genotype likelihoods using ANGSD v.0.931 (Korneliussen et al. 2014).
Briefly, after mapping individual files, we added read groups for each individual, sorted by coordinate and removed duplicates using PicardTools v.2.9.2 (https://broadinstitute.github.io/picard/). 
```
for file in *.bam; \
do \
java -Xmx32g \
-jar picard.jar \
AddOrReplaceReadGroups \
I=$file \
O=./sorted/${file/.bam/.sorted.bam} \
SORT_ORDER=coordinate \
RGID=${file/.bam/} \
RGLB=${file/.bam/} \
RGPL=illumina \
RGPU=machine \
RGSM=${file/.bam/}; \
done

cd sorted/

for file in *.bam; \
do \
java -Xmx32g \
-jar picard.jar \
MarkDuplicates \
I=$file \
O=./deduplicated/${file/.sorted.bam/.dedupli.bam} \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
M=$file.metrics; \
done
```
We then used PicardTools to create a reference dictionary and samtools to index the reference genome.
```
java -Xmx32g \
-jar picard.jar \
CreateSequenceDictionary \
R=/path/to/genome.fasta \
O=/path/to/genome/directory/genome.dict

cd /path/to/genome/directory/

samtools faidx genome.scf.fasta 
```
Then individual bam files were processed using GATK v.3.7.0 function IndelRealigner to locally improve read alignments around indels. This process included a first step to split reads that contain Ns in their cigar string and reassign mapping qualities.
```
mkdir splitted

for file in ./*.bam; \
do \
java -Xmx32g \
-jar GenomeAnalysisTK.jar \
-T SplitNCigarReads \
-R /path/to/genome.fasta \
-I ./$file \
-o ./splitted/${file/.dedupli.bam/.split.bam} \
-rf ReassignOneMappingQuality \
-RMQF 255 \
-RMQT 60 \
-U ALLOW_N_CIGAR_READS; \
done
```
Then realignment around indels was performed.
```
for file in ./*split.bam; \
do \
java -Xmx32g \
-jar /GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /path/to/genome.fasta \
-I ./$file \
-o ./$file.intervals; \
done

for file in ./*split.bam; \
do \
java -Xmx32g \
-jar GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /path/to/genome.fasta \
-I ./$file -o ./${file/split.bam/.realigned.bam} \
-dcov 100 \
-targetIntervals ./$file.intervals; \
done
```
After realignment, all bam files were indexed using samtools.
```
for file in *.bam; \
do \
samtools index $file; \
done
```
ANGSD v.0.931 was then run using as input a filelist listing all bamfile. The -sites option was used to limit the analyses to asubset of sites at fourfold degenerate positions.

ls *realigned.bam > bam.filelist

angsd \
-bam bam.filelist \
-GL 2 \
-doMajorMinor 1 \
-doMaf 1 \
-SNP_pval 2e-6 \
-minMapQ 20 \
-minQ 20 \
-minInd 12 \
-minMaf 0.045 \
-nThreads 12 \
-doGlf 2 \
-out ./genolike \
-sites FFD_ContigPos.txt
```
