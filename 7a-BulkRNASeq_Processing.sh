##Before starting, download raw genome and bulk RNA-Seq files from the NCBI SRA database. The accession numbers are: (add once paper is published)
##In addition, refer to the documentation for each of the following software packages for installation instructions of the versions used in this analysis.
#trimmomatic/0.39
#fastqc/0.12.1
#hisat2/2.2.1
#samtools/1.19.2
#htseq/2.0.3

#First step: trim the reads using trimmomatic

module load trimmomatic/0.39

Loop over all files

or sample in `ls *_R1.fastq.gz`
do
prefix=$(basename $sample "_R1.fastq.gz")

trimmomatic PE -phred33 ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz ${prefix}_1P.fastq ${prefix}_1U.fastq ${prefix}_2P.fastq ${prefix}_2U.fastq ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done

# Generate fastqc files for each paired output of trimmomatic

module load fastqc/0.12.1

fastqc *P.fastq

# Build reference index with hisat2 and align the RNA-Seq reads to the genome

module load hisat2/2.2.1

hisat2_extract_splice_sites.py ICE_V1_polished_MAKER2_AllModels.gtf > splicesites.tsv ## This command creates the slipce site file

hisat2_extract_exons.py ICE_V1_polished_MAKER2_AllModels.gtf > exons.tsv ## This command creates the exon file

hisat2-build -p 8 --ss splicesites.tsv --exon exons.tsv IcePlantAssembly_Polished.fasta IcePlant_index ## This command created the genome indexing

mkdir RNA_Alingment_Dir
cd RNA_Alingment_Dir

## Running Hisat2

for sample in `ls *_1P.fastq`
do

prefix=$(basename $sample "_1P.fastq")

hisat2 -p 8 --rg-id=${prefix} --rg PL:ILLUMINA -x IcePlant_index --dta --rna-strandness RF -1 ${prefix}_1P.fastq -2 ${prefix}_2P.fastq -S ${prefix}.sam

done

module load samtools/1.19.2

##Run Samtools to convert SAM files to BAM

for sample in `ls RNA_Alingment_Dir/*.sam`
do
prefix=$(basename $sample ".sam")

samtools sort -@ 8 -o ${prefix}.bam RNA_Alingment_Dir/${prefix}.sam

done

#Index bam files

for sample in *.bam
do
  samtools index $sample
done

##Run HTSEQ-COUNT to estimate gene expression:

module load htseq/2.0.3

for sample in `ls *.bam`
do
prefix=$(basename $sample ".bam")

htseq-count --format bam --order pos --mode union --stranded reverse --minaqual 1 --type exon --idattr gene_id ${prefix}.bam ICE_V1_polished_MAKER2_AllModels.gtf > ${prefix}_gene.tsv

done

## Joining all files

mkdir FinalFiles

join -j 1 1_D21_8AM_C1_gene.tsv 2_D21_8AM_C2_gene.tsv > temp1.tsv
join -j 1 temp1.tsv 3_D21_8AM_C3_gene.tsv > temp2.tsv
join -j 1 temp2.tsv 4_D21_8AM_S1_gene.tsv > temp3.tsv
join -j 1 temp3.tsv 5_D21_8AM_S2_gene.tsv > temp4.tsv
join -j 1 temp4.tsv 6_D21_8AM_S3_gene.tsv > temp5.tsv
join -j 1 temp5.tsv 7_D21_12PM_C1_gene.tsv > temp6.tsv
join -j 1 temp6.tsv 8_D21_12PM_C2_gene.tsv > temp7.tsv
join -j 1 temp7.tsv 9_D21_12PM_C3_gene.tsv > temp8.tsv
join -j 1 temp8.tsv 10_D21_12PM_S1_gene.tsv > temp9.tsv
join -j 1 temp9.tsv 11_D21_12PM_S2_gene.tsv > temp10.tsv
join -j 1 temp10.tsv 12_D21_12PM_S3_gene.tsv > temp11.tsv
join -j 1 temp11.tsv 13_D21_8PM_C1_gene.tsv > temp12.tsv
join -j 1 temp12.tsv 14_D21_8PM_C2_gene.tsv > temp13.tsv
join -j 1 temp13.tsv 15_D21_8PM_C3_gene.tsv > temp14.tsv
join -j 1 temp14.tsv 16_D21_8PM_S1_gene.tsv > temp15.tsv
join -j 1 temp15.tsv 17_D21_8PM_S2_gene.tsv > temp16.tsv
join -j 1 temp16.tsv 18_D21_8PM_S3_gene.tsv > temp17.tsv
join -j 1 temp17.tsv 19_D21_0AM_C1_gene.tsv > temp18.tsv
join -j 1 temp18.tsv 20_D21_0AM_C2_gene.tsv > temp19.tsv
join -j 1 temp19.tsv 21_D21_0AM_C3_gene.tsv > temp20.tsv
join -j 1 temp20.tsv 22_D21_0AM_S1_gene.tsv > temp21.tsv
join -j 1 temp21.tsv 23_D21_0AM_S2_gene.tsv > temp22.tsv
join -j 1 temp22.tsv 24_D21_0AM_S3_gene.tsv > temp23.tsv
join -j 1 temp23.tsv 25_D21_4PM_C1_gene.tsv > temp24.tsv
join -j 1 temp24.tsv 26_D21_4PM_C2_gene.tsv > temp25.tsv
join -j 1 temp25.tsv 27_D21_4PM_C3_gene.tsv > temp26.tsv
join -j 1 temp26.tsv 28_D21_4PM_S1_gene.tsv > temp27.tsv
join -j 1 temp27.tsv 29_D21_4PM_S2_gene.tsv > temp28.tsv
join -j 1 temp28.tsv 30_D21_4PM_S3_gene.tsv > temp29.tsv
join -j 1 temp29.tsv 31_D21_4AM_C1_gene.tsv > temp30.tsv
join -j 1 temp30.tsv 32_D21_4AM_C2_gene.tsv > temp31.tsv
join -j 1 temp31.tsv 33_D21_4AM_C3_gene.tsv > temp32.tsv
join -j 1 temp32.tsv 34_D21_4AM_S1_gene.tsv > temp33.tsv
join -j 1 temp33.tsv 35_D21_4AM_S2_gene.tsv > temp34.tsv
join -j 1 temp34.tsv 36_D21_4AM_S3_gene.tsv > FinalFiles/gene_read_counts_table_D21.tsv

echo "GeneID 8C1 8C2 8C3 8S1 8S2 8S3 12C1 12C2 12C3 12S1 12S2 12S3 20C1 20C2 20C3 20S1 20S2 20S3 0C1 0C2 0C3 0S1 0S2 0S3 16C1 16C2 16C3 16S1 16S2 16S3 4C1 4C2 4C3 4S1 4S2 4S3" > header.txt
cat header.txt FinalFiles/gene_read_counts_table_D21.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > FinalFiles/gene_read_counts_table_D21_final.tsv


