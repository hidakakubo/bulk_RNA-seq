#!/bin/bash

SECONDS=0

# change working directory
cd /Users/hidaka/Desktop/bulk_kinoshita/mapping


# make sample name file
# サンプルによって作り方が違う
cd data
ls | awk -F'_' '{if(NR%2==0){print $1"_"$2"_"$3}}' > sample_name.txt
cd ..


# fastp
mkdir -p fastp_outputs
for SAMPLE in `cat data/sample_name.txt`
do
fastp \
-i data/${SAMPLE}_R1.fastq.gz \
-I data/${SAMPLE}_R2.fastq.gz \
-3 \
-o fastp_outputs/out_${SAMPLE}_R1.fastq.gz \
-O fastp_outputs/out_${SAMPLE}_R2.fastq.gz \
-h fastp_outputs/sample_${SAMPLE}.html \
-j fastp_outputs/sample_${SAMPLE}.json \
-q 15 \
-n 10 \
-t 1 \
-T 1 \
-w 6
done

echo "fastp finished running!"


# HISAT2でMappng
# リファレンスゲノムのインデックスファイルの作成
hisat2-build HISAT2/refdata-gex-mm10-2020-A_snca.fasta HISAT2/genome_index

# 実行
mkdir -p HISAT2_outputs
for SAMPLE in `cat data/sample_name.txt`
do
hisat2 -x HISAT2/genome_index \
-1 fastp_outputs/out_${SAMPLE}_R1.fastq.gz \
-2 fastp_outputs/out_${SAMPLE}_R2.fastq.gz \
-k 3 \
-p 6 \
|samtools sort -@ 6 -O BAM - > HISAT2_outputs/${SAMPLE}.sort.bam && \
samtools index -@ 6 HISAT2_outputs/${SAMPLE}.sort.bam
done

echo "HISAT2 finished running!"


# stringtieでQuantification
mkdir -p stringtie_outputs
for SAMPLE in `cat data/sample_name.txt`
do
stringtie HISAT2_outputs/${SAMPLE}.sort.bam \
-e \
-p 6 \
-o stringtie_outputs/${SAMPLE}.gtf \
-G HISAT2/refdata-gex-mm10-2020-A_snca.gtf
done

echo "StringTie finished running!"


# make sample name path file
cd stringtie_outputs
ls | awk -F'.' -v P="$(pwd)" '{print $1" "P"/"$0}' > ../../DEG_analysis/sample_name_path.txt
cd ..



duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."