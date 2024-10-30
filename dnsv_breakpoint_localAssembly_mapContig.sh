#!/bin/bash

# Global variables
REGION="$1"
INFILE="$2"
OUTDIR="$3"
GENOME="$4"

# Extract discordant reads (those where both mates map to different locations)
samtools view -b -h \
-F 1294 \
-o ${OUTDIR}/${REGION}_discordant_reads.bam \
${INFILE} \
${REGION}

# Extract split read pairs (reads with CIGAR string containing 'S' or 'H', indicating soft-clipping or hard-clipping)
samtools view -h -f 0x2 -F 0x100 ${INFILE} ${REGION} | \
awk '($1 ~ /^@/ || $6 ~ /[SH]/) {print $1}' > ${OUTDIR}/${REGION}_clipped_read_names.txt

samtools view -h -f 0x2 -F 0x100 ${INFILE} ${REGION} | \
grep -Ff ${OUTDIR}/${REGION}_clipped_read_names.txt | \
samtools view -h -b -o ${OUTDIR}/${REGION}_split_reads.bam -

# EXTRACT READS ADJACENT TO THE SV REGION
contig=$(echo $REGION | cut -d ':' -f1)
positions=$(echo $REGION | cut -d ':' -f2)
start_pos=$(echo $positions | cut -d '-' -f1)
end_pos=$(echo $positions | cut -d '-' -f2)

extend_start_start=$(($start_pos - 5000))
extend_start_end=$start_pos
extend_start_region="${contig}:${extend_start_start}-${extend_start_end}"

extend_end_start=$end_pos
extend_end_end=$(($end_pos + 5000))
extend_end_region="${contig}:${extend_end_start}-${extend_end_end}"

echo $extend_start_region
echo $extend_end_region

samtools view -h -f 0x2 -F 0x100 \
-o ${OUTDIR}/${REGION}_start_adjacent_extended_reads.bam \
${INFILE} \
${extend_start_region}

samtools view -h -f 0x2 -F 0x100 \
-o ${OUTDIR}/${REGION}_end_adjacent_extended_reads.bam \
${INFILE} \
${extend_end_region}

# Extract reads with unmapped pair (one read is mapped, the other is unmapped)
samtools view -b -h \
-f 8 \
-F 4 \
-o ${OUTDIR}/${REGION}_unmapped_pair_reads.bam \
${INFILE} \
${REGION}

# Extract split reads (reads with CIGAR string containing 'S' or 'H', indicating soft-clipping or hard-clipping)
#samtools view -h ${INFILE} ${REGION} | \
#awk '{if($1 ~ /^@/ || $6 ~ /S|H/) print $0}' | \
#samtools view -b -o ${OUTDIR}/split_reads.bam

# Merge all three bams to one
samtools merge -f \
-o ${OUTDIR}/${REGION}_merged_reads.bam \
${OUTDIR}/${REGION}_discordant_reads.bam \
${OUTDIR}/${REGION}_unmapped_pair_reads.bam \
${OUTDIR}/${REGION}_split_reads.bam \
${OUTDIR}/${REGION}_start_adjacent_extended_reads.bam \
${OUTDIR}/${REGION}_end_adjacent_extended_reads.bam


# Indexing for visualization
samtools index ${OUTDIR}/${REGION}_merged_reads.bam

# Sort the merged bam
samtools sort -n \
-o ${OUTDIR}/${REGION}_name_sorted_merged_reads.bam \
${OUTDIR}/${REGION}_merged_reads.bam

# Bam to fastq
samtools fastq \
${OUTDIR}/${REGION}_name_sorted_merged_reads.bam \
-1 ${OUTDIR}/${REGION}_name_sorted_merged_reads_1.fastq \
-2 ${OUTDIR}/${REGION}_name_sorted_merged_reads_2.fastq \
-s ${OUTDIR}/${REGION}_name_sorted_merged_reads_single.fastq \
-0 ${OUTDIR}/${REGION}_name_sorted_merged_reads_non.fastq


# PERFORM LOCAL ASSEMBLY USING VELVET
cd $OUTDIR
/home/hyunwoo/programs/velvet/velvet-1.2.10/velveth ${OUTDIR} 31 -fastq -separate ${OUTDIR}/${REGION}_name_sorted_merged_reads_1.fastq ${OUTDIR}/${REGION}_name_sorted_merged_reads_2.fastq
/home/hyunwoo/programs/velvet/velvet-1.2.10/velvetg ${OUTDIR} -exp_cov auto -cov_cutoff auto

# ALIGN THE ASEMBLED CONTIG ON TO THE REFERENCE
/home/hyunwoo/programs/minimap2/minimap2/minimap2 -t 30 -a ${GENOME} ${OUTDIR}/contigs.fa | samtools sort -o ${REGION}_local_assembly_alignment.bam

samtools index ${OUTDIR}/${REGION}_local_assembly_alignment.bam
