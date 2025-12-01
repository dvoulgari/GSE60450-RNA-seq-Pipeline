#!/bin/bash

# Creating necessary directories
BASE_DIR=$(pwd)
INDEX_DIR="$BASE_DIR/1_Index/index"
FASTQ_DIR="FASTQ_FILES/sra"
FASTQC_DIR="$BASE_DIR/2_FASTQC_FILES"
TRIMMED_DIR="$BASE_DIR/3_TRIMMED_DATA"
FASTQC_TRIMMED_DIR="$BASE_DIR/4_FASTQC_TRIMMED"
ALIGNMENT_DIR="$BASE_DIR/5_Alignment_outputs"

# Directory containing the FASTQ files.
FASTQ_DIR="$BASE_DIR/FASTQ_FILES/sra"

# Files of genome and index. Using absolute paths here for reliability.
GTF_FILE="$BASE_DIR/1_Index/Mus_musculus.GRCm39.114.gtf"
GENOME_FA="$BASE_DIR/1_Index/Mus_musculus.GRCm39.dna.primary_assembly.fa"

# Samples list
SAMPLES=(
"SRR1552444.fastq.gz:Sample44"
"SRR1552445.fastq.gz:Sample45"
"SRR1552446.fastq.gz:Sample46"
"SRR1552447.fastq.gz:Sample47"
"SRR1552448.fastq.gz:Sample48"
"SRR1552449.fastq.gz:Sample49"
"SRR1552450.fastq.gz:Sample50"
"SRR1552451.fastq.gz:Sample51"
"SRR1552452.fastq.gz:Sample52"
"SRR1552453.fastq.gz:Sample53"
"SRR1552454.fastq.gz:Sample54"
"SRR1552455.fastq.gz:Sample55"
)

# --- Function Definitions ---

process_step() {
    local step="$1"

    # Create output directories if they don't exist
    mkdir -p "$FASTQC_DIR" "$TRIMMED_DIR" "$FASTQC_TRIMMED_DIR" "$ALIGNMENT_DIR"

    for sample in "${SAMPLES[@]}"; do
        fastq_file="${sample%:*}"
        sample_name="${sample##*:}"

        echo "Processing $step for $sample_name..."
        case "$step" in
        "fastq_qc")
            echo "Running FASTQC for ${sample_name}..."
            fastqc "$FASTQ_DIR/$fastq_file" -o "$FASTQC_DIR" --threads 16
            ;;

        "trim")
            echo "Running Cutadapt for ${sample_name}..."
            cutadapt -m 1 -o "$TRIMMED_DIR/${sample_name}_trimmed.fastq.gz" \
            "$FASTQ_DIR/$fastq_file"
            ;;

        "fastqc_trimmed")
            echo "Running FASTQC for trimmed data of ${sample_name}..."
            fastqc "$TRIMMED_DIR/${sample_name}_trimmed.fastq.gz" -o "$FASTQC_TRIMMED_DIR" --threads 16
            ;;

        "alignment")
            echo "Running STAR alignment for ${sample_name}..."
            STAR --runThreadN 16 \
            --runMode alignReads \
            --genomeDir "$INDEX_DIR" \
            --sjdbGTFfile "$GTF_FILE" \
            --readFilesCommand gunzip -c \
            --readFilesIn "$TRIMMED_DIR/${sample_name}_trimmed.fastq.gz" \
            --outFileNamePrefix "$ALIGNMENT_DIR/${sample_name}_" \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts
            ;;
        esac
    done
}


# Main Script
# Step 1: Indexing reference genome (only once)
echo "Building index for Mus Musculus GRCm39 genome..."
mkdir -p "$INDEX_DIR"

STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir "$INDEX_DIR" \
--genomeFastaFiles "$GENOME_FA" \
--sjdbGTFfile "$GTF_FILE"

echo "Genome built."

# Step 2: Quality Control 
echo "Creating FASTQC Reports..."
process_step "fastq_qc"

# Step 3: Trim the data 
echo "Trimming the data..."
process_step "trim"

# Step 4: Quality Control on Trimmed Data
echo "Creating FASTQC Reports for the Trimmed Data"
process_step "fastqc_trimmed"

# Step 5: Alignment
echo "Aligning the data to the reference genome..."
process_step "alignment"

echo "All samples processed."