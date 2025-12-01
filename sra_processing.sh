#!/bin/bash

echo -e "\n--- Processing SRA Files ---"

# Define the directory containing the  SRA files
SRA_DIR="/mnt/raid1/despvoulg/Desktop/FASTQ_FILES/sra/"

# Define the output directory for the FASTQ.GZ files
OUTPUT_DIR="/mnt/raid1/despvoulg/Desktop/FASTQ_FILES/sra/"

# Loop through all .sra files in the specified directory
for sra_file in "$SRA_DIR"/*.sra; do
  echo "Processing: $(basename "$sra_file")..."
  # Run fastq-dump for each .sra file
  fastq-dump --gzip --split-3 "$sra_file" --outdir "$OUTPUT_DIR"
  if [ $? -eq 0 ]; then
    echo "Successfully processed $(basename "$sra_file")"
  else
    echo "Error processing $(basename "$sra_file")"
  fi
done

echo "--- SRA File Processing Complete ---"
