#!/bin/bash

#================================================================
#               MEGAHIT METAGENOME ASSEMBLY SCRIPT
#
# This script first subsamples paired-end FASTQ data using seqkit,
# then runs MEGAHIT assembler on the subsampled files.
#
# This version corrects the seqkit command syntax.
#
#================================================================


# --- Step 1: Configuration ---
# --- PLEASE EDIT THESE VARIABLES TO MATCH YOUR SETUP ---

# Set the fraction for subsampling (e.g., 0.5 for 50%, 0.3 for 30%)
FRACTION=0.3

# Directory where your original FASTQ files are located
INPUT_DIR="/share/org/YZWL/yzwl_lixg/tmp/megahit/data"

# Directory to store the new subsampled FASTQ files
# A unique name is generated based on the fraction
SUBSAMPLE_DIR="${INPUT_DIR}/subsampled_${FRACTION}"

# Directory for the final MEGAHIT assembly output
# A unique name is generated based on the fraction
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/tmp/megahit/megahit_output_${FRACTION}"

# Original input filenames
READ1="all_1.fq"
READ2="all_2.fq"

# Subsampled output filenames
READ1_SUB="all_1.sub.fq"
READ2_SUB="all_2.sub.fq"

# Number of CPU threads to use (this command gets all available cores)
THREADS=$(nproc)

# Random seed for sampling to ensure reproducibility and read pairing
SEED=100


# --- Step 2: Subsampling with seqkit ---

echo "====================================================="
echo "=== STEP 1: SUBSAMPLING PAIRED-END READS"
echo "====================================================="
echo "Input Directory: ${INPUT_DIR}"
echo "Output Directory for subsamples: ${SUBSAMPLE_DIR}"
echo "Sampling Fraction: ${FRACTION}"
echo "-----------------------------------------------------"

# Create the output directory for subsampled files if it does not exist
mkdir -p "${SUBSAMPLE_DIR}"

# --- Corrected seqkit commands ---
# We run the command separately for each file but use the same seed (-s)
# to ensure the same lines are sampled, thus preserving pairs.

echo "Subsampling R1 file: ${READ1}"
seqkit sample -p ${FRACTION} -s ${SEED} "${INPUT_DIR}/${READ1}" -o "${SUBSAMPLE_DIR}/${READ1_SUB}"
if [ $? -ne 0 ]; then
    echo "[ERROR] seqkit subsampling failed for R1. Terminating script."
    exit 1
fi

echo "Subsampling R2 file: ${READ2}"
seqkit sample -p ${FRACTION} -s ${SEED} "${INPUT_DIR}/${READ2}" -o "${SUBSAMPLE_DIR}/${READ2_SUB}"
if [ $? -ne 0 ]; then
    echo "[ERROR] seqkit subsampling failed for R2. Terminating script."
    exit 1
fi

echo "Subsampling completed successfully."
echo


# --- Step 3: Assembly with MEGAHIT ---

echo "====================================================="
echo "=== STEP 2: RUNNING MEGAHIT ASSEMBLY"
echo "====================================================="
echo "Using subsampled files from: ${SUBSAMPLE_DIR}"
echo "Final assembly output will be in: ${OUTPUT_DIR}"
echo "Number of threads: ${THREADS}"
echo "-----------------------------------------------------"

# Create the final MEGAHIT output directory
# mkdir -p "${OUTPUT_DIR}"

# Run MEGAHIT using the smaller, subsampled files as input
megahit -1 "${SUBSAMPLE_DIR}/${READ1_SUB}" \
        -2 "${SUBSAMPLE_DIR}/${READ2_SUB}" \
        -o "${OUTPUT_DIR}" \
        -t "${THREADS}"

# Check the exit status of the MEGAHIT command
if [ $? -eq 0 ]; then
    echo
    echo "====================================================="
    echo "MEGAHIT assembly finished successfully!"
    echo "Final contigs file: ${OUTPUT_DIR}/final.contigs.fa"
    echo "====================================================="
else
    echo
    echo "====================================================="
    echo "[ERROR] MEGAHIT assembly failed."
    echo "Check logs in ${OUTPUT_DIR} for details."
    echo "====================================================="
    exit 1
fi
