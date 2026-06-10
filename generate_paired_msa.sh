#!/bin/bash

# This is the master script to generate the input for cofolding models with custom pairedMSA.
# Currently supporting AF3, Protenix, and Boltz.
# Usage: bash generate_paired_msa.sh <protein1 (UniprotID)> <protein2 (UniprotID)> <output_dir>

#TODO: Need to change this after test
cd /n/home01/kibumpark/group_folder/p18.3__ASR_AF/PhyloMSA/src/phylomsa

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: bash generate_paired_msa.sh <protein1(UniprotID)> <protein2(UniprotID)> <output_dir>"
    exit 1
fi

# Assign input arguments to variables
PROT1=$1
PROT2=$2
OUTPUT_DIR=$3

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Fetch sequences for the two proteins using UniProt IDs
for PROT in "$PROT1" "$PROT2"; do
    if python get_data.py --uniprot_id "$PROT" --output_dir "$OUTPUT_DIR"; then
        echo "Sequences fetched successfully for $PROT."
    else
        echo "Error fetching sequences for $PROT."
        exit 1
    fi
done

# Extract the pure sequence strings from the fetched FASTA files
SEQUENCE1=$(grep -v "^>" "$OUTPUT_DIR/${PROT1}.fasta" | tr -d '\n')
SEQUENCE2=$(grep -v "^>" "$OUTPUT_DIR/${PROT2}.fasta" | tr -d '\n')

# Filter based on phylogenetic prior and CD-HIT clustering
python create_filtered_fasta.py \
    --protein1_id "$PROT1" --protein1_seq "$SEQUENCE1" \
    --protein2_id "$PROT2" --protein2_seq "$SEQUENCE2" \
    --output_dir "$OUTPUT_DIR"

# MAFFT alignment of the filtered sequences
for PROT in "$PROT1" "$PROT2"; do
    if [ ! -f "$OUTPUT_DIR/${PROT}_paired_cdhit.fasta" ]; then
        echo "Error: Filtered sequences file not found for $PROT."
        exit 1
    fi
    if mafft --auto "$OUTPUT_DIR/${PROT}_paired_cdhit.fasta" > "$OUTPUT_DIR/${PROT}_aligned.fasta"; then
        echo "MAFFT alignment completed successfully for $PROT."
    else
        echo "Error during MAFFT alignment for $PROT."
        exit 1
    fi
done

# Convert aligned FASTA to a3m format (removing gaps and lowercase letters)
for PROT in "$PROT1" "$PROT2"; do
    if [ ! -f "$OUTPUT_DIR/${PROT}_aligned.fasta" ]; then
        echo "Error: Aligned FASTA file not found for $PROT."
        exit 1
    fi
    perl /n/home01/kibumpark/group_folder/p18.3__ASR_AF/PhyloMSA/scripts/reformat.pl fas a3m "$OUTPUT_DIR/${PROT}_aligned.fasta" "$OUTPUT_DIR/${PROT}_aligned.a3m"
    if [ $? -eq 0 ]; then
        echo "Conversion to a3m format completed successfully for $PROT."
    else
        echo "Error during conversion to a3m format for $PROT."
        exit 1
    fi
done

# Generate input configs using the Python script (with line continuations)
if python generate_inputs.py \
    --job-name "$OUTPUT_DIR/${PROT1}_${PROT2}" \
    --target all \
    --chain "$PROT1" "$SEQUENCE1" "$OUTPUT_DIR/${PROT1}_aligned.a3m" none none \
    --chain "$PROT2" "$SEQUENCE2" "$OUTPUT_DIR/${PROT2}_aligned.a3m" none none; then
    echo "Configurations generated successfully for $PROT1 and $PROT2."
else
    echo "Error generating configurations for $PROT1 and $PROT2."
    exit 1
fi