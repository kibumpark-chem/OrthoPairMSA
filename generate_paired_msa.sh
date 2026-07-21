#!/bin/bash

# Generate cofolding inputs with custom paired MSA (AF3 / Protenix / Boltz).
# Usage:
#   bash generate_paired_msa.sh <protein1> <protein2> <output_dir> [source]
# source: eggnog (default) | orthodb | interpro

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHYLOMSA_SRC="${SCRIPT_DIR}/src/phylomsa"
REFORMAT_PL="${SCRIPT_DIR}/scripts/reformat.pl"

cd "$PHYLOMSA_SRC"

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: bash generate_paired_msa.sh <protein1(UniprotID)> <protein2(UniprotID)> <output_dir> [eggnog|orthodb|interpro]"
    exit 1
fi

PROT1=$1
PROT2=$2
OUTPUT_DIR=$3
SOURCE=${4:-eggnog}

case "$SOURCE" in
    eggnog|orthodb|interpro) ;;
    *)
        echo "Error: source must be one of: eggnog, orthodb, interpro (got: $SOURCE)"
        exit 1
        ;;
esac

echo "Using ortholog source: $SOURCE"

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

for PROT in "$PROT1" "$PROT2"; do
    if python get_data.py --uniprot_id "$PROT" --output_dir "$OUTPUT_DIR" --source "$SOURCE"; then
        echo "Sequences fetched successfully for $PROT."
    else
        echo "Error fetching sequences for $PROT."
        exit 1
    fi
done

SEQUENCE1=$(grep -v "^>" "$OUTPUT_DIR/${PROT1}.fasta" | tr -d '\n')
SEQUENCE2=$(grep -v "^>" "$OUTPUT_DIR/${PROT2}.fasta" | tr -d '\n')

python create_filtered_fasta.py \
    --protein1_id "$PROT1" --protein1_seq "$SEQUENCE1" \
    --protein2_id "$PROT2" --protein2_seq "$SEQUENCE2" \
    --output_dir "$OUTPUT_DIR"

for PROT in "$PROT1" "$PROT2"; do
    if [ ! -f "$OUTPUT_DIR/${PROT}_paired_cdhit.fasta" ]; then
        echo "Error: Filtered sequences file not found for $PROT."
        exit 1
    fi
    # Redirect stderr to a log file: some environments lack a usable /dev/stderr
    if mafft --auto "$OUTPUT_DIR/${PROT}_paired_cdhit.fasta" \
        > "$OUTPUT_DIR/${PROT}_aligned.fasta" \
        2> "$OUTPUT_DIR/${PROT}_mafft.log"; then
        echo "MAFFT alignment completed successfully for $PROT."
    else
        echo "Error during MAFFT alignment for $PROT. See $OUTPUT_DIR/${PROT}_mafft.log"
        exit 1
    fi
done

for PROT in "$PROT1" "$PROT2"; do
    if [ ! -f "$OUTPUT_DIR/${PROT}_aligned.fasta" ]; then
        echo "Error: Aligned FASTA file not found for $PROT."
        exit 1
    fi
    if perl "$REFORMAT_PL" fas a3m "$OUTPUT_DIR/${PROT}_aligned.fasta" "$OUTPUT_DIR/${PROT}_aligned.a3m"; then
        echo "Conversion to a3m format completed successfully for $PROT."
    else
        echo "Error during conversion to a3m format for $PROT."
        exit 1
    fi
done

# Job basename only (AF3 rejects path-like names / UniProt chain ids are remapped to A,B,...)
JOB_NAME="${PROT1}_${PROT2}"
if python generate_inputs.py \
    --job-name "$OUTPUT_DIR/${JOB_NAME}" \
    --target all \
    --chain "$PROT1" "$SEQUENCE1" "$OUTPUT_DIR/${PROT1}_aligned.a3m" none none \
    --chain "$PROT2" "$SEQUENCE2" "$OUTPUT_DIR/${PROT2}_aligned.a3m" none none; then
    echo "Configurations generated successfully for $PROT1 and $PROT2."
    # AF3/Protenix remap UniProt ids -> A,B; print a short mapping note
    echo "Chain mapping: A=${PROT1}, B=${PROT2}"
else
    echo "Error generating configurations for $PROT1 and $PROT2."
    exit 1
fi
