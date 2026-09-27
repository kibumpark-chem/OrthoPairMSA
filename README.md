# OrthoPairMSA

Phylogeny-aware paired MSA generator for protein–protein complexes in co-folding models (AlphaFold 3, Protenix, Boltz).

Given two UniProt IDs, PhyloMSA builds **species-paired** MSAs: for each shared species it keeps one ortholog from each protein, clusters redundant pairs with CD-HIT, aligns each chain, and writes co-folding input configs that point at those custom MSAs.

## Pipeline

```text
UniProt IDs (prot1, prot2)
        │
        ▼
1. Fetch query sequences + ortholog sets
   (eggNOG | OrthoDB | InterPro)  →  *_interpro.fasta / *.tsv
        │
        ▼
2. Pair by species taxid (one rep per species)
   CD-HIT on concatenated pairs (default 90% id)
        │
        ▼
3. MAFFT align each chain  →  *_aligned.fasta
   reformat.pl → *_aligned.a3m
        │
        ▼
4. generate_inputs.py
   → {prot1}_{prot2}_af3.json
   → {prot1}_{prot2}_protenix.json
   → {prot1}_{prot2}_boltz.yaml
```

Ortholog sources:

| Source     | Meaning |
|------------|---------|
| `eggnog`   | UniProt eggNOG cross-ref → stream group members (default) |
| `orthodb`  | UniProt OrthoDB cross-ref → stream group members |
| `interpro` | InterPro family members via InterPro API |

AF3 / Protenix remap UniProt chain ids to single letters (`A`, `B`, …).

## Requirements

- Python ≥ 3.8 with: `biopython`, `pandas`, `numpy`, `ete3`, `requests`, `pyyaml`
- On `PATH`: `mafft`, `cd-hit`, `perl`
- Network access to UniProt (and InterPro if using that source)

No package install is required for the shell workflow — `generate_paired_msa.sh` runs scripts from `src/phylomsa/` directly. Use any env that has the deps above (e.g. `phylomsa` / `phylo`).

First-time `ete3` taxonomy DB (needed for species mapping):

```bash
python -c "from ete3 import NCBITaxa; NCBITaxa().update_taxonomy_database()"
```

## Quick start

```bash
bash generate_paired_msa.sh <prot1_uniprot> <prot2_uniprot> <output_dir> [source]
```

`source` is optional: `eggnog` (default), `orthodb`, or `interpro`.

### Example (yeast eEF2–Hgh1)

```bash
bash generate_paired_msa.sh P32324 P48362 ./eef2_hgh1 eggnog
bash generate_paired_msa.sh P32324 P48362 ./eef2_hgh1 orthodb
```

Then point AF3 at the generated JSON, e.g. `eef2_hgh1/P32324_P48362_af3.json`.

## Main outputs

For proteins `P1` and `P2` in `<output_dir>`:

| File | Description |
|------|-------------|
| `P1.fasta`, `P2.fasta` | Query sequences |
| `P1_interpro.fasta/.tsv`, `P2_...` | Ortholog set (+ taxonomy metadata) |
| `P1_paired_cdhit.fasta`, `P2_...` | Species-paired, CD-HIT–filtered FASTAs |
| `P1_aligned.a3m`, `P2_aligned.a3m` | Per-chain MSAs for co-folding |
| `P1_P2_af3.json` | AlphaFold 3 input (custom MSA) |
| `P1_P2_protenix.json` | Protenix input |
| `P1_P2_boltz.yaml` | Boltz input |

Logs: `*_mafft.log`. Intermediate: `concatenated_for_cdhit.fasta`, `concatenated_cdhit_out.fasta`.

## Advanced / modular use

Steps can be run separately from `src/phylomsa/`:

```bash
# 1) Fetch one protein’s ortholog set
python get_data.py --uniprot_id P32324 --output_dir OUT --source eggnog

# 2) Species-pair + CD-HIT (needs both proteins’ *_interpro.* files)
python create_filtered_fasta.py \
  --protein1_id P32324 --protein1_seq "<SEQ1>" \
  --protein2_id P48362 --protein2_seq "<SEQ2>" \
  --output_dir OUT

# 3) Align + convert to a3m (as in generate_paired_msa.sh)

# 4) Write co-folding configs
python generate_inputs.py \
  --job-name OUT/P32324_P48362 \
  --target all \
  --chain P32324 SEQ1 OUT/P32324_aligned.a3m none none \
  --chain P48362 SEQ2 OUT/P48362_aligned.a3m none none
```

`--target` can be `af3`, `protenix`, `boltz`, or `all`.

## Notes

- Pairing is by **species** (NCBI taxid via ete3), not by reciprocal best BLAST.
- CD-HIT defaults: identity `0.90`, word size `5`, spacer of 40 `X`s between concatenated sequences.
- Prefer `eggnog` / `orthodb` for true ortholog groups; `interpro` is broader (family-level) and can be noisier for pairing.
