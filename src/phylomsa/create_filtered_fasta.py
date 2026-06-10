#!/usr/bin/env python3
"""
Pair two InterPro-derived sequence sets by SPECIES (NCBI species taxid),
cluster pairs with CD-HIT on concatenated sequences, and write paired FASTAs.

Requires:
  - ete3 (NCBITaxa)
  - biopython
  - pandas
  - cd-hit on PATH

Expected files in output_dir:
  {protein_id}_interpro.tsv
  {protein_id}_interpro.fasta
"""

from ete3 import NCBITaxa
import pandas as pd
import os
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import argparse

ncbi = NCBITaxa()

UNIPROT_ACC_RE = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z0-9]{6,10})\b")

def normalize_acc(acc: str) -> str:
    """Normalize accession for matching: drop version (.1) and isoform (-2)."""
    if acc is None:
        return None
    acc = str(acc).strip()
    acc = acc.split()[0]
    acc = acc.split(".")[0]      # e.g. A0A123.1 -> A0A123
    acc = acc.split("-")[0]      # e.g. P12345-2 -> P12345
    return acc

def acc_from_fasta_record(record) -> str:
    """
    Robustly extract accession from FASTA headers.

    Handles:
      - UniProt:  sp|P12345|NAME   => P12345
      - UniProt:  tr|Q9XXXX|NAME   => Q9XXXX
      - Custom:   A0A0....|unreviewed|Some protein => A0A0....
      - Plain:    A0A0.... => A0A0....
    """
    rid = record.id.strip()

    if "|" in rid:
        parts = rid.split("|")
        # UniProt convention
        if parts[0] in ("sp", "tr") and len(parts) >= 2:
            return normalize_acc(parts[1])
        # Custom convention: accession is first token
        return normalize_acc(parts[0])

    # fallback
    return normalize_acc(rid)

def read_interpro_tsv(tsv_path: str, debug=False) -> pd.DataFrame:
    # Your API TSV is 8 columns (no header)
    tsv_cols = ["Accession", "Database", "Name", "Tax ID", "Scientific Name", "Length", "Entry Accession", "Locations"]

    # Use python engine + fixed 8 cols to be more tolerant of weird lines
    df = pd.read_csv(
        tsv_path,
        sep="\t",
        header=None,
        names=tsv_cols,
        usecols=list(range(8)),
        dtype=str,
        engine="python",
        on_bad_lines="skip",
        encoding="utf-8",
    )

    # If someone accidentally wrote a header row, drop it
    if not df.empty and str(df.iloc[0]["Accession"]).lower() == "accession":
        df = df.iloc[1:].copy()

    df["Accession_raw"] = df["Accession"]
    df["Accession_norm"] = df["Accession"].map(normalize_acc)

    df["Tax ID"] = pd.to_numeric(df["Tax ID"], errors="coerce").astype("Int64")
    df["Length"] = pd.to_numeric(df["Length"], errors="coerce").astype("Int64")

    if debug:
        print(f"[DEBUG] TSV rows read: {len(df)}")
        print("[DEBUG] TSV sample accessions:", df["Accession_raw"].head(5).tolist())
        print("[DEBUG] TSV sample norm accessions:", df["Accession_norm"].head(5).tolist())
        print("[DEBUG] TSV Tax ID non-null:", int(df["Tax ID"].notna().sum()))

    return df

def read_interpro_fasta(fasta_path: str, debug=False) -> dict:
    seq_dict = {}
    raw_ids = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        acc = acc_from_fasta_record(record)
        if acc:
            seq_dict[acc] = str(record.seq)
        raw_ids.append(record.id)

    if debug:
        print(f"[DEBUG] FASTA sequences read: {len(raw_ids)}")
        print("[DEBUG] FASTA sample record.id:", raw_ids[:5])
        print("[DEBUG] FASTA sample extracted acc keys:", list(seq_dict.keys())[:5])

    return seq_dict

def load_and_merge_interpro_data(protein_id, output_dir, debug=False):
    tsv_path = os.path.join(output_dir, f"{protein_id}_interpro.tsv")
    fasta_path = os.path.join(output_dir, f"{protein_id}_interpro.fasta")

    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"Missing TSV: {tsv_path}")
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Missing FASTA: {fasta_path}")

    df = read_interpro_tsv(tsv_path, debug=debug)
    seq_dict = read_interpro_fasta(fasta_path, debug=debug)

    # Map sequences using normalized accession
    df["Sequence"] = df["Accession_norm"].map(seq_dict)

    if debug:
        missing_seq = df["Sequence"].isna().sum()
        print(f"[DEBUG] Rows missing Sequence after mapping: {missing_seq} / {len(df)}")
        if missing_seq:
            print("[DEBUG] Example missing accessions (raw):",
                  df.loc[df["Sequence"].isna(), "Accession_raw"].head(10).tolist())
            print("[DEBUG] Example missing accessions (norm):",
                  df.loc[df["Sequence"].isna(), "Accession_norm"].head(10).tolist())

    # Drop unusable
    df = df.dropna(subset=["Sequence", "Tax ID"]).copy()

    # Fill Length if missing
    df.loc[df["Length"].isna(), "Length"] = df.loc[df["Length"].isna(), "Sequence"].str.len().astype("Int64")

    if debug:
        print(f"[DEBUG] Rows after dropna(Sequence, Tax ID): {len(df)}")

    return df

def get_species_map(taxids, debug=False):
    """
    Build a dict: QueryTaxID -> species_taxid using ete3.
    Robust to inputs like 9606.0 etc by numeric coercion before int.
    """
    taxids = pd.Series(list(taxids))
    taxids = pd.to_numeric(taxids, errors="coerce").dropna().astype(int).unique().tolist()

    species_map = {}
    for t in taxids:
        try:
            lineage = ncbi.get_lineage(int(t))
            ranks = ncbi.get_rank(lineage)
            species_tax = None
            for node in lineage:
                if ranks.get(node) == "species":
                    species_tax = node
                    break
            species_map[int(t)] = species_tax
        except Exception:
            species_map[int(t)] = None

    if debug:
        ok = sum(v is not None for v in species_map.values())
        print(f"[DEBUG] Species mapping success: {ok} / {len(species_map)}")

    return species_map

def add_species_taxid(df, debug=False):
    df = df.copy()
    smap = get_species_map(df["Tax ID"].dropna().astype(int).unique().tolist(), debug=debug)
    df["species_taxid"] = df["Tax ID"].astype(int).map(smap).astype("Int64")

    if debug:
        print("[DEBUG] species_taxid non-null:", int(df["species_taxid"].notna().sum()))
        print("[DEBUG] species_taxid sample:", df["species_taxid"].dropna().head(10).tolist())

    return df

def pick_one_rep_per_species(df):
    """Pick one representative per species to avoid many-to-many merge (choose longest)."""
    df = df.dropna(subset=["species_taxid", "Sequence"]).copy()
    df = df.sort_values(["species_taxid", "Length"], ascending=[True, False])
    return df.drop_duplicates(subset=["species_taxid"], keep="first")

def run_cdhit(in_fasta, out_fasta, c=0.90, n=5):
    subprocess.run(
        ["cd-hit", "-i", in_fasta, "-o", out_fasta, "-c", str(c), "-n", str(n), "-M", "0", "-d", "0","-T", "0"],
        check=True
    )

def create_combined_fasta(protein1_id, protein1_seq, protein2_id, protein2_seq, output_dir,
                          spacer_len=40, cdhit_c=0.90, cdhit_n=5, debug=False):

    os.makedirs(output_dir, exist_ok=True)

    prot1 = load_and_merge_interpro_data(protein1_id, output_dir, debug=debug)
    prot2 = load_and_merge_interpro_data(protein2_id, output_dir, debug=debug)

    # If either side has 0 rows, stop with guidance
    if prot1.empty or prot2.empty:
        print(f"prot1 hits: {len(prot1)}")
        print(f"prot2 hits: {len(prot2)}")
        print("No hits loaded. This is almost always TSV<->FASTA accession mismatch or empty files.")
        print("Run with --debug and paste the DEBUG samples (TSV accessions + FASTA ids).")
        return

    prot1 = add_species_taxid(prot1, debug=debug)
    prot2 = add_species_taxid(prot2, debug=debug)

    print(f"prot1 hits: {len(prot1)} | species-mapped: {int(prot1['species_taxid'].notna().sum())}")
    print(f"prot2 hits: {len(prot2)} | species-mapped: {int(prot2['species_taxid'].notna().sum())}")

    common_species = (
        set(prot1["species_taxid"].dropna().astype(int)) &
        set(prot2["species_taxid"].dropna().astype(int))
    )
    print(f"common species: {len(common_species)}")

    if not common_species:
        print("No common species found. Exiting without running CD-HIT.")
        return

    prot1_rep = pick_one_rep_per_species(prot1[prot1["species_taxid"].isin(common_species)])
    prot2_rep = pick_one_rep_per_species(prot2[prot2["species_taxid"].isin(common_species)])

    merged_df = pd.merge(
        prot1_rep[["species_taxid", "Sequence"]].rename(columns={"Sequence": "Sequence_1"}),
        prot2_rep[["species_taxid", "Sequence"]].rename(columns={"Sequence": "Sequence_2"}),
        on="species_taxid",
        how="inner"
    )

    print(f"paired rows before CD-HIT: {len(merged_df)}")
    if merged_df.empty:
        print("No pairs after merge. (Unexpected if common_species > 0).")
        return

    spacer = "X" * int(spacer_len)

    concat_records = []
    for _, row in merged_df.iterrows():
        stax = str(int(row["species_taxid"]))
        concat_seq = str(row["Sequence_1"]) + spacer + str(row["Sequence_2"])
        concat_records.append(SeqRecord(Seq(concat_seq), id=stax, description=""))

    concat_fa = os.path.join(output_dir, "concatenated_for_cdhit.fasta")
    cdhit_out = os.path.join(output_dir, "concatenated_cdhit_out.fasta")
    SeqIO.write(concat_records, concat_fa, "fasta")

    print("Running CD-HIT on concatenated sequences...")
    run_cdhit(concat_fa, cdhit_out, c=cdhit_c, n=cdhit_n)

    surviving = {int(r.id) for r in SeqIO.parse(cdhit_out, "fasta")}
    final_pairs = merged_df[merged_df["species_taxid"].astype(int).isin(surviving)]
    print(f"pairs after CD-HIT: {len(final_pairs)}")

    prot1_out_records = [SeqRecord(Seq(protein1_seq), id=protein1_id, description="query")] + [
        SeqRecord(Seq(row["Sequence_1"]), id=str(int(row["species_taxid"])), description="")
        for _, row in final_pairs.iterrows()
    ]
    prot2_out_records = [SeqRecord(Seq(protein2_seq), id=protein2_id, description="query")] + [
        SeqRecord(Seq(row["Sequence_2"]), id=str(int(row["species_taxid"])), description="")
        for _, row in final_pairs.iterrows()
    ]

    out1 = os.path.join(output_dir, f"{protein1_id}_paired_cdhit.fasta")
    out2 = os.path.join(output_dir, f"{protein2_id}_paired_cdhit.fasta")
    SeqIO.write(prot1_out_records, out1, "fasta")
    SeqIO.write(prot2_out_records, out2, "fasta")

    print(f"Wrote: {out1}")
    print(f"Wrote: {out2}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein1_id", required=True)
    ap.add_argument("--protein1_seq", required=True)
    ap.add_argument("--protein2_id", required=True)
    ap.add_argument("--protein2_seq", required=True)
    ap.add_argument("--output_dir", required=True)

    ap.add_argument("--spacer_len", type=int, default=40)
    ap.add_argument("--cdhit_c", type=float, default=0.90)
    ap.add_argument("--cdhit_n", type=int, default=5)

    ap.add_argument("--debug", action="store_true",
                    help="Print TSV/FASTA samples and mapping diagnostics.")
    ap.add_argument("--update_taxonomy", action="store_true",
                    help="Update ete3 taxonomy DB (downloads data).")

    args = ap.parse_args()

    if args.update_taxonomy:
        print("Updating ete3 taxonomy database...")
        ncbi.update_taxonomy_database()

    create_combined_fasta(
        args.protein1_id, args.protein1_seq,
        args.protein2_id, args.protein2_seq,
        args.output_dir,
        spacer_len=args.spacer_len,
        cdhit_c=args.cdhit_c,
        cdhit_n=args.cdhit_n,
        debug=args.debug
    )

if __name__ == "__main__":
    main()