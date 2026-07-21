"""Fetch ortholog groups via UniProt cross-references (eggNOG / OrthoDB).

Writes InterPro-compatible outputs expected by create_filtered_fasta.py:
  {uniprot_id}_interpro.fasta
  {uniprot_id}_interpro.tsv
"""

from __future__ import annotations

import time
from typing import List, Optional, Tuple

import requests

UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"

SOURCE_DB = {
    "eggnog": "eggNOG",
    "orthodb": "OrthoDB",
}


def get_ortholog_group_ids(uniprot_id: str, source: str) -> List[str]:
    """Return ortholog-group IDs for *source* from a UniProt entry."""
    source = source.lower()
    if source not in SOURCE_DB:
        raise ValueError(f"Unsupported source: {source}. Choose from {list(SOURCE_DB)}")

    db_name = SOURCE_DB[source]
    url = UNIPROT_ENTRY.format(accession=uniprot_id)
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    data = response.json()

    ids = []
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == db_name:
            oid = ref.get("id")
            if oid and oid not in ids:
                ids.append(oid)
    return ids


def pick_group_id(group_ids: List[str], source: str) -> str:
    """Pick a representative OG id when multiple are present."""
    if not group_ids:
        raise ValueError(f"No {source} cross-references found on UniProt entry.")
    if len(group_ids) == 1:
        print(f"One {source} group ID found: {group_ids[0]}")
        return group_ids[0]

    print(f"Multiple {source} group IDs found: {', '.join(group_ids)}")
    # Prefer KOG / COG / ENOG-style ids for eggNOG; otherwise first id.
    if source == "eggnog":
        preferred = [g for g in group_ids if g.startswith(("KOG", "COG", "ENOG", "arCOG"))]
        if preferred:
            print(f"Using preferred eggNOG id: {preferred[0]}")
            return preferred[0]
    print(f"Assuming the first id to be the representative group: {group_ids[0]}")
    return group_ids[0]


def _uniprot_query(source: str, group_id: str) -> str:
    # UniProt REST: xref:eggnog-KOG0469 / xref:orthodb-364892at2759
    db = SOURCE_DB[source].lower()
    return f"xref:{db}-{group_id}"


def stream_fasta(source: str, group_id: str, out_path: str) -> int:
    """Download UniProt FASTA for all members of an OG. Returns sequence count."""
    query = _uniprot_query(source, group_id)
    params = {"query": query, "format": "fasta"}
    print(f"Streaming FASTA from UniProt for {source}={group_id} ...")
    with requests.get(UNIPROT_STREAM, params=params, stream=True, timeout=600) as r:
        r.raise_for_status()
        n = 0
        with open(out_path, "w") as fh:
            for line in r.iter_lines(decode_unicode=True):
                if line is None:
                    continue
                if line.startswith(">"):
                    n += 1
                fh.write(line + "\n")
    print(f"Wrote {n} sequences to {out_path}")
    return n


def stream_tsv(source: str, group_id: str, out_path: str) -> int:
    """Download metadata TSV in InterPro-compatible 8-column layout."""
    query = _uniprot_query(source, group_id)
    fields = "accession,id,protein_name,organism_id,organism_name,length"
    params = {"query": query, "format": "tsv", "fields": fields}
    print(f"Streaming TSV from UniProt for {source}={group_id} ...")

    n = 0
    with requests.get(UNIPROT_STREAM, params=params, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(out_path, "w") as fh:
            header_skipped = False
            for line in r.iter_lines(decode_unicode=True):
                if line is None:
                    continue
                if not header_skipped:
                    header_skipped = True
                    continue  # skip UniProt header
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                accession, _entry_name, protein_name, tax_id, sci_name, length = parts[:6]
                # Match InterPro download_tsv columns:
                # Accession, Database, Name, Tax ID, Scientific Name, Length,
                # Entry Accession, Locations
                fh.write(
                    "\t".join(
                        [
                            accession,
                            SOURCE_DB[source],
                            protein_name,
                            tax_id,
                            sci_name,
                            length,
                            group_id,
                            "",
                        ]
                    )
                    + "\n"
                )
                n += 1
    print(f"Wrote {n} TSV rows to {out_path}")
    return n


def fetch_ortholog_group(
    uniprot_id: str,
    source: str,
    fasta_out_path: str,
    tsv_out_path: str,
    group_id: Optional[str] = None,
) -> Tuple[str, int, int]:
    """
    Resolve OG id (unless provided), download FASTA + TSV via UniProt.

    Returns (group_id, n_fasta, n_tsv).
    """
    source = source.lower()
    if group_id is None:
        group_ids = get_ortholog_group_ids(uniprot_id, source)
        group_id = pick_group_id(group_ids, source)
    else:
        print(f"Using provided {source} group ID: {group_id}")

    n_fasta = stream_fasta(source, group_id, fasta_out_path)
    # Be polite between large UniProt streams
    time.sleep(1)
    n_tsv = stream_tsv(source, group_id, tsv_out_path)

    if n_fasta == 0 or n_tsv == 0:
        raise ValueError(
            f"No UniProt members found for {source} group {group_id} "
            f"(fasta={n_fasta}, tsv={n_tsv})."
        )
    return group_id, n_fasta, n_tsv
