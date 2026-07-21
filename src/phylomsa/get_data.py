import argparse
import os
import subprocess
import sys

from fetch_ortholog_group import fetch_ortholog_group
from fetch_protein_family import filter_redundant_interpro
from fetch_sequences import fetch_data_by_uniprot_id

VALID_SOURCES = ("eggnog", "orthodb", "interpro")


def run_interpro(uniprot_id: str, output_dir: str, interpro_ids: list) -> str:
    family_id = filter_redundant_interpro(interpro_ids)
    if len(family_id) == 0:
        raise ValueError("No valid InterPro family IDs found after filtering.")
    elif len(family_id) == 1:
        print(f"One valid InterPro family ID found: {family_id[0]['id']}")
    else:
        print(f"Multiple valid InterPro family IDs found: {', '.join([id['id'] for id in family_id])}")
        print("Assuming the first id to be the representative family.")
    family_id = family_id[0]["id"]

    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(current_dir))
    fasta_out_path = os.path.join(output_dir, f"{uniprot_id}_interpro.fasta")
    tsv_out_path = os.path.join(output_dir, f"{uniprot_id}_interpro.tsv")
    fasta_script_path = os.path.join(project_root, "scripts", "interpro_api", "download_fasta_via_api.py")
    tsv_script_path = os.path.join(project_root, "scripts", "interpro_api", "download_tsv_via_api.py")

    with open(fasta_out_path, "w") as fasta_file:
        subprocess.run(
            [sys.executable, fasta_script_path, family_id],
            stdout=fasta_file,
            check=True,
        )
    with open(tsv_out_path, "w") as tsv_file:
        subprocess.run(
            [sys.executable, tsv_script_path, family_id],
            stdout=tsv_file,
            check=True,
        )
    return family_id


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fetch query sequence and ortholog-group members for a UniProt ID."
    )
    parser.add_argument("--uniprot_id", required=True, help="UniProt ID of the protein to fetch.")
    parser.add_argument("--output_dir", required=True, help="Directory to save outputs.")
    parser.add_argument(
        "--source",
        default="eggnog",
        choices=VALID_SOURCES,
        help="Ortholog source: eggnog (default), orthodb, or interpro.",
    )
    parser.add_argument(
        "--group_id",
        default=None,
        help="Optional explicit eggNOG/OrthoDB/InterPro group id (skips auto-resolve).",
    )
    args = parser.parse_args()
    source = args.source.lower()

    os.makedirs(args.output_dir, exist_ok=True)

    fetched_data = fetch_data_by_uniprot_id(args.uniprot_id, args.output_dir)
    if not fetched_data["sequence"]:
        print(f"Failed to fetch sequence for {args.uniprot_id}.")
        return 1

    print(f"Successfully fetched sequence for {args.uniprot_id}.")
    if fetched_data.get("interpro_ids"):
        print(f"Associated InterPro IDs: {', '.join(fetched_data['interpro_ids'])}")

    fasta_out_path = os.path.join(args.output_dir, f"{args.uniprot_id}_interpro.fasta")
    tsv_out_path = os.path.join(args.output_dir, f"{args.uniprot_id}_interpro.tsv")

    if source == "interpro":
        if args.group_id:
            # Direct InterPro download bypassing family filtering
            current_dir = os.path.dirname(os.path.abspath(__file__))
            project_root = os.path.dirname(os.path.dirname(current_dir))
            fasta_script_path = os.path.join(
                project_root, "scripts", "interpro_api", "download_fasta_via_api.py"
            )
            tsv_script_path = os.path.join(
                project_root, "scripts", "interpro_api", "download_tsv_via_api.py"
            )
            print(f"Using provided InterPro ID: {args.group_id}")
            with open(fasta_out_path, "w") as fasta_file:
                subprocess.run(
                    [sys.executable, fasta_script_path, args.group_id],
                    stdout=fasta_file,
                    check=True,
                )
            with open(tsv_out_path, "w") as tsv_file:
                subprocess.run(
                    [sys.executable, tsv_script_path, args.group_id],
                    stdout=tsv_file,
                    check=True,
                )
            family_id = args.group_id
        else:
            family_id = run_interpro(args.uniprot_id, args.output_dir, fetched_data["interpro_ids"])
        print(f"InterPro download complete for family {family_id}.")
    else:
        group_id, n_fasta, n_tsv = fetch_ortholog_group(
            uniprot_id=args.uniprot_id,
            source=source,
            fasta_out_path=fasta_out_path,
            tsv_out_path=tsv_out_path,
            group_id=args.group_id,
        )
        print(
            f"{source} download complete for group {group_id} "
            f"(fasta={n_fasta}, tsv={n_tsv})."
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
