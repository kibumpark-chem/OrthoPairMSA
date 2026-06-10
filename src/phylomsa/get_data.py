import argparse
from fetch_sequences import fetch_data_by_uniprot_id
from fetch_protein_family import filter_redundant_interpro
import subprocess
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch sequences and InterPro family data for a given UniProt ID.")
    parser.add_argument("--uniprot_id", required=True, help="UniProt ID of the protein to fetch.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the fetched sequence and metadata.")
    
    args = parser.parse_args()
    
    # Fetch sequence and InterPro IDs using the provided UniProt ID
    fetched_data = fetch_data_by_uniprot_id(args.uniprot_id, args.output_dir)
    
    if fetched_data["sequence"]:
        print(f"Successfully fetched sequence for {args.uniprot_id}.")
        print(f"Associated InterPro IDs: {', '.join(fetched_data['interpro_ids'])}")
    else:
        print(f"Failed to fetch sequence for {args.uniprot_id}.")

    # Fetch interpro family data using the InterPro IDs
    family_id = filter_redundant_interpro(fetched_data["interpro_ids"])
    if len(family_id) == 0:
        raise ValueError("No valid InterPro family IDs found after filtering.")
    elif len(family_id) == 1:
        print(f"One valid InterPro family ID found: {family_id[0]['id']}")
    else:
        print(f"Multiple valid InterPro family IDs found: {', '.join([id['id'] for id in family_id])}")
        print("Assuming the first id to be the representative family.")
    family_id = family_id[0]['id']

    # Run InterPro API script
    # TODO: Need to clean up the hardcoded paths after testing
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(current_dir))
    fasta_out_path = os.path.join(args.output_dir, f"{args.uniprot_id}_interpro.fasta")
    tsv_out_path = os.path.join(args.output_dir, f"{args.uniprot_id}_interpro.tsv")
    fasta_script_path = os.path.join(project_root, "scripts", "interpro_api", "download_fasta_via_api.py")
    with open(fasta_out_path, "w") as fasta_file:
        subprocess.run(
            ["python", fasta_script_path, family_id], 
            stdout=fasta_file, 
            check=True
        )
    tsv_script_path = os.path.join(project_root, "scripts", "interpro_api", "download_tsv_via_api.py")
    with open(tsv_out_path, "w") as tsv_file:
        subprocess.run(
            ["python", tsv_script_path, family_id], 
            stdout=tsv_file, 
            check=True
        )
