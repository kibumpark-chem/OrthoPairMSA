import requests
import json
import os

def fetch_data_by_uniprot_id(uniprot_id, output_dir):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.json"
    response = requests.get(url)
    fetched_data = {
        "uniprot_id": uniprot_id,
        "status_code": response.status_code,
        "sequence": None,
        "interpro_ids": []
    }
    
    if response.status_code == 200:
        data = response.json()
        sequence = data['sequence']['value']
        fetched_data["sequence"] = sequence
        with open(os.path.join(output_dir, f"{uniprot_id}.fasta"), 'w') as f:
            f.write(f">{uniprot_id}\n{sequence}\n")
        print(f"Sequence for {uniprot_id} saved successfully.")

        interpro_ids = []
        for ref in response.json().get("uniProtKBCrossReferences", []):
            if ref["database"] == "InterPro":
                interpro_ids.append(ref["id"])
        fetched_data["interpro_ids"] = interpro_ids
    else:
        print(f"Failed to fetch sequence for {uniprot_id}. Status code: {response.status_code}")

    return fetched_data