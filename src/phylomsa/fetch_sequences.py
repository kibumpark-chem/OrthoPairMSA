import os

import requests


def fetch_data_by_uniprot_id(uniprot_id, output_dir):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url, timeout=60)
    fetched_data = {
        "uniprot_id": uniprot_id,
        "status_code": response.status_code,
        "sequence": None,
        "interpro_ids": [],
        "eggnog_ids": [],
        "orthodb_ids": [],
    }

    if response.status_code == 200:
        data = response.json()
        sequence = data["sequence"]["value"]
        fetched_data["sequence"] = sequence
        os.makedirs(output_dir, exist_ok=True)
        with open(os.path.join(output_dir, f"{uniprot_id}.fasta"), "w") as f:
            f.write(f">{uniprot_id}\n{sequence}\n")
        print(f"Sequence for {uniprot_id} saved successfully.")

        interpro_ids = []
        eggnog_ids = []
        orthodb_ids = []
        for ref in data.get("uniProtKBCrossReferences", []):
            db = ref.get("database")
            rid = ref.get("id")
            if not rid:
                continue
            if db == "InterPro":
                interpro_ids.append(rid)
            elif db == "eggNOG":
                eggnog_ids.append(rid)
            elif db == "OrthoDB":
                orthodb_ids.append(rid)
        fetched_data["interpro_ids"] = interpro_ids
        fetched_data["eggnog_ids"] = eggnog_ids
        fetched_data["orthodb_ids"] = orthodb_ids
    else:
        print(f"Failed to fetch sequence for {uniprot_id}. Status code: {response.status_code}")

    return fetched_data
