
import requests

def fetch_entry_details(ipr_id):
    """Fetches type, name, and hierarchy info from the InterPro API."""
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/{ipr_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        metadata = data.get("metadata", {})
        return {
            "id": ipr_id,
            "name": metadata.get("name"),
            "type": metadata.get("type", "").lower(),
            "hierarchy": metadata.get("hierarchy")  # Dictionary containing parent/child mappings
        }
    return None

def is_parent_of(parent_id, child_id, hierarchy_map):
    """
    Recursively checks if parent_id is an ancestor of child_id 
    using the hierarchy structure returned by the InterPro API.
    """
    if not hierarchy_map:
        return False
    
    # If we are currently looking at the parent node, check if child is in its descendants
    if hierarchy_map.get("accession") == parent_id:
        return check_descendants_for_id(hierarchy_map, child_id)
        
    return False

def check_descendants_for_id(node, target_id):
    """Helper to traverse down the children arrays looking for a specific target ID."""
    children = node.get("children", [])
    for child in children:
        if child.get("accession") == target_id:
            return True
        if check_descendants_for_id(child, target_id):
            return True
    return False

def filter_redundant_interpro(interpro_ids):
    # Fetch metadata for all IDs
    print(f"Analyzing {len(interpro_ids)} InterPro records and pulling hierarchies...")
    entries = []
    for ipr_id in interpro_ids:
        details = fetch_entry_details(ipr_id)
        if details:
            entries.append(details)
            
    ## Separate into families and non-families (superfamilies, domains, repeats, sites)
    families = [e for e in entries if e["type"] == "family"]
    other_types = [e for e in entries if e["type"] != "family"]
    
    # Identify redundant parent families
    redundant_parents = set()
    
    for f1 in families:
        for f2 in families:
            if f1["id"] == f2["id"]:
                continue
            
            # Check if f1 is a parent of f2 using the hierarchy layout
            if is_parent_of(f1["id"], f2["id"], f2["hierarchy"]) or is_parent_of(f1["id"], f2["id"], f1["hierarchy"]):
                redundant_parents.add(f1["id"])

    # Filter out the parents to isolate the most specific child family
    filtered_families = [f for f in families if f["id"] not in redundant_parents]
    
    # Report specifically the surviving Family
    print("THE MOST SPECIFIC EVOLUTIONARY FAMILY:")
    if not filtered_families:
        print("--> No specific family assigned.")
    for f in filtered_families:
        print(f"--> {f['id']}: {f['name']}")        
    return filtered_families
