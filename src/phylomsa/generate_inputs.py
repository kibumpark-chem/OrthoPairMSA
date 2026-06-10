import argparse
import json
import string
import pathlib
import yaml # Required for Boltz output
from enum import Enum
from typing import Union, List

_ALPHABET = string.ascii_uppercase + string.ascii_lowercase + string.digits

class Mode(Enum):
    DEFAULT = "default"
    CUSTOM_MSA = "custom_msa"
    PRECOMPUTED_MSA = "precomputed_msa"

def _get_chain_id(n: int) -> str:
    if n < len(_ALPHABET):
        return _ALPHABET[n]
    raise ValueError("Exceeded maximum number of chain IDs.")

def _check_required_keys(
        required_keys: List[Union[str, tuple]], provided_keys: List[str]
    ) -> None:
    missing_keys = []
    for key in required_keys:
        if isinstance(key, tuple):
            if not any(k in provided_keys for k in key):
                missing_keys.append(key)
        else:
            if key not in provided_keys:
                missing_keys.append(key)
    if missing_keys:
        raise ValueError(f"Missing required keys: {missing_keys}")

def _protein_config(
        idx: int,
        version: int = 1,
        **kwargs
    ) -> dict:
    config = {
        "id": kwargs.get("id", _get_chain_id(idx)), 
    }

    _check_required_keys(["sequence"], list(kwargs.keys()))
    config["sequence"] = kwargs["sequence"]

    if version == 2:
        unpaired = kwargs.get("unpairedMsaPath") or kwargs.get("unpairedMSA")
        if unpaired:
            config["unpairedMsaPath"] = str(unpaired)
            
        paired = kwargs.get("pairedMsaPath") or kwargs.get("pairedMsa")
        if paired:
            config["pairedMsaPath"] = str(paired)
            
        if "templates" in kwargs and kwargs["templates"] is not None:
            config["templatesPath"] = str(kwargs["templates"])

    return config

def _ligand_config(
        idx: int,
        **kwargs
    ) -> dict:
    config = {
        "id": kwargs.get("id", _get_chain_id(idx)),
    }
    if 'smiles' not in kwargs and 'ccdCodes' not in kwargs:
        raise ValueError("Ligand configuration requires 'smiles' or 'ccdCodes'.")

    if "ccdCodes" in kwargs:
        config["ccdCodes"] = kwargs["ccdCodes"]
    if "smiles" in kwargs:
        config["smiles"] = kwargs["smiles"]

    return config

# ---------------------------------------------------------
# TARGET 1: ALPHAFOLD 3
# ---------------------------------------------------------
def create_af3_config(
        output_path: Union[str, pathlib.Path],
        mode: Union[Mode, str],
        job_name: str,
        entities: List[dict],
        model_seeds: List[int] = None
    ) -> None:
    if model_seeds is None:
        model_seeds = [1]

    config = {
        "name": job_name,
        "modelSeeds": model_seeds,
    }

    if isinstance(mode, str):
        mode = Mode(mode)

    if mode == Mode.DEFAULT:
        version = 1
    elif mode in (Mode.CUSTOM_MSA, Mode.PRECOMPUTED_MSA):
        version = 2

    sequences = []
    for idx, entity in enumerate(entities):
        if entity['type'] == 'protein':
            safe_params = entity.get('params', {}).copy()
            safe_params.pop("templates", None) 
            sequence_config = _protein_config(idx, version, **safe_params)
        elif entity['type'] == 'ligand':
            sequence_config = _ligand_config(idx, **entity.get('params', {}))
        else:
            raise NotImplementedError(f"Handling for {entity['type']} not implemented.")
        
        sequences.append({entity['type']: sequence_config})

    config["sequences"] = sequences
    config["dialect"] = "alphafold3"
    config["version"] = version

    with open(output_path, 'w') as config_file:
        json.dump(config, config_file, indent=4)

# ---------------------------------------------------------
# TARGET 2: PROTENIX
# ---------------------------------------------------------
def create_protenix_config(
        output_path: Union[str, pathlib.Path],
        job_name: str,
        entities: List[dict],
        model_seeds: List[int] = None
    ) -> None:
    if model_seeds is None:
        model_seeds = [1]

    config = {
        "name": job_name,
        "modelSeeds": model_seeds,
    }

    sequences = []
    for idx, entity in enumerate(entities):
        if entity['type'] == 'protein':
            sequence_config = _protein_config(idx, version=2, **entity.get('params', {}))
        elif entity['type'] == 'ligand':
            sequence_config = _ligand_config(idx, **entity.get('params', {}))
        else:
            raise NotImplementedError(f"Handling for {entity['type']} not implemented.")
        
        sequences.append({entity['type']: sequence_config})

    config["sequences"] = sequences
    final_output = [config] 

    with open(output_path, 'w') as config_file:
        json.dump(final_output, config_file, indent=4)

# ---------------------------------------------------------
# TARGET 3: BOLTZ
# ---------------------------------------------------------
def create_boltz_config(
        output_path: Union[str, pathlib.Path],
        job_name: str, # Not strictly used in YAML, but kept for signature parity
        entities: List[dict],
        model_seeds: List[int] = None
    ) -> None:
    """
    Creates a configuration file for Boltz-1/Boltz-2 (YAML format).
    """
    config = {
        "version": 1,
        "sequences": []
    }

    for idx, entity in enumerate(entities):
        if entity['type'] == 'protein':
            params = entity.get('params', {})
            seq_config = {
                "id": params.get("id", _get_chain_id(idx)),
                "sequence": params["sequence"]
            }
            
            # Map Unpaired MSA to Boltz's simple 'msa' key
            if "unpairedMsaPath" in params:
                seq_config["msa"] = str(params["unpairedMsaPath"])
                
            # Warn about Paired MSA
            if "pairedMsaPath" in params:
                print(f"Boltz Note: Chain {seq_config['id']} provided a paired MSA path. "
                      "Boltz strictly requires CSV files for custom paired alignments. "
                      "This path was excluded from the YAML.")

            config["sequences"].append({"protein": seq_config})
            
        elif entity['type'] == 'ligand':
            params = entity.get('params', {})
            lig_config = {
                "id": params.get("id", _get_chain_id(idx)),
            }
            
            # Convert list of CCD codes back to a single string for Boltz
            if "ccdCodes" in params and params["ccdCodes"]:
                lig_config["ccd"] = params["ccdCodes"][0]
            if "smiles" in params:
                lig_config["smiles"] = params["smiles"]
                
            config["sequences"].append({"ligand": lig_config})
            
        else:
            raise NotImplementedError(f"Handling for {entity['type']} not implemented.")

    with open(output_path, 'w') as config_file:
        yaml.dump(config, config_file, sort_keys=False, default_flow_style=False)

# ---------------------------------------------------------
# CLI RUNNER
# ---------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate JSON/YAML input files dynamically for AlphaFold 3, Protenix, and Boltz.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("--job-name", required=True, help="Name of the job.")
    parser.add_argument("--seed", type=int, default=1, help="Model seed (default: 1).")
    parser.add_argument("--target", choices=["af3", "protenix", "boltz", "all"], default="all", 
                        help="Which system to generate the JSON for (default: all).")
    
    parser.add_argument(
        "--chain", 
        action="append", 
        nargs="+", 
        required=True,
        help=(
            "Define a protein chain. Pass this argument multiple times for multimers.\n"
            "Format: ID SEQUENCE [UNPAIRED_MSA] [PAIRED_MSA] [TEMPLATE_PATH]\n"
            "Only ID and SEQUENCE are required. Use 'none' to skip a middle argument."
        )
    )

    parser.add_argument(
        "--ligand", 
        action="append", 
        nargs="+", 
        help=(
            "Define a small molecule/ligand.\n"
            "Format: ID TYPE VALUE\n"
            "TYPE must be 'ccd' or 'smiles'. Example: Z ccd ATP or Y smiles CCO"
        )
    )

    args = parser.parse_args()
    entities = []

    # 1. Parse Proteins
    if args.chain:
        for chain_data in args.chain:
            if len(chain_data) < 2:
                parser.error(f"Chain definition missing arguments. Got {len(chain_data)}, expected at least 2: ID SEQUENCE")
            
            params = {
                "id": chain_data[0],
                "sequence": chain_data[1]
            }
            
            if len(chain_data) >= 3 and chain_data[2].lower() != 'none':
                params["unpairedMsaPath"] = chain_data[2]
            if len(chain_data) >= 4 and chain_data[3].lower() != 'none':
                params["pairedMsaPath"] = chain_data[3]
            if len(chain_data) >= 5 and chain_data[4].lower() != 'none':
                params["templates"] = chain_data[4]
                
            entities.append({"type": "protein", "params": params})

    # 2. Parse Ligands
    if args.ligand:
        for lig_data in args.ligand:
            if len(lig_data) != 3:
                parser.error(f"Ligand definition missing arguments. Got {len(lig_data)}, expected exactly 3: ID TYPE VALUE")
            
            lig_id, lig_type, lig_value = lig_data
            params = {"id": lig_id}
            
            if lig_type.lower() == "ccd":
                params["ccdCodes"] = [lig_value] 
            elif lig_type.lower() == "smiles":
                params["smiles"] = lig_value
            else:
                parser.error(f"Invalid ligand type '{lig_type}'. Must be 'ccd' or 'smiles'.")
                
            entities.append({"type": "ligand", "params": params})

    # 3. Execute generators
    if args.target in ["af3", "all"]:
        has_msa = any("unpairedMsaPath" in e["params"] or "pairedMsaPath" in e["params"] for e in entities)
        af3_mode = Mode.CUSTOM_MSA if has_msa else Mode.DEFAULT
        
        out_name = f"{args.job_name}_af3.json"
        create_af3_config(
            output_path=out_name,
            mode=af3_mode,
            job_name=args.job_name,
            entities=entities,
            model_seeds=[args.seed]
        )
        print(f"Generated: {out_name}")

    if args.target in ["protenix", "all"]:
        out_name = f"{args.job_name}_protenix.json"
        create_protenix_config(
            output_path=out_name,
            job_name=args.job_name,
            entities=entities,
            model_seeds=[args.seed]
        )
        print(f"Generated: {out_name}")

    if args.target in ["boltz", "all"]:
        out_name = f"{args.job_name}_boltz.yaml"
        create_boltz_config(
            output_path=out_name,
            job_name=args.job_name,
            entities=entities,
            model_seeds=[args.seed]
        )
        print(f"Generated: {out_name}")

if __name__ == "__main__":
    main()