import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree,NCBITaxa

#TODO: Add sequence clustering and sequence filtering based on clustering results
#TODO: May need to account for more general fasta headers

_ncbi_instance = None

def _get_ncbi():
    global _ncbi_instance
    if _ncbi_instance is None:
        # This only runs the VERY FIRST time get_ncbi() is called
        _ncbi_instance = NCBITaxa()
    return _ncbi_instance

def _get_ranked_lineage(taxid_list, desired_ranks=None):
    ncbi = _get_ncbi()
    unique_ids = list(set(taxid_list))

    all_lineage_data = []
    for taxid in unique_ids:
        try:
            lineage = ncbi.get_lineage(taxid)
            ranks = ncbi.get_rank(lineage)
            names = ncbi.get_taxid_translator(lineage)
            lineage_map = {
                ranks[t]: names[t] 
                for t in lineage 
                if ranks[t] != 'no rank' # Exclude nodes with no official rank
            }
            lineage_map['Tax ID'] = taxid
            all_lineage_data.append(lineage_map)
        except Exception as e:
            all_lineage_data.append({'Tax ID': taxid})
        
    df = pd.DataFrame(all_lineage_data)

    if desired_ranks:
        for rank in desired_ranks:
            if rank not in df.columns:
                df[rank] = None
        df = df[['Tax ID'] + desired_ranks]

    return df

def add_taxonomic_info(data_df, taxid_col='Tax ID'):
    unique_lineage = _get_ranked_lineage(
        data_df[taxid_col].dropna().unique(), 
        desired_ranks=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain']
    )
    return data_df.merge(unique_lineage, on=taxid_col, how='left')

def create_id_to_seq_dict(seq_records):
    # Assuming each record has a unique identifier in the format 'id|other_info'
    id_map = {rec.id.split('|')[0]: str(rec.seq) for rec in seq_records}
    return id_map

def write_sequences_to_fasta(
        seq_df,
        id_map,
        output_file
    ):
    records_list = []
    for i, row in seq_df.iterrows():
        tax_id = row['Tax ID']
        accession = row['Accession']
        domain = row['domain']
        species = row['species']

        records_list.append(SeqRecord(Seq(id_map[accession]), id=f"{tax_id}", description=f"{species}"))
    new_fasta = SeqIO.write(records_list, str(output_file), 'fasta')
    
    return new_fasta