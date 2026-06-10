import pandas as pd
from Bio import SeqIO
from ete3 import NCBITaxa
from .sequences import SequenceRecord
from ..utils.taxonomy import add_taxonomic_info

class BioDataLoader:
    """Handles the parsing of FASTA and TSV files into SequenceRecords."""

    @classmethod
    def load_data(cls, fasta_path: str, tsv_path: str, fetch_lineage: bool = False) -> list[SequenceRecord]:
        records = []
        seq_df = pd.read_csv(tsv_path, sep='\t')
        seq_records = list(SeqIO.parse(fasta_path, 'fasta'))
        id_to_seq = cls.create_id_to_seq_dict(seq_records)

        if fetch_lineage:
            seq_df = add_taxonomic_info(seq_df)

        for _, row in seq_df.iterrows():
            accession = row['Accession']
            if accession in id_to_seq:
                record = SequenceRecord(
                    id=accession,
                    sequence=str(id_to_seq[accession]),
                    taxonomy=row.get('Taxonomy', None),
                    metadata=row.to_dict()
                )
                records.append(record)
        return records
    
    @staticmethod
    def update_taxonomy():
        from ete3 import NCBITaxa
        print("Updating NCBI Taxonomy database... this may take several minutes.")
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()
        print("Update complete.")

    @staticmethod
    def create_id_to_seq_dict(seq_records):
        return {rec.id.split('|')[0]: rec.seq for rec in seq_records}