import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree,NCBITaxa

from dataclasses import dataclass, field
from typing import Dict, Optional

#TODO: Filter by domain, kingdom, etc.
#TODO: Overlap between two protein families (list of records)
#TODO: One representative sequence per species

@dataclass(frozen=True)
class SequenceRecord:
    """
    Represents a biological sequence combined with metadata from a TSV.
    """
    id: str
    sequence: str
    taxonomy: Optional[str] = None
    metadata: Dict = field(default_factory=dict)

    def __len__(self) -> int:
        """Returns the length of the sequence."""
        return len(self.sequence)

    def add_taxonomy(self, taxonomy: str):
        """Adds a taxonomy string to the record."""
        object.setattr(self, 'taxonomy', taxonomy)

    def to_fasta(self) -> str:
        """Returns the record in standard FASTA format."""
        return f">{self.id} {self.taxonomy or ''}\n{self.sequence}"


    
