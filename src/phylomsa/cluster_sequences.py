import argparse
import subprocess
import logging
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)

def _check_dependencies(tools: list[str]) -> None:
    """
    Verifies that required binaries exist in the system PATH.
    """
    for tool in tools:
        if not shutil.which(tool):
            raise FileNotFoundError(f"Dependency '{tool}' not found in PATH.")

def _prepend_sequence_to_msa(
        records: SeqRecord,
        sequence: str,
        seq_id="query_sequence"
    ) -> SeqRecord:
    new_rec = SeqRecord(Seq(sequence), id=seq_id, description="")
    records.insert(0, new_rec)
    return records

def _check_first_sequence_in_msa(
        msa_path: Union[str, Path],
        sequence: str
    ) -> None:
    """
    Checks if the first sequence in the MSA file matches the given sequence.
    """
    with open(msa_path, 'r') as msa_file:
        records = list(SeqIO.parse(msa_file, "fasta"))
        if not records:
            raise ValueError(f"No sequences found in MSA file: {msa_path}")
        first_seq = str(records[0].seq).replace('-', '')
        if first_seq != sequence:
            new_records = _prepend_sequence_to_msa(records, sequence)
            SeqIO.write(new_records, msa_path, "fasta")

def _convert_msa_format(
        input_file: Union[str, Path],
        output_format: str,
        output_file: Union[str, Path, None] = None
    ) -> None:
    """
    Converts MSA files to the required format using an external tool.
    By default, it will save the converted file in place.
    """
    _check_dependencies(['reformat.pl'])

    input_path = Path(input_file)
    input_format = input_path.suffix.lstrip('.')
    if output_file is None:
        output_path = input_path
    else:
        output_path = Path(output_file)

    cmd = [
        'reformat.pl',
        str(input_path),
        str(output_path)
    ]

    logger.debug(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def _run_multi_seq_alignment(
        input_file: Union[str, Path],
        output_file: Union[str, Path],
        **kwargs
    ) -> None:
    """
    Runs MAFFT to perform multiple sequence alignment.
    """
    #TODO: Add MUSCLE support later
    _check_dependencies(['mafft'])

    cmd_string = f"mafft {input_file} > {output_file}"
    subprocess.run(cmd_string, shell=True, check=True)

    logger.debug(f"Running command: {cmd_string}")

def run_alignment(
        query_sequence: str,
        input_file: Union[str, Path],
        output_format: str = 'a3m',
        mafft_params: Optional[dict] = {}
    ) -> None:
    """
    Main function to run the alignment conversion.
    """

    alignment_output_path = Path(input_file).with_suffix('.aligned.fasta')
    a3m_output_path = Path(input_file).with_suffix(f'.{output_format}')

    _check_first_sequence_in_msa(input_file, query_sequence)
    _run_multi_seq_alignment(input_file, alignment_output_path, **mafft_params)
    _convert_msa_format(alignment_output_path, output_format, a3m_output_path)