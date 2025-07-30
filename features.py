import numpy as np

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from enum import Enum, auto
from typing import Optional, List, Tuple, Union


class Feature(Enum):
    SIRNA_TOKEN = auto()
    MRNA_5_FLANK_TOKEN = auto()
    MRNA_3_FLANK_TOKEN = auto()


class mRNA:
    """Represents an mRNA sequence and methods for featurization."""
    def __init__(self, seq_record: SeqRecord) -> None:
        """
        Initializes the mRNA object from a Biopython SeqRecord.

        Args:
            seq_record (SeqRecord): The Biopython record for the mRNA.
        """
        self.seq = seq_record.seq.back_transcribe().upper()
        self.id = seq_record.id
        self.max_mrna_len = 9756

    def reduce_mrna(self, sirna: 'siRNA', extension: int = 500, rnai_size: Optional[int] = None, pad: bool = False) -> None:
        """
        Reduces the mRNA to a window around the siRNA target site.

        Args:
            sirna (siRNA): The siRNA defining the target region.
            extension (int): The number of nucleotides to include on either side.
            rnai_size (Optional[int]): The length of the siRNA match region.
            pad (bool): If True, pads with 'N' to maintain a constant length.
        """
        if rnai_size is None:
            rnai_match = sirna.seq.reverse_complement().back_transcribe().upper()
        else:
            rnai_match = (
                sirna.seq[:rnai_size].reverse_complement().back_transcribe().upper()
            )
        self.max_mrna_size = 2 * extension + len(rnai_match)
        rna_index = self.seq.find(rnai_match)
        mrna_start = rna_index - extension
        mrna_end = rna_index + len(rnai_match) + extension
        sequence_length = len(self.seq)
        self.seq = self.seq[max(mrna_start, 0) : min(mrna_end, sequence_length)]
        if pad:
            for _ in range(mrna_start, 0):
                self.seq = Seq("N") + self.seq
            for _ in range(sequence_length, mrna_end):
                self.seq += Seq("N")

    def nucleotide_tokenization(self, sirna: 'siRNA', split: bool) -> Union[List[int], Tuple[List[int], List[int]]]:
        """
        Converts the mRNA sequence flanking an siRNA into integer tokens.

        Args:
            sirna (siRNA): The reference siRNA to define the flanking regions.
            split (bool): If True, returns the 5' and 3' flanks as a tuple.
                          If False, returns the combined flanking sequence.

        Returns:
            Union[List[int], Tuple[List[int], List[int]]]: A list of integer tokens,
            or a tuple containing two such lists for the 5' and 3' flanks.
        """
        nt_conversion = {"N": 0, "R": 0, "Y": 0, "G": 1, "C": 2, "A": 3, "U": 4, "T": 4}
        self.reduce_mrna(sirna, extension=19, rnai_size=19, pad=True)
        if split:
            return [nt_conversion[nt] for nt in self.seq[:19]], [
                nt_conversion[nt] for nt in self.seq[-19:]
            ]
        return [nt_conversion[nt] for nt in self.seq]


class siRNA:
    """Represents an siRNA sequence and methods for featurization."""
    def __init__(self, seq_record_antisense: SeqRecord) -> None:
        """
        Initializes the siRNA object from an antisense strand SeqRecord.

        Args:
            seq_record_antisense (SeqRecord): The Biopython record for the
                siRNA antisense strand.
        """
        assert (
            len(seq_record_antisense.seq) >= 19
        ), "siRNA is shorter than 19nt which is the minimum"
        self.max_sirna_len = 21
        self.seq = seq_record_antisense.seq.transcribe().upper()[: self.max_sirna_len]
        self.id = seq_record_antisense.id
        self.thermo_row = []

    def nucleotide_tokenization(self):
        """
        Converts the siRNA sequence into a fixed-length list of integer tokens.

        The sequence is padded with 0s to a maximum length of 21.

        Returns:
            List[int]: The list of integer tokens representing the siRNA sequence.
        """
        nt_conversion = {"G": 1, "C": 2, "A": 3, "U": 4}
        nt_list = []
        for index in range(self.max_sirna_len):
            if index < len(self.seq):
                nt_list.append(nt_conversion[self.seq[index]])
            else:
                nt_list.append(0)
        return nt_list
