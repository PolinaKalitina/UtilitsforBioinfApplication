from abc import ABC, abstractmethod
from typing import Dict, Tuple, Union
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

class BiologicalSequence(ABC):

    @abstractmethod
    def __len__(self) -> int:
        pass

    @abstractmethod
    def __getitem__(self, index: Union[int, slice, Tuple[int, str]]) -> str:
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def check_alphabet(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str) -> None:
        self.sequence: str = sequence
        self.complement_dict: Dict[str, str] = {}

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index: Union[int, slice, Tuple[int, str]]) -> str:
        if isinstance(index, int):
            return self.sequence[index]
        elif isinstance(index, slice):
            return self.sequence[index]
        elif isinstance(index, tuple):
            idx, slice_type = index
            if slice_type == "prefix":
                return self.sequence[:idx]
            elif slice_type == "suffix":
                return self.sequence[idx:]

    def __str__(self) -> str:
        return f"Nucleic Acid Sequence: {self.sequence}"

    def check_alphabet(self) -> bool:
        return all(nucleotide in self.complement_dict for nucleotide in self.sequence)

    def complement(self) -> str:
        return ''.join([self.complement_dict[nucl] for nucl in self.sequence])

    def reverse(self) -> str:
        return self.sequence[::-1]

    def reverse_complement(self) -> str:
        return self.complement()[::-1]


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str) -> None:
        super().__init__(sequence)
        self.complement_dict: Dict[str, str] = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def transcribe(self) -> 'RNASequence':
        return RNASequence(self.sequence.replace('T', 'U'))


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str) -> None:
        super().__init__(sequence)
        self.complement_dict: Dict[str, str] = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str) -> None:
        self.sequence: str = sequence

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index: Union[int, Tuple[int, str]]) -> str:
        if isinstance(index, int):
            return self.sequence[index]
        elif isinstance(index, tuple):
            idx, slice_type = index
            if slice_type == "prefix":
                return self.sequence[:idx]
            elif slice_type == "suffix":
                return self.sequence[idx:]

    def __str__(self) -> str:
        return f"Amino Acid Sequence: {self.sequence}"

    def check_alphabet(self) -> bool:
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        return all(aa in amino_acids for aa in self.sequence)

    def molecular_weight(self) -> float:
        aa_weights = {
            'A': 89.09, 'C': 121.16, 'D': 133.10, 'E': 147.13, 'F': 165.19,
            'G': 75.07, 'H': 155.16, 'I': 131.17, 'K': 146.19, 'L': 131.17,
            'M': 149.21, 'N': 132.12, 'P': 115.13, 'Q': 146.15, 'R': 174.20,
            'S': 105.09, 'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19
        }
        return sum(aa_weights[aa] for aa in self.sequence)

def quality_check(quality_scores: List[int], quality_threshold: int) -> bool:
    return min(quality_scores) >= quality_threshold


def length_check(sequence: str, length_bounds: Tuple[int, int]) -> bool:
    return length_bounds[0] <= len(sequence) <= length_bounds[1]


def gc_check(sequence: str, gc_bounds: Tuple[float, float]) -> bool:
    gc_content = 100 * gc_fraction(sequence)
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Tuple[float, float] = (0, 100),
    length_bounds: Tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0,
) -> int:
    # Чтение последовательностей из файла
    sequences = SeqIO.parse(input_fastq, "fastq")

    # Фильтрация
    filtered_sequences = []
    for record in sequences:
        sequence = str(record.seq)
        quality_scores = record.letter_annotations["phred_quality"]

        if (gc_check(sequence, gc_bounds) and
            length_check(sequence, length_bounds) and
            quality_check(quality_scores, quality_threshold)):
            filtered_sequences.append(record)

    # Запись отфильтрованных последовательностей в новый файл
    count = SeqIO.write(filtered_sequences, output_fastq, "fastq")
    print(f"Saved {count} reads to {output_fastq}")



