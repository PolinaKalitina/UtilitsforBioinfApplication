from abc import ABC, abstractmethod
from typing import Dict, Tuple, Union


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
