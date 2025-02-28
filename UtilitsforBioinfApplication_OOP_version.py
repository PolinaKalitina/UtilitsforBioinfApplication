from abc import ABC, abstractmethod


class BiologicalSequence(abs.ABS):

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        self.sequence = sequence
        self.complement_dict = {}

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
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

    def __str__(self):
        return f"Nucleic Acid Sequence: {self.sequence}"

    def check_alphabet(self):
        return all(nucleotide in self.complement_dict for nucleotide in self.sequence)

    def complement(self):
        return ''.join([self.complement_dict[nucl] for nucl in self.sequence])

    def reverse(self):
        return self.sequence[::-1]

    def reverse_complement(self):
        return self.complement()[::-1]


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)
        self.complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def transcribe(self):
        return RNASequence(self.sequence.replace('T', 'U'))


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)
        self.complement_dict = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.sequence[index]
        elif isinstance(index, tuple):
            idx, slice_type = index
            if slice_type == "prefix":
                return self.sequence[:idx]
            elif slice_type == "suffix":
                return self.sequence[idx:]

    def __str__(self):
        return f"Amino Acid Sequence: {self.sequence}"

    def check_alphabet(self):
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        return all(aa in amino_acids for aa in self.sequence)

    def molecular_weight(self):
        aa_weights = {
            'A': 89.09, 'C': 121.16, 'D': 133.10, 'E': 147.13, 'F': 165.19,
            'G': 75.07, 'H': 155.16, 'I': 131.17, 'K': 146.19, 'L': 131.17,
            'M': 149.21, 'N': 132.12, 'P': 115.13, 'Q': 146.15, 'R': 174.20,
            'S': 105.09, 'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19
        }
        return sum(aa_weights[aa] for aa in self.sequence)
