from abc import ABC, abstractmethod
from typing import Dict, Tuple, Union, List
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import argparse
from loguru import logger


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
    if not quality_scores:
        return False
    return min(quality_scores) >= quality_threshold


def length_check(sequence: str, length_bounds: Tuple[int, int]) -> bool:
    if not sequence:
        return False
    return length_bounds[0] <= len(sequence) <= length_bounds[1]


def gc_check(sequence: str, gc_bounds: Tuple[float, float]) -> bool:
    if not sequence:
        return False
    if gc_bounds[0] < 0 or gc_bounds[1] > 100:
        raise ValueError("GC bounds must be between 0 and 100")
    gc_content = 100 * (sequence.count('G') + sequence.count('C')) / len(sequence)
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def parse_args():
    parser = argparse.ArgumentParser(description='FASTQ quality filter')
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTQ file')
    parser.add_argument('--gc-bounds', type=float, nargs=2, default=[0, 100],
                        metavar=('MIN', 'MAX'), help='GC content bounds (0-100%)')
    parser.add_argument('--length-bounds', type=int, nargs=2, default=[0, 2**32],
                        metavar=('MIN', 'MAX'), help='Length bounds')
    parser.add_argument('-q', '--quality', type=int, default=0,
                        help='Minimum quality threshold')
    parser.add_argument('--log', default='fastq_filter.log', help='Log file path')
    return parser.parse_args()


def setup_logging(log_file: str):
    logger.add(log_file, rotation="10 MB", level="INFO")


def is_fastq_file(filename: str) -> bool:
    try:
        with open(filename) as f:
            for i, line in enumerate(f):
                if i >= 4:
                    break
                if i == 0 and not line.startswith('@'):
                    return False
                if i == 2 and not line.startswith('+'):
                    return False
        return True
    except:
        return False


def filter_fastq():
    args = parse_args()
    if not is_fastq_file(args.input):
        logger.error(f"Input file {args.input} is not in FASTQ format")
        raise ValueError("Input file is not in FASTQ format")

    setup_logging(args.log)

    try:
        logger.info(f"Starting filtering with params: {vars(args)}")
        try:
            sequences = list(SeqIO.parse(args.input, "fastq"))
            total_reads = len(sequences)
            logger.info(f"Total reads: {total_reads}")
        except FileNotFoundError:
            logger.error(f"Input file {args.input} not found")
            raise
        except Exception as e:
            logger.error(f"Error reading input file: {str(e)}")
            raise
        
        stats = {
            'total': total_reads,
            'passed': 0,
            'failed_quality': 0,
            'failed_length': 0,
            'failed_gc': 0
        }

        filtered = []
        for record in sequences:
            seq = str(record.seq)
            quals = record.letter_annotations.get("phred_quality", [])

            quality_ok = quality_check(quals, args.quality)
            length_ok = length_check(seq, args.length_bounds)
            gc_ok = gc_check(seq, args.gc_bounds)

            if not quality_ok:
                stats['failed_quality'] += 1
            if not length_ok:
                stats['failed_length'] += 1
            if not gc_ok:
                stats['failed_gc'] += 1

            if quality_ok and length_ok and gc_ok:
                filtered.append(record)
                stats['passed'] += 1

        logger.info("\nFiltering statistics:")
        logger.info(f"Total reads: {stats['total']}")
        logger.info(f"Passed filters: {stats['passed']} ({stats['passed']/stats['total']:.1%})")
        logger.info(f"Failed by quality: {stats['failed_quality']}")
        logger.info(f"Failed by length: {stats['failed_length']}")
        logger.info(f"Failed by GC content: {stats['failed_gc']}")

        if filtered:
            logger.info("\nFirst 5 passed reads:")
            for record in filtered[:5]:
                gc_content = 100 * (record.seq.count('G') + record.seq.count('C')) / len(record.seq)
                logger.info(f"ID: {record.id}, Length: {len(record.seq)}, GC: {gc_content:.1f}%, "
                          f"Min quality: {min(record.letter_annotations.get('phred_quality', []))}")
        else:
            logger.warning("No reads passed all filters!")

        count = SeqIO.write(filtered, args.output, "fastq")
        logger.info(f"\nSaved {count} reads to {args.output}")
        return count

    except Exception as e:
        logger.error(f"Processing failed: {str(e)}")
        raise


if __name__ == "__main__":
    filter_fastq()
