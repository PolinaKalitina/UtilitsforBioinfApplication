from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def quality_check(quality_scores, quality_threshold):
    return min(quality_scores) >= quality_threshold

def length_check(sequence, length_bounds):
    return length_bounds[0] <= len(sequence) <= length_bounds[1]


def gc_check(sequence, gc_bounds):
    gc_content = 100 * gc_fraction(sequence)
    return gc_bounds[0] <= gc_content <= gc_bounds[1]

def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
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