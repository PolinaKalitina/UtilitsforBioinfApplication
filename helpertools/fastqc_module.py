def gc_check(sequence, gc_bounds):
    '''
    Checks if gc-content of the string is within set bounds
    Parametrs
    ---------------------
    seqence : string
        sequence of nucleotides from FASQC-dictionary
    length_bounds : tuple
        too numeric values -  lower and upper bounds for gc-content

    Returns
    ---------------------
    bool

    '''

    gc_count = 0
    gc_count = sum(True for nucl in list(sequence) if nucl == "G" or nucl == "C")
    gc = 0
    gc = gc_count / len(sequence)*100
    return (gc_bounds[1] >= gc and gc >= gc_bounds[0])


def length_check(sequence, length_bounds):
    '''
    Checks if length of the string is within set bounds
    Parametrs
    ---------------------
    seqence : string
        sequence of nucleotides from FASQC-dictionary
    length_bounds : tuple
        too numeric values -  lower and upper bounds for length

    Returns
    ---------------------
    bool

    '''
    return (length_bounds[1] >= len(sequence) and len(sequence) >= length_bounds[0])


def quality_check(quality, quality_threshold):
    '''
    Checks if quality of reads for the sequenceis above set threshold.

    Parametrs
    ---------------------
    quality : string
        sequence of symbols for quality of read for sequence from FASQC-dictionary
    quality_threshold : int
        quality threshold for provided string

    Returns
    ---------------------
    bool
    '''
    q_list = []
    q_score = 0
    for el in quality:
        q_list.append((ord(str(el))) - 33)
    q_score = sum(q_list) / len(q_list)
    return q_score >= quality_threshold


seqs = {}


def fastqc_to_dict(input_fastq):
    '''
    Converts FASTQC-file to dictionary.

    Parametrs
    ---------------------
    input_fastq : file

    Returns
    ---------------------
    seqs : dict
        Dictionary with structure (name : (sequence, quality))
    '''
    with open(input_fastq) as inpt:
        index_name = 1
        index_seq = 2
        index_qual = 4
        name = ''
        sequence = ''
        quality = ''
        for i, line in enumerate(inpt, 1):
            if i == index_name:
                name = line.strip()
                index_name = i + 4
                continue
            if i == index_seq:
                sequence = line.strip()
                index_seq = i + 4
                continue
            if i == index_qual:
                quality = line.strip()
                index_qual = i + 4
                seqs[name] = (sequence, quality)
                name = ''
                sequence = ''
                quality = ''
                continue
            else:
                continue
        return (seqs)


def dict_to_fastqc(filtered_seqs, output_fastq):
    '''
    Converts dictionary to FASTQC-file.

    Parametrs
    ---------------------
    filtered_seqs : dict
        Dictionary with structure (name : (sequence, quality))

    Returns
    ---------------------
    output_fastq : file
    '''
    with open(output_fastq, 'a') as fastq:
        name_plus = ''
        el = '+'
        for name, (sequence, quality) in filtered_seqs.items():
            name_plus = el + name[1:]
            fastq.write(f'{name}\n')
            fastq.write(f'{sequence}\n')
            fastq.write(f'{name_plus}\n')
            fastq.write(f'{quality}\n')
        return (output_fastq)
