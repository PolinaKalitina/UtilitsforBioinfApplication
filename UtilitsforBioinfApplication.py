

def transcr(seq):  # accepts a string - a sequence of nucleotides of any length
    newseq = ''
    for nucl in seq:
        if nucl == 'T':
            newseq += 'U'
        elif nucl == 't':
            newseq += 'u'
        else:
            newseq += nucl
    return (newseq)  # returns a string - a sequence translated to RNA form


def reverse(seq):  # accepts a string - a sequence of nucleotides of any length
    newseq = ''
    newseq = seq[::-1]
    return (newseq)  # returns a string - a reversed sequence


def compl(seq):   # accepts a string - a sequence of nucleotides of any length
    newseq = ''
    dna_comp_dict = {'A': 'T', 'a': 't', 'T': 'A', 't': 'a', 'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c'}
    rna_comp_dict = {'A': 'U', 'a': 'u', 'U': 'A', 'u': 'a', 'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c'}
    res_test = NA_type(seq)
    if res_test == 'DNA':
        newseq = ''.join([dna_comp_dict[nucl] for nucl in seq])
    else:
        newseq = ''.join([rna_comp_dict[nucl] for nucl in seq])
    return (newseq)  # reterns a string - a comlement of given sequence


def rev_compl(seq):  # accepts a string - a sequence of nucleotides of any length
    res1 = reverse(seq)
    res = compl(res1)
    return (res)  # returns a string - a complement of reversed given sequence


def NA_type(seq):  # accepts a string - a sequence of nucleotides of any length
    d_list = ['T', 't']
    r_list = ['U', 'u']
    n_type = ''
    for nucl in seq:
        if nucl in d_list:
            if n_type != 'RNA':
                n_type = 'DNA'
            elif n_type == 'RNA':
                n_type = False
                break
        elif nucl in r_list:
            if n_type != 'DNA':
                n_type = 'RNA'
            elif n_type == 'DNA':
                n_type = False
                break
    return (n_type)  # returns a string - a type of nucleic acid that is a given string, if the type can not be defined, returns None


def is_NA(seq):  # accepts a string - a sequence of symbols of any length
    n_list = ['T', 't', 'U', 'u', 'A', 'a', 'G', 'g', 'C', 'c']
    for nucl in seq:
        if nucl not in n_list:
            return (False)  # returns 'False' if provided string contains non-nucleotide type of characters or is a mix of DNA and RNA, otherwise returns None



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






def run_dna_rna_tools(*args):

    '''
    Operates on NA strings.
    Uses module run_dna_rna_tools from helpertools folder in tha same directory.

    Parametrs
    ---------------------
    *seqs : list or str
        Any number of strings - DNA or RNA sequences.
    op : str
        Operation to be applied to provided sequences: transcribe, reverse, complement or reverse complement

    Returns
    ---------------------
    res : list or str
        if only one sequence was given, resulting sequence returns as a string, more than one return in list
    '''

    *seqs, op = *args
    res = []

    for seq in seqs:
        res_test = ''
        res_test = is_NA(seq)
        if res_test is False:
            print("Not a nucleic acid!")
        res_test = NA_type(seq)
        if res_test is False:
            print("Not a nucleic acid!")
        elif op == "transcribe":
            res.append(transcr(seq))
        elif op == "reverse":
            res.append(reverse(seq))
        elif op == "complement":
            res.append(compl(seq))
        elif op == "reverse_complement":
            res.append(rev_compl(seq))
        else:
            print("Operation specified incorrectly!")
            return ()

    if len(res) == 1:
        return (''.join(res))
    else:
        return (res)


def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    '''
    Filters sequnces in dictionary by gc-content, length and quality.
    Uses module filter_fasqc from helpertools folder in tha same derictory.

    Parametrs
    ---------------------
    input_fastq : file
        FASTQC-sequences id default .fastqc format
    output_fastq : file
        empty file to write input in
        If not provided, input_fastq is altered
    gc_bounds : tuple or int
        Lower and upper bounds for gc-content in provided sequences
        If an int - only upper bound
        If not specified is (0, 100) by default
    length_bounds : tuple or int
        Lower and upper bounds for length of provided sequences
        If an int - only upper bound
        If not specified is (0, 2**32) by default
    quality_threshold : int
        Upper bound for quality of reads for provided sequences
        If not specified is 0 by default

    Returns
    ---------------------
    filtered_seqs : file
    FASTQC-file
    '''
    seqs = fastqc_to_dict(input_fastq)
    gcbounds = [0]
    lbounds = [0]
    filtered_seqs = {}
    if isinstance(gc_bounds, (int, float)):
        gcbounds.append(gc_bounds)
        gc_bounds = gcbounds
    if isinstance(length_bounds, (int, float)):
        lbounds.append(length_bounds)
        length_bounds = lbounds
    for name, (sequence, quality) in list(seqs.items()):
        if gc_check(sequence, gc_bounds) and length_check(sequence, length_bounds) and quality_check(quality, quality_threshold):
            filtered_seqs[name] = (sequence, quality)
    dict_to_fastqc(filtered_seqs, output_fastq)
    return (filtered_seqs)
