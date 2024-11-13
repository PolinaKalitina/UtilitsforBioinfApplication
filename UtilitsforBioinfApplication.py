from helpertools.dna_rna_tools_module import is_NA, NA_type, rev_compl, compl, reverse, transcr
from helpertools.fastqc_module import gc_check, length_check, quality_check, fastqc_to_dict, dict_to_fastqc


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
