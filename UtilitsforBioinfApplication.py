# importing modules for run_dna_rna_tools and filter_fasqc from helpertools folder in tha same derictory
from helpertools.dna_rna_tools_module import is_NA, NA_type, rev_compl, compl, reverse, transcr
from helpertools.fastqc_module import gc_check, length_check, quality_check


def run_dna_rna_tools(*args):  # accepts a list of any number of strings - DNA or RNA sequences.
    # The last string should specify the operation to be applied to provided sequences - transcribe, reverse, complement or reverse complement
    *seqs, op = args  # separating sequences and an operator
    res = []

    for seq in seqs:
        res_test = ''
        res_test = is_NA(seq)  # checking if provided seqence is in correct format
        if res_test is False:
            print("Not a nucleic acid!")
        res_test = NA_type(seq)  # checking the type of NA of provided sequence
        if res_test is False:
            print("Not a nucleic acid!")
        elif op == "transcribe":  # transcribing sequence
            res.append(transcr(seq))
        elif op == "reverse":  # reverting sequence
            res.append(reverse(seq))
        elif op == "complement":  # creating a complement to sequence
            res.append(compl(seq))
        elif op == "reverse_complement":  # creating a complement to a reversed sequence
            res.append(rev_compl(seq))
        else:
            print("Operation specified incorrectly!")
            return ()

    if len(res) == 1:  # if only one sequence was given, resulting sequence returns as a string
        return (''.join(res))
    else:
        return (res)  # if more than one sequence was given, resulting sequences return as a list


seqs = {}

# accepts four arguments:
# 'seqs' - a dictionary of FASTQC-sequences in format {name of the sequence: (sequence, quality)}, name is a string, sequence and its quality are strings in a tuple;
# 'gc_bounds' - a tuple of lower and upper bounds for gc-content in provided sequences or an int - only upper bound. If not specified is (0, 100) by default;
# 'length_bounds' - a tuple of lower and upper bounds for length of provided sequences or an int - only upper bound. If not specified is (0, 2**32) by default;
# 'quality_threshold' - an int - an upper bound for quality of reads for provided sequences. If not specified is 0 by default;


def filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    gcbounds = [0]
    lbounds = [0]
    filtered_seqs = {}
    if isinstance(gc_bounds, int):  # sets lower bound if only upper is provided
        gcbounds.append(gc_bounds)
        gc_bounds = gcbounds
    if isinstance(length_bounds, int):   # sets lower bound if only upper is provided
        lbounds.append(length_bounds)
        length_bounds = lbounds
    for name, (sequence, quality) in list(seqs.items()):  # filters sequnces in dictionary by gc-content, length and quality, adds pairs that do pass in the new dictionary
        if gc_check(sequence, gc_bounds) is True and length_check(sequence, length_bounds) is True and quality_check(quality, quality_threshold) is True:
            filtered_seqs[name] = (sequence, quality)
    return (filtered_seqs)  # returns dictionary of reads within set parameters
