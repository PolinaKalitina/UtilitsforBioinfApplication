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
