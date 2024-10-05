def gc_check(sequence, gc_bounds):  # accepts a string - sequence of nucleotides from FASQC-dictionary - and a tuple of too numeric values -  lower and upper bounds for gc-content in provided string
    gc_count = 0
    gc_count = sum(True for nucl in list(sequence) if nucl == "G" or nucl == "C")
    gc = 0
    gc = gc_count/len(sequence)*100
    if gc_bounds[1] >= gc and gc >= gc_bounds[0]:
        return (True)  # returns 'True' if gc-content of the string is within set bounds


def length_check(sequence, length_bounds):  # accepts a string - sequence of nucleotides from FASQC-dictionary - and a tuple of too numeric values -  lower and upper bounds for length of provided string
    if length_bounds[1] >= len(sequence) and len(sequence) >= length_bounds[0]:
        return (True)  # returns 'True' if length of the string is within set bounds


def quality_check(quality, quality_threshold):  # accepts a string - sequence of symbols for quality of read for sequence from FASQC-dictionary - and a numeric value -  quality threshold for provided string
    q_list = []
    q_score = 0
    for el in quality:
        q_list.append((ord(str(el))) - 33)
    q_score = sum(q_list)/len(q_list)
    if q_score >= quality_threshold:
        return (True)   # returns 'True' if quality of the read is above set threshold
