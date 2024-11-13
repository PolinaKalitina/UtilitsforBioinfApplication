def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    '''
    Converts multiline FASTA-file to ineline.

    Parametrs
    ---------------------
    input_file : file
        File in .fasta format

    Returns
    ---------------------
    output_file : file
        File in .fasta format.
        If file not provided, input_file is modified and returned.

    '''
    with open(input_fasta, 'r') as original_input:
        current_sequence = ''
        for line in original_input:
            if line.startswith('>'):
                current_sequence += '\n' + line
            elif not line.startswith('>'):
                current_sequence += line.strip()
        if output_fasta is not None:
            with open(output_fasta, 'w') as output_file:
                output_file.write(current_sequence)
                return (output_file)
        else:
            with open(input_fasta, 'w+') as output_file:
                output_file.write(current_sequence)
                return (output_file)


def parse_blast_output(input_file, output_file):
    '''
    Selects names of genes with significant alignments: first line of Discription for every Query.

    Parametrs
    ---------------------
    input_file : file
        File in .txt format

    Returns
    ---------------------
    output_file : file
        File in .txt format with discriptions sorted in alphabetical order.

    '''
    with open(input_file) as search_file:
        search_string = 'Sequences producing significant alignments:'
        discription = []
        n = 0
        line_to_read = 0
        line_to_check = ''
        for line in search_file:
            n += 1
            if search_string in line:
                line_to_read = n + 3
            elif n == line_to_read:
                line_to_check = line.strip()
                if line_to_check.endswith('.1'):
                    line = line_to_check.split(']')
                    final_line = line[0]
                    discription.append(final_line + ']' + '\n')
        for line in discription:
            line.strip()
            if line.endswith('.1]\n'):
                new_line = line.split('.. ')[0] + ']' + '\n'
                ind = discription.index(line)
                discription.pop(ind)
                discription.insert(ind, new_line)
        discription_no_double = []
        discription.sort(key=str.lower)
        for line in discription:
            if line not in discription_no_double:
                discription_no_double.append(line)
        print(discription_no_double)
        with open(output_file, 'a') as output_file_to_write:
            for line in discription_no_double:
                output_file_to_write.write(line)
            return (output_file)
