# UtilitsForBioinfApplication 

## Discription
Utilits for bioinformatical application is a script that consists of four main functions, `run_dna_rna_tools` and `filter_fastqc` in `UtilitsForBioinfApplication` script and `convert_multiline_fasta_to_oneline` and `parse_blast_output` in `bio_files_processor` script.

**run_dna_rna_tools** works with DNA or RNA sequences. 
It can:
 - transcribe
 - reverse
 - create a complement for rovided sequence
 - create a complement foror reversed sequence
All operations will keep the register of the characters - see example.

**filter_fastqc** works with FASTQC-files that are transformed into a dictionary of FASTQC-sequences in format `{name of the sequence: (sequence, quality)}`, name is a string, sequence and its quality are strings in a tuple.
Sequences are filtered by their:
 - GC-content
 - length
 - quality of reads

Both functions use modules from helpertools folder to perform.

**convert_multiline_fasta_to_oneline** works with multiline FASTA-files and converts them to oneline.

**parse_blast_output** works with BLAST search output and selects names of proteins with significant alignments.

## Input and output format

For **run_dna_rna_tools**:

**Input:** list of any number of strings - DNA or RNA sequences. The last string should specify the operation to be applied to provided sequences - `transcribe'`, `reverse`, `complement` or `reverse complement`.

*Input example:* `['AtGC', 'aGGGT', 'transcribe']`

**Output:** changed sequences as a list, if only one sequence was given, it will return as a string. Operation will not be returned.

*Output example:* `['AuGC', 'aGGGU']` 

For **filter_fastqc**:

**Input:** five arguments:
`input_fastq` - FASTQC-file in .fastq format;
`output_fastq` - empty file;
`gc_bounds` - a tuple of lower and upper bounds for gc-content in provided sequences or an int - only upper bound. If not specified is `(0, 100)` by default;
`length_bounds` - a tuple of lower and upper bounds for length of provided sequences or an int - only upper bound. If not specified is `(0, 2**32)` by default;
`quality_threshold` - an int - an upper bound for quality of reads for provided sequences. If not specified is `0` by default;

*Input example:* `'example_fastq.fastq', 'output_fastq.fastq', gc_bounds=(0, 100), length_bounds=(0, 35)`

**Output:** FASTQ-file in .fastq format with sequences that passed by GC-content, length and read quality in format `{name of the sequence: (sequence, quality)}`. If no reads passed, empty file will be returned.

*Output example:* `'output_fastq.fastq'`

For **convert_multiline_fasta_to_oneline**:

**Input:** two arguments
`input_fasta` - multiline FASTA file;
`output_fasta` - empty file to write sequences in. Is optional, if not provided, `input_fasta` is transformed.

*Input format example:*
`>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCT`

`GGTGCTGTG`

**Output:** FASTA-file.

*Output format example:* 
`>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG`

For **parse_blast_output**:

**Input:** two arguments:
`input_file` - .txt file with BLAST search results;
`output_file` - empty file to write sequences in.

*Output:** .txt file with alphabetically ordered list of proteins with significant alignments.

*Output format example:* 
`'conjugal transfer protein TraA [Enterobacteriaceae]
 'conjugal transfer protein TraC [Enterobacteriaceae]
 'DinI-like family protein [Escherichia coli]`

## Limitations
**run_dna_rna_tools** can only work with one typa of NA at a time, you can not provide both in one run. It will also not accept mixed sequences with both 'U' and 'T' present and will not accept sequences if operation is not specified.

**parse_blast_output** is not immune to BLAST output modifications. Strings in output may include additional information or syntax defects.

## Developmend
By Polina Kalitina (iduvzavtra@gmail.com), shaking in the corner.
