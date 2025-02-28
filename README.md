# UtilitsForBioinfApplication 

## Discription
Utilities for Bioinformatics Application is a Python package that provides tools for working with biological sequences (DNA, RNA, and amino acids) and filtering FASTQ files. The package consists of the following components:

**Classes for Biological Sequences:**

- `BiologicalSequence`: Abstract base class defining the interface for biological sequences.

- `NucleicAcidSequence:` Base class for nucleic acid sequences (DNA and RNA). Provides methods for:

     - `complement`: Returns the complementary sequence.
  
     - `reverse`: Returns the reversed sequence.
  
     - `reverse_complement`: Returns the reverse complement of the sequence.

- `DNASequence`: Class for DNA sequences with methods for transcription. Adds a method:

     - `transcribe`: Transcribes DNA to RNA.

- `RNASequence`: Class for RNA sequences. Uses RNA-specific complement rules.

- `AminoAcidSequence`: Class for amino acid sequences with methods for calculating molecular weight. Provides a method:

     - `molecular_weight`: Calculates the molecular weight of the amino acid sequence.
      
For all types of sequences you can also check its length, estract an element by index or a slice from input string and check if the string contains non-NA or non-aa characters.

**Functions for FASTQ Filtering:**

`filter_fastq`: Filters sequences in a FASTQ file based on:

   - GC content: Specified as a tuple (min_gc, max_gc) or a single value (upper bound). If not specified is `(0, 100)` by default;

   - Length: Specified as a tuple (min_length, max_length) or a single value (upper bound). If not specified is `(0, 2**32)` by default;

   - Quality: Minimum quality score threshold. If not specified is `0` by default;

## Input and output format

For **UtilsforBioinfApplication_OOP**:

**Input:** sequences provided as strings via object initialization.

*Input example:* 

`dna = DNASequence("ATGC")`

`print(dna.complement())`

**Output:** methods return strings or lists of strings, depending on the operation.

*Output example:* `Output: TACG` 

For **FASTQ_filtr**:

**Input:** five arguments:
`input_fastq` -  Path to the input FASTQ file;
`output_fastq` - Path to the output FASTQ file;
`gc_bounds` - a tuple of lower and upper bounds for gc-content in provided sequences or an int - only upper bound. If not specified is `(0, 100)` by default;
`length_bounds` - a tuple of lower and upper bounds for length of provided sequences or an int - only upper bound. If not specified is `(0, 2**32)` by default;
`quality_threshold` - an int - an upper bound for quality of reads for provided sequences. If not specified is `0` by default;

*Input example:* `'example_fastq.fastq', 'output_fastq.fastq', gc_bounds=(0, 100), length_bounds=(0, 35)`

**Output:** FASTQ-file in .fastq format with sequences that passed by GC-content, length and read quality. If no reads passed, empty file will be returned.

*Output example:* `'output_fastq.fastq'`

## Limitations
**UtilsforBioinfApplication_OOP** can only work with one typa of NA at a time, you can not provide both in one run. It will also not accept mixed sequences with both 'U' and 'T' present.

## Development
By Polina Kalitina (iduvzavtra@gmail.com)

