# **SequenceMaster**

_A powerful and minimalistic bioinfomatical tool to work with sequences of nucleic acids, proteins such as for work with FASTQ files_

<div style='justify-content: center'>
<img src="img/SequenceMaster.png" align='center', width="100%">
</div>

## Key Features

- **Nucleic Acid Analysis:** SequenceMaster offers comprehensive tools for analyzing DNA and RNA sequences. 

- **Protein Analysis:** Explore protein sequences with advanced tools for pI calculation, protein mass calculation, 

- **FASTQ File Manipulation:** Easily filter FASTQ files depending on quality, sequence length and GC content

## Getting Started

To get started with SequenceMaster, visit our [Documentation](https://yourdocumentationlink.com) for detailed usage instructions and examples.

## Installation

To get the tool clone the git repository::

```bash
git clone https://github.com/SuleimanovShakir/SequenceMaster.git && cd SequenceMaster
```

## Usage

### Nucleic acids

As we live in the world of **central dogma of molecular biology**, this function checks is the sequence DNA or RNA. This is an important limit because some of the functions may use only DNA sequence in it's work.

To start working with nucleic acids, just run this command:

```{python}
nucleic_acid_tools(*sequences: str, action: str)
```
- Firstly, you have to enter your DNA or RNA sequences just as 'str' type as positional arguments
- Last arguments is always an action that you prefer to perform to your sequences. Here is a list of possible actions
```
Supported Actions:
- "transcribe": Transcribes DNA sequence to RNA sequence.
- "reverse": Reverses DNA or RNA sequences.
- "complement": Makes the complement sequence. Takes only DNA sequences.
- "reverse_complement": Gives reversed and complement to inputed DNA sequence.
- "make_binary": Recodes DNA or RNA sequence in binary code.First number is always reffered to nucleic acid type: 0 - DNA, 1 - RNA.
```

### Proteins

To start working with proteins, just run this command:
```{python}
protein_tools(*sequences: str, action: str)
```

- Firstly, you have to enter your protein sequence just as 'str' type as positional arguments. You have to use one-letter coding.
- Last arguments is always an action that you prefer to perform to your sequences. Here is a list of possible actions

```
Supported Actions:
- "get_pI": Calculate isoelectric points for each amino acid in the sequence.
- "calculate_aa_freq": Calculate the frequency of each amino acid in a protein sequence.
- "translate_protein_rna": Translate amino acid sequence to RNA, using random codons for each amino acid.
- "three_letter_code": Convert one-letter amino acid sequence to three-letter coding.
- "protein_mass": Calculate the molecular weight of the protein sequence.
```

### FASTQ

Last but not least - `fastq_filter`

To start filtering your FASTQ files, just run this command:
```{python}
fastq_filter(seqs: dict, gc_bound: Union[tuple, int, float] = (0, 100),
                 length_bound: tuple = (0, 2**32), quality_threshold: Union[int, float] = 0) -> dict
```
List of arguments:
- Firstly, you have to enter dictionary  consisting of fastq sequences. The structure is as follows.
Key: string, sequence name. Value: tuple of two strings (sequence and quality).
- gc_bound: tuple of required range of GC percentage (inclusive), num or float if only higher border of the range is needed (exclusive). (0,100) is a defaul range of GC content.
- length_bound: tuple of required range of sequence length (inclusive). (0,2**32) is a defaul range of length.
- quality_threshold: int of lowest level of Q-score (inclusive). 0 is a default Q-score

There is only one action - filter input dictionary by this parameters and give but dictionary with those entries that meet the conditions.

## Examples

### Nucleic acids tools
```{python}
print(nucl_acid_tools('ATG', action='transcribe')) -> 'AUG'
print(nucl_acid_tools('ATG', 'GAT', action='reverse')) -> ['GTA', 'TAG']
print(nucl_acid_tools('ATG', 'tAg', 'AAA', action='complement')) -> ['TAC', 'aTc', 'TTT']
print(nucl_acid_tools('ATG', 'tAg', 'AAA', action='reverse_complement')) -> ['CAT', 'cTa', 'TTT']
print(nucl_acid_tools('ATG', action='make_binary')) -> '0,0,0,0,1,1,0'
```

### Proteins tools
```{python}
print(protein_tools('KKLMN', action='calculate_aa_freq')) -> {'K': 2, 'L': 1, 'M': 1, 'N': 1}
print(protein_tools('KLMN','PRST', action='translate_protein_rna')) -> ['AAGCUAAUGAAC', 'CCAAGGAGCACG']
print(protein_tools('KLMN', 'pRsT', action='three_letter_code')) -> ['Lys-Leu-Met-Asn', 'Pro-Arg-Ser-Thr']
print(protein_tools('KLMN', action='protein_mass')) -> 486.26244
```

### FASTQ filter
```{python}
print(fastq_filter(sequences_dict, (40,60), (0, 200), 25)) -> filtered_sequences_dict
```

## Contacts 
Shakir Suleimanov,
Please, do not hesitate to contact me via [Git-Hub](https://github.com/SuleimanovShakir) or [e-mail](suleymanovef@gmail.com).

