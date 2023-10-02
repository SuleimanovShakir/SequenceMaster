import source.nucleic_acids_dict as dict


# Check if seq is DNA
def is_dna(seq: str) -> bool:
    unique_chars = set(seq)
    nucleotides = set('AaTtGgCc')
    return unique_chars <= nucleotides


# Check if seq is RNA
def is_rna(seq: str) -> bool:
    unique_chars = set(seq)
    nucleotides = set('AaUuGgCc')
    return unique_chars <= nucleotides


# Transcription function
def transcribe(seq: str) -> str:
    if not is_dna(seq):
        raise ValueError("Sequence is not a DNA, input should be DNA")
    else:
        return ''.join([dict.TRANSCRIPTION_DICT[i] for i in seq])


# Reverse function
def reverse(seq: str) -> str:
    if not is_dna(seq) or is_rna(seq):
        raise ValueError("Sequence is not a DNA/RNA, input should be DNA/RNA")
    else:
        return seq[::-1]


# Function for reverse seq
def complement(seq: str) -> str:
    if not is_dna(seq):
        raise ValueError("Sequence is not a DNA, input should be DNA")
    else:
        return ''.join([dict.COMPLMENTARITY_DICT[i] for i in seq])
         

# Function for make reverse and complement seq
def reverse_complement(seq: str) -> str:
    if not is_dna(seq):
        raise ValueError("Sequence is not a DNA, input should be DNA")
    else:
        return complement(reverse(seq))


# Function for rewriting nucleotide sequence to binary code
def make_binary(seq:str) -> tuple:
    if not is_dna(seq) or is_rna(seq):
        raise ValueError("Sequence is not a DNA/RNA, input should be DNA/RNA")
    else:
        binary_seq = ''.join([dict.BINARY_DICT[i] for i in seq])
        if is_dna(seq):
            result = '0,' + binary_seq[0:-1]
            return result
        if is_rna(seq):
            result = '1,' + binary_seq[0:-1]
            return result