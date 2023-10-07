# Additional script for working with fastq files

# Function to count GC content
def gc_content(seq: str) -> float:
    """
    This function calculates percent of GC nucleotides in sequence.

    Arguments:
        seq = sequence of DNA.

    Output:
        percent of GC nucleotides from the whole sequence.
    """
    gc = 0
    for nucleotide in seq:
        if nucleotide == 'G' or nucleotide == 'C':
            gc += 1
    return round(gc/len(seq)*100, 3)


# Function to calculate the length of the sequence
def seq_length(seq: str) -> float:
    """
    This function calculates the length of the sequence.

    Arguments:
        seq = sequence of DNA.

    Output:
        length of the sequence in number of nucleotides
    """
    return len(seq)


# Function to calculate quality score (Q-score)
def quality_score(quality_seq: str) -> float:
    """
    This function calculates the Q-score for the sequence.

    Arguments:
        quality_seq = sequence with the data of Q-score for each nucleotide.

    Output:
        Q-score.
    """
    score = 0
    for symbol in quality_seq:
        score += ord(symbol) - 33
    return score/len(quality_seq)
