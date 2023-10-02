# Function to count GC content
def gc_content(seq: str) -> float:
    """
    This function calculates percent of GC nucleotides in sequence

    Arguments:
        seq = sequence of DNA

    Output:
        percent of GC nucleotides from the whole sequence 
    """
    gc = 0
    for nucleotide in seq:
        if nucleotide == 'G' or nucleotide == 'C':
            gc += 1
    return round(gc/len(seq)*100, 3)






        