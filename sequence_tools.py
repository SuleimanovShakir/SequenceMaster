from source import fastq_read as fr


# Function for filtering FASTQ file. 
def fastq_filter(seqs: dict, gc_bound: tuple = (0, 100), length_bound: tuple = (0, 2**32), quality_threshold: float = 0) -> dict:
    dict = seqs.copy()
    for value in seqs:
        fastq = seqs[value][0]
        if type(gc_bound) == int or type(gc_bound) == float:
            gc_check = fr.gc_content(fastq) >= gc_bound
        else: 
            gc_check = fr.gc_content(fastq) < gc_bound[0] or fr.gc_content(fastq) > gc_bound[1]
        len_check = fr.seq_length(fastq) < length_bound[0] or fr.seq_length(fastq) > length_bound[1]
        quality_check = fr.quality_score(fastq) < quality_threshold
        if gc_check or len_check or quality_check:
            del dict[value]      
    return dict


# Function for working with protein sequences
def main(*args: str):
    """
    Main function to perform various actions on protein sequences.

    Args:
    - *args: Variable number of arguments. The first n-1 arguments should be protein sequences,
             and the last argument should be a string specifying the action to be performed.

    Returns:
    - The result of the specified action on the input protein sequences.

    Raises:
    - ValueError: If the specified action is not supported or if there is an error in the number of sequences.
                  Also raised if the input sequences are not valid protein sequences.

    Supported Actions:
    - "get_pI": Calculate isoelectric points for each amino acid in the sequence.
    - "needleman_wunsch": Perform global alignment of two sequences using the Needleman-Wunsch algorithm.
    - "build_scoring_matrix": Build a scoring matrix for amino acid pairs.
    - "calculate_aa_freq": Calculate the frequency of each amino acid in a protein sequence.
    - "translate_protein_rna": Translate amino acid sequence to RNA, using random codons for each amino acid.
    - "convert_to_3L_code": Convert one-letter amino acid sequence to three-letter coding.
    - "protein_mass": Calculate the molecular weight of the protein sequence.
    """

    action = args[-1]
    action_list = {
        "get_pI": get_pI,
        "needleman_wunsch": needleman_wunsch,
        "build_scoring_matrix": build_scoring_matrix,
        "calculate_aa_freq": calculate_aa_freq,
        "translate_protein_rna": translate_protein_rna,
        "convert_to_3L_code": convert_to_3L_code,
        "protein_mass": protein_mass,
    }

    if action not in action_list:
        raise ValueError(f"No such action: {action}")

    if not (
        action == "needleman_wunsch"
        and len(args) == 3
        or action != "needleman_wunsch"
        and len(args) == 2
    ):
        raise ValueError("Error in number of sequences")

    for sequence in args[:-1]:
        if not all([letter.capitalize() in AMINO_LETTERS for letter in sequence]):
            raise ValueError(f"The sequence is not protein sequence: {sequence}")

    result = action_list[action](*args[:-1])

    return result
