import source.fastq_read as fr
import source.protein_tools as pt
import source.nucleic_acids_tools as nat


# Function for filtering FASTQ file
def fastq_filter(seqs: dict, gc_bound: tuple = (0, 100), length_bound: tuple = (0, 2**32), quality_threshold: int = 0) -> dict:
    """
    This function work with FASTQ files and filters them by
    GC content, length and Q-score.

    Arguments:
    - seqs: dictionary consisting of fastq sequences. The structure is as follows.
    Key: string, sequence name. Value: tuple of two strings (sequence and quality).
    - gc_bound: tuple of required range of GC percentage (inclusive).
    num of float if only higher border of the range is neede (exclusive).
    - length_bound: tuple of required range of sequence length (inclusive).
    - quality_threshold: int of lowest level of Q-score (inclusive).

    Output:
        - dictionary of those samples, which match all conditions.
    """
    result_dict = seqs.copy()
    for value in seqs:
        fastq = seqs[value][0]
        if isinstance(gc_bound, int) or isinstance(gc_bound, float):
            gc_check = fr.gc_content(fastq) >= gc_bound
        else:
            gc_check = fr.gc_content(fastq) < gc_bound[0] or fr.gc_content(fastq) > gc_bound[1]
        len_check = fr.seq_length(fastq) < length_bound[0] or fr.seq_length(fastq) > length_bound[1]
        quality_check = fr.quality_score(fastq) < quality_threshold
        if gc_check or len_check or quality_check:
            del result_dict[value]      
    return result_dict


# Function for working with protein sequences
def protein_tools(*arguments: str):
    """
    Main function to perform various actions on protein sequences.

    Arguments:
    - *arguments: Variable number of arguments. The first n-1 arguments should be protein sequences,
             and the last argument should be a string specifying the action to be performed.

    Output:
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

    action = arguments[-1]
    action_list = {
        "get_pI": pt.get_pI,
        "needleman_wunsch": pt.needleman_wunsch,
        "build_scoring_matrix": pt.build_scoring_matrix,
        "calculate_aa_freq": pt.calculate_aa_freq,
        "translate_protein_rna": pt.translate_protein_rna,
        "convert_to_3L_code": pt.convert_to_3L_code,
        "protein_mass": pt.protein_mass,
    }

    if action not in action_list:
        raise ValueError(f"No such action: {action}")

    if not (
        action == "needleman_wunsch"
        and len(arguments) == 3
        or action != "needleman_wunsch"
        and len(arguments) == 2
    ):
        raise ValueError("Error in number of sequences")

    for sequence in arguments[:-1]:
        if not all([letter.capitalize() in pt.AMINO_LETTERS for letter in sequence]):
            raise ValueError(f"The sequence is not protein sequence: {sequence}")

    result = action_list[action](*arguments[:-1])

    return result


# Function for working with DNA/RNA sequences
def nucl_acid_tools(*arguments: str):
    """
    This function works with DNA or RNA sequences

    Arguments:
    - *arguments: Variable number of arguments. The first n-1 arguments should be DNA or RNA sequences,
    and the last argument should be a string specifying the action to be performed.

    Output:
    - The result of the specified action on the input nucleic acid sequences.

    Raises:
    - ValueError: If the specified action is not supported or 
    if there is an error in the number of sequences.
    Also raised if the input sequences are not valid DNA or RNA sequences.
    
    Supported Actions:
    - "transcribe": Transcribes DNA sequence to RNA sequence.
    - "reverse": Reverses DNA or RNA sequences.
    - "complement": Makes the complement sequence. Takes only DNA sequences.
    - "reverse_complement": Gives reversed and complement to inputed DNA sequence.
    - "make_binary": Recodes DNA or RNA sequence in binary code.
    First number is always reffered to nucleic acid type: 0 - DNA, 1 - RNA.
    """
    list_of_seq = arguments[:-1]
    action = arguments[-1]
    output_seq_list = []
    FUNC_DICT = { # Function dictionary
        'transcribe': nat.transcribe,
        'reverse': nat.reverse,
        'complement': nat.complement,
        'reverse_complement': nat.reverse_complement,
        'make_binary': nat.make_binary
        }

    if action in FUNC_DICT:
        for seq in list_of_seq:
            function = FUNC_DICT[action]
            output_seq_list.append(function(seq))
        if len(output_seq_list) == 1:
            return output_seq_list[0]
        else:
            return output_seq_list
    else:
        raise ValueError("Function is not found")
