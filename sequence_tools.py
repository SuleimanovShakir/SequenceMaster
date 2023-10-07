# Importing modules
from typing import Union
import source.fastq_read as fr
import source.protein_tools as pt
import source.nucleic_acids_tools as nat


# Function for filtering FASTQ file
def fastq_filter(seqs: dict, gc_bound: Union[tuple, int, float] = (0, 100),
                 length_bound: Union[tuple, int, float] = (0, 2**32),
                 quality_threshold: Union[int, float] = 0) -> dict:
    """
    This function work with FASTQ files and filters them by
    GC content, length and Q-score.

    Arguments:
    - seqs (dict): dictionary consisting of fastq sequences. The structure is as follows.
    Key: string, sequence name. Value: tuple of two strings (sequence and quality).
    - gc_bound (tuple, int, float): tuple of required range of GC percentage (inclusive),
    num or float if only higher border of the range is needed (exclusive).
    - length_bound (tuple, int, float): tuple of required range of sequences length (inclusive),
    num or float if only higher border of the range is needed (exclusive).
    - quality_threshold (int): int of lowest level of Q-score (inclusive).

    Output:
    - dictionary of those samples, which match all conditions.
    """
    if not len(seqs) > 0:
        raise ValueError('There are no fastq sequences')
    seqs_type = isinstance(seqs, dict)
    gc_bound_type = isinstance(gc_bound, (tuple, int, float))
    length_bound_type = isinstance(length_bound, (tuple, int, float))
    quality_thr_type = isinstance(quality_threshold, (int, float))
    if not (seqs_type and gc_bound_type and length_bound_type and quality_thr_type):
        raise ValueError('Your arguments are not suitable!')
    result_dict = seqs.copy()
    for value in seqs:
        fastq = seqs[value][0]
        if isinstance(gc_bound, (int, float)):
            gc_check = fr.gc_content(fastq) >= gc_bound
        else:
            gc_check = fr.gc_content(fastq) < gc_bound[0] or fr.gc_content(fastq) > gc_bound[1]
        if isinstance(length_bound, (int, float)):
            len_check = fr.seq_length(fastq) >= length_bound
        else:
            len_check = fr.seq_length(fastq) < length_bound[0] or fr.gc_content(fastq) > length_bound[1]
        quality_check = fr.quality_score(fastq) < quality_threshold
        if gc_check or len_check or quality_check:
            del result_dict[value]
    return result_dict


# Function for working with protein sequences
def protein_tools(*sequences: str, action: str):
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
    - "calculate_aa_freq": Calculate the frequency of each amino acid in a protein sequence.
    - "translate_protein_rna": Translate amino acid sequence to RNA, using random codons for each amino acid.
    - "three_letter_code": Convert one-letter amino acid sequence to three-letter coding.
    - "protein_mass": Calculate the molecular weight of the protein sequence.
    """

    output_sequence = []
    ACTION_LIST = {
        "get_pI": pt.get_pI,
        "calculate_aa_freq": pt.calculate_aa_freq,
        "translate_protein_rna": pt.translate_protein_rna,
        "three_letter_code": pt.three_letter_code,
        "protein_mass": pt.protein_mass,
    }
    # check for the size of sequence if only 1 is entered
    if len(sequences) == 1:
        if len(sequences[0]) == 0:
            raise ValueError("Your sequence is empty")
    # check for actions
    if action not in ACTION_LIST:
        raise ValueError(f"No such action: {action}")
    for seq in sequences:
        function = ACTION_LIST[action]
        output_sequence.append(function(seq))
    if len(output_sequence) <= 1:
        return output_sequence[0]
    return output_sequence


# Function for working with DNA/RNA sequences
def nucl_acid_tools(*sequences: str, action: str):
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
    output_seq_list = []
    FUNC_DICT = { # Function dictionary
        'transcribe': nat.transcribe,
        'reverse': nat.reverse,
        'complement': nat.complement,
        'reverse_complement': nat.reverse_complement,
        'make_binary': nat.make_binary
        }
    # check for the size of sequence if only 1 is entered
    if len(sequences) == 1:
        if len(sequences[0]) == 0:
            raise ValueError("Your sequence is empty")
    # check for actions
    if action not in FUNC_DICT:
        raise ValueError(f"No such action: {action}")
    for seq in sequences:
        function = FUNC_DICT[action]
        output_seq_list.append(function(seq))
    if len(output_seq_list) <= 1:
        return output_seq_list[0]
    return output_seq_list
