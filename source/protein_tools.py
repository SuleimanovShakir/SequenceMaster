# importing necessary modules
import source.protein_dict as pd
from random import choice


AMINO_LETTERS = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]


# Function to determine is the sequence is a protein or not
def is_protein(seq: str) -> bool:
    """
    This function checks if the sequence is a protein or not

    Arguments:
        seq (str): A sequence of aminoacids

    Output:
        returns True or False
    """
    unique_chars = set(seq)
    aminoacids = set(pd.AA_MONOISOTOPIC_MASS_DICT.keys())
    return unique_chars <= aminoacids


# Function to get pI for each aa
def get_pI(
    sequence: str,
    pI_values: dict = None,
) -> str:
    """
    Gives isoelectric point value for each aminoacid individually

    Args:
    - sequence (str): sequence for which to calculate isoelectric point
    - pI_values (dict): acid dissociation constants for each aminoacid

    Return:
    - str: string, which contains:
            - an original sequence,
            - list of tuple pairs of aminoacid and corresponding isoelectric point,
    """

    if pI_values is None:
        # Default pKa_values if not provided
        pI_values = pd.AA_pI

    aminoacid_pIs = []

    # Calculate pI for each amino acid in the sequence while preserving case
    analysed_aa = []
    for aa in sequence:
        aa_upper = aa.upper()
        if aa_upper not in analysed_aa:
            if aa_upper in pI_values:
                pI = pI_values[aa_upper]
                analysed_aa.append(aa_upper)
                if aa.isupper():
                    aminoacid_pIs.append((aa_upper, pI))
                else:
                    aminoacid_pIs.append((aa, pI))
        else:
            continue

    return f"Sequence: {sequence}. Isoelectric point of each aminoacid: {aminoacid_pIs}"


# Function to build scoring matrix for needleman_wunsch function
def build_scoring_matrix(
    match_score: int,
    mismatch_score: int,
    amino_acid_alphabet: str = None,
) -> dict:
    """
    Build a default scoring matrix, if not provided in needleman-wunsch function parameter

    Args:
    - match_score (int): integer value of a matching score of aminoacids
    - mismatch_score (int): integer value of a mismatching score of aminoacids
    - amino_acid_alphabet (str): upper case amino acid alphabet

    Returns:
    - a dictionary of dictionaries representing a scoring matrix for aminoacids paris. Key of a dictionary is an aminoacid and its value is a dictionary of scores
    """

    if amino_acid_alphabet is None:
        # Default pKa_values if not provided
        amino_acid_alphabet = "ACDEFGHIKLMNPQRSTVWY"

    scoring_matrix = {}

    for aa1 in amino_acid_alphabet:
        scoring_matrix[aa1] = {}
        for aa2 in amino_acid_alphabet:
            scoring_matrix[aa1][aa2] = (
                match_score if aa1.upper() == aa2.upper() else mismatch_score
            )

    return scoring_matrix


# Function to perform alignment based on needleman_wunsch algorithm
def needleman_wunsch(
    seq1: str,
    seq2: str,
    scoring_matrix: dict = None,
    gap_penalty: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1,
) -> str:
    """
    Uses Needleman-Wunsch algorithm to make a global alignment of two sequences.

    Args:
    - seq1 (str): first aminoacid sequence for alignment
    - seq2 (str): second aminoacid sequence for alignment
    - scoring_matrix (dict): A dictionary representing a scoring matrix for amino acid pairs
      If not provided, a default scoring_matrix is generated based on match and mismatch scores
    - gap_penalty (int): integer va;ue of a penalty score for introducing a gap in the alignment
    - match_score (int): integer value of a matching score for matching aminoacids
    - mismatch_score (int): integer value of a mismatching score for mismatched aminoacids

    Returns:
    - string: a string containing the aligned sequences (str), the aligned score (int)
    """
    if scoring_matrix is None:
        # Default scoring matrix if not provided
        scoring_matrix = build_scoring_matrix(match_score, mismatch_score)

    seq1_upper = seq1.upper()  # Convert seq1 to uppercase
    seq2_upper = seq2.upper()  # Convert seq2 to uppercase

    m, n = len(seq1_upper), len(seq2_upper)

    # Initialize matrices
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[""] * (n + 1) for _ in range(m + 1)]

    # Fill in the scoring matrix and traceback matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i - 1][j - 1] + scoring_matrix.get(seq1_upper[i - 1], {}).get(
                seq2_upper[j - 1], mismatch_score
            )
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty

            dp[i][j] = max(match, delete, insert)

            if dp[i][j] == match:
                traceback[i][j] = "D"  # Diagonal (indicates a match/mismatch)
            elif dp[i][j] == delete:
                traceback[i][j] = "U"  # Up (indicates a gap in seq2)
            else:
                traceback[i][j] = "L"  # Left (indicates a gap in seq1)

    # Traceback to find the aligned sequences while preserving case
    aligned_seq1, aligned_seq2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == "D":
            # check original case of amionacid in seq1
            if seq1[i - 1].isupper():
                aligned_seq1.append(seq1_upper[i - 1])
            else:
                aligned_seq1.append(seq1[i - 1])

            # check original case of amionacid in seq2
            if seq2[j - 1].isupper():
                aligned_seq2.append(seq2_upper[j - 1])
            else:
                aligned_seq2.append(seq2[j - 1])

            i -= 1
            j -= 1
        elif i > 0 and traceback[i][j] == "U":
            # check original case of amionacid in seq1
            if seq1[i - 1].isupper():
                aligned_seq1.append(seq1_upper[i - 1])
            else:
                aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")

            i -= 1
        else:
            aligned_seq1.append("-")
            # check original case of amionacid in seq2
            if seq2[j - 1].isupper():
                aligned_seq2.append(seq2_upper[j - 1])
            else:
                aligned_seq2.append(seq2[j - 1])

            j -= 1

    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return f"{aligned_seq1}, {aligned_seq2}, final score: {dp[m][n]}"


# Function to calculate frequency of unique aminoacid in the sequence
def calculate_aa_freq(sequences: str) -> dict:
    """
    Calculates the frequency of each amino acid in a protein sequence or sequences.

    :param sequences: protein sequence or sequences
    :type sequences: str or list of str
    :return: dictionary with the frequency of each amino acid
    :rtype: dict
    """

    # Creating a dictionary with aminoacid frequencies:
    amino_acid_frequency = {}

    for amino_acid in sequences:
        # If the aminoacid has been already in:
        if amino_acid in amino_acid_frequency:
            amino_acid_frequency[amino_acid] += 1
        # If the aminoacid hasn't been already in:
        else:
            amino_acid_frequency[amino_acid] = 1

    return amino_acid_frequency


# Function to convert one-letter protein sequence to three-letter protein sequence
def three_letter_code(seq: str) -> str:
    """
    This function takes one letter aminoacids sequence and convert's it to three leter coding

    Arguments:
        seq (str): A sequence of aminoacids

    Output:
        same sequence but in three-letter coding
    """
    seq = seq.upper()
    if is_protein(seq):
        sequence = "".join(pd.AA_ONE_TO_THREE_LETTER.get(aa) for aa in seq)
        return sequence[:-1]
    else:
        raise ValueError("Sequence is not a protein, input should be protein")


# Function to calculate protein mass
def protein_mass(seq: str) -> float:
    """
    This function takes aminoacids sequence and counts it's summary molecular weight using monoisotopic masses

    Arguments:
        seq (str): A sequence of aminoacids

    Output:
        returns molecular weight
    """
    seq = seq.upper()
    if is_protein(seq):
        mass = sum(pd.AA_MONOISOTOPIC_MASS_DICT.get(aa) for aa in seq)
        return mass
    else:
        raise ValueError("Sequence is not a protein, input should be protein")


# Function to translate Protein to RNA
def translate_protein_rna(seq: str) -> str:
    """
    This function takes  aminoacid sequence and translates in to the RNA.
    As most of the aminoacids are coded with several different codons,
    this function will take a random codon of the set for such aminoacids.

    Arguments:
        seq (str): A sequence of RNA molecule

    Output:
        returns sequence of aminoacids
    """
    seq = seq.upper()
    if is_protein(seq):
        rna = ""
        for aa in seq:
            codon = choice(pd.AA_CODON_DICT.get(aa))
            rna += codon
        return rna
    else:
        raise ValueError("Sequence is not a protein, input should be a protein")
