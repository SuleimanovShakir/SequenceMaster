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

    if action not in FUNC_DICT:
        raise ValueError(f"No such action: {action}")
    for seq in sequences:
        function = FUNC_DICT[action]
        output_seq_list.append(function(seq))
    if len(output_seq_list) <= 1:
        return output_seq_list[0]
    return output_seq_list


EXAMPLE_FASTQ = {
    # 'name' : ('sequence', 'quality')
    '@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
    '@SRX079804:1:SRR292678:1:1101:47176:47176': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'),
    '@SRX079804:1:SRR292678:1:1101:149302:149302': ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA', '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A'),
    '@SRX079804:1:SRR292678:1:1101:170868:170868': ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT', 'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB'),
    '@SRX079804:1:SRR292678:1:1101:171075:171075': ('CATTATAGTAATACGGAAGATGACTTGCTGTTATCATTACAGCTCCATCGCATGAATAATTCTCTAATATAGTTGTCAT', 'HGHHHHGFHHHHFHHEHHHHFGEHFGFGGGHHEEGHHEEHBHHFGDDECEGGGEFGF<FGGIIGEBGDFFFGFFGGFGF'),
    '@SRX079804:1:SRR292678:1:1101:175500:175500': ('GACGCCGTGGCTGCACTATTTGAGGCACCTGTCCTCGAAGGGAAGTTCATCTCGACGCGTGTCACTATGACATGAATG', 'GGGGGFFCFEEEFFDGFBGGGA5DG@5DDCBDDE=GFADDFF5BE49<<<BDD?CE<A<8:59;@C.C9CECBAC=DE'),
    '@SRX079804:1:SRR292678:1:1101:190136:190136': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT', 'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE'),
    '@SRX079804:1:SRR292678:1:1101:190845:190845': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC', 'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE'),
    '@SRX079804:1:SRR292678:1:1101:198993:198993': ('AGTTATTTATGCATCATTCTCATGTATGAGCCAACAAGATAGTACAAGTTTTATTGCTATGAGTTCAGTACAACA', '<<<=;@B??@<>@><48876EADEG6B<A@*;398@.=BB<7:>.BB@.?+98204<:<>@?A=@EFEFFFEEFB'),
    '@SRX079804:1:SRR292678:1:1101:204480:204480': ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG', '<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')
    }  
