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


def fastq_to_dict(input_path: str) -> dict:
    """
    This function converts given FASTQ file to dictionary

    Arguments:
        input_path = PATH to your file on computer

    Output:
        dictionary containing FASTQ file as: key = ID, value = tuple(sequence, quality)
    """
    with open(input_path, 'r') as file:
        lines = []
        for line in file:
            line = line.strip('\n')
            lines.append(line)
    output_dict = {}
    size = len(lines)
    fastq = lines.copy()
    for value in range(0, size, 4):
        id = fastq[value]
        seq = fastq[value+1]
        quality = fastq[value+3]
        output_dict[id] = (seq, quality)
    return output_dict


def dict_to_fastq(input_dict: dict, output_path: str = None):
    """
    This function converts dictionary to FASTQ

    Arguments:
        output_path = PATH to your file in which you want to save your processed file

    Output:
        dictionary containing FASTQ file as: key = ID, value = tuple(sequence, quality)
    """
    final_list = []
    for i in range(0, len(input_dict)):
        keys = list(input_dict.keys())[i]
        final_list.append(keys)
        third_line = keys.replace('@', '+')
        values = list(input_dict.values())[i]
        final_list.append(values[0])
        final_list.append(third_line)
        final_list.append(values[1])
    with open(output_path, 'w') as output_file:
        for line in final_list:
            output_file.write(line+'\n')
        print('File is written!')
