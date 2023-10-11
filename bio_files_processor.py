# Import modules
import os


# Function to convert multiline fatsa to oneline
def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """
    This function converts multiline FASTA to oneline FASTA

    Arguments (positional):
    - input_path (str): full path to the file that you want to work with
    - output_filename (str): enter just a name of the file, don't add extention
    """
    # Chech if folder exist and create outout_path if not given
    input_folder = input_fasta.rsplit('/', 1)[0]
    input_name = input_fasta.rsplit('/', 1)[1]
    is_exist = os.path.exists(f'{input_folder}/results/')
    if not is_exist:
        os.makedirs(f'{input_folder}/results/')
    if output_fasta is None:
        output_path = f'{input_folder}/results/{input_name}'
    else:
        output_path = f'{input_folder}/results/{output_fasta}.fasta'
    with open(input_fasta, 'r') as file:
        lines = []
        for line in file:
            line = line.strip('')
            lines.append(line)
        output_file = []
    with open(output_path, 'w') as output_file:
        first_line = lines[0]
        output_file.write(first_line+'')
        for element in lines[1:]:
            if element.startswith('>'):
                output_file.write('\n'+element+'')
            else:
                element = element.strip('\n')
                output_file.write(element+'')
        print('File is written!')


# Function to change starting position in fasta sequence
def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: str):
    """
    This function shifts sequence to n positions from the start. But nothing is deleted!
    Part of sequence which is now goes before shift is located in the end of the new sequence

    Arguments (positional):
    - input_path (str): full path to the file that you want to work with
    - shift (int): the number of nucleotides by which the starting position in the file must be shifted
    - output_filename (str): enter just a name of the file, don't add extention
    """
    if not isinstance(shift, int):
        raise ValueError('Shift arguments has to be numeric')
    input_folder = input_fasta.rsplit('/', 1)[0]
    input_name = input_fasta.rsplit('/', 1)[1]
    is_exist = os.path.exists(f'{input_folder}/results/')
    if not is_exist:
        os.makedirs(f'{input_folder}/results/')
    if output_fasta is None:
        output_path = f'{input_folder}/results/{input_name}'
    else:
        output_path = f'{input_folder}/results/{output_fasta}.fasta'
    with open(input_fasta, 'r') as file:
        lines = []
        for line in file:
            line = line.strip('')
            lines.append(line)
        output_file = []
    with open(output_path, 'w') as output_file:
        for element in lines:
            if element.startswith('>'):
                output_file.write(element+'\n')
            else:
                element = element.strip('\n')
                shifted_element = element[shift:]+element[:shift]
                output_file.write(shifted_element+'\n')
        print('File is written!')

GBK = '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Python/HW6_Files/example_data/example_gbk.gbk'

with open(GBK, 'r') as file:
    lines = file.readlines()
    print(lines[40:60])