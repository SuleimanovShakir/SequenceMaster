# Import modules
import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """
    This function converts multiline FASTA to oneline FASTA
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
