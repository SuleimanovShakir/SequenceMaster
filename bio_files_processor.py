def convert_multiline_fasta_to_oneline(input_path: str, output_path: str = None):
    """
    This function converts multiline FASTA to oneline FASTA
    """
    with open(input_path, 'r') as file:
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
