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
def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: str = None):
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


# Function to read GBK file to FASTA
def select_genes_from_gbk_to_fasta(input_gbk: str, output_fasta: str = None,
                                   *, genes_of_interest: list,
                                   n_before: int = 1, n_after: int = 1):
    """
    This function uses gbk database files to extract FASTA sequence of proteins.
    Uses only those genes that have a name in gbk database. User enters name of 
    interested genes and how much genes before and after this gene of interest
    he want to extract with proteins sequences for these genes.

    Arguments (positional):
    - input_gbk (str): full path to the file that you want to work with
    - output_fasta (str): enter just a name of the file, don't add extention!

    Arguments (keyword):
    - genes_of_interest (list): list of gene names user is interested in
    - n_before (int): how many of genes before particular genes you want to extract
    - n_after (int): how many of genes after particular genes you want to extract
    """
    # Check if arguments are in rigth type
    if not (isinstance(genes_of_interest, list) and isinstance(n_before, int) and isinstance(n_after, int)):
        raise ValueError('Your arguments are not in suitable type')
    # Check if genes list is empty
    if len(genes_of_interest) == 0:
        raise ValueError('You have entered an empty list of genes')
    # Chech if PATH for input file is given
    if input_gbk is None:
        raise ValueError("You didn't enter any PATH to file")
    # Chech if folder exist and create outout_path if not given
    input_folder = input_gbk.rsplit('/', 1)[0]
    input_name = input_gbk.rsplit('/', 1)[1]
    is_exist = os.path.exists(f'{input_folder}/gbk_fasta_resuls/')
    if not is_exist:
        os.makedirs(f'{input_folder}/gbk_fasta_resuls/')
    if output_fasta is None:
        output_path = f'{input_folder}/gbk_fasta_resuls/{input_name}.fasta'
    else:
        output_path = f'{input_folder}/gbk_fasta_resuls/{output_fasta}.fasta'
    # Reading GBK
    with open(GBK, 'r') as file:
        # Convert GBK to list
        lines = []
        for line in file:
            line = line.strip(' ')
            line = line.strip('\n')
            lines.append(line)
    # Create list for individual CDS. It will be lists
    # of [CDS, gene, translation] in total list
    individual_cds = []
    for element in lines:
        # Looking for block starting with CDS
        if element.startswith('CDS'):
            cds_index = lines.index(element)
            element = element.replace('             ', ';')
            element = element.replace('..', ',')
            cds_gene_translation = []
            cds_gene_translation.append(element)
            # Specially a little bit more if '/=gene' would be further by accident
            for i in range(cds_index, cds_index+5):
                if lines[i].startswith('/gene='):
                    gene = lines[i]
                    break
                else:
                    gene = element
            cds_gene_translation.append(gene)
            # Specially a little bit more if '/=translation' would be further by accident
            for i in range(cds_index, cds_index+12): 
                if lines[i].startswith('/translation='):
                    seq = lines[i]
            cds_gene_translation.append(seq)
            individual_cds.append(cds_gene_translation)
    # Create output list. It will consist only from gene:translation pair
    output_list = []
    for cds in individual_cds:
        for sample in genes_of_interest:
            gene = f'/gene="{sample}"'
            # Check if gene of interest in CDS
            if gene == cds[1]:
                # Define defore and after in indices
                index_cds = individual_cds.index(cds)
                before = (index_cds - n_before) if (index_cds - n_before)>= 0 else 0
                after = (index_cds + n_after) if (index_cds + n_after) <= len(individual_cds) else len(individual_cds)
                area_of_interest = individual_cds[before: after+1]
                # Add all gene:seq pairs to output_list
                for element in area_of_interest:
                    gene = element[1]
                    sequence = element[2]
                    output_list.append(gene)
                    output_list.append(sequence)
    # Write the final list to the FASTA file
    with open(output_path, 'w') as output_file:
        for element in output_list:
            if element.startswith('/gene='):
                output_file.write(element.lstrip('/gene=').replace('"', '')+'\n')
            if element.startswith('/translation='):
                output_file.write(element.strip('/translation=').replace('"', '')+'\n')
            if element.startswith('CDS'):
                output_file.write(element+'\n')
    print('Your file is written')

