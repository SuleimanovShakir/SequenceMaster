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
    This function uses gbk database files to extract FASTA sequence of gene
    Uses only those genes that have a name in gbk
    """
    # Check if arguments are in rigth type
    if not (isinstance(genes_of_interest, list) and isinstance(n_before, int) and isinstance(n_after, int)):
        raise ValueError('Your arguments are not in suitable type')
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
        # Create file with all genes
        all_genes = list(filter(lambda element: '/gene=' in element, lines))
        # Create dictionary {key = gene, value = translation}
        gbk_dict = {}
        for element in lines:
            if element.startswith('/gene='):
                if not element.startswith(all_genes[-1]):
                    for i in range(lines.index(element), lines.index(element)+9):
                        if lines[i].startswith('/translation='):
                            seq = lines[i]
                            gbk_dict[element.strip('/')] = seq.strip('/')
                else:
                    for i in range(lines.index(element), lines.index(element)+9):
                        if lines[i].startswith('/translation='):
                            seq = lines[i]
                            gbk_dict[element.strip('/')] = seq.strip('/')
                    break
    # Create final list
    keys_gbk_dict = list(gbk_dict.keys())
    final_list = []
    # Iterate by each gene in gene of interest which is given by user
    for gene in genes_of_interest:
        gene = f'gene="{gene}"'
        if gene in keys_gbk_dict:
            index_gene = keys_gbk_dict.index(gene)
            before = (index_gene - n_before) if (index_gene - n_before)>= 0 else 0
            after = (index_gene + n_after) if (index_gene + n_after) <= len(keys_gbk_dict) else len(keys_gbk_dict)
            list_of_int_genes = keys_gbk_dict[before: after+1]
            # Create local list for genes around one concrete gene from the list 
            # of gene of interest which is given by user
            gene_seq_list = []
            # Add pair (gene, seq) to local list of genes
            for name in list_of_int_genes:
                sequence = gbk_dict[name]
                gene_seq_list.append(name)
                gene_seq_list.append(sequence)
        final_list.append(gene_seq_list)
    # Write the final list to the FASTA file
    with open(output_path, 'w') as output_file:
        for element in lines:
            if element.startswith('>'):
                output_file.write(element+'\n')
        for pairs in final_list:
            for sample in pairs:
                if sample.startswith('gene='):
                    output_file.write('>'+sample+'\n')
                if sample.startswith('translation='):
                    sample = sample.strip('translation="')
                    output_file.write(sample+'\n')
    print('Your file is written')
