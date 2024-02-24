import os

def folder_parser(input: str = None, output: str = None,
                  folder_name: str = None, format: str = None) -> str:
    """
    This function takes as input a folder path, processes folder name
    and creates output_path.

    Arguments (positional):
    - input (str): The path to the input folder containing the files
      to be processed. If not provided, the current working directory is used.
    - output (str): Name of the output file.
    - folder_name (str): A custom name for the folder to be created
      for the output files.
    - format (str): The format for the saved file. 

    Output:
    - str: output path
    """
    input_folder = input.rsplit('/', 1)[0]
    input_name = input.rsplit('/', 1)[1]
    is_exist = os.path.exists(f'{input_folder}/{folder_name}/')
    if not is_exist:
        os.makedirs(f'{input_folder}/{folder_name}/')
    if output is None:
        output_path = f'{input_folder}/{folder_name}/{input_name}'
    else:
        output_path = f'{input_folder}/{folder_name}/{output}.{format}'
    return output_path
