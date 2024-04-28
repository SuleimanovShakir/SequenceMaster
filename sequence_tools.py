# Importing modules
import os
import sys
import datetime
import io
from dataclasses import dataclass
from typing import Union

import requests
from dotenv import load_dotenv
from bs4 import BeautifulSoup
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def fastq_filter(input_path: str = None, output_filename: str = None, *,
                 gc_bound: Union[tuple, int, float] = (0, 100),
                 length_bound: Union[tuple, int, float] = (0, 2**32),
                 quality_threshold: Union[int, float] = 0) -> None:
    """
    This function work with FASTQ files and filters them by
    GC content, length and Q-score.

    Arguments (positional):
    - input_path (str): full path to the file that you want to work with
    - output_filename (str): enter just a name of the file, don't add extention

    Arguments (keyword):
    - gc_bound (tuple, int, float): tuple of required range of GC percentage (inclusive),
    num or float if only higher border of the range is needed (exclusive).
    - length_bound (tuple, int, float): tuple of required range of sequences length (inclusive),
    num or float if only higher border of the range is needed (exclusive).
    - quality_threshold (int): int of lowest level of Q-score (inclusive).

    Output:
    - list of BioSeq records. Write file to .fastq
    """
    # Chech if PATH for input file is given
    if input_path is None:
        raise ValueError("You didn't enter any PATH to file")
    # Chech if folder exist and create outout_path if not given
    input_folder = input_path.rsplit('/', 1)[0]
    input_name = input_path.rsplit('/', 1)[1]
    is_exist = os.path.exists(f'{input_folder}/fastq_filtrator_resuls/')
    if not is_exist:
        os.makedirs(f'{input_folder}/fastq_filtrator_resuls/')
    if output_filename is None:
        output_path = f'{input_folder}/fastq_filtrator_resuls/{input_name}'
    else:
        output_path = f'{input_folder}/fastq_filtrator_resuls/{output_filename}.fastq'
    # Create dict from FASTQ
    seqs = list(SeqIO.parse(input_path, "fastq"))
    # Check that this dict is not empty
    if len(seqs) <= 0:
        raise ValueError('There are no fastq sequences')
    # Check if all given argumets have relevant type
    gc_bound_type = isinstance(gc_bound, (tuple, int, float))
    length_bound_type = isinstance(length_bound, (tuple, int, float))
    quality_thr_type = isinstance(quality_threshold, (int, float))
    if not (gc_bound_type and length_bound_type and quality_thr_type):
        raise ValueError('Your arguments are not suitable!')
    # Create filtered list
    filtered_fastq = []

    for line in seqs:
        if isinstance(gc_bound, (int, float)):
            gc_check = gc_fraction(line)*100 >= gc_bound
        else:
            gc_check = gc_fraction(line)*100 < gc_bound[0] or gc_fraction(line)*100 > gc_bound[1]
        if isinstance(length_bound, (int, float)):
            len_check = len(line) >= length_bound
        else:
            len_check = len(line) < length_bound[0] or len(line) > length_bound[1]
        quality_check = sum(line.letter_annotations['phred_quality'])/len(line.letter_annotations['phred_quality']) < quality_threshold
        if not (gc_check or len_check or quality_check):
            filtered_fastq.append(line)

    # Write  filtered data into new .fastq file
    SeqIO.write(filtered_fastq, output_path, 'fastq')


class WrongSequence(ValueError):
    '''
    Class for wrong sequence. Custom error.
    '''
    pass


class BiologicalSequence(str):
    '''
    Class for biological sequences.

    Attributes:
    ----------
        sequence (str): sequence of biological sequence

    Methods:
    ----------
    slice():
        function to take slices in sequences

    alphabet_checking():
        function to check if the sequence consist of only allowed nucleotides
    '''
    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def slice(self, start_index, stop_index):
        '''
        Function to take slices in sequences
        '''
        return self.sequence[start_index:stop_index]

    def alphabet_checking(self):
        '''
        Function to check if the sequence consist of only allowed nucleotides
        '''
        if not set(self.sequence) <= set(type(self).dictionary.keys()):
            raise WrongSequence('Wrong sequence')
        return True


class NucleicAcidSequence(BiologicalSequence):
    '''
    Class to represent nucleic acid sequences

    Attributes:
    ----------
    sequence (str): sequence of nucleic acid sequences

    Methods:
    -------
    complement():
        returns the complement of the sequence

    gc_content():
        returns the GC content of the sequence  

    '''
    def __init__(self, sequence):
        super().__init__(sequence)
        if not self.alphabet_checking():
            del self.sequence
            raise WrongSequence('You have entered a wrong sequence')
        self.gc_cont = None

    def complement(self):
        '''
        Return the complement of the sequence

        Returns:
        -------
        NucleicAcidSequence: complement of the sequence
        '''
        return type(self)(''.join([type(self).dictionary[i] for i in self.sequence]))

    def gc_content(self):
        '''
        Returns the GC content of the sequence

        Returns:
        -------
        float: GC content of the sequence
        '''
        gc = 0
        for nucleotide in self.sequence:
            if nucleotide in set('CGcg'):
                gc += 1
        self.gc_cont = round(gc / len(self.sequence), 3)
        return self.gc_cont


class DNAsequence(NucleicAcidSequence):
    '''
    Class for DNA sequences

    Attributes:
    -----------
    sequence: str
        DNA sequence

    Methods:
    --------
    transcribe():
        Returns the RNA sequence of the DNA sequence
    '''
    dictionary = {
        'A': 'T',
        'G': 'C',
        'T': 'A',
        'C': 'G',
        'a': 't',
        'g': 'c',
        't': 'a',
        'c': 'g',
}
    trans_dict = {
        'A': 'U',
        'G': 'C',
        'T': 'A',
        'C': 'G',
        'a': 'u',
        'g': 'c',
        't': 'a',
        'c': 'g',
}
    def transcribe(self):
        return RNAsequence(''.join([self.trans_dict[i] for i in self.sequence]))


class RNAsequence(NucleicAcidSequence):
    '''
    Class for RNA sequences
    Attributes:
        sequence: string of RNA sequence

    '''
    dictionary = {
        'A': 'U',
        'G': 'C',
        'U': 'A',
        'C': 'G',
        'a': 'u',
        'g': 'c',
        'u': 'a',
        'c': 'g',
}


class AminoAcidSequence(BiologicalSequence):
    '''
    Class for aminoacid sequences
    Attributes:
    ----------
        sequence: string of aminoacid sequence

    Methods:
    ----------
    protein_mass():
        function to calculate the mass of the protein

    '''
    dictionary = {
        "A": 71.03711,
        "C": 103.00919,
        "D": 115.02694,
        "E": 129.04259,
        "F": 147.06841,
        "G": 57.02146,
        "H": 137.05891,
        "I": 113.08406,
        "K": 128.09496,
        "L": 113.08406,
        "M": 131.04049,
        "N": 114.04293,
        "P": 97.05276,
        "Q": 128.05858,
        "R": 156.10111,
        "S": 87.03203,
        "T": 101.04768,
        "V": 99.06841,
        "W": 186.07931,
        "Y": 163.06333,
    }

    def __init__(self, sequence):
        super().__init__(sequence)
        if not self.alphabet_checking():
            del self.sequence
            raise WrongSequence('You have entered a wrong sequence')

    def protein_mass(self):
        '''
        Fucntion to calculate the mass of the protein.

        return (float): The mass of the protein.
        '''
        mass = sum(self.dictionary.get(aa) for aa in self.sequence)
        return mass


#Add TOKEN
load_dotenv()
TOKEN = os.getenv('TG_API_TOKEN')


def reformat_time(time_delta):
    '''
    Function to rewrite the time delta in a more readable format.

    :param time_delta: The time delta to be rewritten.
    '''
    days = time_delta.days
    if days > 1:
        sec = time_delta.seconds
        hours = sec // 3600
        sec = sec - (hours * 3600)
        minutes = sec // 60
        seconds = sec - (minutes * 60)
        return f'{days} days, {int(hours):02}:{int(minutes):02}:{int(seconds):02}'
    return str(time_delta)


def telegram_logger(chat_id):
    '''
    A decorator that sends a log message to a specific Telegram bot.

    :param chat_id: The ID of the Telegram chat.
    '''
    def actual_decorator(func):
        def wrapper(*args, **kwargs):
            #Redirect stdout and stderror to IO files
            temp_stdout = io.StringIO()
            temp_stderr = io.StringIO()
            sys.stdout = temp_stdout
            sys.stderr = temp_stderr
            
            #Try to run the function
            try:
                #Start time counting
                start_time = datetime.datetime.now()

                #Call the function
                func(*args, **kwargs)

                #End time counting
                end_time = datetime.datetime.now()

                #Calculate time taken
                time_taken = end_time - start_time

                #Create a message
                message = f'{func.__name__} has finished in {reformat_time(time_taken)}'

                #Send message to telegram
                requests.get(f'https://api.telegram.org/bot{TOKEN}/sendMessage?text={message}&chat_id={chat_id}')

                #Create a file
                file = {}
                file["document"] = (f"{func.__name__}_success.log", temp_stdout.getvalue() + temp_stderr.getvalue())

                #Send file to telegram
                requests.post(f'https://api.telegram.org/bot{TOKEN}/sendDocument', params={'chat_id': chat_id}, files=file)
            
            #Catch any possible error
            except Exception as exception:
                
                print(exception, file=temp_stderr)
                #Create a message
                message = f'Lol, your {func.__name__} has failed to run, because of the {type(exception).__name__}: {exception}'

                #Send message to telegram
                requests.get(f'https://api.telegram.org/bot{TOKEN}/sendMessage?text={message}&chat_id={chat_id}')

                #Create a file
                file = {}
                file["document"] = (f"{func.__name__}_failed.log", temp_stdout.getvalue() + temp_stderr.getvalue())

                #Send file to telegram
                requests.post(f'https://api.telegram.org/bot{TOKEN}/sendDocument', params={'chat_id': chat_id}, files=file)
            
            #Reassign stdout and stderr
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
        return wrapper
    return actual_decorator


@dataclass
class GenscanOutput:
    '''
    Dataclass that stores the output of the genscan parser
    '''
    status: int
    cds_list: list
    intron_list: list
    exon_list: list

    def __repr__(self):
        return f'GenscanOutput run with status code {self.status}. {len(self.cds_list)} cds were found. {len(self.intron_list)} intorns were found. {len(self.exon_list)} exons were found'


def run_genscan(sequence=None, sequence_file=None, organism="Vertebrate", exon_cutoff=1.00, sequence_name=""):
    '''
    Run the Genscan gene prediction tool and parse the output.

    This function sends a request to the Genscan web server with the given parameters,
    parses the output using BeautifulSoup, and returns a GenscanOutput object containing
    the predicted CDS, intron, and exon sequences.

    :param sequence: The DNA sequence to analyze (either sequence or sequence_file must be provided).
    :param sequence_file: A file containing the DNA sequence to analyze (either sequence or sequence_file must be provided).
    :param organism: The taxa to use for the analysis (default: "Vertebrate").
    :param exon_cutoff: The minimum exon probability (default: 1.00).
    :param sequence_name: The name of the sequence (default: "").
    :return: A GenscanOutput object containing the predicted request status, CDS, intron, and exon sequences.
    '''

    #URL of the genscan web server
    genscan_url = 'http://argonaute.mit.edu/cgi-bin/genscanw_py.cgi'

    # file={'files': open(sequence_file,'r')}

    # Read sequence from file if provided
    if sequence_file:
        with open(sequence_file, "r") as file:
            sequence = file.read()

    # Payload - arguments of the post request
    payload = {
        '-o': organism,
        '-e': exon_cutoff,
        '-n': sequence_name,
        '-p': "Predicted peptides only",
        '-s': sequence,
    }

    # Send request
    genscan_run = requests.post(genscan_url, data=payload)

    # Create soup
    soup = BeautifulSoup(genscan_run.text, 'html.parser')

    # Find only part containing information
    record = [i for i in str(soup.find_all('pre')).split('\n') if i != '']

    # Find start and end of exon records
    exon_start = record.index('Predicted genes/exons:')
    exon_end = record.index('Suboptimal exons with probability &gt; 1.000')

    # Create empty dictionary for exons
    exons = {}

    # Fill dict with exons
    for line in record[exon_start: exon_end]:
        if line.startswith(' '):
            splitted_line = [i for i in line.split(' ') if i != '']
            exons[splitted_line[0]] = (splitted_line[1], int(splitted_line[3]), int(splitted_line[4]))

    # Create empty dictionary for introns
    introns = {}

    # IDs of exons
    exon_ids = list(exons.keys())

    # Fill dict with introns
    for i in range(1, len(exons) - 1):
        introns[f'{exon_ids[i-1]}:{exon_ids[i]}'] = ('Intronic region', exons[exon_ids[i-1]][2] + 1, exons[exon_ids[i]][1] - 1)

    # Find start of CDS records
    cds_start = record.index('Predicted peptide sequence(s):')

    # Define names of ids and indices
    cds_names = []
    for line in record[cds_start:-1]:
        if line.startswith('&gt'):
            name = line.split('|')[1]
            ids = record.index(line)
            cds_names.append((name, ids))

    # Create last technical CDS
    cds_names.append(('Last_technical_CDS', len(record)-1))

    # Create empty CDS dictionary
    cds = {}

    # Create a dictionary with CDS names and sequences
    for i in range(0, len(cds_names)-1):
        cds_rec = ''
        for k in range(cds_names[i][1]+1, cds_names[i+1][1]):
            #print(record[k+1])
            cds_rec += ''.join(record[k])
        cds[cds_names[i][0]] = cds_rec

    return GenscanOutput(genscan_run.status_code, cds, introns, exons)