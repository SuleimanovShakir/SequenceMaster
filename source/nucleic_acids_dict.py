# Create dictionary for transcription function
TRANSCRIPTION_DICT = {
    'A': 'A',
    'G': 'G',
    'T': 'U',
    'C': 'C',
    'a': 'a',
    'g': 'g',
    't': 'u',
    'c': 'c',
}


# Create dictionary for complimentarity function
COMPLMENTARITY_DICT = {
    'A': 'T',
    'G': 'C',
    'T': 'A',
    'C': 'G',
    'a': 't',
    'g': 'c',
    't': 'a',
    'c': 'g',
}


# Binary code dictionary
BINARY_DICT = {
    'A': '0,0,',
    'T': '0,1,',
    'G': '1,0,',
    'C': '1,1,',
    'U': '0,1,',
}
