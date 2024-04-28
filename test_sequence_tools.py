import pytest

from sequence_tools import RNAsequence, AminoAcidSequence, DNAsequence, GenscanOutput, run_genscan


def test_rna_sequence_gc_content():
    '''
    Function to test the gc_content function.
    '''
    target = 0.4
    result = RNAsequence('AUGCGCAUUU').gc_content()
    assert target == result


def test_rna_sequence_slice():
    '''
    Function to test the slice function.
    '''
    target = RNAsequence('AGUCGA')
    result = RNAsequence('AUAGUCGAUU').slice(2,8)
    assert target == result


def test_amino_acid_sequence_mass():
    '''
    Function to test the protein_mass function.
    '''
    target = 946.5
    result = AminoAcidSequence('KSTARAGCTA').protein_mass()
    assert target == pytest.approx(result, abs=0.1)


def test_dna_sequence_transcribe():
    '''
    Function to test the transcribe function.
    '''
    target = RNAsequence('UACGUACG')
    result = DNAsequence('ATGCATGC').transcribe()
    assert target == result


def test_wrong_sequence_error():
    '''
    Function to test the custom wrong sequence error.
    '''
    inp = 'ATGCU'
    with pytest.raises(ValueError):
        DNAsequence(inp)


def test_sequence_length():
    '''
    Function to test the length function.
    '''
    target = 8
    result = len(RNAsequence('AUGCGUUU'))
    assert target == result


def test_run_genscan_output_type():
    '''
    Function to test the run_genscan function's output type.
    '''
    target = type(GenscanOutput(status=None, cds_list=None, intron_list=None, exon_list=None))
    result = type(run_genscan('example_data/transferrin.fasta'))
    assert target == result


def test_run_genscan_output_status():
    '''
    Function to test the run_genscan function's output status.
    '''
    target = 200
    result = run_genscan('example_data/transferrin.fasta').status
    assert target == result
