import pytest

from sequence_tools import RNAsequence, AminoAcidSequence, DNAsequence, GenscanOutput, run_genscan, fastq_filter


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


def test_fastq_filter_writing_file():
    '''
    Function to test the fastq_filter function's writing file.
    '''
    fastq_filter('example_data/example_fastq.fastq', 'filtered_fastq', gc_bound=(40,60), length_bound=(88, 200), quality_threshold=15)
    with open('example_data/fastq_filtrator_results/filtered_fastq.fastq', 'r') as file:
        result = file.read()
    target = '@SRX079804:1:SRR292678:1:1101:24563:24563 1:N:0:1 BH:failed\nATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG\n+\nBFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D\n'
    assert target == result