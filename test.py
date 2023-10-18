# Import functions
from bio_files_processor import convert_multiline_fasta_to_oneline
from bio_files_processor import change_fasta_start_pos
from bio_files_processor import select_genes_from_gbk_to_fasta
from sequence_tools import fastq_filter

# Create input paths
INPUT_FASTA = 'example_data/example_multiline_fasta.fasta'
INPUT_ONELINE_FASTA = 'example_data/example_oneline_fasta.fasta'
INPUT_GBK = 'example_data/example_gbk.gbk'
INPUT_FASTQ = 'example_data/example_fastq.fastq'

# Run functions
convert_multiline_fasta_to_oneline(INPUT_FASTA, 'oneline_fasta')
change_fasta_start_pos(INPUT_ONELINE_FASTA, 3, 'shifted_fasta')
select_genes_from_gbk_to_fasta(INPUT_GBK, 'genes_of_interest', genes_of_interest=['mngA', 'araE', 'trpA'], n_before=4, n_after=5)
fastq_filter(INPUT_FASTQ, 'filtered_fastq', gc_bound=(40,60), length_bound=(0, 200), quality_threshold=25) 